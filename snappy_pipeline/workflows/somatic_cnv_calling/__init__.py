# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_cnv_calling`` step

This step allows for the detection of CNV events for cancer samples from targeted sequenced (e.g.,
exomes or large panels) or whole genome sequencing.
The wrapped tools start from the aligned reads (thus off ``ngs_mapping``) and generate CNV calls for somatic variants.

The wrapped tools implement different strategies.  Some work "reference free" and just use the
somatic BAM files for their input, some work in "matched cancer normal mode" and need the cancer
and normal BAM files, others again cancer BAM files, and additionally a
set of non-cancer BAM files for their background (the panel of normals).

Some tools may also use germline & somatic variants to estimate allele-specific copy number changes,
and resolve loss-of-heterozygocity. In this case, the small variants need to be computed separately from the ``somatic_variant_calling`` step.

Finally, some tools can use external estimation of tumor purity and ploidy.
This estimation can be either be provided in the sample sheet, or computed from the sequencing data by tools.

==========
Step Input
==========

Gene somatic CNV calling for targeted sequencing starts off the aligned reads, i.e.,
``ngs_mapping``.

Tools that use panel of normals can obtain their input in two different ways:

- A static file, from another cohort or from public datasets.
  In this case, the user is responsible to make sure that the data & methods used to create the panel are compatible to the cohort's.
- The ``panel_of_normals`` step.
  The panel will be created if necessary, using the same conditions that for the cohort (genome release, exome kit assignment, ...)

When requested, the optional germline and somatic small variant calls are created using a modified version of the ``somatic_variant_calling`` step.
The ``somatic__cnv_calling`` step generates the small variants (TODO: how exactly) and stores them (TODO: where exactly).

Likewise, purity estimations can be automatically computed by the ``somatic__cnv_calling`` step,
to supplement or replace the estimations that may be provided in the samplesheet.

===========
Step Output
===========

TODO: The whole section of output needs revision. Main question is: what is the best format to encode CNAs?

There is no widely used standard to report copy number alterations.
In absence of a better solution, all CNV tools implemented in somatic pipeline output the segmentation table loosely following the `DNAcopy format <https://bioconductor.org/packages/devel/bioc/manuals/DNAcopy/man/DNAcopy.pdf>`_.`
The copy number call may or may not be present, and the chromosome number is replaced by its name.
The segmentation output is in file ``output/<mapper>.<cnv caller>.<lib name>/out/<mapper>.<cnv caller>.<lib name>_dnacopy.seg``.

::

    output/
    +-- bwa.cnvkit.P001-N1-DNA1-WES1
    |   |-- out
    |   |   |-- bwa.cnvkitP001-N1-DNA1-WES1_dnacopy.seg
            [...]

Note that tool ``cnvetti`` doesn't follow the snappy convention above:
the tool name is followed by an underscore & the action, where the action is one of ``coverage``, ``segment`` and ``postprocess``.
For example, the output directory would contain a directory named ``bwa.cnvetti_coverage.P002-T1-DNA1-WES1``.

.. note:: Tool-Specific Output

    Each tool produces its own set of outputs, generally not in standard format.
    Some of these files are linked from ``work`` to ``output``, but not necessarily all of them.
    Some tools (for example ``cnvkit``) also produces a report, with tables and figures.


=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_targeted_seq_cnv_calling.rst

=====================================
Available Somatic Targeted CNV Caller
=====================================

- ``cnvkit`` (for both WGS & WES)
- ``sequenza`` (only WES)
- ``purecn`` (only WES)
- ``Control-FREEC`` (only WGS - this tools might not be supported)

================================
Logic of the step for ``cnvkit``
================================

--------
Overview
--------

``cnvkit`` was designed to call CNV on whole exome data. It has the concept of _targets_ (the regions enriched by the exome kit),
and the _antitargets_ (those regions outside of enrichment).
The coverage of _targets_ and _antitargets_ are expected to be very different,
but there is still information to be gained in the _antitarget_ regions,
albeit at a much lower resolution than for _target_ regions.

``cnvkit`` was later used with some success on whole genome data.
WGS data was defined as _target_ regions covering the whole genome, with empty _antitarget_ regions.

------------------------
Sample-independent files
------------------------

``cvnkit`` allows the user to define _accessible_ regions (_via_ the ``access`` bed file).
This excludes repeats, low complexity or PAR regions, that cannot be properly mapped, and therefore used for CNV calling.

For exome data, the _target_ regions are supposed to be well curated, so they are not affected by the _access_ regions.
The _antitarget_ regions, however, are only defined within _accessible_ regions.
For WGS data, the _antitarget_ regions are empty, and the _target_ regions are set to the _accessible_ regions, when present.
Even in the absence of user-defined _accessible_ regions, the _target_ and _antitarget_ regions will not contain long ``N`` sequences.

Finally, the pipeline builds separates ``bed`` files for _target_ and _antitarget_ regions, for each exome kit present in the cohort,
and for WGS data if there is any.

---------
Reference
---------

The ``cnvkit`` authors recommend to use a panel of normals to normalize the coverage over bins.
This is usually created by running the ``panel_of_normals`` step.
The ``somatic_cnv_calling`` step will create a reference (panel of normals) if requested.
Otherwise, it is possible to use references created for different cohorts, but the user
must ensure that the data & methods used for the current cohort and to create the reference are compatible.
In particular, the exome enrichment kit must be identical, and the sex of the donors should be
similar (not to use a female-only reference for a male cohort, for example).

If there are not enough normal samples to create such a reference, the corresponding normal sample
can be used, in a normal/tumor pair setting similar to the somatic small variant calling situation.

In case no normals are available at all, a flat prior can be used.

------------
Calling CNVs
------------

The _target_ and _antitarget_ ``bed`` files created in the earlier sub-steps are used as input,
based on the exome kit (or WGS status).

The coverage is computed for the tumor sample, and normalised using the reference.
As seen previously, the reference can be either exome kit-based, or sample-specific.

The normalised coverage is the segmented, and copy numbers are called, optionally using
small variants and/or purity estimates.

If B-allele fractions are used, the pipeline will create the small variants, only for samples
with a corresponding normal.
If purity is used, the user can choose to override the values in the sample sheet (when present)
with the output of the tool of her choice.
"""

import os
import os.path
import re
import typing

from biomedsheets.models import BioEntity, BioSample, TestSample, NGSLibrary
from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import OutputFiles, Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import SomaticCnvCalling as SomaticCnvCallingConfigModel
from .model import Sex, LibraryKitDefinition, PanelOfNormalsOrigin

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the somatic_targeted_seq_cnv_calling step
DEFAULT_CONFIG = SomaticCnvCallingConfigModel.default_config_yaml_string()

#: JSON key for "isCancer"
KEY_IS_CANCER = "isCancer"

#: Value for "libraryType" is whole exome sequencing
VALUE_WES = "WES"

#: Value for "libraryType" is panel sequencing
VALUE_PANEL = "Panel-seq"

#: Values for targeted sequencing
VALUES_TARGETED_SEQ = (VALUE_WES, VALUE_PANEL)

#: Standard key/extension values for BCF files
BCF_KEY_EXTS = (
    ("bcf", ".bcf"),
    ("bcf_md5", ".bcf.md5"),
    ("bcf_csi", ".bcf.csi"),
    ("bcf_csi_md5", ".bcf.csi.md5"),
)


class SomaticCnvCallingStepPart(BaseStepPart):
    """Shared code for all caller classes in somatic_targeted_seq_cnv_calling"""

    def __init__(self, parent: "SomaticCnvCallingWorkflow"):
        super().__init__(parent)

    def _get_sample_sex(self, library_name: str) -> Sex:
        if self.config.sex == Sex.MALE or self.config.sex == Sex.FEMALE:
            sample_sex = self.config.sex
        elif self.config.sex == Sex.SAMPLESHEET and library_name in self.parent.sex:
            sample_sex = self.parent.sex[library_name]
        else:
            sample_sex = Sex.UNKNOWN
        return sample_sex

    @staticmethod
    @dictify
    def _get_log_file_from_prefix(prefix: str) -> typing.Iterator[typing.Dict[str, str]]:
        key_ext = (
            ("log", ".log"),
            ("sh", ".sh"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class CnvKitStepPart(SomaticCnvCallingStepPart):
    """Perform somatic targeted CNV calling using cnvkit"""

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "access",
        "target",
        "antitarget",
        "coverage",
        "reference",
        "flat_reference_panel",
        "flat_reference_wgs",
        "fix",
        "segment",
        "call",
        "bintest",
        "plot/diagram",
        "plot/scatter",
        "report/metrics",
        "report/segmetrics",
    )

    # Overwrite defaults
    default_resource_usage = ResourceUsage(threads=1, time="03:59:59", memory="7680M")  # 4h

    def __init__(self, parent: SomaticCnvCallingStepPart):
        super().__init__(parent)

    def get_input_files(self, action: str) -> typing.Callable:
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action.replace("/", "_")))

    def get_params(self, action: str) -> typing.Callable:
        """Return parameters input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_params_{}".format(action.replace("/", "_")))

    def get_output_files(self, action: str) -> typing.Callable:
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        f = getattr(self, "_get_output_files_{}".format(action.replace("/", "_")))
        return f()

    def get_log_file(self, action: str) -> typing.Dict[str, str]:
        """Return log files, dependent on rule"""
        # Validate action
        self._validate_action(action)
        base_name = os.path.join("work", f"{{mapper}}.{self.name}.{{library_name}}", "log")
        # Access, target & antitarget steps are cohort-wide, the others are library-dependent
        if action in ("access",):
            prefix = f"work/{self.name}/log/{action}"
        elif action in ("target", "antitarget"):
            prefix = f"work/{self.name}/log/{action}" + ".{panel_name}"
        elif action in ("coverage",):
            prefix = os.path.join(base_name, action + ".{region}")
        elif action in (
            "reference",
            "fix",
            "segment",
            "call",
            "bintest",
            "report/metrics",
            "report/segmetrics",
        ):
            prefix = os.path.join(base_name, action.replace("/", "_"))
        elif action in ("plot/diagram", "plot/scatter"):
            prefix = os.path.join(base_name, action.replace("/", "_") + ".{contig_name}")
        elif action == "flat_reference_panel":
            prefix = f"work/{{mapper}}.{self.name}/log/reference.{{panel_name}}"
        elif action == "flat_reference_wgs":
            prefix = f"work/{{mapper}}.{self.name}/log/reference"
        return SomaticCnvCallingStepPart._get_log_file_from_prefix(prefix)

    def get_result_files(self, library_name: str, mapper: str) -> typing.List[str]:
        """Files to symlink to output"""
        base_name = f"{mapper}.{self.name}.{library_name}"
        result_files = []
        # Tumor samples
        if library_name in self.parent.normal_library:
            # Main results
            prefix = os.path.join("output", base_name, "out", base_name)
            for suffix in ("cnr", "segments.cns", "cns", "bintest.cnr"):
                result_files.append(prefix + "." + suffix)
            # Log files
            prefix = os.path.join("output", base_name, "log")
            for ext in ("log", "conda_list.txt", "conda_info.txt", "sh"):
                result_files.append(os.path.join(prefix, f"coverage.target.{ext}"))
                result_files.append(os.path.join(prefix, f"coverage.antitarget.{ext}"))
            for suffix in ("fix", "segment", "call", "bintest"):
                for ext in ("log", "conda_list.txt", "conda_info.txt", "sh"):
                    result_files.append(prefix + "/" + suffix + "." + ext)
            # Log of reference is no panel of normals
            # if not self.config[self.name]["panel_of_normals"]["enabled"]:
            #     normal_library = self.parent.normal_library[library_name]
            #     prefix = os.path.join("output", f"{mapper}.{self.name}.{normal_library}", "log", f"{mapper}.{self.name}.{normal_library}.reference")
            #     for ext in ("log", "conda_list.txt", "conda_info.txt", "sh"):
            #         result_files.append(prefix + "." + ext)
            # Reports
            if "reports" in self.config[self.name]:
                prefix = os.path.join("output", base_name, "report", base_name)
                for report in ("metrics", "segmetrics"):
                    if report in self.config[self.name]["reports"]:
                        result_files.append(prefix + "." + report + ".tsv")
            # Plots (per chromosome)
            if "plots" in self.config[self.name]:
                prefix = os.path.join("output", base_name, "plot")
                for plot in ("diagram", "scatter"):
                    if plot in self.config[self.name]["plots"]:
                        for contig in self.parent.contigs:
                            result_files.append(os.path.join(prefix, plot, contig + ".png"))
        # else:  # Normal samples
        #     prefix = os.path.join("output", base_name, "log", "reference")
        #     for ext in ("log", "conda_list.txt", "conda_info.txt", "sh"):
        #         result_files.append(prefix + "." + ext)
        return result_files

    # ----- Access --------------------------------------------------------------------------------

    def _get_input_files_access(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return None

    def _get_params_access(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        params = {"reference": self.w_config.static_data_config.reference.path}
        params["min_gap_size"] = self.config[self.name]["access"]["min_gap_size"]
        access = self.config[self.name]["access"].get("exclude", None)
        if access:
            params["access"] = access

    def _get_output_files_access(self) -> typing.Dict[str, str]:
        return {"access": f"work/{self.name}/out/access.bed"}

    # ----- Target --------------------------------------------------------------------------------

    def _get_input_files_target(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        for panel in self.config.path_target_interval_list_mapping:
            if panel.name == wildcards.panel_name:
                return {"region": panel.path}

    def _get_params_target(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "split": self.config[self.name]["target"]["split"],
            "avg_size": self.config[self.name]["target"]["avg_size"],
        }

    def _get_output_files_target(self) -> typing.Dict[str, str]:
        return {"region": f"work/{self.name}/out/{{panel_name}}_target.bed"}

    # ----- Antitarget ----------------------------------------------------------------------------

    def _get_input_files_antitarget(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # No antitarget for WGS
        return {
            "target": f"work/{self.name}/out/{wildcards.panel_name}_target.bed",
            "access": f"work/{self.name}/out/access.bed",
        }

    def _get_params_antitarget(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "avg_size": self.config[self.name]["antitarget"]["avg_size"],
            "min_size": self.config[self.name]["antitarget"]["min_size"],
        }

    def _get_output_files_antitarget(self) -> typing.Dict[str, str]:
        return {"region": f"work/{self.name}/out/{{panel_name}}_antitarget.bed"}

    # ----- Coverage ------------------------------------------------------------------------------

    def _get_input_files_coverage(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # BAM/BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam.bai"),
        }

        # Region (target or antitarget) file
        panel = self.parent.libraryKit.get(wildcards.library_name, None)
        if panel is None:
            input_files["region"] = f"work/{self.name}/out/access.bed"
        else:
            input_files["region"] = f"work/{self.name}/out/{panel.name}_{wildcards.region}.bed"
        return input_files

    def _get_params_coverage(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "fasta": self.w_config.static_data_config.reference.path,
            "count": self.config[self.name]["coverage"]["count"],
            "min_mapq": self.config[self.name]["coverage"]["min_mapq"],
            "processes": self.default_resource_usage.threads,
        }

    def _get_output_files_coverage(self) -> typing.Dict[str, str]:
        return {"coverage": f"work/{{mapper}}.{self.name}.{{library_name}}/out/{{region}}.cnn"}

    # ----- Reference -----------------------------------------------------------------------------

    def _get_input_files_reference(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        """Builds reference from the paired normal, or flat prior in absence of normal"""
        input_files = {}
        normal_library = self.parent.normal_library.get(wildcards.library_name, None)
        input_files["normals"] = [
            f"work/{wildcards.mapper}.{self.name}.{normal_library}/out/target.cnn",
            f"work/{wildcards.mapper}.{self.name}.{normal_library}/out/antitarget.cnn",
        ]
        return input_files

    def _get_params_reference(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        params = {
            "fasta": self.w_config.static_data_config.reference.path,
            "cluster": self.config[self.name]["reference"]["cluster"],
            "min_cluster_size": self.config[self.name]["reference"]["min_cluster_size"],
            "male_reference": self.config[self.name]["use_male_reference"],
            "no_gc": self.config[self.name]["reference"]["no_gc"],
            "no_edge": self.config[self.name]["reference"]["no_edge"],
            "no_rmask": self.config[self.name]["reference"]["no_rmask"],
        }
        sample_sex = self._get_sample_sex(wildcards.library_name)
        if sample_sex == Sex.MALE or sample_sex == Sex.FEMALE:
            params["sample_sex"] = str(sample_sex)
        return

    def _get_output_files_reference(self) -> typing.Dict[str, str]:
        """TODO: flat prior reference should be library-independent"""
        return {"reference": f"work/{{mapper}}.{self.name}.{{library_name}}/out/reference.cnn"}

    def _get_input_files_flat_reference_panel(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        """Builds reference from the paired normal, or flat prior in absence of normal"""
        input_files = {}
        panel = self.parent.libraryKit.get(wildcards.library_name, None)
        if panel is None:  # WGS, target is access, no antitarget
            input_files["target"] = f"work/{self.name}/out/access.bed"
        else:  # WES, both target & antitarget
            input_files["target"] = f"work/{self.name}/out/{panel.name}_target.bed"
            input_files["antitarget"] = f"work/{self.name}/out/{panel.name}_antitarget.bed"
        return input_files

    def _get_input_files_flat_reference_panel(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return self._get_input_files_flat_reference_panel(wildcards)

    def _get_params_flat_reference_panel(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return self._get_params_reference(wildcards)

    def _get_output_files_flat_reference_panel(self) -> typing.Dict[str, str]:
        """TODO: flat prior reference should be library-independent"""
        return {"reference": f"work/{{mapper}}.{self.name}/out/reference.{{panel_name}}.cnn"}

    def _get_input_files_flat_reference_wgs(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return self._get_input_files_flat_reference_panel(wildcards)

    def _get_params_flat_reference_wgs(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return self._get_params_reference(wildcards)

    def _get_output_files_flat_reference_wgs(self) -> typing.Dict[str, str]:
        """TODO: flat prior reference should be library-independent"""
        return {"reference": f"work/{{mapper}}.{self.name}/out/reference.cnn"}

    # ----- Fix -----------------------------------------------------------------------------------

    def _get_input_files_fix(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # Coverage on targets
        input_files = {
            "target": f"work/{wildcards.mapper}.{self.name}.{wildcards.library_name}/out/target.cnn"
        }
        # Coverage on antitargets when present (absent for WGS)
        panel = self.parent.libraryKit.get(wildcards.library_name, None)
        if panel is not None:  # WGS - no antitarget
            input_files["antitarget"] = (
                f"work/{wildcards.mapper}.{self.name}.{wildcards.library_name}/out/antitarget.cnn"
            )
        # Get reference from panel of normals if available, otherwise from normal or flat when no normal
        if not self.config[self.name]["panel_of_normals"]["enabled"]:  # Paired normal or flat
            normal_library = self.parent.normal_library.get(wildcards.library_name, None)
            if normal_library:
                input_files["reference"] = (
                    f"work/{{mapper}}.{self.name}.{normal_library}/out/reference.cnn"
                )
            else:
                if panel:
                    input_files["reference"] = (
                        f"work/{{mapper}}.{self.name}/out/reference.{panel.name}.cnn"
                    )
                else:
                    input_files["reference"] = f"work/{{mapper}}.{self.name}/out/reference.cnn"
        elif (
            self.config[self.name]["panel_of_normals"]["origin"]
            == PanelOfNormalsOrigin.PREVIOUS_STEP
        ):  # Panel_of_normals step
            input_files["reference"] = self.parent._get_panel_of_normals_path(self.name, panel)
        else:
            input_files["reference"] = self.config[self.name]["panel_of_normals"][
                "path_panel_of_normals"
            ]
        return input_files

    def _get_params_fix(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "sample_id": wildcards.library_name,
            "cluster": self.config[self.name]["fix"]["cluster"],
            "no_gc": self.config[self.name]["fix"]["no_gc"],
            "no_edge": self.config[self.name]["fix"]["no_edge"],
            "no_rmask": self.config[self.name]["fix"]["no_rmask"],
        }

    def _get_output_files_fix(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"coverage": os.path.join("work", base_name, "out", base_name + ".cnr")}

    # ----- Segment -------------------------------------------------------------------------------

    def _get_input_files_segment(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # Coverage
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        input_files = {"coverage": f"work/{base_name}/out/{base_name}.cnr"}
        # Segmentation using SNVs if requested and available (normal must be present)
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            input_files["variants"] = f"work/{base_name}/out/{base_name}.vcf.gz"
        return input_files

    def _get_params_segment(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # Segmentation parameters
        params = {
            "method": self.config[self.name]["segment"]["method"],
            "threshold": self.config[self.name]["segment"]["threshold"],
            "drop_low_coverage": self.config[self.name]["segment"]["drop_low_coverage"],
            "drop_outliers": self.config[self.name]["segment"]["drop_outliers"],
        }
        if self.config[self.name]["segment"]["method"] == "cbs":
            params["smooth_cbs"] = self.config[self.name]["segment"]["smooth_cbs"]
        params["processes"] = self.default_resource_usage.threads
        # Normal & tumor sample ids if SNVs
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            params["sample_id"] = wildcards.library_name
            params["normal_id"] = self.normal_library[wildcards.library_name]
            params["min_variant_depth"] = self.config[self.name]["segment"]["min_variant_depth"]
            params["zygocity_freq"] = self.config[self.name]["segment"]["zygocity_freq"]
        return params

    def _get_output_files_segment(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {
            "segments": os.path.join("work", base_name, "out", base_name + ".segments.cns"),
            "dataframe": os.path.join("work", base_name, "out", "dataframe.rds"),
        }

    # ----- Call ----------------------------------------------------------------------------------

    def _get_input_files_call(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # Segmentation
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        input_files = {"segments": f"work/{base_name}/out/{base_name}.segments.cns"}
        # SNVs if requested and available (normal must be present)
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            input_files["variants"] = f"work/{base_name}/out/{base_name}.vcf.gz"
        # Purity from the tool if requested and not from the samplesheet
        if (
            self.config[self.name]["purity"]["enabled"] and self.config[self.name]["purity"]["tool"]
        ):  # Need purity, and can use tool to obain it
            if (
                self.config[self.name]["purity"]["ignore_samplesheet"]
                or wildcards.library_name not in self.parent.purity
            ):
                # Don't use samplesheet
                input_files["purity"] = (
                    f"work/{base_name}/out/{wildcards.mapper}.{self.config.purity.tool}.txt"
                )
        return input_files

    def _get_params_call(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        # Call parameters
        params = {
            "method": self.config[self.name]["call"]["method"],
            "thresholds": self.config[self.name]["call"]["thresholds"],
            "filter": self.config[self.name]["call"]["filter"],
            "drop_low_coverage": self.config[self.name]["call"]["drop_low_coverage"],
            "male_reference": self.config[self.name]["use_male_reference"],
        }
        # If center_at defined, use it, otherwise use the center method
        center = self.config[self.name]["call"].get("center_at", None)
        if center is not None:
            params["center_at"] = center
        else:
            params["center"] = self.config[self.name]["call"].get("center", "None")
        # Normal & tumor sample ids if SNVs
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            params["sample_id"] = wildcards.library_name
            params["normal_id"] = self.normal_library[wildcards.library_name]
        # Sample sex if known, otherwise guessed by the tool
        sample_sex = self._get_sample_sex(wildcards.library_name)
        if sample_sex == Sex.MALE or sample_sex == Sex.FEMALE:
            params["sample_sex"] = sample_sex
        # If requested, purity from samplesheet or from default if no tool
        if self.config[self.name]["purity"]["enabled"]:
            purity = self.parent.purity.get(
                wildcards.library_name, self.config.purity.default_purity
            )
            if purity is not None and not self.config[self.name]["purity"]["ignore_samplesheet"]:
                params["purity"] = purity
                if self.config.default_ploidy:
                    params["ploidy"] = self.config.default_ploidy
        return params

    def _get_output_files_call(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"calls": os.path.join("work", base_name, "out", base_name + ".cns")}

    # ----- Bintest -------------------------------------------------------------------------------

    def _get_input_files_bintest(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        return {
            "coverage": f"work/{base_name}/out/{base_name}.cnr",
            "segments": f"work/{base_name}/out/{base_name}.segments.cns",
        }

    def _get_params_bintest(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "alpha": self.config[self.name]["bintest"]["alpha"],
            "target": self.config[self.name]["bintest"]["target"],
        }

    def _get_output_files_bintest(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"coverage": os.path.join("work", base_name, "out", base_name + ".bintest.cnr")}

    # ----- Plots --------------------------------------------------------------------------------

    def _get_input_files_plot_diagram(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        input_files = {
            "coverage": f"work/{base_name}/out/{base_name}.cnr",
            "segments": f"work/{base_name}/out/{base_name}.segments.cns",
        }
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            input_files["variants"] = f"work/{base_name}/out/{base_name}.vcf.gz"
        return input_files

    def _get_params_plot_diagram(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "threshold": self.config[self.name]["plots"]["diagram"]["threshold"],
            "min_probes": self.config[self.name]["plots"]["diagram"]["min_probes"],
            "no_shift_xy": self.config[self.name]["plots"]["diagram"]["no_shift_xy"],
        }

    def _get_output_files_plot_diagram(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"figure": os.path.join("work", base_name, "plot", "diagram", "{contig_name}.pdf")}

    def _get_input_files_plot_scatter(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        input_files = {
            "coverage": f"work/{base_name}/out/{base_name}.cnr",
            "segments": f"work/{base_name}/out/{base_name}.segments.cns",
        }
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            input_files["variants"] = f"work/{base_name}/out/{base_name}.vcf.gz"
        return input_files

    def _get_params_plot_scatter(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        params = {
            "chromosome": wildcards.contig_name,
            "antitarget_marker": self.config[self.name]["plots"]["scatter"]["antitarget_marker"],
            "by_bin": self.config[self.name]["plots"]["scatter"]["by_bin"],
            "segment_color": self.config[self.name]["plots"]["scatter"]["segment_color"],
            "trend": self.config[self.name]["plots"]["scatter"]["trend"],
            "y_max": self.config[self.name]["plots"]["scatter"]["y_max"],
            "y_min": self.config[self.name]["plots"]["scatter"]["y_min"],
            "fig_size": self.config[self.name]["plots"]["scatter"]["fig_size"],
            "sample_id": wildcards.library_name,
        }
        variants = self.config[self.name].get("variants", None)
        if variants and wildcards.library_name in self.normal_library:
            params["normal_id"] = self.normal_library[wildcards.library_name]
            params["min_variant_depth"] = self.config[self.name]["plots"]["scatter"][
                "min_variant_depth"
            ]
            params["zygocity_freq"] = self.config[self.name]["plots"]["scatter"]["zygocity_freq"]
        return params

    def _get_output_files_plot_scatter(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"figure": os.path.join("work", base_name, "plot", "scatter", "{contig_name}.pdf")}

    # ----- Metrics (metrics & segmetrics) --------------------------------------------------------

    def _get_input_files_report_metrics(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        return {
            "coverage": f"work/{base_name}/out/{base_name}.cnr",
            "segments": f"work/{base_name}/out/{base_name}.segments.cns",
        }

    def _get_params_report_metrics(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {"drop_low_coverage": self.config[self.name]["reports"]["drop_low_coverage"]}

    def _get_output_files_report_metrics(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"report": os.path.join("work", base_name, "report", base_name + ".metrics.tsv")}

    def _get_input_files_report_segmetrics(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        base_name = f"{wildcards.mapper}.{self.name}.{wildcards.library_name}"
        return {
            "coverage": f"work/{base_name}/out/{base_name}.cnr",
            "segments": f"work/{base_name}/out/{base_name}.segments.cns",
        }

    def _get_params_report_segmetrics(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        return {
            "drop_low_coverage": self.config[self.name]["reports"]["drop_low_coverage"],
            "stats": (
                "mean",
                "median",
                "mode",
                "t-test",
                "stdev",
                "sem",
                "mad",
                "mse",
                "iqr",
                "bivar",
                "ci",
                "pi",
            ),
            "alpha": self.config[self.name]["reports"]["alpha"],
            "bootstrap": self.config[self.name]["reports"]["bootstrap"],
        }

    def _get_output_files_report_segmetrics(self) -> typing.Dict[str, str]:
        base_name = f"{{mapper}}.{self.name}.{{library_name}}"
        return {"report": os.path.join("work", base_name, "report", base_name + ".segmetrics.tsv")}


class SomaticCnvCallingWorkflow(BaseStep):
    """Perform somatic targeted sequencing CNV calling"""

    #: Workflow name
    name = "somatic_cnv_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=False)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticCnvCallingConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                CnvKitStepPart,
                # ControlfreecStepPart,
                # SequenzaStepPart,
                # PureCNStepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        self.registered_pons = self._optionally_register_pon()

        # Collect extra information per library
        self.normal_library = self._get_normal_library()
        self.libraryKit = self._get_panel_information()
        self.sex = self._get_sex()
        self.purity = self._get_purity()

    def get_result_files(self) -> OutputFiles:
        fns = []
        for seq_type, tools in self.config.tools:
            for library in self._get_libraries():
                if library.extra_infos.get("libraryType", "").lower() != seq_type:
                    continue
                test_sample = library.test_sample
                if test_sample.extra_infos.get("extractionType", "") != "DNA":
                    continue
                bio_sample = test_sample.bio_sample
                is_tumor = bio_sample.extra_infos.get("isTumor", True)
                if is_tumor:
                    for tool in tools:
                        f = self.substep_getattr(tool, "get_result_files")
                        for mapper in self.w_config.step_config["ngs_mapping"]["tools"]["dna"]:
                            for fn in f(library.name, mapper):
                                fns.append(fn)
        return OutputFiles(fns)

    def _get_libraries(self) -> typing.Iterator[NGSLibrary]:
        for sheet in self.shortcut_sheets:
            for donor in sheet.sheet.bio_entities.values():
                for bio_sample in donor.bio_samples.values():
                    for test_sample in bio_sample.test_samples.values():
                        for library in test_sample.ngs_libraries.values():
                            yield library

    def _get_normal_library(self) -> typing.Dict[str, str]:
        normal_for_donor = {}
        for library in self._get_libraries():
            test_sample = library.test_sample
            if test_sample.extra_infos.get("extractionType", "") != "DNA":
                continue
            bio_sample = test_sample.bio_sample
            is_tumor = bio_sample.extra_infos.get("isTumor", None)
            if is_tumor is None:
                raise ValueError(f"Missing 'isTumor' value for library '{library.name}'")
            if is_tumor:
                continue
            donor = bio_sample.bio_entity
            if donor.name in normal_for_donor:
                raise ValueError(f"Multiple normals for donor '{donor.name}'")
            normal_for_donor[donor.name] = library.name

        normal_library = {}
        for library in self._get_libraries():
            test_sample = library.test_sample
            if test_sample.extra_infos.get("extractionType", "") != "DNA":
                continue
            bio_sample = test_sample.bio_sample
            donor = bio_sample.bio_entity
            if bio_sample.extra_infos.get("isTumor", True):
                normal_library[library.name] = normal_for_donor[donor.name]
        return normal_library

    def _optionally_register_pon(self) -> typing.Dict[str, str]:
        """
        Register all possible combination of panel of normals:
        - WGS PON for all configured WGS tools which require/can use it
        - WES PON for all configured WES tools which require/can use it, one for each enrichment kit

        Note that there is no need to specify the genome release,
        because the panel_of_normals step used here MUST be in the same project,
        so it has the same configuration, and only one genome release is allowed per configuration.
        """
        registered_pons = list()
        for tool in self.config.tools.wgs:
            pon_name = f"wgs.{tool}"
            if pon_name in registered_pons:
                continue
            if self.config[tool].get("panel_of_normals", None) and self.config[
                tool
            ].panel_of_normals.get("path_panel_of_normals_step", None):
                self.register_sub_workflow(
                    "panel_of_normals",
                    self.config[tool].panel_of_normals.path_panel_of_normals_step,
                    pon_name,
                )
                registered_pons.append(pon_name)
        for tool in self.config.tools.wes:
            for panel in self.config.path_target_interval_list_mapping:
                pon_name = f"wes.{tool}.{panel.name}"
                if pon_name in registered_pons:
                    continue
                if self.config[tool].get("panel_of_normals", None) and self.config[
                    tool
                ].panel_of_normals.get("path_panel_of_normals_step", None):
                    self.register_sub_workflow(
                        "panel_of_normals",
                        self.config[tool].panel_of_normals.path_panel_of_normals_step,
                        pon_name,
                    )
                    registered_pons.append(pon_name)
        return registered_pons

    def _get_panel_information(self) -> typing.Dict[str, str]:
        # Set default panel
        default = None
        for panel in self.config.path_target_interval_list_mapping:
            if panel.name == "__default__":
                default = panel
                break

        # Extract library pattern (the "libraryKit" column in samplesheet)
        # On output:
        # - the panel name and panel path if libraryKit is present & known
        # - the default panel path if libraryKit is undefined or not found
        # - None for WGS
        # - ValueError if libraryType is missing or unknown (not WES nor WGS)
        libraryKit = {}
        for library in self._get_libraries():
            test_sample = library.test_sample
            if test_sample.extra_infos.get("extractionType", "") != "DNA":
                continue

            libraryType = library.extra_infos.get("libraryType", None)
            if libraryType is None:
                raise ValueError(f"Missing library type for library '{library.name}'")
            elif libraryType == "WES":
                if library.extra_infos.get("libraryKit", None):
                    for panel in self.config.path_target_interval_list_mapping:
                        if re.match(panel.pattern, library.extra_infos.get("libraryKit")):
                            libraryKit[library.name] = panel
                            break
                    if library.name not in libraryKit:
                        libraryKit[library.name] = default
                else:
                    libraryKit[library.name] = default
                if libraryKit[library.name] is None:
                    raise ValueError(f"Undefined panel for library '{library.name}")
            elif libraryType == "WGS":
                libraryKit[library.name] = None
            else:
                raise ValueError(
                    f"Unknown library type '{libraryType}' for library '{library.name}'"
                )

        return libraryKit

    def _get_purity(self) -> typing.Dict[str, str]:
        """Returns the purity value from the 'purity' library extra_infos. Missing otherwise"""
        purity = {}
        for library in self._get_libraries():
            p = library.extra_infos.get("purity", None)
            if p:
                try:
                    p = float(p)
                    if 0 <= p and p <= 1:
                        purity[library.name] = p
                except:
                    pass
        return purity

    def _get_sex(self) -> typing.Dict[str, Sex]:
        sex = {}
        for library in self._get_libraries():
            donor = library.test_sample.bio_sample.bio_entity
            donor_sex = donor.extra_infos.get("sex", None)
            if donor_sex == "male":
                donor_sex = Sex.MALE
            elif donor_sex == "female":
                donor_sex = Sex.FEMALE
            else:
                donor_sex = Sex.UNKNOWN
            sex[library.name] = donor_sex
        return sex

    def _get_panel_of_normals_path(self, tool: str, panel: LibraryKitDefinition | None) -> str:
        pon_path = None
        assert self.config[tool]["panel_of_normals"][
            "enabled"
        ], f"Panel of normals not enabled for '{tool}'"
        assert (
            self.config[tool]["panel_of_normals"]["origin"] == PanelOfNormalsOrigin.PREVIOUS_STEP
        ), f"'{tool}' panel of normals not from previous step"
        if panel is None:
            pon_id = f"wgs.{tool}"
        else:
            pon_id = f"wes.{tool}.{panel.name}"
        assert pon_id in self.registered_pons, f"Requested panel '{pon_id}' not registered"
        pon = self.parent.sub_workflows[pon_id]
        pon_path = pon(f"output/{{mapper}}.{tool}/out/{panel.name}.ext")
        return pon_path
