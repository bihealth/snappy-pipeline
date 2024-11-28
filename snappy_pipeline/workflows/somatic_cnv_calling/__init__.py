# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_cnv_calling`` step

This step allows for the detection of CNV events for cancer samples from targeted sequenced (e.g.,
exomes or large panels) or whole genome sequencing. Panel sequencing is not implemented yet, it might be in a later release.
The wrapped tools start from the aligned reads (thus off ``ngs_mapping``) and generate CNV calls for somatic variants.

The wrapped tools implement different strategies.  Some work "reference free" and just use the
somatic BAM files for their input, some work in "matched cancer normal mode" and need the cancer
and normal BAM files, and finally others use a set of non-cancer BAM files for their background (the panel of normals).

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

- A static file or files, from another cohort or from public datasets.
  In this case, the user is responsible to make sure that the data & methods used to create the panel are compatible to the cohort's.
- The ``panel_of_normals`` step.
  The panel will be created if necessary, using data from normal samples from the cohort.

When requested, the optional germline and somatic small variant calls are created in the ``somatic_variant_calling`` step.
Once again, it is the responsability of the user to make sure that variants creeated in that way are suitable for CNV calling.

Likewise, purity estimations can be automatically computed by the ``somatic__cnv_calling`` step,
when estimations are not provided in the samplesheet or configuration file.

===========
Step Output
===========

TODO: The whole section of output needs revision. Main question is: what is the best format to encode CNAs?
``vcf`` is an possibility, the benefits are a (more or less) well-defined format, but the major drawback is
that (as far as I know), most CNV analysis tools (from ``R`` in particular) don't recognize this format (for CNAs).

Currently, the only implemented tool is ``cnvkit``. Therefore, the ``cnvkit`` output is left as produced by the software.

------
cnvkit
------

The structure of the output is:

::

    output/
    +-- bwa.cnvkit.P001-N1-DNA1-WES1
    |   |-- out
    |   |   |-- bwa.cnvkitP001-N1-DNA1-WES1.<extension>
            [...]

There are 4 main outputs:

- The ratios (extension ``.cnr``) contains the ratio of expected coverage between tumor and the reference
  in each bin, or logarithmic scale.
  This can be used to examine the data or experiment with different segmentation algorithms.
- The segments (extension ``.segments.cns``) contains the output of the segmentation. A single log2 ratio value is
  attributed to each segment.
  The segmentation covers most of the part of the genome accessible to mapping.
- The calls (extension ``calls.cns``) contains only the non-diploid segments, called after thresholding.
- The results of differential coverage tests by bins (extension ``bintest.cns``). Only significant tests are listed.

Reports & plots are also available on user's request, found in ``report`` and ``plot`` sub-directories.


=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_cnv_calling.rst

=====================================
Available Somatic Targeted CNV Caller
=====================================

- ``cnvkit`` (for both WGS & WES)

===========================================
Description of the step for ``cnvkit`` tool
===========================================

--------
Overview
--------

``cnvkit`` was designed to call CNV on whole exome data. It has the concept of _targets_ (the regions enriched by the exome kit),
and the _antitargets_ (those regions accessible for CNV calling, but outside of enrichment).
The coverage of _targets_ and _antitargets_ are expected to be very different,
but there is still information to be gained in the _antitarget_ regions,
albeit at a much lower resolution than for _target_ regions.

``cnvkit`` was later used with some success on whole genome data.
WGS data was defined as _target_ regions covering the whole genome, with empty _antitarget_ regions.

.. tip:: For Agilent kits, the ``cnvkit`` authors recommend to use baits as targets. Baits are in the ``*_Covered.bed`` file.

---------------------------------
Regions accessible to CNV calling
---------------------------------

``cnvkit`` needs to know about the regions of the genome accessible to CNV calling.
Typically, regions masked with ``N`` are excluded, but a user may also want to exclude
repeats, segmental duplications, or low complexity regions.

There are multiple ways to get this information:

1. The user can provide a ``bed`` file detailing accessible regions, using the ``path_access`` option in the ``access`` part of the configuration.
2. The user can specifically exclude regions using the ``exclude`` option is the ``access`` part of the configuration.
   In this case, the pipeline will create an accessible regions file from the whole genome and the excluded parts.
3. Otherwise, the pipeline creates the accessible regions file from the reference genome sequence, by removing
   only parts masked with ``N``.

-----------------------
The reference coverage
-----------------------

``cnvkit`` build a reference coverage to compensate locii effects when assessing coverage changes in tumor samples.
Because this reference is best constructed from multiple normal samples, the pipeline implements it under the
``panel_of_normals`` section.
The pipeline offers 4 different ways to built this reference:

1. ``cohort``: the reference is taken from the pipeline's ``panel_of_normals`` step.
   If it isn't there, it will be created, according to its configuration (which must be present).
   The ``cnvkit`` authors suggest that `10 to 20 normal samples <https://www.biostars.org/p/296422/>`_
   are sufficient to build a good reference. When there are more samples, the selected one should be taken from those with average file size.
2. ``file``: the reference is taken from another panel of normals, possibly from another cohort or from public data.
   Beside the reference coverage itself, the target and antitarget bed files must also be provided.
   Note that it is the user's responsability to make sure that the panel of normal is suitable for the cohort.
3. ``paired``: the reference is built only from one normal sample, paired with the tumor.
   It is _not_ recommended by the ``cnvkit`` authors, but can be beneficial in some circumstances
   (for example experimental designs with treated cell lines).
4. ``flat``: a _flat_ reference is computed, which discards locus-specific effects.
   It should only be used when there are no normals nor suitable panel, or for benchmarking purposes.

------------------
WGS-specific notes
------------------

In WGS mode, the _antitarget_ regions (defined as accessible regions without target regions) are empty.

But _target_ regions need to be carefully selected: the size of the bins used to identify CNAs must be
selected so that discovery of focal events is possible.
When the reference is taken from the current cohort, or from another panel of normals, the target regions
are taken from either panel of normals, and no other consideration is necessary.
But otherwise, it must be computed provided by the user (configuration option ``avg_size`` in the ``target`` section),
or from the available data.

In the latter case, the algorithm is as follows:

- When no normal sample data is available, then ``cvnkit`` default value is used.
- Otherwise, the value is computed by ``cnvkit``'s ``autobin`` module.
  It is run in ``wgs`` mode if the access regions cover the complete genome
  (no access file provided, no exclude files, and the ``min_gap_size`` parameter not set),
  else it is run in ``amplicon`` mode.
  Note that the latter should only be used when user-defined accessible regions are quite restricted
  (limited to protein-coding exons devoid of repeats or low complexity regions, for example).

-----------------------------
Other notes on implementation
-----------------------------

.. note:: CNA calling on panel data is not implemented yet, even though ``cnvkit`` allows it in principle.

.. note:: The current pipeline tries to replicate the behaviour of the ``batch`` module of ``cnvkit``,
   while keeping the flexibility to diverge from it.
   In particular, the possibility of obtaining the reference coverage from the paired normal is implemented.

.. note:: The current implementation doesn't allow to mix multiple exome enrichment kits.
   Future versions will hopefully lift this restriction. However, mixing WES, WGS & possibly panel data is
   more challenging, and is not on the roadmap for future improvements.

"""

import os
import os.path
import re
from typing import Callable, NamedTuple, Any

from biomedsheets.models import NGSLibrary
from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from biomedsheets.io_tsv.base import LIBRARY_TO_EXTRACTION, EXTRACTION_TYPE_DNA
from snakemake.io import OutputFiles, Wildcards, InputFiles

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from snappy_pipeline.models.common import Sex, SexOrigin, SexValue
from snappy_pipeline.models.cnvkit import SegmentationMethod as CnvkitSegmentationMethod

from .model import SomaticCnvCalling as SomaticCnvCallingConfigModel
from .model import Cnvkit as CnvkitConfig
from .model import PanelOfNormalsOrigin, PurityOrigin, VariantOrigin

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the somatic_targeted_seq_cnv_calling step
DEFAULT_CONFIG = SomaticCnvCallingConfigModel.default_config_yaml_string()


class SomaticCnvCallingStepPart(BaseStepPart):
    """Shared code for all caller classes in somatic_targeted_seq_cnv_calling"""

    def __init__(self, parent: "SomaticCnvCallingWorkflow"):
        super().__init__(parent)


class CnvKitStepPart(SomaticCnvCallingStepPart):
    """Perform somatic targeted CNV calling using cnvkit"""

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "access",
        "autobin",
        "target",
        "antitarget",
        "coverage",
        "reference",
        "fix",
        "segment",
        "call",
        "bintest",
        "scatter",
        "metrics",
        "genemetrics",
        "segmetrics",
    )

    # Overwrite defaults
    default_resource_usage = ResourceUsage(threads=1, time="03:59:59", memory="7680M")  # 4h

    def __init__(self, parent: SomaticCnvCallingStepPart):
        super().__init__(parent)

        self.is_wgs = (
            any([libraryKit is None for libraryKit in self.parent.tumors.keys()])
            and self.name in self.config.tools.wgs
        )
        self.is_wes = (
            any([libraryKit is not None for libraryKit in self.parent.tumors.keys()])
            and self.name in self.config.tools.wes
        )
        assert not (self.is_wgs and self.is_wes), "WES & WGS are mixed"

        if self.is_wgs or self.is_wes:
            assert (
                len(self.parent.tumors) == 1
            ), "Current cnvkit tool implementation can't handle multiple library types or kits"

            self.libraryKit = list(self.parent.tumors.keys())[0]
            self.tumors = {x.library.name: x for x in self.parent.tumors[self.libraryKit]}

            self.cfg: CnvkitConfig = self.config.get(self.name)
            self.pon_source = self.cfg.panel_of_normals.source

            self._set_cnvkit_pipeline_logic()

            self.path_baits = self._get_path_baits()

            if (
                self.cfg.somatic_purity_ploidy_estimate.enabled
                and self.cfg.somatic_purity_ploidy_estimate.source == PurityOrigin.SAMPLESHEET
            ):
                assert not any(
                    [x.purity is None for x in self.tumors.values()]
                ), "Missing purity value from samplesheet"

            self.base_out = "work/{mapper}.cnvkit/out/cnvkit."
            self.base_out_lib = (
                "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}."
            )

    def _set_cnvkit_pipeline_logic(self):
        """
        Creates instance variables to choose path in cnvkit pipeline

        Access: regions accessible for CNV calling (unmasked)
            path_access or when missing build from genome reference + optional list of excluded region

        Target: regions of good coverage
            From baits (WES) or accessible regions (WGS) + estimate of target size from config or autobin step

        Antitarget: regions of low coverage
            antitarget = access - target, only WES, otherwise empty

        Reference:
            Flat: based on targets & antitargets only
            Cohort: from panel_of_normals step
            File: from another cohort or public data (reference + target + antitarget [WES only])
            Paired: reference built from the target & antitarget coverage of one normal sample only (paired with the tumor)

        Therefore, a reference must be created for flat & paired choices (one reference per normal sample in the latter case).
        The logic to create the reference is (panel of normal is pon):
        - access created if path_access is missing or average target size must be estimated
        - average target size estimated if value not in config and dataset is WGS
        - target created always
        - antitarget created when dataset is WES
        """
        self.paired = self.pon_source == PanelOfNormalsOrigin.PAIRED
        self.build_ref = self.paired or self.pon_source == PanelOfNormalsOrigin.FLAT
        self.compute_avg_target_size = (
            self.is_wgs and self.paired and self.cfg.target.avg_size is None
        )
        self.create_access = self.build_ref and (not self.cfg.path_access)
        self.plain_access = (
            not self.cfg.path_access
            and len(self.cfg.access.exclude) == 0
            and self.cfg.access.min_gap_size is None
        )

        self.variants_from_cohort = (
            self.cfg.somatic_variant_calling.enabled
            and self.cfg.somatic_variant_calling.source == VariantOrigin.COHORT
        )
        self.variants_from_file = (
            self.cfg.somatic_variant_calling.enabled
            and self.cfg.somatic_variant_calling.source == VariantOrigin.FILE
        )

    def _get_sample_sex(self, library_name: str | None) -> SexValue | None:
        if self.cfg.sample_sex.source == SexOrigin.SAMPLESHEET and library_name:
            sample_sex = self.tumors[library_name].sex
        elif self.cfg.sample_sex.source == SexOrigin.CONFIG:
            sample_sex = self.cfg.sample_sex.default
        else:
            sample_sex = None
        return sample_sex

    def _get_path_baits(self) -> str | None:
        if not self.is_wes:
            return None
        default = None
        for item in self.cfg.path_target_interval_list_mapping:
            if item.name == self.libraryKit:
                return item.path
            elif item.name == "__default__":
                default = item.path
        if default is None:
            raise ValueError(f"Missing library kit definition for {self.libraryKit}")
        return default

    def get_input_files(self, action: str) -> Callable:
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action.replace("/", "_")))

    def get_args(self, action: str) -> Callable:
        """Return parameters input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_args_{}".format(action.replace("/", "_")))

    @dictify
    def get_output_files(self, action: str):
        """
        Return output paths, dependent on rule

        It is important to take good care of wildcards, because
        when a paired reference is used on WGS without setting the avg target size,
        the output of autobin and target are built for the normal library.
        So in this case, library_name stands for the normal library, rather than
        for the tumor.
        """
        self._validate_action(action)

        base_report_lib = (
            "work/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}."
        )

        output_files = {}
        match action:
            case "access":
                output_files = {"access": self.base_out + "access.bed"}
            case "autobin":
                output_files = {"result": self.base_out_lib + "autobin.txt"}
            case "target":
                if self.compute_avg_target_size and self.paired:
                    output_files = {"target": self.base_out_lib + "target.bed"}
                else:
                    output_files = {"target": self.base_out + "target.bed"}
            case "antitarget":
                output_files = {"antitarget": self.base_out + "antitarget.bed"}
            case "coverage":
                output_files = {"coverage": self.base_out_lib + "{region,(target|antitarget)}.cnn"}
            case "reference":
                if self.paired:
                    output_files = {"reference": self.base_out_lib + "reference.cnn"}
                else:
                    output_files = {"reference": self.base_out + "reference.cnn"}
            case "fix":
                output_files = {"ratios": self.base_out_lib + "cnr"}
            case "segment":
                output_files = {
                    "segments": self.base_out_lib + "segments.cns",
                    "dataframe": self.base_out_lib + "rds",
                }
            case "call":
                output_files = {"calls": self.base_out_lib + "calls.cns"}
            case "bintest":
                output_files = {"tests": self.base_out_lib + "bintest.cns"}
            case "metrics":
                output_files = {"report": base_report_lib + "metrics.tsv"}
            case "segmetrics":
                output_files = {"report": base_report_lib + "segmetrics.tsv"}
            case "genemetrics":
                output_files = {"report": base_report_lib + "genemetrics.tsv"}
            case "scatter":
                output_files = {
                    "plot": "work/{mapper}.cnvkit.{library_name}/plot/{mapper}.cnvkit.{library_name}.scatter.{contig_name}.jpeg"
                }

        for k, v in output_files.items():
            yield k, v
            yield k + "_md5", v + ".md5"

    @dictify
    def get_log_file(self, action):
        """Return panel of normal files"""
        # Validate action
        self._validate_action(action)

        base_log = "work/{mapper}.cnvkit/log/cnvkit."
        base_log_lib = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}."

        if action in ("access", "antitarget"):
            tpl = base_log + action
        elif action in (
            "autobin",
            "fix",
            "segment",
            "call",
            "bintest",
            "metrics",
            "segmetrics",
            "genemetrics",
        ):
            tpl = base_log_lib + action
        elif action == "target":
            if self.compute_avg_target_size and self.paired:
                tpl = base_log_lib + "target"
            else:
                tpl = base_log + "target"
        elif action == "reference":
            if self.paired:
                tpl = base_log_lib + "reference"
            else:
                tpl = base_log + "reference"
        elif action == "coverage":
            tpl = base_log_lib + "{region,(target|antitarget)}coverage"
        elif action in ("scatter",):
            tpl = base_log_lib + action + ".{contig_name}"
        else:
            raise ValueError(f"Logs of action '{action}' not implemented yet")

        for key, ext in (
            ("conda_list", ".conda_list.txt"),
            ("conda_info", ".conda_info.txt"),
            ("log", ".log"),
            ("sh", ".sh"),
        ):
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    def get_result_files(self, library_name: str, mapper: str) -> list[str]:
        """Files to symlink to output"""
        base_out_lib = (
            "output/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}."
        ).format(mapper=mapper, library_name=library_name)

        base_report_lib = (
            "output/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}."
        ).format(mapper=mapper, library_name=library_name)

        base_plot_lib = (
            "output/{mapper}.cnvkit.{library_name}/plot/{mapper}.cnvkit.{library_name}."
        ).format(mapper=mapper, library_name=library_name)

        result_files = []

        for suffix in ("cnr", "segments.cns", "calls.cns", "bintest.cns"):
            result_files.append(base_out_lib + suffix)

        actions_to_log = ("fix", "segment", "call", "bintest")
        for action in actions_to_log:
            result_files += [
                path.replace("work", "output", 1).format(mapper=mapper, library_name=library_name)
                for path in filter(
                    lambda p: not p.endswith(".md5"), self.get_log_file(action).values()
                )
            ]

        # Logs of metrics not linked
        for report in ("metrics", "segmetrics", "genemetrics"):
            if self.cfg.get(report).get("enabled"):
                result_files.append(base_report_lib + report + ".tsv")

        # Logs of plots not links
        # TODO: Mouse date: only chromosomes 1 to 19
        chrs = ["all"] + list(map(str, range(1, 23))) + ["X"]
        if (
            self.cfg.sample_sex.source != SexOrigin.CONFIG
            or self.cfg.sample_sex.default == SexValue.FEMALE
        ):
            chrs.append("Y")

        for plot in ("scatter",):
            if self.cfg.get(plot).get("enabled"):
                for chr in chrs:
                    result_files.append(base_plot_lib + f"{plot}.{chr}.jpeg")

        result_files += [x + ".md5" for x in result_files]
        return result_files

    # ----- Access --------------------------------------------------------------------------------

    def _get_input_files_access(self, wildcards: Wildcards) -> dict[str, str]:
        assert self.create_access, "Should not build access, already available"
        return {}

    def _get_args_access(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        """
        Arguments used to compute accessible regions for mapping

        When accessible regions are needed to compute average target size
        (WGS without average target size set in the config)
        then accessible region must cover the full genome (except masked).
        Otherwise, access is built with excluded regions.
        This happens when the average target size is set in the config in WGS,
        or for WES.
        """
        assert self.create_access, "Should not build access, already available"
        return dict(input) | {
            "reference": self.w_config.static_data_config.reference.path,
            "min-gap-size": self.cfg.access.min_gap_size,
            "exclude": self.cfg.access.exclude,
        }

    # ----- Autobin (never used directly, only to compute target size in WGS settings) ------------

    def _get_input_files_autobin(self, wildcards: Wildcards) -> dict[str, str]:
        """
        Input files used to get a good estimate of the average target size

        This is only used for WGS data when the average target size isn't set in the config.
        The access must be computed over the whole genome (no exclude files)
        """
        assert wildcards["library_name"] not in self.tumors, "Autobin always computed on normals"
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam".format(**wildcards)
        input_files = {"bams": [ngs_mapping(tpl)]}
        if self.create_access:
            if self.plain_access:
                input_files["access"] = self.base_out.format(**wildcards) + "access.bed"
            else:
                input_files["target"] = self.base_out.format(**wildcards) + "access.bed"
        return input_files

    def _get_args_autobin(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        assert (
            self.compute_avg_target_size
        ), "Trying to estimate average target size for non-WGS samples"
        args = dict(input) | {"bp-per-bin": 50000}
        if self.plain_access:
            args["method"] = "wgs"
        else:
            args["method"] = "amplicon"
            if "target" not in args:
                args["target"] = self.cfg.path_access
        return args

    # ----- Target --------------------------------------------------------------------------------

    def _get_input_files_target(self, wildcards: Wildcards) -> dict[str, str]:
        """Input files to compute the target regions

        For WES, no input files, it comes from the baits (in arguments) or
        the pon, a previously computed file or the baits (no reference needed)

        For WGS, target is access, with avg size from the config, or 5000 when
        no normal is available (flat prior) or autobin-computed avg size when paired.
        In the latter case, the access must be computed from whole genome
        (no exclude, no min_avg_size)
        """
        assert self.build_ref, "Should not build targets, already available"
        input_files = {}
        if self.is_wgs:
            if self.create_access:
                input_files["interval"] = self.base_out.format(**wildcards) + "access.bed"
            if self.compute_avg_target_size:
                input_files["avg-size"] = self.base_out_lib.format(**wildcards) + "autobin.txt"
        return input_files

    def _get_args_target(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        assert self.build_ref, "Should not build targets, already available"
        if self.is_wes:
            args = {
                "avg-size": self.cfg.target.avg_size,
                "split": self.cfg.target.split,
                "interval": self.path_baits,
            }
        else:
            assert self.is_wgs, "Panel not implemented yet"
            args = dict(input) | {"split": self.cfg.target.split}
            if args.get("avg-size", None) is not None:
                args["avg-size"] = self._read_autobin_output(args["avg-size"])
            elif self.cfg.target.avg_size is not None:
                args["avg-size"] = self.cfg.target.avg_size
            else:
                args["avg-size"] = 5000
        if self.w_config.static_data_config.get("features", None):
            args["annotate"] = self.w_config.static_data_config.features.path
            args["short-names"] = self.cfg.target.short_names
        return args

    # ----- Antitarget ----------------------------------------------------------------------------

    def _get_input_files_antitarget(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {"target": self.base_out.format(**wildcards) + "target.bed"}
        if self.create_access:
            input_files["access"] = self.base_out.format(**wildcards) + "access.bed"
        return input_files

    def _get_args_antitarget(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "avg-size": self.cfg.antitarget.avg_size,
            "min-size": self.cfg.antitarget.min_size,
        }
        if "access" not in args:
            args["access"] = self.cfg.path_access
        return args

    # ----- Coverage ------------------------------------------------------------------------------

    def _get_input_files_coverage(self, wildcards: Wildcards) -> dict[str, str]:
        """
        Compute coverage of region (either target or antitarget)

        Except when region provided with file, the region is computed by the pipeline,
        and must be inculded with the inputs (possibly from the panel_of_normals step).
        For WGS paired, the target regions are sample-dependent, because the optimal
        average target size is sample-dependent (via the rough normal sample coverage).
        In that case, the target regions must be taken from the normal sample, to
        avoid requesting to build targets from the tumor sample.
        """
        # BAM/BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {"bam": ngs_mapping(base_path + ".bam")}

        # Region (target or antitarget) file
        if self.build_ref:
            if self.compute_avg_target_size:
                tpl = self.base_out_lib + "{region}.bed"
                if wildcards["library_name"] in self.tumors:
                    input_files["intervals"] = tpl.format(
                        mapper=wildcards["mapper"],
                        library_name=self.parent.matched_normal[wildcards["library_name"]],
                        region=wildcards["region"],
                    )
                else:
                    input_files["intervals"] = tpl.format(**wildcards)
            else:
                input_files["intervals"] = self.base_out.format(**wildcards) + "{region}.bed"
        elif self.pon_source == PanelOfNormalsOrigin.COHORT:
            panel_of_normals = self.parent.sub_workflows["panel_of_normals_cnvkit"]
            base_path = "output/{mapper}.cnvkit/out/cnvkit.{region}.bed"
            input_files["intervals"] = panel_of_normals(base_path)

        return input_files

    def _get_args_coverage(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "reference": self.w_config.static_data_config.reference.path,
            "min-mapq": self.cfg.coverage.min_mapq,
            "count": self.cfg.coverage.count,
        }
        if "intervals" not in args:
            intervals = self.cfg.panel_of_normals.get("path_{region}".format(**wildcards), "")
            assert intervals != "", "Missing path to {region}".format(**wildcards)
            args["intervals"] = intervals
        return args

    # ----- Reference (flat or pairwise) ----------------------------------------------------------

    def _get_input_files_reference(self, wildcards: Wildcards) -> dict[str, str]:
        """Builds reference from the paired normal, or flat prior in absence of normal"""
        assert self.build_ref, "Should not build reference"
        input_files = {}
        if self.paired:
            input_files["normals"] = [self.base_out_lib.format(**wildcards) + "target.cnn"]
            if self.is_wes:
                input_files["normals"].append(
                    self.base_out_lib.format(**wildcards) + "antitarget.cnn"
                )
        elif self.pon_source == PanelOfNormalsOrigin.FLAT:
            input_files["target"] = self.base_out.format(**wildcards) + "target.bed"
            if self.is_wes:
                input_files["antitarget"] = self.base_out.format(**wildcards) + "antitarget.bed"
        return input_files

    def _get_args_reference(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        assert self.build_ref, "Should not build reference"
        args = dict(input) | {
            "reference": self.w_config.static_data_config.reference.path,
            "cluster": self.cfg.cluster,
            "no-gc": not self.cfg.gc,
            "no-rmask": not self.cfg.rmask,
            "no-edge": not self.cfg.get("edge", self.is_wes),
        }
        if self.cfg.cluster:
            args["min-cluster-size"] = self.cfg.min_cluster_size
        if self.cfg.diploid_parx_genome:
            args["diploid-parx-genome"] = self.cfg.diploid_parx_genome
        sample_sex = self._get_sample_sex(wildcards.get("library_name", None))
        if sample_sex is not None:
            args["sample-sex"] = str(sample_sex)
            if sample_sex == SexValue.MALE and self.paired:
                args["male-reference"] = True
            else:
                args["male-reference"] = self.cfg.male_reference
        else:
            args["male-reference"] = self.cfg.male_reference
        return args

    # ----- Fix -----------------------------------------------------------------------------------

    def _get_input_files_fix(self, wildcards: Wildcards) -> dict[str, str]:
        # Coverage on targets & optionally on antitargets
        input_files = {"target": self.base_out_lib.format(**wildcards) + "target.cnn"}
        if self.is_wes:
            input_files["antitarget"] = self.base_out_lib.format(**wildcards) + "antitarget.cnn"
        if self.paired:
            tpl = "{mapper}.cnvkit.{normal_library}".format(
                mapper=wildcards["mapper"],
                normal_library=self.parent.matched_normal[wildcards["library_name"]],
            )
            input_files["reference"] = os.path.join("work", tpl, "out", tpl + ".reference.cnn")
        elif self.pon_source == PanelOfNormalsOrigin.FLAT:
            input_files["reference"] = self.base_out.format(**wildcards) + "reference.cnn"
        elif self.pon_source == PanelOfNormalsOrigin.COHORT:
            panel_of_normals = self.parent.sub_workflows["panel_of_normals_cnvkit"]
            base_path = "output/{mapper}.cnvkit/out/cnvkit.panel_of_normals.cnn"
            input_files["reference"] = panel_of_normals(base_path)
        return input_files

    def _get_args_fix(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "cluster": self.cfg.cluster,
            "no-gc": not self.cfg.gc,
            "no-rmask": not self.cfg.rmask,
            "no-edge": not self.cfg.get("edge", self.is_wes),
        }
        args["sample-id"] = wildcards.library_name
        if self.cfg.diploid_parx_genome:
            args["diploid-parx-genome"] = self.cfg.diploid_parx_genome
        if "reference" not in args:
            args["reference"] = self.cfg.panel_of_normals.path_panel_of_normals
        return args

    # ----- Variant-related convenience functions -------------------------------------------------

    def _variants_from_cohort_input(self) -> str:
        variants = self.parent.sub_workflows["somatic_variant_calling_cnvkit"]
        tpl = f"{{mapper}}.{self.cfg.somatic_variant_calling.tool}.{{library_name}}"
        base_path = os.path.join("output", tpl, "out", tpl + ".vcf.gz")
        return variants(base_path)

    def _variants_args(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "min-variant-depth": self.cfg.somatic_variant_calling.min_variant_depth,
            "sample-id": wildcards.library_name,
            "normal-id": self.parent.matched_normal[wildcards.library_name],
        }
        if self.cfg.somatic_variant_calling.zygocity_freq is not None:
            args["zygicity-freq"] = self.cfg.somatic_variant_calling.zygocity_freq
        return args

    # ----- Segment -------------------------------------------------------------------------------

    def _get_input_files_segment(self, wildcards: Wildcards) -> dict[str, str]:
        # Coverage
        input_files = {"ratios": self.base_out_lib.format(**wildcards) + "cnr"}
        # Segmentation using SNVs from cohort
        if self.variants_from_cohort:
            input_files["variants"] = self._variants_from_cohort_input().format(**wildcards)
        return input_files

    def _get_args_segment(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        # Segmentation parameters
        args = dict(input) | {
            "method": self.cfg.segment.method,
            "threshold": self.cfg.segment.threshold,
            "drop-outliers": self.cfg.segment.drop_outliers,
            "drop-low-coverage": self.cfg.drop_low_coverage,
        }
        if self.cfg.segment.method == CnvkitSegmentationMethod.CBS:
            args["smooth-cbs"] = self.cfg.segment.smooth_cbs
        if self.cfg.somatic_variant_calling.enabled:
            args |= self._variants_args(wildcards, input)
            if "variants" not in args:
                args["variants"] = self.cfg.somatic_variant_calling.path_somatic_variant_calling
        return args

    # ----- Call ----------------------------------------------------------------------------------

    def _get_input_files_call(self, wildcards: Wildcards) -> dict[str, str]:
        # Segmentation
        input_files = {"segments": self.base_out_lib.format(**wildcards) + "segments.cns"}
        # Segmentation using SNVs from cohort
        if self.variants_from_cohort:
            input_files["variants"] = self._variants_from_cohort_input().format(**wildcards)
        # Purity from the tool
        if (
            self.cfg.somatic_purity_ploidy_estimate.enabled
            and self.cfg.somatic_purity_ploidy_estimate.source == PurityOrigin.COHORT
        ):
            purity = self.parent.sub_workflows["somatic_purity_ploidy_estimate_cnvkit"]
            tpl = f"{{mapper}}.{self.cfg.somatic_purity_ploidy_estimate.tool}.{{library_name}}"
            base_path = os.path.join("output", tpl, "out", tpl + ".txt")
            input_files["purity_file"] = purity(base_path).format(**wildcards)
        return input_files

    def _get_args_call(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        # Call parameters
        args = dict(input) | {
            "method": self.cfg.call.method,
            "thresholds": self.cfg.call.thresholds,
            "drop-low-coverage": self.cfg.drop_low_coverage,
            "male-reference": self.cfg.male_reference,
        }
        if self.cfg.call.center_at is not None:
            args["center-at"] = self.cfg.call.center_at
        else:
            if self.cfg.call.center is not None:
                args["center"] = self.cfg.call.center
        if self.cfg.diploid_parx_genome:
            args["diploid-parx-genome"] = self.cfg.diploid_parx_genome
        if self.cfg.somatic_variant_calling.enabled:
            args |= self._variants_args(wildcards, input)
            if "variants" not in args:
                args["variants"] = self.cfg.somatic_variant_calling.path_somatic_variant_calling
        # Sample sex if known, otherwise guessed by the tool
        sample_sex = self._get_sample_sex(wildcards.library_name)
        if sample_sex is not None:
            args["sample-sex"] = str(sample_sex)
            if sample_sex == SexValue.MALE and self.paired:
                args["male-reference"] = True
        # If requested, purity from samplesheet or from default
        if self.cfg.somatic_purity_ploidy_estimate.enabled:
            if args.get("purity_file", None) is not None:
                (purity, ploidy) = self._read_purity_ploidy_output(args["purity_file"])
            elif self.cfg.somatic_purity_ploidy_estimate.source == PurityOrigin.SAMPLESHEET:
                purity = self.tumors[wildcards.library_name].purity
                ploidy = self.tumors[wildcards.library_name].ploidy
            elif self.cfg.somatic_purity_ploidy_estimate.source == PurityOrigin.CONFIG:
                purity = self.cfg.purity.purity
                ploidy = self.cfg.purity.ploidy
            args["purity"] = purity
            args["ploidy"] = ploidy
        return args

    # ----- Bintest -------------------------------------------------------------------------------

    def _get_input_files_bintest(self, wildcards: Wildcards) -> dict[str, str]:
        return {
            "ratios": self.base_out_lib.format(**wildcards) + "cnr",
            "segments": self.base_out_lib.format(**wildcards) + "segments.cns",
        }

    def _get_args_bintest(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        return dict(input) | {
            "alpha": self.cfg.bintest.alpha,
            "target": self.cfg.bintest.target,
        }

    # ----- Plots --------------------------------------------------------------------------------

    def _get_input_files_scatter(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {
            "ratios": self.base_out_lib.format(**wildcards) + "cnr",
            "segments": self.base_out_lib.format(**wildcards) + "segments.cns",
        }
        if self.variants_from_cohort:
            input_files["variants"] = self._variants_from_cohort_input().format(**wildcards)
        return input_files

    def _get_args_scatter(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "antitarget-marker": self.cfg.scatter.antitarget_marker,
            "by-bin": self.cfg.scatter.by_bin,
            "segment-color": self.cfg.scatter.segment_color,
            "trend": self.cfg.scatter.trend,
            "fig-size": self.cfg.scatter.fig_size,
            "width": self.cfg.scatter.width,
        }
        if self.cfg.scatter.y_min is not None:
            args["y-min"] = self.cfg.scatter.y_min
        if self.cfg.scatter.y_min is not None:
            args["y-min"] = self.cfg.scatter.y_min
        if wildcards["contig_name"] != "all":
            args["chromosome"] = wildcards["contig_name"]
        if self.cfg.somatic_variant_calling.enabled:
            args |= self._variants_args(wildcards, input)
            if "variants" not in args:
                args["variants"] = self.cfg.somatic_variant_calling.path_somatic_variant_calling
        args["title"] = f"{wildcards['library_name']} - {wildcards['contig_name']}"
        return args

    # ----- Metrics (metrics & segmetrics) --------------------------------------------------------

    def _get_input_files_metrics(self, wildcards: Wildcards) -> dict[str, str]:
        return {
            "ratios": self.base_out_lib.format(**wildcards) + "cnr",
            "segments": self.base_out_lib.format(**wildcards) + "segments.cns",
        }

    def _get_args_metrics(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        return dict(input) | {"drop-low-coverage": self.cfg.drop_low_coverage}

    def _get_input_files_segmetrics(self, wildcards: Wildcards) -> dict[str, str]:
        return {
            "ratios": self.base_out_lib.format(**wildcards) + "cnr",
            "segments": self.base_out_lib.format(**wildcards) + "segments.cns",
        }

    def _get_args_segmetrics(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        return dict(input) | {
            "drop-low-coverage": self.cfg.drop_low_coverage,
            "alpha": self.cfg.segmetrics.alpha,
            "bootstrap": self.cfg.segmetrics.bootstrap,
            "smooth-bootstrap": self.cfg.segmetrics.smooth_bootstrap,
            "stats": self.cfg.segmetrics.stats,
        }

    def _get_input_files_genemetrics(self, wildcards: Wildcards) -> dict[str, str]:
        return {
            "ratios": self.base_out_lib.format(**wildcards) + "cnr",
            "segments": self.base_out_lib.format(**wildcards) + "segments.cns",
        }

    def _get_args_genemetrics(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "drop-low-coverage": self.cfg.drop_low_coverage,
            "male-reference": self.cfg.male_reference,
            "threshold": self.cfg.genemetrics.threshold,
            "min-probes": self.cfg.genemetrics.min_probes,
            "alpha": self.cfg.genemetrics.alpha,
            "bootstrap": self.cfg.genemetrics.bootstrap,
            "stats": [x.replace("t-test", "ttest") for x in self.cfg.genemetrics.stats],
        }
        if self.cfg.diploid_parx_genome:
            args["diploid-parx-genome"] = self.cfg.diploid_parx_genome
        sample_sex = self._get_sample_sex(wildcards.library_name)
        if sample_sex is not None:
            args["sample-sex"] = str(sample_sex)
            if sample_sex == SexValue.MALE and self.paired:
                args["male-reference"] = True
        return args

    # ----- Read small files to put values in parameters

    def _read_autobin_output(self, filename: str) -> int:
        nb = r"([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)"
        pattern = re.compile("^Target:[ \t]+" + nb + "[ \t]+" + nb + "$")
        with open(filename) as f:
            for line in f:
                m = pattern.match(line)
                if m:
                    return int(float(m.groups()[4]))
        return -1

    def _read_purity_ploidy_output(self, filename: str) -> tuple[float, float]:
        # TODO: Tool-dependent parsing of purity/ploidy file
        nb = r"([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)"
        pattern = re.compile("^Purity/ploidy:[ \t]+" + nb + "[ \t]+" + nb + "$")
        with open(filename) as f:
            for line in f:
                m = pattern.match(line)
                if m:
                    return (float(m.groups()[1]), float(m.groups()[4]))
        return (-1.0, -1.0)


class LibraryInfo(NamedTuple):
    library: NGSLibrary
    donor: str
    is_tumor: bool
    libraryType: str
    libraryKit: str | None
    sex: Sex | None
    purity: float | None
    ploidy: float = 2


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
        # Collect extra information per library
        self.valid_dna_libraries = {}
        for sheet in self.shortcut_sheets:
            self.valid_dna_libraries |= SomaticCnvCallingWorkflow._get_dna_libraries(sheet)

        # All tumor samples, by libraryKit, with None for WGS
        self.tumors = SomaticCnvCallingWorkflow._split_by(
            SomaticCnvCallingWorkflow._filter_by(
                self.valid_dna_libraries.values(), "is_tumor", lambda x: x
            ),
            "libraryKit",
        )

        self.matched_normal = self._match_normals(self.valid_dna_libraries)

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
        for subworkflow in (
            "panel_of_normals",
            "somatic_variant_calling",
            "somatic_purity_ploidy_estimate",
        ):
            self._optionally_register_subworkflow(subworkflow)

    def get_result_files(self) -> OutputFiles:
        fns = []

        for tool in self.config.tools.wgs:
            for mapper in self.w_config.step_config["ngs_mapping"].tools.dna:
                for library in self.tumors.get(None):
                    fns += self.sub_steps.get(tool).get_result_files(library.library.name, mapper)

        for tool in self.config.tools.wes:
            for mapper in self.w_config.step_config["ngs_mapping"].tools.dna:
                for libraryKit in self.tumors.keys():
                    if libraryKit is None:
                        continue
                    for library in self.tumors.get(libraryKit):
                        fns += self.sub_steps.get(tool).get_result_files(
                            library.library.name, mapper
                        )

        return OutputFiles(fns)

    def _match_normals(self, valid_dna_libraries: list[LibraryInfo]) -> dict[str, str]:
        normals = SomaticCnvCallingWorkflow._split_by(
            SomaticCnvCallingWorkflow._filter_by(
                valid_dna_libraries.values(), "is_tumor", lambda x: not x
            ),
            "libraryKit",
        )

        # Pairing between tumor & normals (must share the same libraryKit)
        matched_normal = {
            sample.library.name: None for samples in self.tumors.values() for sample in samples
        }
        for libraryKit, samples in self.tumors.items():
            if libraryKit in normals:
                normals_by_donor = SomaticCnvCallingWorkflow._split_by(normals[libraryKit], "donor")
                for sample in samples:
                    donor = sample.donor
                    normal = normals_by_donor.get(donor, [])
                    assert (
                        len(normal) < 2
                    ), f"Muliple valid donor samples for tumor library {sample.library.name}"
                    if normal:
                        assert (
                            normal[0].sex == sample.sex
                        ), f"Normal & tumor samples {normal[0].library.name} & {sample.library.name} from donor {donor} have different sex"
                        normal_library = normal[0].library
                        matched_normal[sample.library.name] = normal_library.name
        return matched_normal

    def _optionally_register_subworkflow(self, subworkflow: str):
        for tool in set(self.config.tools.wgs + self.config.tools.wes):
            assert self.config.get(tool) is not None, f"Requested tool '{tool}' not configured"
            cfg = self.config.get(tool)
            subworkflow_config = cfg.get(subworkflow)
            if (
                subworkflow_config
                and subworkflow_config.enabled
                and str(subworkflow_config.source) == "cohort"
            ):
                self.register_sub_workflow(
                    subworkflow,
                    subworkflow_config.get(f"path_{subworkflow}"),
                    f"{subworkflow}_{tool}",
                )

    @staticmethod
    def _get_dna_libraries(sheet) -> dict[str, LibraryInfo]:
        allowed_library_types = [
            k for k, v in LIBRARY_TO_EXTRACTION.items() if v == EXTRACTION_TYPE_DNA
        ]

        valid_dna_libraries = {}
        for donor in sheet.sheet.bio_entities.values():
            sex: SexValue = donor.extra_infos.get("sex", None)
            for bio_sample in donor.bio_samples.values():
                is_tumor = bio_sample.extra_infos.get("isTumor", None)
                assert (
                    is_tumor is not None
                ), f"Missing 'isTumor' value for sample '{donor.name}-{bio_sample.name}'"
                if is_tumor:
                    purity = bio_sample.extra_infos.get("purity", None)
                    ploidy = bio_sample.extra_infos.get("ploidy", 2)
                else:
                    purity = 0
                    ploidy = 2
                for test_sample in bio_sample.test_samples.values():
                    if (
                        test_sample.extra_infos.get("extractionType", "").upper()
                        != EXTRACTION_TYPE_DNA
                    ):
                        continue
                    for library in test_sample.ngs_libraries.values():
                        assert (
                            library.name not in valid_dna_libraries
                        ), f"Duplicate entry for library {library.name}"
                        libraryType = library.extra_infos.get("libraryType", None)
                        assert (
                            libraryType is not None
                        ), f"Missing library type for library '{library.name}'"
                        if libraryType.upper() not in allowed_library_types:
                            continue
                        libraryKit = None
                        if libraryType.upper() == "WES" or libraryType.upper() == "Panel":
                            libraryKit = library.extra_infos.get("libraryKit", None)
                            assert (
                                libraryKit is not None
                            ), f"Missing library kit for library '{library.name}'"
                        valid_dna_libraries[library.name] = LibraryInfo(
                            library,
                            donor.name,
                            is_tumor,
                            libraryType,
                            libraryKit,
                            sex,
                            purity,
                            ploidy,
                        )

        return valid_dna_libraries

    @staticmethod
    def _split_by(
        valid_dna_libraries: list[LibraryInfo], i: str = "library"
    ) -> dict[Any, list[LibraryInfo]]:
        split = {}
        for entry in valid_dna_libraries:
            index = getattr(entry, i)
            if isinstance(index, (int, float, complex, bool)):
                index = str(index)
            if index not in split:
                split[index] = []
            split[index].append(entry)
        return split

    @staticmethod
    def _filter_by(
        valid_dna_libraries: list[LibraryInfo],
        i: str = "library",
        f: Callable[[Any], bool] = lambda x: True,
    ) -> list[LibraryInfo]:
        filtered = []
        for entry in valid_dna_libraries:
            index = getattr(entry, i)
            if f(index):
                filtered.append(entry)
        return filtered