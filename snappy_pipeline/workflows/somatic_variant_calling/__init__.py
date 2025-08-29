# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_calling`` step

The ``somatic_variant_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs somatic variant calling.  The result are variant files
with somatic variants (bgzip-ed and indexed VCF files).

Usually, the somatic variant calling step is followed by the ``somatic_variant_annotation`` step.

==========
Step Input
==========

The somatic variant calling step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk`` and each read mapper
``mapper`` that the library has been aligned with, and the variant caller ``var_caller``, the
pipeline step will create a directory ``output/{mapper}.{var_caller}.{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

Two ``vcf`` files are produced:

- ``{mapper}.{var_caller}.{lib_name}.vcf.gz`` which contains only the variants that have passed all filters, or that were protected, and
- ``{mapper}.{var_caller}.{lib_name}.full.vcf.gz`` which contains all variants, with the reason for rejection in the ``FILTER`` column.

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.mutect2.P001-N1-DNA1-WES1-4
    |   `-- out
    |       |-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz
    |       |-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz.tbi
    |       |-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz.md5
    |       `-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz.tbi.md5
    [...]

Generally, the callers are set up to return many variants to avoid missing clinically important ones.
They are likely to contain low-quality variants and for some callers variants flagged as being non-somatic.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_calling.rst

=================================
Available Somatic Variant Callers
=================================

The following somatic variant callers are currently available

- ``mutect2`` is the recommended caller
- ``mutect`` & ``scalpel`` are deprecated
- the joint variant callers ``bcftools``, ``platypus``, ``gatk`` & ``varscan`` are unsupported.

==========================
Efficient usage of mutect2
==========================

The recommended caller ``mutect2`` is implemented using a parallelisation mechanism to reduce execution time.
During parallelisation, the genome is split into small chunks, and parallel jobs perform the somatic variant calling in their region only.
All results are then assembled to generate the final output files.

There is at least one chunk by contig defined in the reference genome, with the chunk size upper limit given by the ``window_length`` configuration option.
So large chromosomes can be split into several chunks, while smaller contigs are left in one chunk.
Even for large chunk size, this parallelisation can create hundreds of jobs when the reference genome contains many contigs
(unplaced or unlocalized contigs, viral sequences, decoys, HLA alleles, ...).
Somatic variant calling is generally meaningless for many of these small contigs.
It is possible to configure the somatic variant calling to avoid the contigs irrelevant for downstream analysis, for example:

.. code-block:: yaml

  mutect2:
    ignore_chroms: ['*_random', 'chrUn_*', '*_decoy', "EBV", "HPV*", "HBV", "HCV-*", "HIV-*", "HTLV-1", "CMV", "KSHV", "MCV", "SV40"] # GRCh38.d1.vd1
    window_length: 300000000   # OK for WES, too large for WGS
    keep_tmpdir: onerror       # Facilitates debugging
    ...


=======
Reports
=======

Currently, no reports are generated.
"""

import os
import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand, Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import SomaticVariantCalling as SomaticVariantCallingConfigModel
from .model import TumorNormalMode as TumorNormalMode

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

EXT_MATCHED = {
    "mutect": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
        "txt": ".txt",
        "txt_md5": ".txt.md5",
        "wig": ".wig",
        "wig_md5": ".wig.md5",
    },
    "scalpel": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
        "tar": ".tar.gz",
        "tar_md5": ".tar.gz.md5",
    },
    "mutect2": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
    },
}

#: Available somatic variant callers
SOMATIC_VARIANT_CALLERS = {"mutect2"}

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = SomaticVariantCallingConfigModel.default_config_yaml_string()


class SomaticVariantCallingStepPart(BaseStepPart):
    """Base class for somatic variant calling step parts

    Variant calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/out/"
            "{{mapper}}.{var_caller}.{{tumor_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_input_files(self, action: str):
        """Return generic input function.

        :param action: Action (i.e., step) in the workflow, examples: 'run', 'filter',
        'contamination'.
        :type action: str

        :return: Returns input function based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        """Helper wrapper function"""
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        tumor_base_path = ("output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}").format(
            **wildcards
        )
        input_files = {
            "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
            "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
        }

        normal_library = self.get_normal_lib_name(wildcards)
        if normal_library:
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=normal_library, **wildcards
                )
            )
            input_files.update(
                {
                    "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                    "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                }
            )

        return input_files

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair.get(wildcards.tumor_library, None)
        return pair.normal_sample.dna_ngs_library.name if pair else None

    def get_tumor_lib_name(self, wildcards):
        """Return name of tumor library"""
        pair = self.tumor_ngs_library_to_sample_pair.get(wildcards.tumor_library, None)
        return pair.tumor_sample.dna_ngs_library.name if pair else wildcards.tumor_library

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        # Validate action
        self._validate_action(action)
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, var_caller=[self.name], ext=EXT_VALUES))
        )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)

        prefix = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/log/"
            "{{mapper}}.{var_caller}.{{tumor_library}}"
        ).format(var_caller=self.__class__.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class Mutect2StepPart(SomaticVariantCallingStepPart):
    """Somatic variant calling with Mutect2"""

    #: Step name
    name = "mutect2"

    #: Class available actions
    actions = [
        "scatter",
        "run",
        "gather",
        "filter",
    ]  # "contamination", "pileup_normal", "pileup_tumor")

    #: Class resource usage dictionary. Key: action (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "scatter": ResourceUsage(
            threads=1,
            time="00:02:00",
            memory="1000M",
        ),
        "run": ResourceUsage(
            threads=1,
            time="5-00:00:00",
            memory="8000M",
        ),
        "gather": ResourceUsage(
            threads=1,
            time="03:59:00",
            memory="32768M",
        ),
        "filter": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="15872M",
        ),
        "contamination": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="7680M",
        ),
        "pileup_normal": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="8000M",
        ),
        "pileup_tumor": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="8000M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        run_resource_usage = self.resource_usage_dict["run"]
        self.resource_usage_dict["run"] = ResourceUsage(
            threads=self.config.mutect2.num_cores or run_resource_usage.threads,
            time=run_resource_usage.time,
            memory=run_resource_usage.memory,
        )

    def check_config(self):
        if self.name not in self.config.tools:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        # Return requested function
        return getattr(self, "_get_input_files_{}".format(action))

    def get_params(self, action):
        self._validate_action(action)
        if action == "scatter":
            ignore_chroms = list(
                set(
                    self.w_config.get("ignore_chroms", [])
                    + self.config.get("ignore_chroms", [])
                    + self.config.get(self.name).get("ignore_chroms", [])
                )
            )
            return {
                "ignore_chroms": sorted(list(ignore_chroms)),
                "padding": self.config.mutect2.padding,
            }
        return {}

    def _get_input_files_scatter(self, wildcards):
        return {"fai": self.w_config.static_data_config.reference.path + ".fai"}

    def _get_input_files_run(self, wildcards):
        """Get input files for rule ``run``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'run', BAM and BAI files.
        """
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        tumor_base_path = ("output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}").format(
            **wildcards
        )

        scatteritem_base_path = "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}/mutect2par/scatter/{scatteritem}".format(
            **wildcards
        )

        input_files = {
            "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
            "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            "region": scatteritem_base_path + ".region.bed",
        }

        # Adjustment tumor_only mode

        tumor_normal_mode = self.config.get(self.name, {}).get("tumor_normal_mode")
        if tumor_normal_mode is None:
            raise ValueError(
                "'tumor_normal_mode' not defined in step configuration for %s" % self.name
            )

        if tumor_normal_mode == TumorNormalMode.AUTOMATIC:
            normal_library = self.get_normal_lib_name(wildcards)
            if normal_library:
                # Use paired mode if normal is present
                normal_base_path = (
                    "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                        normal_library=normal_library, **wildcards
                    )
                )
                input_files.update(
                    {
                        "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                        "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                    }
                )
            # else: tumor_only — do not add normal BAMs
        elif tumor_normal_mode == TumorNormalMode.PAIRED:
            normal_library = self.get_normal_lib_name(wildcards)
            if not normal_library:
                raise ValueError("Normal sample required but not found.")
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=normal_library, **wildcards
                )
            )
            input_files.update(
                {
                    "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                    "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                }
            )
        elif tumor_normal_mode == TumorNormalMode.TUMOR_ONLY:
            # No normal BAMs to include
            pass
        else:
            raise ValueError(f"Unsupported tumor_normal_mode: {tumor_normal_mode}")

        input_files["reference"] = self.w_config.static_data_config.reference.path
        return input_files

    def _get_input_files_gather(self, wildcards):
        gather = self.parent.workflow.globals.get("gather")
        gather = getattr(gather, self.name)
        scatteritem_base_path = "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}/mutect2par/run/{{scatteritem}}".format(
            **wildcards
        )
        input_files = {
            "vcf": scatteritem_base_path + ".raw.vcf.gz",
            "stats": scatteritem_base_path + ".raw.vcf.stats",
            "f1r2": scatteritem_base_path + ".raw.f1r2.tar.gz",
        }

        return dict(map(lambda item: (item[0], gather(item[1])), input_files.items()))

    def _get_input_files_filter(self, wildcards):
        """Get input files for rule ``filter``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'filter'.
        """
        base_path = (
            "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}".format(
                **wildcards
            )
        )
        input_files = {
            "raw": base_path + ".raw.vcf.gz",
            "stats": base_path + ".raw.vcf.stats",
            "orientation": base_path + ".raw.read_orientation_model.tar.gz",
            "reference": self.w_config.static_data_config.reference.path,
        }
        if self.get_normal_lib_name(wildcards):
            if "contamination" in self.actions:
                input_files["table"] = base_path + ".contamination.tbl"
                input_files["segments"] = base_path + ".segments.tbl"
        return input_files

    def _get_input_files_pileup_normal(self, wildcards):
        """Get input files for rule ``pileup_normal``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'pileup_normal', BAM and BAI files.
        """
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        base_path = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
            normal_library=self.get_normal_lib_name(wildcards), **wildcards
        )
        return {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam"),
            "reference": self.w_config.static_data_config.reference.path,
            "common_variants": self.config.mutect2.common_variants,
        }

    def _get_input_files_pileup_tumor(self, wildcards):
        """Get input files for rule ``pileup_tumor``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'pileup_tumor', BAM and BAI files.
        """
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}".format(
            **wildcards
        )
        return {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam"),
            "reference": self.w_config.static_data_config.reference.path,
            "common_variants": self.config.mutect2.common_variants,
        }

    def _get_input_files_contamination(self, wildcards: Wildcards):
        """Get input files for rule ``contamination``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'contamination', Normal and Tumor
        pileup files.
        """
        base_path = (
            "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}".format(
                **wildcards
            )
        )
        return {
            "normal": base_path + ".normal.pileup",
            "tumor": base_path + ".tumor.pileup",
            "reference": self.w_config.static_data_config.reference.path,
        }

    def get_output_files(self, action):
        """Get output files for Mutect2 rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns dictionary with expected output files based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Initialise variables
        exts = {}
        output_files = {}

        # Validate action
        self._validate_action(action)

        # Set expected extensions and basepath based on action
        base_path_out = self.base_path_out

        if action == "scatter":
            scatter = self.parent.workflow.globals.get("scatter")
            scatter = getattr(scatter, self.name)
            template = "work/{{{{mapper}}}}.{var_caller}.{{{{tumor_library}}}}/out/{{{{mapper}}}}.{var_caller}.{{{{tumor_library}}}}/{var_caller}par/scatter/{{scatteritem}}.region.bed".format(
                var_caller=self.name
            )
            return {"regions": scatter(template)}

        if action == "run":
            base_path_out = "work/{{mapper}}.{var_caller}.{{tumor_library}}/out/{{mapper}}.{var_caller}.{{tumor_library}}/{var_caller}par/run/{{scatteritem}}{ext}"
            exts = {
                "raw": ".raw.vcf.gz",
                "raw_md5": ".raw.vcf.gz.md5",
                "raw_tbi": ".raw.vcf.gz.tbi",
                "raw_tbi_md5": ".raw.vcf.gz.tbi.md5",
                "stats": ".raw.vcf.stats",
                "stats_md5": ".raw.vcf.stats.md5",
                "f1r2": ".raw.f1r2.tar.gz",
                "f1r2_md5": ".raw.f1r2.tar.gz.md5",
            }
        if action == "gather":
            exts = {
                "vcf": ".raw.vcf.gz",
                "vcf_md5": ".raw.vcf.gz.md5",
                "vcf_tbi": ".raw.vcf.gz.tbi",
                "vcf_tbi_md5": ".raw.vcf.gz.tbi.md5",
                "stats": ".raw.vcf.stats",
                "stats_md5": ".raw.vcf.stats.md5",
                "orientation": ".raw.read_orientation_model.tar.gz",
                "orientation_md5": ".raw.read_orientation_model.tar.gz.md5",
            }
        if action == "filter":
            exts = {
                "full_vcf": ".full.vcf.gz",
                "full_vcf_md5": ".full.vcf.gz.md5",
                "full_vcf_tbi": ".full.vcf.gz.tbi",
                "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
                "vcf": ".vcf.gz",
                "vcf_md5": ".vcf.gz.md5",
                "vcf_tbi": ".vcf.gz.tbi",
                "vcf_tbi_md5": ".vcf.gz.tbi.md5",
            }
        if action == "contamination":
            exts = {
                "table": ".contamination.tbl",
                "table_md5": ".contamination.tbl.md5",
                "segments": ".segments.tbl",
                "segments_md5": ".segments.tbl.md5",
            }
        if action == "pileup_normal":
            exts = {"pileup": ".normal.pileup", "pileup_md5": ".normal.pileup.md5"}
        if action == "pileup_tumor":
            exts = {"pileup": ".tumor.pileup", "pileup_md5": ".tumor.pileup.md5"}

        # Define output dictionary
        for k, v in exts.items():
            output_files[k] = base_path_out.format(var_caller=self.name, ext=v)
        return output_files

    def get_log_file(self, action):
        """Get log files for Mutect2 rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns dictionary with expected log files based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Initialise variables
        postfix = ""
        log_files = {}
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )

        # Validate action
        self._validate_action(action)

        # Set expected format based on action
        if action != "gather":
            if action == "run":
                postfix = ".{{scatteritem}}"
            else:
                postfix = "." + action

        prefix = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/log/"
            "{{mapper}}.{var_caller}.{{tumor_library}}{postfix}"
        ).format(var_caller=self.name, postfix=postfix)

        # Define output dictionary
        for key, ext in key_ext:
            log_files[key] = prefix + ext
            log_files[key + "_md5"] = prefix + ext + ".md5"
        return log_files

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return self.resource_usage_dict.get(action)


class SomaticVariantCallingWorkflow(BaseStep):
    """Perform somatic variant calling"""

    #: Workflow name
    name = "somatic_variant_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
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
            config_model_class=SomaticVariantCallingConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                Mutect2StepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

        if "mutect2" in self.config.tools:
            if self.config.mutect2.common_variants:
                self.sub_steps["mutect2"].actions.extend(
                    ["contamination", "pileup_normal", "pileup_tumor"]
                )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.{caller}.{tumor_library.name}"
        for caller in set(self.config.tools) & set(SOMATIC_VARIANT_CALLERS):
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller=caller,
                ext=EXT_MATCHED[caller].values() if caller in EXT_MATCHED else EXT_VALUES,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller=caller,
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for bio_entity in sheet.sheet.bio_entities.values():
                for bio_sample in bio_entity.bio_samples.values():
                    if not bio_sample.extra_infos.get("isTumor", False):
                        continue
                    for test_sample in bio_sample.test_samples.values():
                        extraction_type = test_sample.extra_infos.get("extractionType", "unknown")
                        if extraction_type.lower() != "dna":
                            if extraction_type == "unknown":
                                msg = "INFO: sample {} has missing extraction type, ignored"
                                print(msg.format(test_sample.name), file=sys.stderr)
                            continue
                        for ngs_library in test_sample.ngs_libraries.values():
                            yield from expand(tpl, tumor_library=[ngs_library], **kwargs)
