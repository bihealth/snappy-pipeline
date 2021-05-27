# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_filtration`` step

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_filtration.rst

=========
Important
=========

Because the EB Filter step is so time consuming, the data going
can be heavily prefiltered! (e.g. using Jannovar with the offExome flag).

TODO: document filter, for now see the eb_filter wrapper!

=======
Concept
=======
All variants are annotated with the dkfz-bias-filter to remove sequencing
and PCR artifacts. The variants annotatated with EBFilter are variable, i.e.
only variants that have the PASS flag set because we assume only those will
be kept.

We borrowed the general workflow from variant_filtration, i.e. working with
pre-defined filter sets and exon/region lists.

========
Workflow
========

* 1. Do the filtering genome wide (this file needs to be there, always)
- dkfz-ebfilter-filterset1-genomewide

* 2. optionally, subset to regions defined in bed file, which return
- dkfz-ebfilter-filterset1-regions1

and so on for filterset1 to n

filterset1:
filter bPcr, bSeq flags from dkfz-bias-filter

filterset2:
additionally filter variants with EBscore < x, x is configurable
"""

from collections import OrderedDict
import os
import random
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_annotation
step_config:
  somatic_variant_filtration:
    drmaa_snippet: ''  # default, you can override by step below
    path_somatic_variant_annotation: ../somatic_variant_annotation
    path_ngs_mapping: ../ngs_mapping
    tools_ngs_mapping: null
    tools_somatic_variant_calling: null
    filter_sets:
    # no_filter: no_filters    # implicit, always defined
      dkfz_only: ''  # empty
      dkfz_and_ebfilter:
        ebfilter_threshold: 2.4
      dkfz_and_ebfilter_and_oxog:
        vaf_threshold: 0.08
        coverage_threshold: 5
      dkfz_and_oxog:
        vaf_threshold: 0.08
        coverage_threshold: 5
    exon_lists: {}
    # genome_wide: null         # implicit, always defined
    # ensembl74: path/to/ensembl47.bed
    eb_filter:
      shuffle_seed: 1
      panel_of_normals_size: 25
      min_mapq: 20
      min_baseq: 15
      # Parallelization configuration
      drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
      window_length: 10000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_drmaa: true           # use drmaa for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2    # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
"""


class SomaticVariantFiltrationStepPart(BaseStepPart):
    """Shared code for all tools in somatic_variant_filtration"""

    def __init__(self, parent):
        super().__init__(parent)
        self.log_path = (
            r"work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            r"dkfz_bias_filter.{tumor_library,[^\.]+}/log/snakemake.dkfz_bias_filter.log"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor

    @dictify
    def _get_log_file(self, action):
        """Return path to log file for the given action"""
        assert action in self.actions, "Invalid action"
        if action == "write_panel":
            return (
                "work/{mapper}.eb_filter.panel_of_normals/log/"
                "{mapper}.eb_filter.panel_of_normals.log"
            )
        else:
            name_pattern = self.token
            key_ext = (
                ("log", ".log"),
                ("conda_info", ".conda_info.txt"),
                ("conda_list", ".conda_list.txt"),
            )
            for key, ext in key_ext:
                yield key, os.path.join("work", name_pattern, "log", name_pattern + ext)

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_params(self, _action):
        """Return arguments to pass down."""

        def params_function(wildcards):
            if wildcards.tumor_library not in self.donors:
                return {
                    "tumor_library": wildcards.tumor_library,
                    "normal_library": self.get_normal_lib_name(wildcards),
                }
            else:
                return {}

        return params_function


class DkfzBiasFilterStepPart(SomaticVariantFiltrationStepPart):
    """Flag variants with the DKFZ bias filter"""

    name = "dkfz_bias_filter"

    def __init__(self, parent):
        super().__init__(parent)
        self.actions = ("run",)
        self.token = (
            "{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.{tumor_library}"
        )

    @dictify
    def get_input_files(self, action):
        """Return path to jannovar-annotated vcf input file"""
        assert action == "run"
        # VCF file and index
        tpl = (
            "output/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf.{tumor_library}/out/"
            "{mapper}.{var_caller}.jannovar_annotate_somatic_vcf.{tumor_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        variant_annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        for key, ext in key_ext.items():
            yield key, variant_annotation(tpl + ext)
        # BAM file and index
        tpl = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}"
        key_ext = {"bam": ".bam", "bai": ".bam.bai"}
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for key, ext in key_ext.items():
            yield key, ngs_mapping(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        assert action == "run"
        prefix = (
            r"work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            r"dkfz_bias_filter.{tumor_library,[^\.]+}/out/{mapper}.{var_caller}."
            r"jannovar_annotate_somatic_vcf.dkfz_bias_filter.{tumor_library}"
        )
        key_ext = {
            "vcf": ".vcf.gz",
            "tbi": ".vcf.gz.tbi",
            "vcf_md5": ".vcf.gz.md5",
            "tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in key_ext.items():
            yield key, prefix + ext

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["somatic_variant_filtration_dkfz_bias_filter_run"] = {
            "mem": 3 * 1024,
            "time": "72:00",
            "ntasks": 1,
        }


class EbFilterStepPart(SomaticVariantFiltrationStepPart):
    """Flag variants with EBFilter"""

    name = "eb_filter"

    def __init__(self, parent):
        super().__init__(parent)
        self.actions = ("run", "write_panel")
        self.token = (
            "{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}"
        )

    def get_input_files(self, action):
        assert action in self.actions
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_run(self, wildcards):
        # VCF file and index
        tpl = (
            "work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.{tumor_library}/out/{mapper}.{var_caller}."
            "jannovar_annotate_somatic_vcf.dkfz_bias_filter."
            "{tumor_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, tpl.format(**wildcards) + ext
        # BAM file and index
        tpl = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}"
        key_ext = {"bam": ".bam", "bai": ".bam.bai"}
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for key, ext in key_ext.items():
            yield key, ngs_mapping(tpl.format(**wildcards) + ext)
        # Panel of normals TXT file
        yield "txt", self._get_output_files_write_panel()["txt"].format(**wildcards)

    def _get_input_files_write_panel(self, wildcards):
        bam_paths = self._get_panel_of_normal_bams(wildcards)
        return {"bam": bam_paths, "bai": [p + ".bai" for p in bam_paths]}

    def get_output_files(self, action):
        """Return output files for the filtration"""
        assert action in self.actions
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_run(self):
        prefix = (
            r"work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            r"dkfz_bias_filter.eb_filter.{tumor_library,[^\.]+}/out/"
            r"{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            r"dkfz_bias_filter.eb_filter.{tumor_library}"
        )
        key_ext = {
            "vcf": ".vcf.gz",
            "tbi": ".vcf.gz.tbi",
            "vcf_md5": ".vcf.gz.md5",
            "tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in key_ext.items():
            yield key, prefix + ext

    @dictify
    def _get_output_files_write_panel(self):
        # TODO: add the actual normal sample here?!
        yield "txt", (
            "work/{mapper}.eb_filter.panel_of_normals/out/{mapper}.eb_filter."
            "panel_of_normals.txt"
        )

    def write_panel_of_normals_file(self, wildcards):
        """Write out file with paths to panels-of-normal"""
        output_path = self.get_output_files("write_panel")["txt"].format(**wildcards)
        with open(output_path, "wt") as outf:
            for bam_path in self._get_panel_of_normal_bams(wildcards):
                print(bam_path, file=outf)

    @listify
    def _get_panel_of_normal_bams(self, wildcards):
        """Return list of "panel of normal" BAM files."""
        libraries = []
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    if not bio_sample.extra_infos["isTumor"]:
                        libraries.append(bio_sample.dna_ngs_library.name)
        libraries.sort()
        random.seed(self.config["eb_filter"]["shuffle_seed"])
        lib_count = self.config["eb_filter"]["panel_of_normals_size"]
        random.shuffle(libraries)
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}"
        for library in libraries[:lib_count]:
            yield ngs_mapping(tpl.format(normal_library=library, **wildcards) + ".bam")

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["somatic_variant_filtration_eb_filter_run"] = {
            "mem": 8 * 1024,
            "time": "144:00",
            "ntasks": 1,
        }


class ApplyFiltersStepPartBase(SomaticVariantFiltrationStepPart):
    """Base class for the different filters."""

    name = None

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = (
            "{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}.{exon_list}"
        )
        self.base_path_out = os.path.join("work", name_pattern, "out", name_pattern + "{ext}")
        self.path_log = os.path.join("work", name_pattern, "log", name_pattern + ".log")

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_filtration_{}_run".format(self.name)] = {
            "mem": int(3.75 * 1024 * 2),
            "time": "01:00",
            "ntasks": 2,
        }


class ApplyFiltersStepPart(ApplyFiltersStepPartBase):
    """Apply the configured filters."""

    name = "apply_filters"

    def get_args(self, action):
        def args_function(wildcards):
            result = {
                "normal_sample": self.get_normal_lib_name(wildcards),
                "tumor_sample": wildcards.tumor_library,
            }
            return result

        assert action == "run"
        return args_function

    @dictify
    def get_input_files(self, action):
        assert action == "run", "Unsupported actions"
        tpl = (
            "work/{mapper}.{var_caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}/out/{mapper}.{var_caller}."
            "jannovar_annotate_somatic_vcf.dkfz_bias_filter.eb_filter."
            "{tumor_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, tpl + ext

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{step}", self.name).replace(
                "{exon_list}", "genome_wide"
            ).replace("{ext}", ext)

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log.replace("{step}", self.name).replace("{exon_list}", "genome_wide")


class FilterToExonsStepPart(ApplyFiltersStepPartBase):
    """Apply the configured filters."""

    name = "filter_to_exons"

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                yield key, self.base_path_out.format(
                    step="apply_filters",
                    mapper=wildcards.mapper,
                    var_caller=wildcards.var_caller,
                    filter_set=wildcards.filter_set,
                    exon_list="genome_wide",
                    ext=ext,
                )

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{step}", "filter_to_exons").replace("{ext}", ext)

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log.replace("{step}", self.name)


class SomaticVariantFiltrationWorkflow(BaseStep):
    """Perform somatic variant filtration"""

    name = "somatic_variant_filtration"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow,
            config,
            cluster_config,
            config_lookup_paths,
            config_paths,
            workdir,
            (SomaticVariantAnnotationWorkflow, SomaticVariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                DkfzBiasFilterStepPart,
                EbFilterStepPart,
                ApplyFiltersStepPart,
                FilterToExonsStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow(
            "somatic_variant_annotation", self.config["path_somatic_variant_annotation"]
        )
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_somatic_variant_calling"]:
            self.config["tools_somatic_variant_calling"] = self.w_config["step_config"][
                "somatic_variant_calling"
            ]["tools"]

    @listify
    def get_result_files(self):
        """Return list of result files
        Process all primary DNA libraries and perform pairwise calling for tumor/normal pairs
        """
        callers = set(self.config["tools_somatic_variant_calling"])
        name_pattern = (
            "{mapper}.{caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library.name}."
            "{filter_set}.{exon_list}"
        )
        filter_sets = ["no_filter"]
        filter_sets += self.config["filter_sets"].keys()
        exon_lists = ["genome_wide"]
        exon_lists += list(self.config["exon_lists"].keys())
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            filter_set=filter_sets,
            exon_list=exon_lists,
            ext=EXT_VALUES,
        )
        # TODO: filtration for joint calling not implemented yet

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl, tumor_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "somatic_variant_filtration", "path_somatic_variant_annotation"),
            "Path to variant calling not configured but required for somatic variant annotation",
        )
