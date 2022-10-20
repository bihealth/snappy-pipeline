# -*- coding: utf-8 -*-
"""Implementation of the ``helper_gcnv_model_wgs`` step

The ``helper_gcnv_model_wgs`` step takes as the input the results of the ``ngs_mapping`` step
(aligned germline reads) and builds a model that can be used by GATK4 gCNV. Important: the workflow
assumes that all samples in the cohort use the same library kit and all are WGS.

==========
Step Input
==========

The step uses Snakemake sub workflows for the result of the ``ngs_mapping``
(aligned reads BAM files).

===========
Step Output
===========

All donors will be used to generate the two parts of the required gCNV model, specifically:
``ploidy-model`` and ``cnv_calls-model``. Both are required to execute gCNV in CASE mode.

For example, the relevant directories might look as follows:

::

    work/
    +-- bwa.gcnv_contig_ploidy.default
        `-- out
            `-- bwa.gcnv_contig_ploidy.default
                |-- SAMPLE_0
                |   |-- contig_ploidy.tsv
                |   |-- global_read_depth.tsv
                |   |-- mu_psi_s_log__.tsv
                |   |-- sample_name.txt
                |   `-- std_psi_s_log__.tsv
                |-- [...]
                `-- bwa.gcnv_contig_ploidy.default
                    `-- ploidy-model
                        |-- contig_ploidy_prior.tsv
                        |-- gcnvkernel_version.json
                        |-- interval_list.tsv
                        |-- mu_mean_bias_j_lowerbound__.tsv
                        |-- mu_psi_j_log__.tsv
                        |-- ploidy_config.json
                        |-- std_mean_bias_j_lowerbound__.tsv
                        `-- std_psi_j_log__.tsv
    +-- bwa.gcnv_call_cnvs.default.***_of_***
        `-- out
            `-- bwa.gcnv_call_cnvs.default.***_of_***
                |-- cnv_calls-calls
                |   |-- SAMPLE_0
                |       `-- [...]
                |   |-- [...]
                |-- cnv_calls-model
                |  |-- denoising_config.json
                |  |-- gcnvkernel_version.json
                |  |-- interval_list.tsv
                |  |-- log_q_tau_tk.tsv
                |  |-- mu_W_tu.tsv
                |  |-- mu_ard_u_log__.tsv
                |  |-- mu_log_mean_bias_t.tsv
                |  |-- mu_psi_t_log__.tsv
                |  |-- std_W_tu.tsv
                |  |-- std_ard_u_log__.tsv
                |  |-- std_log_mean_bias_t.tsv
                |  `-- std_psi_t_log__.tsv
                `-- cnv_calls-tracking
                    `-- [...]

====================
Global Configuration
====================

- At the moment, no global configuration is used.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_helper_gcnv_model_wgs.rst

"""

import os

import attr
from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand, glob_wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, ResourceUsage
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.wgs_cnv_calling import GcnvWgsStepPart

#: Default configuration for the helper_gcnv_model_wgs schema
DEFAULT_CONFIG = r"""
# Default configuration helper_gcnv_model_wgs
step_config:
  helper_gcnv_model_wgs:
    path_ngs_mapping: ../ngs_mapping  # REQUIRED

    gcnv:
      # Path to BED file with uniquely mappable regions.
      path_uniquely_mapable_bed: null  # REQUIRED
"""


class BuildGcnvModelStepPart(GcnvWgsStepPart):
    """Build model for GATK4 gCNV"""

    #: Step name
    name = "gcnv"

    #: Class available actions
    actions = (
        "preprocess_intervals",
        "annotate_gc",
        "filter_intervals",
        "scatter_intervals",
        "coverage",
        "contig_ploidy",
        "scatter_intervals",
        "call_cnvs_cohort_mode",
        "post_germline_calls",
        "post_germline_calls_cohort_mode",
        "merge_cohort_vcfs",
    )

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "high_resource": ResourceUsage(
            threads=16,
            time="2-00:00:00",
            memory="46080M",
        ),
        "default": ResourceUsage(
            threads=1,
            time="04:00:00",
            memory="7680M",
        ),
    }

    @dictify
    def _get_input_files_call_cnvs_cohort_mode(self, wildcards):
        """Yield input files for ``call_cnvs`` in COHORT mode.

        Overwrites original method to set `library_kit` to 'default'. The general assumption is
        that all samples selected for the cohort have the same library kit.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa').
        :type wildcards: snakemake.io.Wildcards
        """
        path_pattern = (
            "work/{name_pattern}/out/{name_pattern}/temp_{{shard}}/scattered.interval_list"
        )
        name_pattern = "{mapper}.gcnv_scatter_intervals.default"
        yield "interval_list_shard", path_pattern.format(name_pattern=name_pattern)
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            path_pattern = "{mapper}.gcnv_coverage.{library_name}".format(
                mapper=wildcards.mapper, library_name=lib
            )
            tsvs.append(
                "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                    name_pattern=path_pattern, ext=ext
                )
            )
        yield ext, tsvs
        ext = "ploidy"
        path_pattern = "{mapper}.gcnv_contig_ploidy.default".format(**wildcards)
        yield ext, "work/{name_pattern}/out/{name_pattern}/.done".format(name_pattern=path_pattern)
        key = "intervals"
        path_pattern = "gcnv_annotate_gc.default"
        yield key, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=path_pattern, ext="tsv"
        )

    @dictify
    def _get_input_files_post_germline_calls_cohort_mode(self, wildcards, checkpoints):
        checkpoint = checkpoints.build_gcnv_model_scatter_intervals
        library_kit = "default"
        scatter_out = checkpoint.get(library_kit=library_kit, **wildcards).output[0]
        shards = list(
            map(
                os.path.basename,
                glob_wildcards(os.path.join(scatter_out, "temp_{shard}/{file}")).shard,
            )
        )
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield "calls", [
            "work/{name_pattern}.{shard}/out/{name_pattern}.{shard}/.done".format(
                name_pattern=name_pattern, shard=shard
            )
            for shard in shards
        ]
        ext = "ploidy"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield ext, "work/{name_pattern}/out/{name_pattern}/.done".format(name_pattern=name_pattern)

    def get_cnv_model_result_files(self, _unused=None):
        """Get gCNV model results.

        :return: Returns list of result files for the gCNV build model workflow.
        """
        # Initialise variables
        name_pattern_contig_ploidy = (
            "work/{mapper}.gcnv_contig_ploidy.{library_kit}/out/"
            "{mapper}.gcnv_contig_ploidy.{library_kit}/.done"
        )
        name_pattern_filter_intervals = (
            "work/{mapper}.gcnv_filter_intervals.{library_kit}/out/"
            "{mapper}.gcnv_filter_intervals.{library_kit}.interval_list"
        )
        chosen_kits = ["default"]  # as defined in ``GcnvWgsStepPart._build_ngs_library_to_kit()``

        yield from expand(
            name_pattern_contig_ploidy,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            library_kit=chosen_kits,
        )
        yield from expand(
            name_pattern_filter_intervals,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            library_kit=chosen_kits,
        )

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        result = super().get_resource_usage(action)

        def get_memory(wildcards, input=None, threads=None, attempt=None):
            _, _, _ = wildcards, input, threads  # unused but cannot be renamed
            return f"{attempt * 4 * 1024 + 16 * 1024}M"

        if action == "filter_intervals":
            result = attr.evolve(
                result,
                memory=get_memory,
            )

        return result


class HelperBuildWgsGcnvModelWorkflow(BaseStep):
    """Perform gCNV model building for WGS samples"""

    #: Workflow name
    name = "helper_gcnv_model_wgs"

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((BuildGcnvModelStepPart,))
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the gCNV build model workflow."""
        ext_values = (".ratio.tsv", ".interval.vcf.gz", ".vcf.gz")
        donors = {donor.name for donor in self.all_donors() if donor.dna_ngs_library}
        name_pattern = "{mapper}.gcnv_post_germline_calls.{index.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("work", name_pattern, "out", name_pattern + "{ext}"),
            donors,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            ext=ext_values,
        )

    def _yield_result_files(self, tpl, donors, **kwargs):
        """Build output paths from path template and extension list.

        Will only yield the result files for pedigrees where the index is in ``donors``.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if pedigree.index.name in donors:
                    yield from expand(tpl, index=[pedigree.index], **kwargs)

    @listify
    def all_donors(self, include_background=True):
        """Get all donors.

        :param include_background: Boolean flag to defined if background should be included or not.
        Default: True, i.e., background will be included.

        :return: Returns list of all donors in sample sheet.
        """
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        self.ensure_w_config(
            ("step_config", "helper_gcnv_model_wgs", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for gCNV model building.",
        )
