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
    `-- default
        `-- out
            |-- default.contig_ploidy
            |   |-- SAMPLE_0
            |   |   |-- contig_ploidy.tsv
            |   |   |-- global_read_depth.tsv
            |   |   |-- mu_psi_s_log__.tsv
            |   |   |-- sample_name.txt
            |   |   `-- std_psi_s_log__.tsv
            |   |-- [...]
            |   `-- ploidy-model
            |       |-- contig_ploidy_prior.tsv
            |       |-- gcnvkernel_version.json
            |       |-- interval_list.tsv
            |       |-- mu_mean_bias_j_lowerbound__.tsv
            |       |-- mu_psi_j_log__.tsv
            |       |-- ploidy_config.json
            |       |-- std_mean_bias_j_lowerbound__.tsv
            |       `-- std_psi_j_log__.tsv
            `-- default.***_of_***.call_cnvs
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
from typing import Any

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import glob_wildcards
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import BuildGcnvModelStepPart
from snappy_wrappers.resource_usage import ResourceUsage

from .model import HelperGcnvModelWgs as HelperGcnvModelWgsConfigModel

#: Default configuration for the helper_gcnv_model_wgs schema


class BuildGcnvWgsModelStepPart(BuildGcnvModelStepPart):
    """Build model for GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        # No mapping given as WGS, we will use the "default" one for all.
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library:
                yield donor.dna_ngs_library.name, "default"

    @dictify
    def _get_input_files_call_cnvs(self, wildcards):
        """Yield input files for ``call_cnvs`` in COHORT mode.

        Overwrites original method to set `library_kit` to 'default'. The general assumption is
        that all samples selected for the cohort have the same library kit.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa').
        :type wildcards: snakemake.io.Wildcards
        """
        yield (
            "interval_list_shard",
            "work/default/out/default.scatter_intervals/temp_{shard}/scattered.interval_list",
        )
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            tsvs.append(f"work/{lib}/out/{lib}.coverage.tsv")
        yield "tsv", tsvs
        yield "ploidy", "work/default/out/default.contig_ploidy/.done"
        yield "intervals", "work/default/out/default.annotate_gc.tsv"

    @dictify
    def _get_input_files_post_germline_calls(self, wildcards, checkpoints):
        checkpoint = checkpoints.build_gcnv_model_scatter_intervals
        library_kit = "default"
        scatter_out = checkpoint.get(library_kit=library_kit, **wildcards).output[0]
        shards = list(
            map(
                os.path.basename,
                glob_wildcards(os.path.join(scatter_out, "temp_{shard}/{file}")).shard,
            )
        )
        yield (
            "calls",
            [f"work/{library_kit}/out/{library_kit}.{shard}.call_cnvs/.done" for shard in shards],
        )
        yield "ploidy", f"work/{library_kit}/out/{library_kit}.contig_ploidy/.done"

    def _get_params_preprocess_intervals(self, wildcards: Wildcards) -> dict[str, Any]:
        gcnv_config = self.parent.get_task_config("helper_gcnv_model_wgs").gcnv
        return {
            "path_uniquely_mapable_bed": gcnv_config.path_uniquely_mapable_bed,
        }

    _get_params_coverage = _get_params_preprocess_intervals

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        result = super().get_resource_usage(action)

        def get_memory(wildcards, input=None, threads=None, attempt=1):
            _, _, _ = wildcards, input, threads  # unused but cannot be renamed
            return f"{attempt * 4 * 1024 + 16 * 1024}MB"

        if action == "filter_intervals":
            import dataclasses

            result = dataclasses.replace(
                result,
                mem=get_memory,
            )

        return result


class HelperBuildWgsGcnvModelWorkflow(BaseStep):
    """Perform gCNV model building for WGS samples"""

    #: Workflow name
    name = "helper_gcnv_model_wgs"
    config_model_class = HelperGcnvModelWgsConfigModel
    produces = [DataSignature(DataType.MODELS, frozenset({"gcnv"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local helper gCNV WGS model output paths for downstream consumers."""
        _ = kwargs
        return {
            "ploidy_done": "output/default/out/default.contig_ploidy/.done",
            "calls_done": "output/default/out/default.{shard}.call_cnvs/.done",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                BuildGcnvWgsModelStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the gCNV build model workflow"""
        yield from self.sub_steps["gcnv"].get_result_files()

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
