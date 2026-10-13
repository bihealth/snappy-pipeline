# -*- coding: utf-8 -*-
"""GATK gCNV model for WGS data

The step builds a model for GATK4 gCNV from the aligned germline reads of the whole cohort.
All libraries must be whole genomes of the same library kit.
The kit is named ``default``.

The model has two parts, the contig ploidy model and the CNV calls model (``cnv_calls-model``).
Together they let :ref:`step_sv_calling_wgs` run gCNV in case mode.
The step has no ``tool`` field.

Outputs
-------

The model directories stay in ``work/`` and are not published:

- ``default/out/default.contig_ploidy/`` holds the ploidy model (``ploidy-model``) and the ploidy calls per sample.
- ``default/out/default.{shard}_of_{n}.call_cnvs/`` holds, for each shard of the genome, the CNV calls model (``cnv_calls-model``), the calls per sample and the tracking data.

The keys ``ploidy_done`` and ``calls_done`` for downstream tasks are the ``.done`` files in these directories.
The ``gcnv.precomputed_model_paths`` of the ``sv_calling_wgs`` task name these model directories.

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``.
They are the calls on the cohort libraries themselves.

- ``{library_name}/out/{library_name}.post_germline_calls.vcf.gz``
- ``{library_name}/out/{library_name}.post_germline_calls.interval.vcf.gz``
- ``{library_name}/out/{library_name}.post_germline_calls.ratio.tsv``
"""

import dataclasses
import os

from biomedsheets.shortcuts import GermlineCaseSheet
from snakemake.io import glob_wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import BuildGcnvModelStepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_common import build_ngs_library_to_kit

from .model import HelperGcnvModelWgs as HelperGcnvModelWgsConfigModel

#: Default configuration for the helper_gcnv_model_wgs schema


def _filter_intervals_memory(wildcards, input=None, threads=None, attempt=1):
    _, _, _ = wildcards, input, threads  # unused but cannot be renamed
    return f"{attempt * 4 * 1024 + 16 * 1024}MB"


class BuildGcnvWgsModelStepPart(BuildGcnvModelStepPart):
    """Build model for GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        # No mapping given as WGS, we will use the "default" one for all.
        self.ngs_library_to_kit = build_ngs_library_to_kit(
            self.parent.primary_dna_libraries(include_background=True)
        )

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

    resource_usage = {
        **BuildGcnvModelStepPart.resource_usage,
        "filter_intervals": dataclasses.replace(
            BuildGcnvModelStepPart.default_resource_usage, mem=_filter_intervals_memory
        ),
    }


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
            "ploidy_done": "work/default/out/default.contig_ploidy/.done",
            "calls_done": "work/default/out/default.{shard}.call_cnvs/.done",
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
