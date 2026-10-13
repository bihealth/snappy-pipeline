# -*- coding: utf-8 -*-
"""GATK gCNV model for targeted sequencing data

The step builds a model for GATK4 gCNV from the aligned germline reads of the whole cohort.
It builds one model per library kit.
``gcnv.path_target_interval_list_mapping`` maps the ``libraryKit`` of each library to a kit name and its target intervals.
Without a mapping, all libraries belong to the kit ``default``.

The model has two parts, the contig ploidy model and the CNV calls model (``cnv_calls-model``).
Together they let :ref:`step_sv_calling_targeted` run gCNV in case mode.
The step has no ``tool`` field.

Outputs
-------

The model directories stay in ``work/`` and are not published:

- ``{library_kit}/out/{library_kit}.contig_ploidy/`` holds the ploidy model (``ploidy-model``) and the ploidy calls per sample.
- ``{library_kit}/out/{library_kit}.{shard}_of_{n}.call_cnvs/`` holds, for each shard of the genome, the CNV calls model (``cnv_calls-model``), the calls per sample and the tracking data.

The keys ``ploidy_done`` and ``calls_done`` for downstream tasks are the ``.done`` files in these directories.
The ``gcnv.precomputed_model_paths`` of the ``sv_calling_targeted`` task name these model directories.

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``.
They are the calls on the cohort libraries themselves.

- ``{library_name}/out/{library_name}.post_germline_calls.vcf.gz``
- ``{library_name}/out/{library_name}.post_germline_calls.interval.vcf.gz``
- ``{library_name}/out/{library_name}.post_germline_calls.ratio.tsv``
"""

import os

from biomedsheets.shortcuts import GermlineCaseSheet
from snakemake.io import glob_wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import BuildGcnvModelStepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_common import build_ngs_library_to_kit

from .model import HelperGcnvModelTargeted as HelperGcnvModelTargetedConfigModel

#: Default configuration for the helper_gcnv_model_targeted schema


class BuildGcnvTargetSeqModelStepPart(BuildGcnvModelStepPart):
    """Build model for GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        gcnv_config = self.parent.get_task_config("helper_gcnv_model_targeted").gcnv
        self.ngs_library_to_kit = build_ngs_library_to_kit(
            self.parent.primary_dna_libraries(include_background=True),
            gcnv_config.path_target_interval_list_mapping,
        )

    @dictify
    def _get_input_files_post_germline_calls(self, wildcards, checkpoints):
        checkpoint = checkpoints.build_gcnv_model_scatter_intervals
        library_kit = self.ngs_library_to_kit.get(wildcards.library_name)
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


class HelperBuildTargetSeqGcnvModelWorkflow(BaseStep):
    """Perform gCNV model building for WES samples by library kit"""

    #: Workflow name
    name = "helper_gcnv_model_targeted"
    config_model_class = HelperGcnvModelTargetedConfigModel
    produces = [DataSignature(DataType.MODELS, frozenset({"gcnv"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local helper gCNV model output paths for downstream consumers."""
        kit = kwargs.get("library_kit", "{library_kit}")
        return {
            "ploidy_done": f"work/{kit}/out/{kit}.contig_ploidy/.done",
            "calls_done": f"work/{kit}/out/{kit}.{{shard}}.call_cnvs/.done",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                BuildGcnvTargetSeqModelStepPart,
            )
        )

        # Build mapping from NGS DNA library to library kit
        self.ngs_library_to_kit = self.sub_steps["gcnv"].ngs_library_to_kit
