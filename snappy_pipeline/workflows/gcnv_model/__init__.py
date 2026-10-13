# -*- coding: utf-8 -*-
"""GATK gCNV model

The step builds a model for GATK4 gCNV in cohort mode from the aligned germline reads of the whole cohort, background data sets included.
The model has two parts, the contig ploidy model (``ploidy-model``) and the CNV calls model (``cnv_calls-model``) of each shard of the genome.
The ``gcnv`` tool of :ref:`step_sv_calling_targeted` and :ref:`step_sv_calling_wgs` runs gCNV in case mode with this model.
These tasks name the model task in ``depends_on.gcnv_model``.
The step has no ``tool`` field.

The field ``intervals`` selects how the genome is binned:

- ``targeted`` bins the target regions of each library kit, for exomes and panels.
  ``gcnv.path_target_interval_list_mapping`` maps the ``libraryKit`` of each library to a kit name and its target intervals, and the step builds one model per kit.
  Libraries without a matching kit are left out.
- ``wgs`` bins the whole genome.
  All libraries must be whole genomes of the same library kit.
  The kit is named ``wgs``.

Outputs
-------

The model directories stay in ``work/`` and are not published:

- ``{library_kit}/out/{library_kit}.contig_ploidy/`` holds the ploidy model (``ploidy-model``) and the ploidy calls per sample.
- ``{library_kit}/out/{library_kit}.{shard}.call_cnvs/`` holds, for each shard of the genome, the CNV calls model (``cnv_calls-model``), the calls per sample and the tracking data.
- ``{library_kit}/out/{library_kit}.scatter_intervals/`` holds the interval lists of the shards, one directory ``temp_{shard}`` each.
  The number of shards is known after the checkpoint ``scatter_intervals``.

The keys for downstream tasks are ``interval_list`` (the binned intervals of the kit), ``ploidy`` and ``calls`` (the ``.done`` files in the model directories) and ``shards`` (the directory of the shard interval lists).

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

from .model import GcnvModel as GcnvModelConfigModel
from .model import Intervals


def _filter_intervals_memory(wildcards, input=None, threads=None, attempt=1):
    _, _, _ = wildcards, input, threads  # unused but cannot be renamed
    return f"{attempt * 4 * 1024 + 16 * 1024}MB"


class GcnvModelStepPart(BuildGcnvModelStepPart):
    """Build model for GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        self.ngs_library_to_kit = build_ngs_library_to_kit(
            self.parent.primary_dna_libraries(include_background=True), self.config
        )

    @property
    def resource_usage(self):
        """Whole genomes need more memory to filter the intervals."""
        usage = dict(BuildGcnvModelStepPart.resource_usage)
        if self.config.intervals == Intervals.wgs:
            usage["filter_intervals"] = dataclasses.replace(
                BuildGcnvModelStepPart.default_resource_usage, mem=_filter_intervals_memory
            )
        return usage

    @dictify
    def _get_input_files_post_germline_calls(self, wildcards, checkpoints):
        checkpoint = checkpoints.gcnv_model_scatter_intervals
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


class GcnvModelWorkflow(BaseStep):
    """Perform gCNV model building by library kit"""

    #: Workflow name
    name = "gcnv_model"
    config_model_class = GcnvModelConfigModel
    produces = [DataSignature(DataType.MODELS, frozenset({"gcnv"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def task_produces(cls, config, upstream):
        """The model of the task's intervals, ``targeted`` or ``wgs``."""
        return (DataSignature(DataType.MODELS, frozenset({"gcnv", str(config.intervals)})),)

    @classmethod
    def get_output_paths(
        cls, config, signature=None, library_kit="{library_kit}", shard="{shard}"
    ) -> dict[str, str]:
        """Return the local model paths for the case mode of the callers."""
        prefix = f"work/{library_kit}/out/{library_kit}"
        return {
            "interval_list": f"{prefix}.interval_list",
            "ploidy": f"{prefix}.contig_ploidy/.done",
            "calls": f"{prefix}.{shard}.call_cnvs/.done",
            "shards": f"{prefix}.scatter_intervals",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((WritePedigreeStepPart, GcnvModelStepPart))
