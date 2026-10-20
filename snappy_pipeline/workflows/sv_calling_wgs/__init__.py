"""Germline structural variant calling for WGS data

The step calls structural variants in whole genome samples from the aligned reads.
The ``tool`` field selects the caller:

- ``delly2`` calls the variants per library, merges the sites of a pedigree and genotypes them in every member.
  The regions in ``depends_on.excluded_regions`` are optional and excluded.
- ``manta`` calls the variants of all members of a pedigree together.
- ``popdel`` profiles each library and then calls the deletions of all libraries together.
- ``gcnv`` runs GATK gCNV in case mode, with the model of the task in ``depends_on.gcnv_model``.
  That task is a :ref:`step_gcnv_model` task with ``intervals: wgs``, and its intervals are the ones of the calls.
- ``melt`` calls mobile element insertions, with the genes of the task in ``depends_on.gene_regions``.

The task gives one result per pedigree, named after its index library.
``skip_libraries`` of a tool skips the whole pedigree of a library.
The step gives no results if the mapper of the alignments is not a DNA mapper.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
All tools publish the same files.

- ``{library_name}/out/{library_name}.vcf.gz``
- ``{library_name}/out/{library_name}.vcf.gz.tbi``

Here ``library_name`` is the index library of the pedigree.
Downstream tasks read the key ``vcf``.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart
from snappy_pipeline.workflows.common.sv_calling import (
    VCF_OUTPUTS,
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_wrappers.tools.genome_windows import yield_regions

from .model import SvCallingWgs as SvCallingWgsConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the sv_calling_wgs step


def escape_dots_dashes(s: str) -> str:
    """Escape dots and dashes with double-underscore constructs."""
    return s.replace("_", "__under__").replace(".", "__dot__").replace("-", "__hyphen__")


class PopDelStepPart(
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    BaseStepPart,
):
    """WGS SV identification using PopDel.

    PopDel requires whole-cohort wide calling.  We implement cohort-wise calling to speed
    things up a bit more.
    """

    name = "popdel"
    actions = ("profile", "call", "concat_calls", "reorder_vcf")

    temp_outputs = {"call": ("vcf", "vcf_tbi")}

    intermediate_outputs = {"profile": ("profile",), "concat_calls": ("vcf", "vcf_tbi")}

    default_resource_usage = ResourceUsage(threads=2, runtime="4d", mem="24GB")

    @dictify
    def _get_input_files_profile(self, wildcards):
        """Return input files for "call" action"""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield from self.parent.reference_inputs().items()

    @dictify
    def _get_output_files_profile(self):
        prefix = "work/{library_name}/out/{library_name}"
        yield "profile", f"{prefix}.profile"

    @dictify
    def _get_input_files_call(self, wildcards):
        dna = self.parent.primary_dna_libraries(include_background=True)
        paths = [
            f"work/{library_name}/out/{library_name}.profile"
            for library_name in dna["library_name"]
        ]
        yield "profile", paths

    @dictify
    def _get_output_files_call(self):
        infix = self._get_log_file_infix_call()
        yield "vcf", f"work/{infix}/out/{infix}.call.vcf.gz"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.call.vcf.gz.tbi"

    def _get_log_file_infix_call(self):
        return "{chrom}-{begin}-{end}"

    @dictify
    def _get_input_files_concat_calls(self, wildcards):
        window_size = self.config.popdel.window_size
        padding = self.config.popdel.max_sv_size
        vcfs = []
        with open(self._get_fai_path(), "rt") as fai_file:
            for r in yield_regions(
                fai_file,
                window_size=window_size,
                padding=padding,
                ignore_chroms=self._get_ignore_chroms(),
            ):
                if r.begin == 0:
                    r.begin = 1
                chrom = escape_dots_dashes(r.chrom)
                infix = f"{chrom}-{r.begin}-{r.end}"
                vcfs.append(f"work/{infix}/out/{infix}.call.vcf.gz")
        yield "vcf", vcfs
        yield "vcf_tbi", [f"{vcf}.tbi" for vcf in vcfs]

    def _get_fai_path(self):
        return self.parent.get_upstream_paths("reference").fai

    def _get_ignore_chroms(self):
        return self.config.ignore_chroms

    @dictify
    def _get_output_files_concat_calls(self):
        infix = self._get_log_file_infix_concat_calls()
        yield "vcf", f"work/{infix}/out/{infix}.concat_calls.vcf.gz"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.concat_calls.vcf.gz.tbi"

    def _get_log_file_infix_concat_calls(self):
        return "cohort"

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        yield "vcf", "work/cohort/out/cohort.concat_calls.vcf.gz"
        yield "vcf_tbi", "work/cohort/out/cohort.concat_calls.vcf.gz.tbi"

    def _get_output_files_reorder_vcf(self):
        return dict(VCF_OUTPUTS)

    def _get_params_call(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "chrom": wildcards.chrom,
            "begin": wildcards.begin,
            "end": wildcards.end,
        }

    def _get_params_reorder_vcf(self, wildcards: Wildcards) -> dict[str, Any]:
        """Used in Snakefile to rule ``sv_calling_wgs_popdel_reorder_vcf``"""
        return {
            "ped_members": " ".join(
                self.parent.pedigrees(include_background=True)[wildcards.library_name]
            )
        }


class Sniffles2StepPart(BaseStepPart):
    """WGS SV identification using Sniffles 2"""

    name = "sniffles2"
    actions = ("bam_to_snf", "snf_to_vcf")

    default_resource_usage = ResourceUsage(threads=2, runtime="2h", mem="4GB")

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/sniffles2.{index_ngs_library}/out/sniffles2.{index_ngs_library}{ext}"
        )

    @dictify
    def _get_input_files_bam_to_snf(self, wildcards):
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    @dictify
    def _get_output_files_bam_to_snf(self):
        infix = "sniffles2_bam_to_snf.{library_name}"
        yield "snf", f"work/{infix}/out/{infix}.snf"

    @dictify
    def _get_input_files_snf_to_vcf(self, wildcards):
        snfs = []
        for library_name in self.parent.pedigrees(include_background=True)[
            wildcards.index_ngs_library
        ]:
            infix = f"sniffles2_bam_to_snf.{library_name}"
            snfs.append(f"work/{infix}/out/{infix}.snf")
        yield "snf", snfs

    @dictify
    def _get_output_files_snf_to_vcf(self):
        infix = "sniffles2.{index_ngs_library}"
        yield "snf", f"work/{infix}/out/{infix}.snf"

    # FIXME: missing get_result_files


class SvCallingWgsWorkflow(BaseStep):
    """Perform (germline) WGS SV calling"""

    name = "sv_calling_wgs"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"}))]
    #: The pedigrees of background data sets are called, too
    uses_background = True

    config_model_class = SvCallingWgsConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local WGS SV output paths for downstream consumers."""
        library_name = kwargs.get("library_name", "{library_name}")
        return {"vcf": f"work/{library_name}/out/{library_name}.vcf.gz"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        selected_tool = self.config.tool
        match selected_tool:
            case Tool.delly2:
                selected_sub_step = Delly2StepPart
            case Tool.manta:
                selected_sub_step = MantaStepPart
            case Tool.popdel:
                selected_sub_step = PopDelStepPart
            case Tool.gcnv:
                selected_sub_step = RunGcnvStepPart
            case Tool.melt:
                selected_sub_step = MeltStepPart
            # case Tool.sniffles2:
            #     selected_sub_step = Sniffles2StepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {selected_tool}")
        # Register only the selected tool step class.
        self.register_sub_step_classes((selected_sub_step, WritePedigreeStepPart))
