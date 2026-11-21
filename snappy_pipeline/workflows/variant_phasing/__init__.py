"""Variant phasing with GATK or WhatsHap

Phases germline variants of pedigrees, or the combined somatic and germline variants of tumor libraries.

The ``tool`` field selects the phasing:

- ``gatk_phase_by_transmission`` runs PhaseByTransmission on the annotated variants of ``depends_on.variants``. It changes the genotypes of some variants, which is a problem when phasing *de novo* variants.
- ``gatk_read_backed_phasing`` runs ReadBackedPhasing on ``depends_on.variants`` with the alignments of the trio (``depends_on.alignments``), scattered over ``gatk_read_backed_phasing.num_jobs`` chunks of the genome. It is not fully reliable.
- ``gatk_read_backed_phasing_tumor`` runs ReadBackedPhasing on the combined somatic and germline variants of each tumor DNA library (``depends_on.combined_variants``, from ``combine_variants``) with its alignments, as pVACtools describes it for the proximal variants of neoepitope prediction.
- ``whatshap_tumor`` phases the same input with WhatsHap, a maintained read-based phasing tool, and writes the phase into the ``HP`` tag of ReadBackedPhasing (``--tag=HP``), which pVACseq reads; the genotypes stay unphased.
  WhatsHap assumes diploid data, so loss of heterozygosity, rearrangements and subclonal variants of a tumor limit it, as they limit ReadBackedPhasing.

The pedigree tools produce one result per child in a pedigree whose parents are both sequenced, the tumor tool one per tumor DNA library.

One task runs one tool.
To use both, add a ``gatk_read_backed_phasing`` task whose ``depends_on.variants`` names the ``gatk_phase_by_transmission`` task.
Only results of the GATK variant callers can be phased.
It is unclear how useful the phasing is for clinical studies.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

Here ``library_name`` is the library of the child, or the tumor library.
The output has the tags of the input plus ``phased``.
- with ``calls_report``: ``{library_name}/report/calls/{library_name}.calls.tsv``, the calls of the VCF as a table (one row per call), and the same table as a datavzrd site below ``report/calls/html/``, which ``snappy report`` shows.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

from functools import cached_property
from typing import Any

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.common.calls_report import CallsReportStepPart
from snappy_pipeline.workflows.common.scatter import RegionScatter
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import Tool
from .model import VariantPhasing as VariantPhasingConfigModel

#: The phased VCF of an entity, without ``.vcf.gz``
_BASE_PATH_OUT = "work/{library_name}/out/{library_name}"

#: Output keys of a phased VCF
_OUTPUTS = tuple(vcf_outputs(_BASE_PATH_OUT))


class WriteTrioPedigreeStepPart(WritePedigreeStepPart):
    """Write the trio PED of a child: the child and its sequenced parents, in the child's family (its ``related_group``)"""

    #: Step name
    name = "write_trio_pedigree"

    #: The child names the PED
    entity_wildcard = "library_name"

    def _get_input_files_run(self, wildcards):
        return []

    def _members(self, entity: str):
        """The rows of the child and of its sequenced parents"""
        libraries = self.parent.build_library_dataframe(include_background=True)
        child = self.parent.library_row(entity)
        names = [
            name for name in (entity, child["father_name"], child["mother_name"]) if name != "0"
        ]
        return libraries.set_index("library_name", drop=False).loc[names]

    def _family(self, entity: str) -> str:
        return self.parent.library_row(entity)["related_group"]


class VariantPhasingBaseStep(BaseStepPart):
    """Base class of the phasing tools"""

    #: Class available actions
    actions = ("run",)

    def _yield_vcf(self, wildcards):
        """Yield the VCF to phase and its index

        Annotated variants are stored under the library of the pedigree's index, phased variants
        under the library of the child.
        """
        library_name = wildcards.library_name
        if not self.parent.variants_are_phased:
            library_name = self.parent.pedigree_of(library_name, include_background=True)
        variants = self.parent.get_upstream_paths("variants", library_name=library_name)
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}


class PhaseByTransmissionStepPart(VariantPhasingBaseStep):
    """Phasing by transmission."""

    #: Name of the step in the pipeline.
    name = "gatk_phase_by_transmission"
    wrappers = {"run": "gatk_phase_by_transmission"}

    result_keys = _OUTPUTS

    @dictify
    def _get_input_files_run(self, wildcards):
        # Pedigree file required for PhaseByTransmission.
        yield (
            "ped",
            f"work/write_pedigree.{wildcards.library_name}/out/{wildcards.library_name}.ped",
        )
        yield from self._yield_vcf(wildcards)
        yield from self.parent.reference_inputs().items()

    def _get_output_files_run(self):
        return vcf_outputs(_BASE_PATH_OUT)

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        config = self.config.gatk_phase_by_transmission
        return {
            "de_novo_prior": config.de_novo_prior,
            "pedigree_validation_type": config.pedigree_validation_type,
            "father_allele_first": config.father_allele_first,
        }

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem=f"{14 * 1024}MB")


class ReadBackedPhasingStepPart(RegionScatter, VariantPhasingBaseStep):
    """Read-backed phasing, scattered over chunks of the genome and gathered per child."""

    #: Name of the step in the pipeline.
    name = "gatk_read_backed_phasing"

    #: Class available actions
    actions = ("scatter", "run", "gather")
    wrappers = RegionScatter.scatter_wrapper | {
        "gather": "bcftools/concat",
        "run": "gatk_read_backed_phasing",
    }

    extra_arguments_action = "run"

    result_keys = {"gather": _OUTPUTS}

    temp_outputs = {"scatter": ("regions",), "run": ("vcf", "vcf_tbi")}

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="4h", mem="14336MB"),
        "gather": ResourceUsage(threads=1, runtime="4h", mem="2048MB"),
    }

    #: The phased variants of a child in one part of the genome
    chunk_tpl = "work/{library_name}/chunks/{scatteritem}"

    def scatter_ignore_contigs(self) -> list[str]:
        return self.config.ignore_contigs

    def _yield_bams(self, wildcards):
        """Yield the alignments of the trio"""
        libraries = self.parent.build_library_dataframe(include_background=True)
        child = libraries[libraries["library_name"] == wildcards.library_name].iloc[0]
        if "0" in (child["father_name"], child["mother_name"]):
            return
        trio_libs = [child["library_name"], child["father_name"], child["mother_name"]]
        alignments = [
            self.parent.get_upstream_paths("alignments", library_name=lib) for lib in trio_libs
        ]
        yield "bam", [aln.bam for aln in alignments]
        yield "bai", [aln.bai for aln in alignments]

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from self._yield_bams(wildcards)
        yield from self._yield_vcf(wildcards)
        yield from self.parent.reference_inputs().items()
        yield "region", self.region(wildcards)

    def _get_input_files_gather(self, wildcards):
        return {
            "vcf": self.gathered(self.chunk_tpl + ".vcf.gz", wildcards),
            "vcf_tbi": self.gathered(self.chunk_tpl + ".vcf.gz.tbi", wildcards),
        }

    def _get_params_run(self, wildcards):
        config = self.config.gatk_read_backed_phasing
        return {
            "phase_quality_threshold": config.phase_quality_threshold,
            "window_size": config.window_size,
            "max_variants_per_window": config.max_variants_per_window,
        }

    def _get_output_files_run(self):
        return {
            "vcf": self.chunk_tpl + ".vcf.gz",
            "vcf_tbi": self.chunk_tpl + ".vcf.gz.tbi",
        }

    def _get_output_files_gather(self):
        return vcf_outputs(_BASE_PATH_OUT)

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.chunks/{scatteritem}.log"}

    def _get_log_file_gather(self):
        return super()._get_log_file_run()


class TumorPhasingStepPart(VariantPhasingBaseStep):
    """Phasing of a tumor library's combined somatic and germline variants with its alignments."""

    result_keys = _OUTPUTS

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem="32GB")

    def _get_params_run(self, wildcards):
        return {}

    def _get_input_files_run(self, wildcards):
        library_name = wildcards.library_name
        variants = self.parent.get_upstream_paths("combined_variants", library_name=library_name)
        alignments = self.parent.get_upstream_paths("alignments", library_name=library_name)
        return {
            "vcf": variants.vcf,
            "vcf_tbi": variants.vcf_tbi,
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
        }

    def _get_output_files_run(self):
        return vcf_outputs(_BASE_PATH_OUT)


class TumorReadBackedPhasingStepPart(TumorPhasingStepPart):
    """GATK 3 ReadBackedPhasing of a tumor library"""

    name = "gatk_read_backed_phasing_tumor"
    wrappers = {"run": "gatk_read_backed_phasing_tumor"}


class TumorWhatshapStepPart(TumorPhasingStepPart):
    """WhatsHap phasing of a tumor library"""

    name = "whatshap_tumor"
    wrappers = {"run": "whatshap/phase"}


#: The tools that phase tumor libraries instead of pedigrees
_TUMOR_TOOLS = (Tool.gatk_read_backed_phasing_tumor, Tool.whatshap_tumor)


class VariantPhasingWorkflow(BaseStep):
    """Perform (small) variant phasing"""

    name = "variant_phasing"
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "phased"})),
        DataSignature(DataType.VARIANTS, frozenset({"combined", "phased"})),
    ]
    config_model_class = VariantPhasingConfigModel

    tool_parts = {
        Tool.gatk_phase_by_transmission: PhaseByTransmissionStepPart,
        Tool.gatk_read_backed_phasing: ReadBackedPhasingStepPart,
        Tool.gatk_read_backed_phasing_tumor: TumorReadBackedPhasingStepPart,
        Tool.whatshap_tumor: TumorWhatshapStepPart,
    }
    helper_parts = (WriteTrioPedigreeStepPart,)

    @classmethod
    def task_produces(cls, config, upstream):
        """The variants this task reads, tagged ``phased``."""
        field = "combined_variants" if config.tool in _TUMOR_TOOLS else "variants"
        return (upstream[field].with_tags("phased"),)

    @classmethod
    def task_library_requirement(cls, config):
        """The tumor tool phases tumor DNA libraries."""
        return (
            "role == 'tumor' and extraction_type == 'dna'" if config.tool in _TUMOR_TOOLS else None
        )

    @classmethod
    def default_group_by(cls, config):
        """The tumor tools phase each tumor library; the pedigree tools declare the grouping of their variants task."""
        return "library" if config.tool in _TUMOR_TOOLS else None

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local phased-variant output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        paths = vcf_outputs(_BASE_PATH_OUT.format(library_name=lib))
        return {key: paths[key] for key in ("vcf", "vcf_tbi")}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        if self.config.calls_report:
            self.register_sub_step_classes((CallsReportStepPart,))

    @cached_property
    def variants_are_phased(self) -> bool:
        """Whether ``depends_on.variants`` is a phasing task, whose results are per child."""
        return "phased" in self.upstream_signature("variants").tags

    @cached_property
    def output_entities(self) -> list[str]:
        """The children whose parents are both sequenced (their DNA libraries), or the tumor
        libraries."""
        if self.config.tool in _TUMOR_TOOLS:
            return super().output_entities
        return self.children_with_both_parents()
