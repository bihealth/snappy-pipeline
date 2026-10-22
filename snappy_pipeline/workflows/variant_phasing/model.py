import enum
from typing import Annotated, Literal

from pydantic import Field, model_validator

from snappy_pipeline.models import (
    EnumField,
    SnappyModel,
    SnappyStepModel,
    ToolModel,
    validators,
)
from snappy_pipeline.workflows.abstract.protocol import (
    COHORT_LEVEL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation.model import ExpectedAnnotatedVariants


class Tool(enum.StrEnum):
    gatk_phase_by_transmission = "gatk_phase_by_transmission"
    gatk_read_backed_phasing = "gatk_read_backed_phasing"
    gatk_read_backed_phasing_tumor = "gatk_read_backed_phasing_tumor"
    whatshap_tumor = "whatshap_tumor"


#: The ``depends_on`` fields each tool reads
TOOL_DEPENDENCIES = {
    Tool.gatk_phase_by_transmission: ("variants",),
    Tool.gatk_read_backed_phasing: ("variants",),
    Tool.gatk_read_backed_phasing_tumor: ("combined_variants",),
    Tool.whatshap_tumor: ("combined_variants",),
}


class GatkReadBackedPhasing(ToolModel):
    phase_quality_threshold: float = 20.0
    """quality threshold for phasing"""

    num_jobs: int = 24
    """number of chunks the genome is split into, each phased in its own job"""

    window_size: Annotated[int, Field(ge=1)] = 10000
    """size of the windows (bases) in which the variants are counted"""

    max_variants_per_window: int = 100
    """windows with more variants are not phased"""


class GatkPhaseByTransmission(ToolModel):
    de_novo_prior: float = 1e-8
    """use 1e-6 when interested in phasing de novos"""

    pedigree_validation_type: Literal["ALL", "STRICT", "SILENT"] = "SILENT"
    """how strictly GATK validates the pedigree against the VCF (``--pedigreeValidationType``)"""

    father_allele_first: bool = True
    """write the paternal allele first in phased genotypes (``--FatherAlleleFirst``)"""


class ExpectedPhasedVariants(SnappyModel):
    """Consumer-driven contract: expected output keys from variant_phasing."""

    vcf: str
    vcf_tbi: str


class ExpectedCombinedVariants(SnappyModel):
    """Consumer-driven contract: the combined somatic and germline variants of a tumor library."""

    vcf: str
    vcf_tbi: str


class VariantPhasingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ]
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", ("annotated", "phased")})),
        COHORT_LEVEL,
        ExpectedPathSchema(ExpectedAnnotatedVariants),
    ] = ""
    """Annotated germline variants, or the phased variants of a ``gatk_phase_by_transmission`` task
    (pedigree tools)"""

    combined_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"combined"})),
        ExpectedPathSchema(ExpectedCombinedVariants),
    ] = ""
    """The combined somatic and germline variants of the tumor libraries
    (``gatk_read_backed_phasing_tumor``, ``whatshap_tumor``)"""

    reference: Reference


class VariantPhasing(SnappyStepModel):
    depends_on: VariantPhasingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    ignore_contigs: list[str] = ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"]
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""

    gatk_read_backed_phasing: GatkReadBackedPhasing = GatkReadBackedPhasing()

    gatk_phase_by_transmission: GatkPhaseByTransmission = GatkPhaseByTransmission()

    gatk_read_backed_phasing_tumor: ToolModel = ToolModel()

    whatshap_tumor: ToolModel = ToolModel()

    @model_validator(mode="after")
    def ensure_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
