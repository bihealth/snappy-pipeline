import enum
from typing import Annotated

from snappy_pipeline.models import (
    EnumField,
    SnappyModel,
    SnappyStepModel,
)
from snappy_pipeline.workflows.abstract.protocol import (
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


class GatkReadBackedPhasing(SnappyModel):
    phase_quality_threshold: float = 20.0
    """quality threshold for phasing"""

    num_jobs: int = 24
    """number of chunks the genome is split into, each phased in its own job"""


class GatkPhaseByTransmission(SnappyModel):
    de_novo_prior: float = 1e-8
    """use 1e-6 when interested in phasing de novos"""


class ExpectedPhasedVariants(SnappyModel):
    """Consumer-driven contract: expected output keys from variant_phasing."""

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
        ExpectedPathSchema(ExpectedAnnotatedVariants),
    ]
    """Annotated germline variants, or the phased variants of a ``gatk_phase_by_transmission`` task"""

    reference: Reference


class VariantPhasing(SnappyStepModel):
    depends_on: VariantPhasingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    ignore_chroms: list[str] = ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"]
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""

    gatk_read_backed_phasing: GatkReadBackedPhasing = GatkReadBackedPhasing()

    gatk_phase_by_transmission: GatkPhaseByTransmission = GatkPhaseByTransmission()
