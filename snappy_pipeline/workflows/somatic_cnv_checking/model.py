import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class ExpectedSomaticCnvCalls(SnappyModel):
    """Consumer-driven contract: expected output keys from somatic CNV caller steps."""

    vcf: str


class CnvAssayType(enum.StrEnum):
    WES = "WES"
    WGS = "WGS"


class SomaticCnvCheckingDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"
    cnv_calling: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedSomaticCnvCalls),
    ] = "cnv_calling"


class SomaticCnvChecking(SnappyStepModel):
    depends_on: SomaticCnvCheckingDependsOn = Field(default_factory=SomaticCnvCheckingDependsOn)

    cnv_assay_type: CnvAssayType | None = None
    """
    Empty: no CNV,
    WES for somatic_targeted_seq_snv_calling step,
    WGS for somatic_wgs_cnv_calling step
    """

    excluded_regions: str = ""
    """Bed file of regions to be excluded"""

    max_depth: int = 10000
    """Max depth for pileups"""

    min_depth: int = 20
    """Minimum depth for reference and alternative alleles to consider variant"""

    min_baf: Annotated[float, Field(0.4, ge=0, le=0.5)]
    """Maximum BAF to consider variant as heterozygous (between 0 & 1/2)"""
