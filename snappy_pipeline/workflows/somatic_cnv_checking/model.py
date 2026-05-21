import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class CnvAssayType(enum.StrEnum):
    WES = "WES"
    WGS = "WGS"


class SomaticCnvCheckingDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"
    cnv_calling: str = "cnv_calling"


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
