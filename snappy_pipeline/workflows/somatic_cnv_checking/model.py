import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyStepModel


class CnvAssayType(enum.Enum):
    WES = "WES"
    WGS = "WGS"


class SomaticCnvChecking(SnappyStepModel):
    path_ngs_mapping: str

    path_cnv_calling: Annotated[str, Field(examples=["../somatic_targeted_seq_cnv_calling"])] = ""

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

    @model_validator(mode="after")
    def ensure_cnv_assay_type_is_specified(self: Self) -> Self:
        if self.path_cnv_calling and not self.cnv_assay_type:
            raise ValueError("CNV assay type must be specified")
        return self
