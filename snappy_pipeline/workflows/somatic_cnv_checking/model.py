from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel


class SomaticCnvChecking(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    path_cnv_calling: Annotated[str, Field(examples=["../somatic_cnv_calling"])] = ""

    excluded_regions: str = ""
    """Bed file of regions to be excluded"""

    max_depth: int = 10000
    """Max depth for pileups"""

    min_depth: int = 20
    """Minimum depth for reference and alternative alleles to consider variant"""

    min_baf: Annotated[float, Field(0.4, ge=0, le=0.5)]
    """Maximum BAF to consider variant as heterozygous (between 0 & 1/2)"""
