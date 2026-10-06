import enum

from snappy_pipeline.models import SnappyStepModel


class CallingTools(enum.StrEnum):
    MUTECT2 = "mutect2"
    GATKHC = "gatk_hc"


class VariantOrigin(enum.StrEnum):
    SOMATIC = "somatic"
    GERMLINE = "germline"
    ANY = "any"


class AnyVariantCalling(SnappyStepModel):
    variant_origin: VariantOrigin = VariantOrigin.SOMATIC

    tools: list[CallingTools] = []

    path_ngs_mapping: str = "../ngs_mapping"
    """Path to ngs_mapping step"""

    tools_ngs_mapping: list[str] = []
    """Default: use those defined in ngs_mapping step"""
