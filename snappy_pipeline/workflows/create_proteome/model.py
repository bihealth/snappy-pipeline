import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel


class ToolNgsMapping(enum.StrEnum):
    BWA = "bwa"
    BWA_MEM2 = "bwa_mem2"


class ToolVariantCalling(enum.StrEnum):
    MUTECT2 = "mutect2"
    GATK4_HC = "gatk4_hc"


class ToolVariantAnnotation(enum.StrEnum):
    VEP = "vep"


class InputVariantType(enum.StrEnum):
    CALLING = "calling"
    ANNOTATION = "annotation"
    FILTRATION = "filtration"


class CreateProteome(SnappyStepModel):
    tool_ngs_mapping: ToolNgsMapping = ToolNgsMapping.BWA

    variant_type: InputVariantType = InputVariantType.CALLING
    path_variant: str = "../germlinec_variant_calling"
    tool_variant_calling: ToolVariantCalling = ToolVariantCalling.GATK4_HC
    tool_variant_annotation: ToolVariantAnnotation | None = None
    is_filtered: bool = False

    add_reference: bool = False
    """Should the process also compute proteome based on unmutated sequence"""
    path_proteome: Annotated[
        str | None,
        Field(examples=["gencode.v43.pc_translations.fa.gz", "Homo_sapiens.GRCh38.pep.all.fa.gz"]),
    ] = None
    """Path to curated unmutated proteome"""
