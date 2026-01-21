import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyStepModel, ToggleModel, validators
from snappy_pipeline.models.gatk import GATK
from snappy_pipeline.models.parallel import Parallel


class Tool(enum.StrEnum):
    mutect2 = "mutect2"


# Adjustment tumor_only mode
class TumorNormalMode(enum.StrEnum):
    AUTOMATIC = "automatic"
    PAIRED = "paired"
    TUMOR_ONLY = "tumor_only"
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class Contamination(ToggleModel, GATK):
    common_variants: str = ""
    """Common germline variants for contamination estimation"""

    pileup: GATK = GATK()
    """Parameters for GetPileupSummaries used on the tumor bam (& normalk if present)"""

    @model_validator(mode="after")
    def ensure_common_variant_when_enabled(self):
        if self.enabled and not self.common_variants:
            raise ValueError("Common variants must be present when contamination is enabled")
        return self


class Mutect2(Parallel, GATK):
    # Sadly a type of
    # `FilePath | None = None`
    # still applies `FilePath` validation on `None`, which errors
    panel_of_normals: str | None = ""
    """Set path to panel of normals vcf if required"""

    germline_resource: str | None = ""
    """Germline variants resource (same as panel of normals)"""

    contamination: Contamination
    """Estimation of contamination using GetPileupSummaries & CalculateContamination"""

    filtration: GATK = GATK()
    """Additional arguments for filtration"""

    padding: int = 5000
    """Padding around intervals for scatter/gather"""

    tumor_normal_mode: TumorNormalMode = TumorNormalMode.AUTOMATIC
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class SomaticVariantCalling(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [], min_length=1)]
    """List of tools"""

    path_ngs_mapping: str = "../ngs_mapping"
    """Path to ngs_mapping"""

    ignore_chroms: Annotated[
        list[str],
        Field(examples=["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"]),
    ] = []
    """Patterns of contig names to ignore"""

    mutect2: Mutect2 | None = None
    """Configuration for MuTect 2"""
