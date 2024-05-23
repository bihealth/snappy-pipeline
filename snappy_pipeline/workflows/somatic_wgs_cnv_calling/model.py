import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from models import SnappyStepModel, EnumField, SnappyModel
from models.cnvkit import Cnvkit


class Tool(enum.Enum):
    cnvetti = "cnvetti"
    control_freec = "control_freec"


class Canvas(SnappyModel):
    path_reference: str
    path_filter_bed: str
    path_genome_folder: str


class CnvettiPreset(SnappyModel):
    window_length: int
    count_kind: str
    segmentation: str
    normalization: str


class Cnvetti(SnappyModel):
    window_length: int | None = None
    count_kind: str | None = None
    segmentation: str | None = None
    normalization: str | None = None
    presets: dict[str, CnvettiPreset] = {
        "deep_wgs": CnvettiPreset(
            {
                "window_length": 200,
                "count_kind": "Coverage",
                "segmentation": "HaarSeg",
                "normalization": "MedianGcBinned",
            }
        )
    }
    preset: str = "deep_wgs"


class ControlFreecConvert(SnappyModel):
    org_obj: str = "org.Hs.eg.db::org.Hs.eg.db"

    tx_obj: str = "TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene"

    bs_obj: str = "BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5"


class ControlFreec(SnappyModel):
    path_chrlenfile: str

    path_mappability: str

    path_mappability_enabled: bool = False

    window_size: int = -1
    """set to a value >=0 you want a specific fixed window size"""

    convert: ControlFreecConvert


# If defaults need to be overwritten, subclass the model and override the defaults
class CnvkitWgs(Cnvkit):
    pass


class SomaticWgsCnvCalling(SnappyStepModel):
    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]

    path_somatic_variant_calling: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    somatic_variant_calling_tool: str

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.cnvetti])]

    canvas: Canvas | None = None

    cnvetti: Cnvetti | None = None

    control_freec: ControlFreec | None = None

    cnvkit: CnvkitWgs | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
