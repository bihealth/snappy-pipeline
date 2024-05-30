import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.cnvkit import Cnvkit


class Tool(enum.StrEnum):
    cnvetti = "cnvetti"
    control_freec = "control_freec"


class Canvas(SnappyModel):
    path_reference: str
    """Path to Canvas reference file"""

    path_filter_bed: str
    """Path to Canvas filter BED file"""

    path_genome_folder: str
    """Path to Canvas genome folder"""


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
            **{
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
    """Path to ControlFreec ChrLenFile"""

    path_mappability: str
    """Path to ControlFreec mappability file"""

    path_mappability_enabled: bool = False

    window_size: int = -1
    """set to a value >=0 you want a specific fixed window size"""

    convert: ControlFreecConvert


# If defaults need to be overwritten, subclass the model and override the defaults
class CnvkitWgs(Cnvkit):
    pass


class SomaticWgsCnvCalling(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"

    path_somatic_variant_calling: Annotated[str, Field(examples=["../somatic_variant_calling"])]

    somatic_variant_calling_tool: str

    tools: Annotated[list[Tool], EnumField(Tool, [Tool.cnvetti], min_length=1)]

    canvas: Canvas | None = None

    cnvetti: Cnvetti | None = None

    control_freec: ControlFreec | None = None

    cnvkit: CnvkitWgs | None = None
