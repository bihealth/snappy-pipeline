import enum

from typing import Annotated
from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel


class LibraryKitEntry(SnappyModel):
    """
    Mapping from enrichment kit to target region BED file, for either computing per--target
    region coverage or selecting targeted exons.

    The following will match both the stock IDT library kit and the ones
    with spike-ins seen fromr Yale genomics.  The path above would be
    mapped to the name "default".
      - name: IDT_xGen_V1_0
        pattern: "xGen Exome Research Panel V1\\.0*"
        path: "path/to/targets.bed"
    """

    name: Annotated[str, Field(examples=["IDT_xGen_V1_0"])]

    pattern: Annotated[str, Field(examples=["xGen Exome Research Panel V1\\.0*"])]

    path: Annotated[str, Field(examples=["path/to/targets.bed"])]


class LibraryKit(SnappyModel):
    path_target_interval_list_mapping: list[LibraryKitEntry] = []
    """Connects sample-based library kit in sample sheets with corresponding bed files"""


class SexValue(enum.StrEnum):
    MALE = "male"
    FEMALE = "female"


class SexOrigin(enum.StrEnum):
    AUTOMATIC = "auto"
    SAMPLESHEET = "samplesheet"
    CONFIG = "config"


class Sex(SnappyModel):
    source: SexOrigin = SexOrigin.AUTOMATIC
    default: SexValue | None = None

    @model_validator(mode="after")
    def ensure_default_value(self):
        if self.source == SexOrigin.CONFIG and not self.default:
            raise ValueError("Undefined default sex value in configuration file")
        return self
