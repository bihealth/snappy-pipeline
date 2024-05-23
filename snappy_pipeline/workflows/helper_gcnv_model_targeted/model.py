import os
from typing import Annotated

from pydantic import Field

from models import SnappyStepModel, SnappyModel


class TargetIntervalEntry(SnappyModel):
    """
    The following will match both the stock IDT library kit and the ones
    with spike-ins seen fromr Yale genomics.  The path above would be
    mapped to the name "default".
      - name: IDT_xGen_V1_0
        pattern: "xGen Exome Research Panel V1\\.0*"
        path: "path/to/targets.bed"
    """

    name: Annotated[str, Field(examples=["IDT_xGen_V1_0"])]

    pattern: Annotated[str, Field(examples=["xGen Exome Research Panel V1\\.0*"])]

    path: Annotated[os.PathLike, Field(examples=["path/to/targets.bed"])]


class Gcnv(SnappyModel):
    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""

    path_target_interval_list_mapping: list[TargetIntervalEntry] = []


class HelperGcnvModelTargeted(SnappyStepModel):
    path_ngs_mapping: str

    gcnv: Gcnv
