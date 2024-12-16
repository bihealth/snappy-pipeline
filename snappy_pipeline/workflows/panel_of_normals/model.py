import enum
from typing import Annotated, Literal

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, Parallel, validators
from snappy_pipeline.models.common import LibraryKitEntry, Sex
from snappy_pipeline.models.cnvkit import CnvkitToReference as CnvkitGeneric
from snappy_pipeline.models.mutect2 import Mutect2 as Mutect2Generic


class Tool(enum.StrEnum):
    mutect2 = "mutect2"
    cnvkit = "cnvkit"
    purecn = "purecn"
    access = "access"


class Mutect2(Parallel, Mutect2Generic):
    path_normals_list: str = ""
    """Optional file listing libraries to include in panel"""

    java_options: str = " -Xmx16g"
    """Optional java run-time options"""


class CnvKit(CnvkitGeneric):
    path_normals_list: str | None = None
    """Optional file listing libraries to include in panel"""

    path_target_interval_list_mapping: list[LibraryKitEntry] = []
    """
    Bed files of enrichment kit sequences (MergedProbes for Agilent SureSelect),
    recommended by PureCN author
    """

    sample_sex: Sex = Sex()
    """Sets the sex of all normals used in the panel"""

    path_access: str | None = None
    """Overrides access when not None"""

    ignore_chroms: Annotated[
        list[str],
        Field(
            examples=[
                "chrM",
                "MT",
                "*_random",
                "chrUn_*",
                "GL*",
            ]
        ),
    ] = []
    """Additional contigs to ignore, specific to the tool"""


class GenomeName(enum.StrEnum):
    hg18 = "hg18"
    hg19 = "hg19"
    hg38 = "hg38"
    mm9 = "mm9"
    mm10 = "mm10"
    rn4 = "rn4"
    rn5 = "rn5"
    rn6 = "rn6"
    canFam3 = "canFam3"


class PureCn(SnappyModel):
    path_normals_list: str = ""
    """Optional file listing libraries to include in panel"""

    # targets_definition: list[LibraryKitEntry] = []
    """
    Bed files of enrichment kit sequences (MergedProbes for Agilent SureSelect),
    recommended by PureCN author
    """

    path_genomicsDB: str
    """Mutect2 genomicsDB created during panel_of_normals"""

    genome_name: Annotated[
        GenomeName | Literal["unknown"],
        EnumField(GenomeName, json_schema_extra={"options": {"unknown"}}),
    ] = "unknown"

    path_target_interval_list_mapping: list[LibraryKitEntry] = []
    """For filename only..."""

    mappability: str = ""
    """
    GRCh38:
    /fast/work/groups/cubi/projects/biotools/static_data/app_support/PureCN/hg38/mappability.bw
    """

    reptiming: str = ""
    """Nothing for GRCh38"""

    seed: int = 1234567


class PanelOfNormals(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.mutect2], min_length=1)]

    path_ngs_mapping: str = "../ngs_mapping"

    ignore_chroms: Annotated[
        list[str],
        Field(
            examples=[
                "NC_007605",
                "hs37d5",
                "chrEBV",
                "HLA-*",
                "GL000220.*",
                "chrEBV",
                "HPV*",
                "CMV",
                "HBV",
                "HCV-*",
                "HIV-*",
                "KSHV",
                "HTLV-1",
                "MCV",
                "*_decoy",
                "chrUn_GL00220*",
                "SV40",
            ]
        ),
    ] = []
    """Patterns of contig names to ignore"""

    mutect2: Mutect2 | None = None

    cnvkit: CnvKit | None = None

    purecn: PureCn | None = None
