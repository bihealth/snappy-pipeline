import enum
from typing import Annotated, Literal

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.cnvkit import PanelOfNormals as CnvKit
from snappy_pipeline.models.gatk import GATK
from snappy_pipeline.models.parallel import Parallel


class Tool(enum.StrEnum):
    mutect2 = "mutect2"
    cnvkit = "cnvkit"
    purecn = "purecn"
    access = "access"


class GenomicsDb(GATK):
    java_options: str = "-Xms16g -Xmx32g"
    """
    GenomicsDBImport needs lots of memory.
    (mind the caveats in https://gatk.broadinstitute.org/hc/en-us/articles/35967560824859-GenomicsDBImport)
    """


class Mutect2(Parallel, GATK):
    path_normals_list: str = ""
    """Optional file listing libraries to include in panel"""

    germline_resource: str
    """Germline variants resource (same as panel of normals)"""

    padding: int = 5000
    """Padding around intervals for scatter/gather"""

    genomicsdb: GenomicsDb = GenomicsDb()
    """
    Parameters for the creation of the (transient?) genomics database.
    The 2nd part of panel of normals creation (CreateSomaticPanelOfNormals) is not under user control
    """


class Access(SnappyModel):
    """Creates access file for cnvkit, based on genomic sequence & excluded regions (optionally)"""

    exclude: list[str] = []
    """[access] Bed file of regions to exclude (mappability, blacklisted, ...)"""

    min_gap_size: int = 0
    """[access] Minimum gap size between accessible sequence regions (0: use default value)"""


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

    path_bait_regions: str
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

    enrichment_kit_name: str = "unknown"
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
                ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"],
                [
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
                ],
            ]
        ),
    ] = []
    """Patterns of contig names to ignore"""

    mutect2: Mutect2 | None = None

    cnvkit: CnvKit | None = None

    access: Access = Access()

    purecn: PureCn | None = None
