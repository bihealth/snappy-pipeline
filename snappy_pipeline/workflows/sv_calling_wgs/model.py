import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyStepModel, SnappyModel, EnumField
from snappy_pipeline.models.gcnv import PrecomputedModelEntry


class DnaTool(enum.StrEnum):
    delly2 = "delly2"
    manta = "manta"
    popdel = "popdel"
    gcnv = "gcnv"
    melt = "melt"


class DnaLongTool(enum.StrEnum):
    sniffles2 = "sniffles2"
    # These seem to be unused:
    # sniffles = "sniffles"
    # pb_honey_spots = "pb_honey_spots"


class Tools(SnappyModel):
    dna: Annotated[list[DnaTool], EnumField(DnaTool, [DnaTool.delly2])]
    dna_long: Annotated[list[DnaLongTool], EnumField(DnaLongTool, [])]


class Gcnv(SnappyModel):
    path_par_intervals: str
    """Path to interval block list with PAR region for contig calling."""

    precomputed_model_paths: list[PrecomputedModelEntry] = []
    """
    Path to gCNV model - will execute analysis in CASE MODE.
    Example:
      - library: "Agilent SureSelect Human All Exon V6"  # Kit name, match in path_target_interval_list_mapping
        contig_ploidy: /path/to/ploidy-model         # Output from `DetermineGermlineContigPloidy`
        model_pattern: /path/to/model_*              # Output from `GermlineCNVCaller`
    """

    path_uniquely_mapable_bed: str  # "mapable" is a typo in the original code
    """Path to BED file with uniquely mappable regions."""

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Delly2(SnappyModel):
    path_exclude_tsv: str | None = None

    map_qual: int = 1

    geno_qual: int = 5

    qual_tra: int = 20

    mad_cutoff: int = 9

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Manta(SnappyModel):
    num_threads: int = 16

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Melt(SnappyModel):
    me_refs_infix: str = "1KGP_Hg19"
    me_types: list[str] = ["ALU", "LINE1", "SVA"]
    jar_file: str
    genes_file: str = "add_bed_files/1KGP_Hg19 / hg19.genes.bed"
    """adjust, e.g., Hg38/Hg38.genes.bed"""

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Popdel(SnappyModel):
    window_size: int = 10000000

    max_sv_size: int = 20000
    """== padding"""

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Sniffles2(SnappyModel):
    tandem_repeats: Annotated[
        str,
        Field(
            examples=[
                "/fast/groups/cubi/work/projects/biotools/sniffles2/trf/GRCh37/human_hs37d5.trf.bed"
            ]
        ),
    ]

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class SvCallingWgs(SnappyStepModel):
    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]

    tools: Tools

    delly2: Delly2 | None = None

    manta: Manta | None = None

    popdel: Popdel | None = None

    gcnv: Gcnv | None = None

    melt: Melt | None = None

    sniffles2: Sniffles2 | None = None

    ignore_chroms: list[str] = ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "chrUn_*"]

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for data_type in ("dna", "dna_long"):
            tool_list = getattr(self.tools, data_type)
            for tool in tool_list:
                if not getattr(self, tool):
                    raise ValueError(f"Tool {tool} not configured")
        return self
