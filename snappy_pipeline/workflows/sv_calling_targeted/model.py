import enum
from typing import Annotated, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.gcnv import PrecomputedModelEntry, TargetIntervalEntry


class Tool(enum.StrEnum):
    gcnv = "gcnv"
    delly2 = "delly2"
    manta = "manta"
    melt = "melt"


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


class SvCallingTargeted(SnappyStepModel):
    path_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]

    tools: Annotated[
        list[Tool], EnumField(Tool, [Tool.gcnv, Tool.delly2, Tool.manta], min_length=1)
    ]

    path_target_interval_list_mapping: list[TargetIntervalEntry]
    """
    The following allows to define one or more set of target intervals.  This is only used by gcnv.
    Example:
     - name: "Agilent SureSelect Human All Exon V6"
       pattern: "Agilent SureSelect Human All Exon V6.*"
       path: "path/to/targets.bed"
    """

    gcnv: Gcnv | None = None

    dell2: Delly2 | None = None

    manta: Manta | None = None

    melt: Melt | None = None

    @model_validator(mode="after")
    def ensure_tools_are_configured(self: Self) -> Self:
        for tool in self.tools:
            if not getattr(self, str(tool)):
                raise ValueError(f"Tool {tool} not configured")
        return self
