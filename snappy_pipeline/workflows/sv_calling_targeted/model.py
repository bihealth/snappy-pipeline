import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.gcnv import PrecomputedModelEntry, TargetIntervalEntry


class Tool(enum.StrEnum):
    gcnv = "gcnv"
    delly2 = "delly2"
    manta = "manta"
    melt = "melt"

# Note: GCNV config for WGS & targeted has differences
class Gcnv(SnappyModel):
    # TODO: these two will only be used if the wrappers accessing this key are actually run by the workflow
    # path_par_intervals: str = ""
    # """Path to interval block list with PAR region for contig calling."""

    # path_uniquely_mapable_bed: str
    # """path to BED file with uniquely mappable regions."""

    path_target_interval_list_mapping: list[TargetIntervalEntry]
    """
    The following allows to define one or more set of target intervals.  This is only used by gcnv.
    Example:
     - name: "Agilent SureSelect Human All Exon V6"
       pattern: "Agilent SureSelect Human All Exon V6.*"
       path: "path/to/targets.bed"
    """

    # Not sure if this can/should be empty be default; can we set a required flag?
    # TODO: this can probably not be empty
    # TODO: the 'library' fields need to also exist in path_target_interval_list_mapping (as name)
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

# TODO: only a singular Delly2 config for germline WGS & targeted is needed
class Delly2(SnappyModel):
    path_exclude_tsv: str | None = Field(None, examples=[
        "/data/cephfs-1/work/projects/cubit/20.05/static_data/app_support/Delly/0.8.3/human.hg19.excl.tsv",
        "/data/cephfs-1/work/projects/cubit/20.05/static_data/app_support/Delly/0.8.3/human.hg38.excl.tsv"
    ])

    map_qual: int = 1

    geno_qual: int = 5

    qual_tra: int = 20

    mad_cutoff: int = 9

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """

# TODO: only a singular manta config for WGS & targeted is needed
class Manta(SnappyModel):
    num_threads: int = 16

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """

# TODO: only a singular melt config for WGS & targeted is needed
class Melt(SnappyModel):
    me_refs_infix: str = Field(examples=["1KGP_Hg19"])

    me_types: list[str] = ["ALU", "LINE1", "SVA"]

    jar_file: str = Field(examples=["MELT.jar"])

    genes_file: str = Field(examples=["add_bed_files/1KGP_Hg19/hg19.genes.bed"])
    """adjust, e.g., Hg38/Hg38.genes.bed"""

    me_refs_path: str = Field(examples=["me_refs"])
    """Path to MELT reference files"""

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class SvCallingTargeted(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"

    tools: Annotated[
        list[Tool], EnumField(Tool, [Tool.gcnv, Tool.delly2, Tool.manta], min_length=1)
    ]

    gcnv: Gcnv | None = None

    delly2: Delly2 | None = None

    manta: Manta | None = None

    melt: Melt | None = None
