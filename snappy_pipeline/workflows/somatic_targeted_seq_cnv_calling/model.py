import enum
from typing import Annotated, Any, Literal

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    sequenza = "sequenza"
    copywriter = "copywriter"
    cnvetti_on_target = "cnvetti_on_target"
    cnvetti_off_target = "cnvetti_off_target"
    purecn = "purecn"


class SequenzaExtraArgs(SnappyModel):
    hom: float = 0.9
    """Threshold to select homozygous positions"""

    het: float = 0.25
    """Threshold to select heterozygous positions"""

    qlimit: float = 20
    """Minimum nucleotide quality score for inclusion in the counts"""

    qformat: str = "sanger"
    """Quality format, options are "sanger" or "illumina". This will add an offset of 33 or 64 respectively to the qlimit value"""


class SequenzaExtractExtraArgs(SnappyModel):
    gamma: int = 60
    """scarHRD value"""

    kmin: int = 50
    """scarHRD value"""


class SequenzaFitExtraArgs(SnappyModel):
    N_ratio_filter: int = Field(10, alias="N.ratio.filter")
    N_BAF_filter: int = Field(1, alias="N.BAF.filter")
    segment_filter: int = Field(3000000, alias="segment.filter")
    mufreq_treshold: float = Field(0.1, alias="mufreq.threshold")
    ratio_priority: bool = Field(False, alias="ratio_priority")
    ploidy: list[float] = [
        1.0,
        1.1,
        1.2,
        1.3,
        1.4,
        1.5,
        1.6,
        1.7,
        1.8,
        1.9,
        2.0,
        2.1,
        2.2,
        2.3,
        2.4,
        2.5,
        2.6,
        2.7,
        2.8,
        2.9,
        3.0,
        3.1,
        3.2,
        3.3,
        3.4,
        3.5,
        3.6,
        3.7,
        3.8,
        3.9,
        4.0,
        4.1,
        4.2,
        4.3,
        4.4,
        4.5,
        4.6,
        4.7,
        4.8,
        4.9,
        5.0,
        5.1,
        5.2,
        5.3,
        5.4,
        5.5,
    ]


class Sequenza(SnappyModel):
    length: int = 50
    assembly: str = "hg19"
    """Must be hg38 for GRCh38. See copynumber for complete list (augmented with hg38)"""

    extra_args: SequenzaExtraArgs | dict[str, Any] = {}
    """Extra arguments for sequenza bam2seqz"""

    ignore_chroms: list[str] = [
        "X",
        "Y",
        "MT",
        "NC_007605. hs37d5",
        "chrEBV",
        "*_decoy",
        "HLA-*",
        "GL000220.*",
    ]
    """patterns of chromosome names to ignore"""

    extra_args_extract: SequenzaExtractExtraArgs | dict[str, Any] = SequenzaExtractExtraArgs()
    """Valid arguments: see ?sequenza::sequenza.extract in R"""

    extra_args_fit: SequenzaFitExtraArgs | dict[str, Any] = SequenzaFitExtraArgs()
    """Valid arguments: see ?sequenza::sequenza.fit in R"""


class CopyWriter(SnappyModel):
    path_target_regions: str
    """Path to target regions"""

    bin_size: int = 20000  # TODO: make actually configurable

    plot_genes: str
    """Path to civic annotation"""

    genome: str = "hg19"
    """Could be hg38 (consider setting prefix to 'chr' when using GRCh38.v1)"""

    features: str = "EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75"

    prefix: str = ""

    nThread: int = 8


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
    genome_name: Annotated[
        GenomeName | Literal["unknown"],
        EnumField(GenomeName, json_schema_extra={"options": {"unknown"}}),
    ] = "unknown"
    """Must be one from hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6, canFam3"""

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
    extra_commands: dict[str, Any] = {
        "model": "betabin",
        "fun-segmentation": "PSCBS",
        "post-optimize": "",
    }
    """Recommended extra arguments for PureCN, extra_commands: {} to clear them all"""

    path_container: Annotated[
        str, Field(examples=["../panel_of_normals/work/containers/out/purecn.simg"])
    ]
    """
    A PureCN panel of normals is required,
    with the container, the intervals & the PON rds file
    """

    path_intervals: Annotated[
        str,
        Field(
            examples=[
                "../panel_of_normals/output/purecn/out/<enrichement_kit_name>_<genome_name>.list"
            ]
        ),
    ]

    path_panel_of_normals: Annotated[
        str,
        Field(
            examples=["../panel_of_normals/output/bwa.purecn/out/bwa.purecn.panel_of_normals.rds"]
        ),
    ]
    """Path to the PureCN panel of normal"""

    path_mapping_bias: Annotated[
        str,
        Field(examples=["../panel_of_normals/output/bwa.purecn/out/bwa.purecn.mapping_bias.rds"]),
    ]
    """Path to the PureCN mapping bias file"""

    somatic_variant_caller: str = "mutect2"
    """
    IMPORTANT NOTE:
    Mutect2 must be called with "--genotype-germline-sites true --genotype-pon-sites true
    """

    path_somatic_variants: Annotated[str, Field(examples=["../somatic_variant_calling_for_purecn"])]


class CnvettiOnTarget(SnappyModel):
    path_target_regions: str


class CnvettiOffTarget(SnappyModel):
    path_target_regions: str

    window_length: int = 20000


class SomaticTargetedSeqCnvCalling(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.purecn], min_length=1)]
    path_ngs_mapping: str = "../ngs_mapping"

    sequenza: Sequenza | None = None
    copywriter: CopyWriter | None = None
    purecn: PureCn | None = None
    cnvetti_on_target: CnvettiOnTarget | None = None
    cnvetti_off_target: CnvettiOffTarget | None = None
