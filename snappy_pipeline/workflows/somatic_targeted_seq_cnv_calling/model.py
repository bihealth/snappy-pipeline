import enum
from typing import Annotated, Any, Literal

from pydantic import ConfigDict, Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.cnvkit import Cnvkit


class Tool(enum.StrEnum):
    cnvkit = "cnvkit"
    sequenza = "sequenza"
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
    model_config = ConfigDict(populate_by_name=True)

    N_ratio_filter: int = Field(10, alias="N.ratio.filter")
    N_BAF_filter: int = Field(1, alias="N.BAF.filter")
    segment_filter: int = Field(3000000, alias="segment.filter")
    mufreq_treshold: float = Field(0.1, alias="mufreq.treshold")
    ratio_priority: bool = Field(False, alias="ratio.priority")
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

    extra_args: Annotated[SequenzaExtraArgs | dict[str, Any], Field(union_mode="left_to_right")] = (
        SequenzaExtraArgs()
    )
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

    extra_args_extract: Annotated[
        SequenzaExtractExtraArgs | dict[str, Any], Field(union_mode="left_to_right")
    ] = SequenzaExtractExtraArgs()
    """Valid arguments: see ?sequenza::sequenza.extract in R"""

    extra_args_fit: Annotated[
        SequenzaFitExtraArgs | dict[str, Any], Field(union_mode="left_to_right")
    ] = SequenzaFitExtraArgs()
    """Valid arguments: see ?sequenza::sequenza.fit in R"""


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
    """Path to the PureCN apptainer/singularity container image"""

    somatic_variant_caller: str = "mutect2"
    """
    IMPORTANT NOTE:
    Mutect2 must be called with "--genotype-germline-sites true --genotype-pon-sites true
    """


class SomaticTargetedSeqCnvCallingDependsOn(SnappyModel):
    somatic_variants: str = "somatic_variants"
    ngs_mapping: str = "ngs_mapping"
    panel_of_normals: str = ""
    """
    Required when ``tool: cnvkit`` or ``tool: purecn``.
    Must name the upstream ``panel_of_normals`` task that produced the matching PON
    (e.g. ``panel_of_normals_cnvkit`` or ``panel_of_normals_purecn``).
    Snakemake tracks the PON outputs as proper input files via this dependency.
    """


class SomaticTargetedSeqCnvCalling(SnappyStepModel):
    depends_on: SomaticTargetedSeqCnvCallingDependsOn = Field(
        default_factory=SomaticTargetedSeqCnvCallingDependsOn
    )

    tool: Annotated[Tool, EnumField(Tool, default=Tool.cnvkit)]

    cnvkit: Cnvkit | None = None
    sequenza: Sequenza | None = None
    purecn: PureCn | None = None

    @model_validator(mode="after")
    def validate_panel_of_normals_dependency(self) -> "SomaticTargetedSeqCnvCalling":
        """Enforce explicit panel-of-normals dependency for cnvkit and purecn.

        Both tools require a pre-built panel of normals.  The dependency must be
        expressed via ``depends_on.panel_of_normals`` so that Snakemake can track
        the PON outputs as proper input files rather than bare config paths.
        """
        if self.tool in (Tool.cnvkit, Tool.purecn):
            if not self.depends_on.panel_of_normals:
                raise ValueError(
                    f"depends_on.panel_of_normals must be set when tool='{self.tool}'; "
                    "name the upstream panel_of_normals task that produced the matching PON "
                    f"(e.g. 'panel_of_normals_{self.tool}')"
                )
        return self
