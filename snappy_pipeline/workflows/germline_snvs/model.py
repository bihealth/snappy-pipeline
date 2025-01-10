import enum
import re

from typing import Annotated, Any

from pydantic import Field, model_validator, AliasGenerator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators

VCF_TAG_PATTERN = re.compile(r"^((INFO|FORMAT|FMT)/)?([A-Za-z_][0-9A-Za-z_.]*|1000G)$")
ANNOTATION_VCF_TAG_PATTERN = re.compile(
    (
        r"^([\^=~\+\.-]|\.\+)?"
        r"(((INFO|FORMAT|FMT)/)?([A-Za-z_][0-9A-Za-z_.]*|1000G))"
        r"(:=(((INFO|FORMAT|FMT)/)?([A-Za-z_][0-9A-Za-z_.]*|1000G)))?"
    )
)


class Tool(enum.StrEnum):
    bcftools = "bcftools"


def BcftoolsModelAliasGenerator(option: str, exceptions: dict[str, str] = {}) -> str:
    if option in exceptions:
        return exceptions[option]
    return option.replace("_", "_")


class BcftoolsModel(SnappyModel):
    model_config = {
        "alias_generator": AliasGenerator(
            serialization_alias=lambda x: BcftoolsModelAliasGenerator(x)
        )
    }

    extra_args: Annotated[
        dict[str, Any],
        Field(
            examples=[
                {"'max-depth'": 5000},
                {"columns": ["FORMAT/DP", "-INFO/SCBZ"]},
                {"'illumina1.3+'": True},
            ],
            alias="extra_args",
        ),
    ] = {}
    """
    Placeholder for arguments available to the tool, but not present in the model.
    The arguments must be input as 'option: value'. Values should not contain the apostrophe (') character.

    The user is responsible to check that the arguments make sense
    (for example that the maximum allowed value is greater that the minimum allowed value).

    To turn off a flag set in the model, the user should add the flag to the extra arguments and set it to False.
    """


# class BcftoolsBamFlag(enum.StrEnum):
#     PAIRED = enum.auto()
#     PROPER_PAIR = enum.auto()
#     UNMAP = enum.auto()
#     MUNMAP = enum.auto()
#     REVERSE = enum.auto()
#     MREVERSE = enum.auto()
#     READ1 = enum.auto()
#     READ2 = enum.auto()
#     SECONDARY = enum.auto()
#     QCFAIL = enum.auto()
#     DUP = enum.auto()
#     SUPPLEMENTARY = enum.auto()
#
#
# def check_enum_name(x: Any) -> list[int]:
#     if isinstance(x, BcftoolsBamFlags):
#         return [x.value]
#     if isinstance(x, str):
#         if x not in BcftoolsBamFlags.__members__:
#             raise ValueError(f"Illegal flag '{x}'")
#         return [BcftoolsBamFlags[x].value]
#     if isinstance(x, int):
#         all_flags = 0
#         l = []
#         for v in BcftoolsBamFlags:
#             all_flags |= v.value
#             if x & v != 0:
#                 l.append(v.value)
#         if x & all_flags != 0:
#             raise ValueError(f"Unexpected flags in '{x}'")
#         return list(set(l))
#     else:
#         raise ValueError(f"Illegal object '{x}'")
#
#
# def convert_enum_names(x: list[Any] | BcftoolsBamFlags | str | int) -> list[int]:
#     if isinstance(x, list):
#         l = []
#         for v in x:
#             l += check_enum_name(v)
#         l = list(set(l))
#     else:
#         l = check_enum_name(x)
#     return l


class BcftoolsBamFlag(enum.StrEnum):
    PAIRED = "PAIRED"
    PROPER_PAIR = "PROPER_PAIR"
    UNMAP = "UNMAP"
    MUNMAP = "MUNMAP"
    REVERSE = "REVERSE"
    MREVERSE = "MREVERSE"
    READ1 = "READ1"
    READ2 = "READ2"
    SECONDARY = "SECONDARY"
    QCFAIL = "QCFAIL"
    DUP = "DUP"
    SUPPLEMENTARY = "SUPPLEMENTARY"


class Pileup(BcftoolsModel):
    max_depth: int = 5000
    """Max raw per-file depth; avoids excessive memory usage (default 5000 greater than bcftools default 250)"""

    min_MQ: int = 35
    """Skip alignments with mapQ smaller than INT [0]"""
    min_BQ: int = 20
    """Skip bases with baseQ/BAQ smaller than INT [1]"""

    # skip_all_set: list[BcftoolsBamFlag | int] | int = []
    # """Skip reads with all of the bits set"""
    # skip_any_set: list[BcftoolsBamFlag | int] | int = [
    #     BcftoolsBamFlag.UNMAP,
    #     BcftoolsBamFlag.SECONDARY,
    #     BcftoolsBamFlag.QCFAIL,
    #     BcftoolsBamFlag.DUP
    # ]
    # """Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP]"""
    # skip_all_unset: list[BcftoolsBamFlag | int] | int = []
    # """Skip reads with all of the bits unset"""
    # skip_any_unset: list[BcftoolsBamFlag | int] | int = []
    # ""Skip reads with any of the bits unset"""

    # validate_skip_all_set = field_validator("skip_all_set", mode="before")(convert_enum_names)
    # validate_skip_any_set = field_validator("skip_any_set", mode="before")(convert_enum_names)
    # validate_skip_all_unset = field_validator("skip_all_unset", mode="before")(convert_enum_names)
    # validate_skip_any_unset = field_validator("skip_any_unset", mode="before")(convert_enum_names)

    skip_all_set: list[BcftoolsBamFlag] = []
    """Skip reads with all of the bits set"""
    skip_any_set: list[BcftoolsBamFlag] = [
        BcftoolsBamFlag.UNMAP,
        BcftoolsBamFlag.SECONDARY,
        BcftoolsBamFlag.QCFAIL,
        BcftoolsBamFlag.DUP,
        BcftoolsBamFlag.SUPPLEMENTARY,
    ]
    """Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY]"""
    skip_all_unset: list[BcftoolsBamFlag] = []
    """Skip reads with all of the bits unset"""
    skip_any_unset: list[BcftoolsBamFlag] = [BcftoolsBamFlag.PAIRED, BcftoolsBamFlag.PROPER_PAIR]
    """Skip reads with any of the bits unset [PAIRED,PROPER_PAIR]"""


class BcftoolsCaller(enum.StrEnum):
    CONSENSUS = "consensus"
    MULTIALLELIC = "multiallelic"


class BcftoolsPloidy(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    X = "X"
    Y = "Y"
    HAPLOID = "1"
    DIPLOID = "2"


class Call(BcftoolsModel):
    caller: BcftoolsCaller = BcftoolsCaller.CONSENSUS
    """The calling method"""
    variants_only: bool = True
    """Output variant sites only"""

    ploidy: BcftoolsPloidy = BcftoolsPloidy.DIPLOID
    """Predefined ploidy [2]"""

    # @field_serializer("caller")
    # def update_caller_value(x: str) -> str:
    #     return x + "-caller"


class Filter(BcftoolsModel):
    min_depth: int = 50
    """Minimum depth to consider a variant"""
    min_baf: float = 0.45
    """Minimum B-allele fraction (always between 0 & 0.5)"""
    min_qual: float = 20.0
    """Minimum genotyping quality"""


class Annotation(BcftoolsModel):
    path_annotation: Annotated[
        str,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/GATK_Best_Practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/GATK_Best_Practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/DRAGEN/hg38/hg38_1000G_phase1.snps.high_confidence.vcf.gz",
            ]
        ),
    ]
    """
    Path to the annotation file. This must be a vcf/bcf file, tab-delimited ones are not supported.
    """

    columns: list[str] = []

    mark_sites: str | None = None
    """Annotate sites which are present ("+") or absent ("-") in the annotation file with a new INFO/TAG flag"""

    header_line: str | None = None
    """Header line which should be appended to the VCF header"""

    @model_validator(mode="after")
    def ensure_valid_mark_sites_tag(self):
        if self.mark_sites:
            if not (self.mark_sites.startswith("+") or self.mark_sites.startswith("-")):
                raise ValueError(
                    f"Mark sites tag '{self.mark_sites} must be prefixed with '+' (mark present sites) or '-' (mark absent sites)"
                )
            if not VCF_TAG_PATTERN.match(self.mark_sites[1:]):
                raise ValueError(f"Illegal Mark sites tag '{self.mark_sites[1:]}'")
        return self

    @model_validator(mode="after")
    def ensure_valid_tags(self):
        for x in self.columns:
            m = ANNOTATION_VCF_TAG_PATTERN.match(x)
            if not m:
                raise ValueError(f"Illegal annotation tag pattern '{x}'")
            g = m.groups()
            if not VCF_TAG_PATTERN.match(g[1]):
                raise ValueError(f"Illegal annotation tag pattern '{x}'")
            if g[5] is not None and not VCF_TAG_PATTERN.match(g[6]):
                raise ValueError(f"Illegal annotation tag pattern '{x}'")
        return self


class BcfTools(SnappyModel):
    skip_indels: bool = True
    """Do not perform indel calling"""

    path_regions: Annotated[
        str | None,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/GATK_Best_Practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/GATK_Best_Practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/DRAGEN/hg38/hg38_1000G_phase1.snps.high_confidence.vcf.gz",
            ]
        ),
    ] = None
    """
    Path to the regions where to perform germline variant calling. Done everywhere when missing.
    """

    ignore_chroms: Annotated[
        list[str],
        Field(
            examples=[
                ["MT", "chrM"],
                ["MT", "GL000*"],
                ["chrM", "chrUn_*", "*_random"],
            ]
        ),
    ] = []
    """Contigs to be ignored during processing. Joined with ``ignored_chroms`` from upper level"""

    pileup: Pileup = Pileup()
    call: Call = Call()
    filter: Filter = Filter()
    annotate: list[Annotation] = []


class GermlineSnvs(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.bcftools], min_length=1)]

    tool_ngs_mapping: str = "bwa"
    """
    Configuration with read mapper and path to mapping output.
    Will use this for generating coverages over autosomes & sex chromosomes using bcftools.
    """

    path_ngs_mapping: str = "../ngs_mapping"

    ignore_chroms: Annotated[
        list[str],
        Field(
            examples=[
                ["HLA*", "*decoy", "*alt", "EBV"],
                ["GL000220.*", "hs37d5", "NC_007605"],
                [
                    "chrUn_GL000220*",
                    "*_decoy",
                    "chrEBV",
                    "CMV",
                    "HBV",
                    "HCV-*",
                    "HIV-*",
                    "KSHV",
                    "HTLV-1",
                    "MCV",
                    "SV40",
                    "HPV*",
                ],
            ]
        ),
    ] = []
    """Contigs to be ignored during processing"""

    bcftools: BcfTools = BcfTools()
