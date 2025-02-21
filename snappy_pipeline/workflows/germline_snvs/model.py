import enum

from typing import Annotated

from pydantic import Field, AfterValidator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators

from snappy_pipeline.models.bcftools import (
    BcftoolsModel,
    BcftoolsBamFlagMultipleTypes,
    BcftoolsBamFlag,
    BcftoolsCaller,
    BcftoolsPloidy,
    transform_to_flag,
    ensure_valid_tags,
    ensure_valid_mark_sites_tag,
)


class Tool(enum.StrEnum):
    bcftools = "bcftools"


class Pileup(BcftoolsModel):
    max_depth: Annotated[int, Field(ge=0)] = 5000
    """Max raw per-file depth; avoids excessive memory usage (default 5000 greater than bcftools default 250)"""

    min_MQ: Annotated[int, Field(ge=0)] = 35
    """Skip alignments with mapQ smaller than INT [0]"""
    min_BQ: Annotated[int, Field(ge=0)] = 20
    """Skip bases with baseQ/BAQ smaller than INT [1]"""

    skip_all_set: Annotated[
        BcftoolsBamFlagMultipleTypes,
        AfterValidator(transform_to_flag),
    ] = BcftoolsBamFlag.NONE
    """Skip reads with all of the bits set"""
    skip_any_set: Annotated[
        BcftoolsBamFlagMultipleTypes,
        AfterValidator(transform_to_flag),
    ] = (
        BcftoolsBamFlag.UNMAP
        | BcftoolsBamFlag.MUNMAP
        | BcftoolsBamFlag.SECONDARY
        | BcftoolsBamFlag.QCFAIL
        | BcftoolsBamFlag.DUP
        | BcftoolsBamFlag.SUPPLEMENTARY
    )
    """Skip reads with any of the bits set [UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY]"""
    skip_all_unset: Annotated[
        BcftoolsBamFlagMultipleTypes,
        AfterValidator(transform_to_flag),
    ] = BcftoolsBamFlag.NONE
    """Skip reads with all of the bits unset"""
    skip_any_unset: Annotated[
        BcftoolsBamFlagMultipleTypes,
        AfterValidator(transform_to_flag),
    ] = [BcftoolsBamFlag.PAIRED, BcftoolsBamFlag.PROPER_PAIR]


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
    min_depth: Annotated[int, Field(ge=0)] = 50
    """Minimum depth to consider a variant"""
    min_baf: Annotated[float, Field(ge=0, le=0.5)] = 0.45
    """Minimum B-allele fraction (always between 0 & 0.5)"""
    min_qual: Annotated[float, Field(ge=0)] = 20.0
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

    columns: Annotated[list[str], AfterValidator(ensure_valid_tags)] = []
    """INFO & FORMAT tags to collect from annotation file"""

    mark_sites: Annotated[str | None, AfterValidator(ensure_valid_mark_sites_tag)] = None
    """Annotate sites which are present ("+") or absent ("-") in the annotation file with a new INFO/TAG flag"""

    header_line: str | None = None
    """Header line which should be appended to the VCF header"""


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
