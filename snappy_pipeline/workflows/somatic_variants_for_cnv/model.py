import enum

from typing import Annotated

from pydantic import Field, AfterValidator, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators

from snappy_pipeline.models.bcftools import (
    BcftoolsModel,
    BcftoolsBamFlag,
    transform_to_flag,
    ensure_valid_tags,
    ensure_valid_mark_sites_tag,
)


class Tool(enum.StrEnum):
    bcftools = "bcftools"


class Pileup(BcftoolsModel):
    max_depth: int = 5000
    """Max raw per-file depth; avoids excessive memory usage (default 5000 greater than bcftools default 250)"""

    min_MQ: int = 35
    """Skip alignments with mapQ smaller than INT [35]"""
    min_BQ: int = 20
    """Skip bases with baseQ/BAQ smaller than INT [20]"""

    full_BAQ: bool = True
    """Apply BAQ everywhere, not just in problematic regions"""
    redo_BAQ: bool = True
    """Recalculate BAQ on the fly, ignore existing BQs"""

    skip_all_set: Annotated[
        str | int | BcftoolsBamFlag | list[str] | list[int] | list[BcftoolsBamFlag],
        AfterValidator(transform_to_flag),
    ] = BcftoolsBamFlag.NONE
    """Skip reads with all of the bits set"""
    skip_any_set: Annotated[
        str | int | BcftoolsBamFlag | list[str] | list[int] | list[BcftoolsBamFlag],
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
        str | int | BcftoolsBamFlag | list[str] | list[int] | list[BcftoolsBamFlag],
        AfterValidator(transform_to_flag),
    ] = BcftoolsBamFlag.NONE
    """Skip reads with all of the bits unset"""
    skip_any_unset: Annotated[
        str | int | BcftoolsBamFlag | list[str] | list[int] | list[BcftoolsBamFlag],
        AfterValidator(transform_to_flag),
    ] = [BcftoolsBamFlag.PAIRED, BcftoolsBamFlag.PROPER_PAIR]
    """Skip reads with any of the bits unset [PAIRED,PROPER_PAIR]"""


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

    regions_overlap: Annotated[int, Field(ge=0, le=2)] = 2
    """Include if POS in the region (0), record overlaps (1), variant overlaps (2) [2]"""

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
    annotate: list[Annotation] = []


class SomaticVariantsForCnv(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.bcftools], min_length=1)]

    tool_ngs_mapping: str = "bwa"
    """
    Configuration with read mapper and path to mapping output.
    Will use this for generating coverages over autosomes & sex chromosomes using bcftools.
    """

    path_ngs_mapping: str = "../ngs_mapping"

    tool_germline_snvs: str = "bcftools"

    path_germline_snvs: str = "../germline_snvs"

    tool_somatic_variant_calling: str = "mutect2"

    path_somatic_variant_calling: str = "../somatic_variant_calling"

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
