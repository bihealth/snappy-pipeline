import enum
from typing import Annotated

from pydantic import Field, AfterValidator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
from snappy_pipeline.models.common import Sex, LibraryKitEntry
from snappy_pipeline.models.bcftools import BcftoolsBamFlag, transform_to_flag


class Tool(enum.StrEnum):
    ascat = "ascat"


class AlleleCounter(SnappyModel):
    loci_prefix: Annotated[
        str,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/ASCAT/G1000_loci_WGS_hg38/G1000_loci_WGS_hg38_chr",
            ]
        ),
    ]
    """Path to the locii & prefix"""

    allele_prefix: Annotated[
        str | None,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/ASCAT/G1000_loci_WGS_hg38/G1000_alleles_WGS_hg38_chr",
            ]
        ),
    ]
    """Path to the alleles & prefix"""

    path_probloci_file: str | None = None
    """A file (chromosome <tab> position; no header) containing specific loci to ignore"""

    minCounts: int = 10
    """Minimum depth required in the normal for a SNP to be considered"""

    min_BQ: int = 20
    """Minimum base quality required for a read to be counted"""

    min_MQ: int = 35
    """Minimum mapping quality required for a read to be counted"""

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
    """
    Exclude reads with any on these flags:

    - 2048 0x800 SUPPLEMENTARY supplementary alignment
    - 1024 0x400 DUP           PCR or optical duplicate
    -  512 0x200 QCFAIL        not passing filters, such as platform/vendor quality controls
    -  256 0x100 SECONDARY     secondary alignment
    -    8 0x008 MUNMAP        next segment in the template unmapped
    -    4 0x004 UNMAP         segment unmapped
    """
    skip_any_unset: Annotated[
        str | int | BcftoolsBamFlag | list[str] | list[int] | list[BcftoolsBamFlag],
        AfterValidator(transform_to_flag),
    ] = BcftoolsBamFlag.PAIRED | BcftoolsBamFlag.PROPER_PAIR
    """
    Exclude reads without any on these flags:

    -    2 0x002 PROPER_PAIR each segment properly aligned according to the aligner
    -    1 0x001 PAIRED      template having multiple segments in sequencing
    """

    max_depth: int = 8000
    """Max depth to used in ``bcftools`` (not part of the original implementation)"""


class GenomeVersionAscat(enum.StrEnum):
    GRCh37 = "hg19"
    GRCh38 = "hg38"


class AscatPlatform(enum.StrEnum):
    AffySNP6 = "AffySNP6"
    Custom10k = "Custom10k"
    IlluminaASA = "IlluminaASA"
    IlluminaGSAv3 = "IlluminaGSAv3"
    Illumina109k = "Illumina109k"
    IlluminaCytoSNP = "IlluminaCytoSNP"
    IlluminaCytoSNP850k = "IlluminaCytoSNP850k"
    Illumina610k = "Illumina610k"
    Illumina660k = "Illumina660k"
    Illumina700k = "Illumina700k"
    Illumina1M = "Illumina1M"
    Illumina25M = "Illumina2.5M"
    IlluminaOmni15 = "IlluminaOmni15"
    IlluminaGDACyto8 = "IlluminaGDACyto-8"
    Affy10k = "Affy10k"
    Affy100k = "Affy100k"
    Affy250k_sty = "Affy250k_sty"
    Affy250k_nsp = "Affy250k_nsp"
    AffyOncoScan = "AffyOncoScan"
    AffyCytoScanHD = "AffyCytoScanHD"
    HumanCNV370quad = "HumanCNV370quad"
    HumanCore12 = "HumanCore12"
    HumanCoreExome24 = "HumanCoreExome24"
    HumanOmniExpress12 = "HumanOmniExpress12"
    IlluminaOmniExpressExome = "IlluminaOmniExpressExome"
    WGS_hg_50X = "WGS_hg_50X"


class AscatAdvancedParams(SnappyModel):
    penalty: float = 70.0
    """Segmentation penalty"""

    gamma: float = 1.0
    """Compactation of LogR profiles. Must be set to 1 for WGS, WES & Panel"""

    min_purity: float = 0.1
    """Minimum boundary of the purity solution search space"""

    max_purity: float = 1.05
    """Maximum boundary of the purity solution search space"""

    min_ploidy: float = 1.5
    """Minimum boundary of the ploidy solution search space"""

    max_ploidy: float = 5.5
    """Maximum boundary of the ploidy solution search space"""

    rho_manual: float | str = "NA"
    """When present, overrides ASCAT optimization and supply rho parameter"""

    psi_manual: float | str = "NA"
    """When present, overrides ASCAT optimization and supply psi parameter"""


class Ascat(SnappyModel):
    allele_counter: AlleleCounter

    genomeVersion: GenomeVersionAscat

    seed: int = 1234567

    path_gc_content: Annotated[
        str,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/ASCAT/GC_G1000_WGS_hg38/GC_G1000_hg38.txt",
            ]
        ),
    ]
    """Path to the GC correction file"""

    path_reptiming: Annotated[
        str,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/ASCAT/RT_G1000_WGS_hg38/RT_G1000_hg38.txt",
            ]
        ),
    ]
    """Path to the correction file for DNA replication timing"""

    y_limit: float = 5
    """Determes the size of the y axis in the nonrounded plot and ASCAT profile"""

    advanced: AscatAdvancedParams = AscatAdvancedParams()

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
    """
    Contigs to be ignored during processing. Joined with ``ignored_chroms`` from upper level.
    Note that ``ascat`` only considers chromosomes 1 to 22 & X anyway.
    """


class SomaticCnvCalling(SnappyStepModel, validators.ToolsMixin):
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.ascat], min_length=1)]

    tool_ngs_mapping: str = "bwa"
    """
    Configuration with read mapper and path to mapping output.
    Will use this for generating a pileup using samtools for obtaining the b allele fraction and computing coverage.
    """

    path_ngs_mapping: str = "../ngs_mapping"

    sex: Sex

    path_target_interval_list_mapping: list[LibraryKitEntry] = []

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

    ascat: Ascat | None = None
