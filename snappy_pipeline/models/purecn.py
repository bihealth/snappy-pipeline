import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel


# Parameters for each action & those shared between actions
# param_table = {
#     "shared": {
#         "genome": Enum,
#         "seed": int,
#     },
#     "IntervalFile": {
#         "off_target": bool,
#         "average_target_width": float,
#         "min_target_width": float,
#         "small_targets": Enum,
#         "off_target_seqlevels": Enum,
#         "min_mappability": list,
#         "average_reptiming_width": float,
#     },
#     "Coverage": {
#         "keep_duplicates": bool,
#         "remove_mapq0": bool,
#         "skip_gc_norm": bool,
#     },
#     "NormalDB": {
#         "genomicsdb_af_field": str,
#         "min_normals_position_specific_fit": float,
#     },
#     "PureCN": {
#         "sex": Enum,
#         "min_af": float,
#         "error": float,
#         "base_quality_offset": int,
#         "min_supporting_reads": int,
#         "db_info_flag": str,
#         "popaf_info_field": str,
#         "cosmic_cnt_info_field": str,
#         "min_cosmic_cnt": int,
#         "interval_padding": int,
#         "min_total_counts": int,
#         "min_fraction_offtarget": float,
#         "fun_segmentation": Enum,
#         "alpha": float,
#         "undo_sd": str,
#         "changpoints_penalty": int,
#         "additional_cmd_args": str,
#         "max_segments": int,
#         "min_logr_sdev": float,
#         "min_purity": float,
#         "max_purity": float,
#         "min_ploidy": float,
#         "max_ploidy": float,
#         "max_copy_number": int,
#         "post_optimize": bool,
#         "bootstrap_n": int,
#         "speedup_heuristics": int,
#         "model_homozygous": bool,
#         "model": Enum,
#         "max_non_clonal": float,
#         "max_homozygous_loss": list,
#     },
# }


class Genome(enum.StrEnum):
    HG18 = "hg18"
    HG19 = "hg19"
    HG38 = "hg38"
    MM9 = "mm9"
    MM10 = "mm10"
    RN4 = "rn4"
    RN5 = "rn5"
    RN6 = "rn6"
    CANFAM3 = "canFam3"


class SmallTargets(enum.StrEnum):
    RESIZE = "resize"
    DROP = "drop"


class OffTargetSeqLevels(enum.StrEnum):
    TARGETED = "targeted"
    ALL = "all"


class FilterMethod(enum.StrEnum):
    AMPDEL = "ampdel"
    CN = "cn"
    CI = "ci"
    SEM = "sem"


class CallingMethod(enum.StrEnum):
    THRESHOLD = "threshold"
    CLONAL = "clonal"
    NONE = "none"


class IntervalFile(SnappyModel):
    off_target: bool = False
    """Include off-target regions"""
    average_target_width: int = 400
    """Split large targets to approximately that size"""
    min_target_width: int = 100
    """Either resize or drop targets smaller than specified"""
    small_targets: SmallTargets = SmallTargets.RESIZE
    """Either 'resize' or 'drop' small targets"""
    average_off_target_width: int = 200000
    """Bin off-target regions to approximately that size"""
    off_target_seqlevels: OffTargetSeqLevels = OffTargetSeqLevels.TARGETED
    """Controls how to deal with chromosomes/contigs not found in baits"""
    mappability: Annotated[
        str,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/PureCN/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/PureCN/hg19/wgEncodeCrgMapabilityAlign75mer.bigWig",
            ]
        ),
    ] = ""
    """``rtracklayer``-parsable file with mappability scores in 1st metadata column"""
    min_mappability: tuple[float, float, float] = (0.6, 0.1, 0.7)
    """Minimum mappability for on-target, off-target and chrY regions"""
    reptiming: Annotated[
        str,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/PureCN/hg19/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
                "<no example for hg38>",
            ]
        ),
    ] = ""
    """``rtracklayer``-parsable file with replication timing scores in 1st metadata column"""
    average_reptiming_width: int = 100000
    """Average the replication timing data into bins of the specified size"""
    exclude: str | None = None
    """File parsable by rtracklayer specifying baits that should be excluded from baits file"""


class Coverage(SnappyModel):
    keep_duplicates: bool = False
    """SCount reads marked as duplicates"""
    remove_mapq0: bool = False
    """Not count reads marked with mapping quality 0"""
    skip_gc_norm: bool = False
    """Skips GC-normalization"""


class NormalDB(SnappyModel):
    genomicsdb_af_field: str = "AF"
    """Info field name where the allelic fraction is stored"""
    min_normals_position_specific_fit: float = 10.0
    """Only change if you know what you are doing"""


class PureCNBase(SnappyModel):
    genome: Genome
    """Genome version. One of hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6, canFam3"""
    seed: int | None = None
    """Seed for random number generator"""


class PureCNPon(PureCNBase):
    intervals: IntervalFile = IntervalFile()
    normaldb: NormalDB = NormalDB()
    coverage: Coverage = Coverage()


class Variant(SnappyModel):
    min_af: float = 0.03
    """minimum allelic fraction"""
    snp_blacklist: str | None = None
    """File parsable by rtracklayer that defines blacklisted regions"""
    error: float = 0.001
    """Estimated default sequencing error rate for artifact filtering. Can be overriden by base quality scores"""
    base_quality_offset: int = 1
    """Subtracts the specified value from the base quality score"""
    min_supporting_reads: int | None = None
    """Instead of calculating the min. number of supporting reads, use specified one"""
    db_info_flag: str = "DB"
    """VCF INFO flag indicating presence in common germline databases"""
    popaf_info_field: str = "POP_AF"
    """VCF INFO field providing population allele frequency"""
    cosmic_cnt_info_field: str = "Cosmic.CNT"
    """VCF INFO field providing counts in the Cosmic database"""
    cosmic_vcf_file: str | None = None
    """Adds a Cosmic.CNT INFO annotation using a Cosmic VCF. Added for convenience, we recommend adding annotations upstream"""
    min_cosmic_cnt: int = 6
    """Min number of COSMIC hits"""
    interval_padding: int = 50
    """Keep variants in the flanking region of specified size"""


class IntervalFilter(SnappyModel):
    min_total_counts: int = 100
    """Keep only intervals with at least that many counts in both tumor and (tanget) normal"""
    min_fraction_offtarget: float = 0.05
    """Ignore off-target internals when only the specified fraction of all intervals are off-target intervals"""


class SegmentationMethod(enum.StrEnum):
    CBS = "CBS"
    PSCBS = "PSCBS"
    GATK4 = "GATK4"
    HCLUST = "Hclust"


class Segmentation(SnappyModel):
    enabled: bool = True
    method: SegmentationMethod = SegmentationMethod.CBS
    alpha: float = 0.005
    """Significance of breakpoints"""
    undo_sd: str | None = None
    """DNAcopy undo.SD argument. If None, tries to find a sensible default"""
    changepoints_penalty: float | None = None
    """GATK4 ModelSegments --number-of-changepoints-penalty-factor argument. If NULL, tries to find a sensible default"""
    additional_cmd_args: str = ""
    """Used in GATK4 segmentation function to add additional ModelSegments arguments"""
    max_segments: int = 300
    """Flag noisy samples with many segments"""
    min_logr_sdev: float = 0.15
    """Set minimum log-ratio standard deviation to this value. Useful when uncorrected biases exceed the log-ratio noise"""

    seg_file: str | None = None
    """External segmentation file (from cnvkit, for example)"""

    @model_validator(mode="after")
    def ensure_args_gatk4(self):
        if self.changepoints_penalty or self.additional_cmd_args:
            if self.method != SegmentationMethod.GATK4:
                raise ValueError(
                    "Segmentation method 'GATK4' must be selected when parameters 'changepoints_penalty' or 'additional_cmd_args' are set"
                )
        return self

    @model_validator(mode="after")
    def ensure_segmentation(self):
        if self.enabled and self.seg_file is not None:
            raise ValueError("Segmentation cannot be enabled when a segmentation file is provided")
        if not self.enabled and not self.seg_file:
            raise ValueError("Segmentation must be either enabled or provided using 'seg_file'")
        return self


class Model(enum.StrEnum):
    BETA = "beta"
    BETABIN = "betabin"


class PureCN(PureCNBase):
    min_mapq: int = 0
    """Minimum mapping quality score (phred scale 0-60) to count a read for coverage depth."""
    min_purity: float = 0.15
    """Minimum considered purity"""
    max_purity: float = 0.95
    """Maximum considered purity"""
    min_ploidy: float = 1.4
    """Minimum considered loidy"""
    max_ploidy: float = 6.0
    """Maximum considered ploidy"""
    max_copy_number: int = 7
    """Maximum allele-specific integer copy number"""
    post_optimize: bool = False
    """Post-optimization"""
    bootstrap_n: int = 0
    """Number of bootstrap replicates"""
    speedup_heuristics: float = 2.0
    """Tries to avoid spending computation time on unlikely local optima"""
    homozygous_model: bool = False
    """Model homozygous variants in very pure samples. Should be 'model_homozygous', but model_* doesn't play well with pytest"""
    fit_model: Model = Model.BETA
    """Model used to fit variants. Either beta or betabin. Should be 'model', but model_* doesn't play well with pytest"""
    log_ratio_calibration: float = 0.1
    """Parameter defining the extend to which log-ratios might be miscalibrated"""
    max_non_clonal: float = 0.2
    """Maximum genomic fraction assigned to a subclonal copy number state"""
    max_homozygous_loss: tuple[float, float] = (0.05, 10000000.0)
    """Maximum genomic fraction assigned to a complete loss and maximum size of a loss in bp"""

    log_ratio_file: str | None = None
    """External log2 copy number ratio file"""

    # TODO: allow PureCN to merge all tumors from the same donor
    additional_tumors: list[str] = []
    """tumor coverages from additional biopsies from the SAME patient, GC-normalized"""

    interval_filter: IntervalFilter = IntervalFilter()
    segmentation: Segmentation = Segmentation()
