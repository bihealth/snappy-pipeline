import enum
from typing import Self

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel


class SexOrigin(enum.StrEnum):
    AUTOMATIC = "auto"
    """Sex determined from the data"""
    SAMPLESHEET = "samplesheet"
    """Donor sex obtained from sample sheet"""
    CONFIG = "config"
    """Donor sex obtained from the configuration (all donors have the same sex)"""


class SexValue(enum.StrEnum):
    MALE = "male"
    FEMALE = "female"


class Sex(SnappyModel):
    source: SexOrigin = SexOrigin.AUTOMATIC

    sample_sex: SexValue | None = None

    @model_validator(mode="after")
    def ensure_valid_sex_value(self):
        if self.source == SexOrigin.CONFIG and self.sample_sex is None:
            raise ValueError("No definition of donors' sex from the configuration")
        return self


class SegmentationMethod(enum.StrEnum):
    cbs = "cbs"
    flasso = "flasso"
    haar = "haar"
    hmm = "hmm"
    hmm_tumor = "hmm-tumor"
    hmm_germline = "hmm-germline"
    none = "none"


class CenterMethod(enum.StrEnum):
    mean = "mean"
    median = "median"
    mode = "mode"
    biweight = "biweight"


class FilterMethod(enum.StrEnum):
    ampdel = "ampdel"
    cn = "cn"
    ci = "ci"
    sem = "sem"


class CallingMethod(enum.StrEnum):
    threshold = "threshold"
    clonal = "clonal"
    none = ""


class Access(SnappyModel):
    exclude: list[str] = []
    """Regions accessible to mapping"""
    min_gap_size: int = 5000
    """Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together."""


class Target(SnappyModel):
    path_baits: str | None = None
    """Path to baits file (Agilent Covered), unset for WGS data"""
    split: bool = True
    """Split large tiled intervals into smaller, consecutive targets."""
    avg_size: float = 800 / 3
    """Average size of split target bins (results are approximate)"""
    short_names: bool = False
    """Reduce multi-accession bait labels to be short and consistent"""


class Antitarget(SnappyModel):
    avg_size: float = 150000
    """Average size of split antitarget bins (results are approximate)"""
    min_size: float | None = None
    """Minimum size of antitarget bins (smaller regions are dropped). When missing, 1/16 avg size"""


class Coverage(SnappyModel):
    count: bool = False
    """Get read depths by counting read midpoints within each bin."""
    min_mapq: int = 0
    """Minimum mapping quality score (phred scale 0-60) to count a read for coverage depth."""


class Fix(SnappyModel):
    smoothing_window_fraction: float | None = None
    """Smoothing window fraction for rolling median bias smoothing. Defaults to 1/sqrt(len(data))"""


class Segment(SnappyModel):
    method: SegmentationMethod = SegmentationMethod.CBS
    """Segmentation method, or 'NONE' for chromosome arm-level averages as segments"""
    threshold: float = 0.0001
    """Significance threshold (p-value or FDR, depending on method) to accept breakpoints during segmentation. For HMM methods, this is the smoothing window size."""
    drop_outliers: int = 10
    """Drop outlier bins more than this many multiples of the 95th quantile away from the average within a rolling window. Set to 0 for no outlier filtering."""
    smooth_cbs: bool = False

    @model_validator(mode="after")
    def ensure_smooth_for_cbs_only(self) -> Self:
        if self.smooth_cbs and self.method != SegmentationMethod.CBS:
            raise ValueError("'smooth_cbs' option can be used only with 'CBS' segmentation method")
        return self


class Call(SnappyModel):
    method: CallingMethod | None = None
    """Calling method."""
    thresholds: list[float] = [-1.1, -0.25, 0.2, 0.7]
    """Hard thresholds for calling each integer copy number, separated by commas"""
    center: CenterMethod = CenterMethod.MEDIAN
    """Re-center the log2 ratio values using this estimator of the center or average value. ('median' if no argument given.)"""
    center_at: float | None = None
    """Subtract a constant number from all log2 ratios. For "manual" re-centering."""
    filter: FilterMethod | None = None
    """Merge segments flagged by the specified filter(s) with the adjacent segment(s)."""

    @model_validator(mode="after")
    def avoid_center_center_at_conflict(self) -> Self:
        if self.center is not None and self.center_at is not None:
            raise ValueError("'call' options 'center' and 'center_at' cannot be used together")
        return self


class Bintest(SnappyModel):
    alpha: float = 0.005
    """Significance threhold."""
    target: bool = False
    """Test target bins only; ignore off-target bins."""


class Plot(SnappyModel):
    enabled: bool = False


class PlotDiagram(Plot):
    chromosome: str | None = None
    """Chromosome to display (full genome when missing)"""
    threshold: float = 0.5
    """Copy number change threshold to label genes."""
    min_probes: int = 3
    """Minimum number of covered probes to label a gene."""
    no_shift_xy: bool = False


class PlotScatter(Plot):
    path_range_list: str | None = None
    """File listing the chromosomal ranges to display, as BED, interval list or 'chr:start-end' text"""
    gene: str | None = None
    """Name of gene or genes (comma-separated) to display."""
    width: int = 1000000
    """Width of margin to show around the selected gene(s)"""
    antitarget_marker: str = "o"
    """Plot antitargets using this symbol when plotting in a selected chromosomal region."""
    by_bin: bool = False
    """Plot data x-coordinates by bin indices instead of genomic coordinates."""
    segment_color: str = "darkorange"
    """Plot segment lines in this color. Value can be any string accepted by matplotlib."""
    trend: bool = False
    """Draw a smoothed local trendline on the scatter plot."""
    y_max: float | None = None
    """y-axis upper limit."""
    y_min: float | None = None
    """y-axis lower limit."""
    fig_size: tuple[float, float] = (6.4, 4.8)
    """Width and height of the plot in inches."""

    @model_validator(mode="after")
    def ensure_range_list_with_gene(self) -> Self:
        if self.gene is not None and not self.path_range_list:
            raise ValueError("'gene' option requires a valid range list")
        return self


class Report(SnappyModel):
    enabled: bool = True


class ReportSegmetrics(Report):
    alpha: float = 0.05
    """Level to estimate confidence and prediction intervals; use with --ci and --pi."""
    bootstrap: int = 100
    """Number of bootstrap iterations to estimate confidence interval; use with --ci."""
    smooth_bootstrap: bool = False
    """Apply Gaussian noise to bootstrap samples, a.k.a. smoothed bootstrap, to estimate confidence interval"""


class ReportGenemetrics(Report):
    alpha: float = 0.05
    """Level to estimate confidence and prediction intervals; use with --ci and --pi."""
    bootstrap: int = 100
    """Number of bootstrap iterations to estimate confidence interval; use with --ci."""
    threshold: float = 0.2
    """Copy number change threshold to report a gene gain/loss"""
    min_probes: int = 3
    """Minimum number of covered probes to report a gain/loss"""


class Report(enum.StrEnum):
    GENEMETRICS = "genemetrics"
    SEGMETRICS = "segmetrics"


class CnvkitToReference(SnappyModel):
    # Substep-secific parameters
    access: Access
    target: Target
    antitarget: Antitarget

    coverage: Coverage

    metrics: Report
    segmetrics: ReportSegmetrics
    genemetrics: ReportGenemetrics

    # Generic parameters (used in different substeps & must agree)
    male_reference: bool = False
    """Create/use male reference (for shifting chrX & chrY)"""
    diploid_parx_genome: str | None = None
    """Considers the given human genome's PAR of chromosome X as autosomal. Example: 'grch38'"""

    cluster: bool = False
    """Calculate and store summary stats for clustered subsets of the normal samples with similar coverage profiles."""
    min_cluster_size: int = 4
    """Minimum cluster size to keep in reference profiles."""

    gc: bool = False
    """Skip GC correction."""
    edge: bool = None
    """Skip edge correction. Automatic selection when None (True for WGS & Panel, False for WES)"""
    rmask: bool = False
    """Skip RepeatMasker correction."""

    drop_low_coverage: bool = False
    """Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples."""

    @model_validator(mode="after")
    def ensure_males_for_reference(self):
        if (
            self.male_reference
            and self.sex.source == SexOrigin.CONFIG
            and self.sex.sample_sex == SexValue.FEMALE
        ):
            raise ValueError("Male reference requested for female cohort")
        return self


class Cnvkit(CnvkitToReference):
    fix: Fix
    segment: Segment
    call: Call
    bintest: Bintest

    diagram: PlotDiagram
    scatter: PlotScatter

    min_variant_depth: int = 20
    """Minimum read depth for a SNV to be displayed in the b-allele frequency plot."""
    zygocity_freq: float = 0.25
    """Ignore VCF's genotypes (GT field) and instead infer zygosity from allele frequencies."""
