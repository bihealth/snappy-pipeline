import enum
from typing import Self

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel


# Parameters for each action & those shared between actions
# param_table = {
#     "shared": {
#         "short_names": bool,
#         "drop_low_coverage": bool,
#         "male_reference": bool,
#         "sample_sex": Enum,
#         "zigocity_freq": float,
#         "min_variant_depth": int,
#         "diploid_parx_genome": str,
#         "normal_id": str,
#         "sample_id": str,
#         "cluster": bool,
#     },
#     "access": {"min_gap_size": int, "exclude": list},
#     "antitarget": {"avg_size": int, "min_size": float},
#     "autobin": {
#         "method": Enum,
#         "bp_per_bin": float,
#         "antitarget_max_size": int,
#         "antitarget_min_size": int,
#         "target_max_size": int,
#         "target_min_size": int,
#     },
#     "bintest": {"target": bool, "alpha": float},
#     "call": {
#         "center": Enum,
#         "filter": Enum,
#         "method": Enum,
#         "center_at": float,
#         "purity": float,
#         "ploidy": float,
#         "thresholds": list,
#     },
#     "coverage": {"count": bool, "min_mapq": int},
#     "fix": {"smoothing_window_fraction": float},
#     "genemetrics": {"alpha": float, "threshold": float, "bootstrap": int},
#     "metrics": {},
#     "reference": {"min_cluster_size": int},
#     "segment": {
#         "smooth_cbs": bool,
#         "method": Enum,
#         "threshold": float,
#         "drop_outliers": int,
#     },
#     "segmetrics": {
#         "alpha": float,
#         "threshold": float,
#         "bootstrap": int,
#         "min_probes": int,
#     },
#     "target": {"split": bool, "avg_size": float},
# }


class SegmentationMethod(enum.StrEnum):
    CBS = "cbs"
    FLASSO = "flasso"
    HAAR = "haar"
    HMM = "hmm"
    HMM_TUMOR = "hmm-tumor"
    HMM_GERMLINE = "hmm-germline"
    NONE = "none"


class CenterMethod(enum.StrEnum):
    MEAN = "mean"
    MEDIAN = "median"
    MODE = "mode"
    BIWEIGHT = "biweight"


class FilterMethod(enum.StrEnum):
    AMPDEL = "ampdel"
    CN = "cn"
    CI = "ci"
    SEM = "sem"


class CallingMethod(enum.StrEnum):
    THRESHOLD = "threshold"
    CLONAL = "clonal"
    NONE = "none"


class Access(SnappyModel):
    """
    In WGS mode, the _target_ regions are set to the accessible regions in the genome.
    These accessible regions can be provided by the user, or computed by the `access`
    module. In the latter case, the optimal bin size is computed by the `autobin` module
    unless this value is provided by the user.
    `autobin` uses the `wgs` method _only_ if the list of excluded region is empty and if
    the `min_gap_size` parameter remains unassigned. If any of these conditions is not met,
    or if a files of accessible regions is provided by the user, then the `amplicon` method
    is used.
    It is recommended to leave the excluded regions empty and not set the `min_gap_size`
    parameter for WGS data, unless the accessible regions are much reduced (for example excluding
    all intergenic regions, repeats, low complexity, ...)
    """
    
    exclude: list[str] = []
    """Regions accessible to mapping"""
    min_gap_size: int | None = None
    """Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together."""


class Target(SnappyModel):
    split: bool = True
    """Split large tiled intervals into smaller, consecutive targets."""
    avg_size: float | None = None
    """
    Average size of split target bins (results are approximate).

    When the parameter is left unassigned, the cnvkit default is used for WES data,
    and an optimal value is computed for WGS data, if there is data for normal control(s).
    """
    short_names: bool = True
    """
    Reduce multi-accession bait labels to be short and consistent.

    Only valid when a gff/gtf features file is defined in the static part of the configuration.
    """


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
    """Segmentation method, or 'none' for chromosome arm-level averages as segments"""
    threshold: float
    """
    Significance threshold (p-value or FDR, depending on method) to accept breakpoints during segmentation.

    For HMM methods, this is the smoothing window size.
    """
    drop_outliers: int = 10
    """Drop outlier bins more than this many multiples of the 95th quantile away from the average within a rolling window. Set to 0 for no outlier filtering."""
    smooth_cbs: bool = False

    @model_validator(mode="after")
    def ensure_smooth_for_cbs_only(self) -> Self:
        if self.smooth_cbs and self.method != SegmentationMethod.CBS:
            raise ValueError("'smooth_cbs' option can be used only with 'CBS' segmentation method")
        return self


class Call(SnappyModel):
    method: CallingMethod = CallingMethod.THRESHOLD
    """Calling method."""
    thresholds: list[float] = [-1.1, -0.25, 0.2, 0.7]
    """Hard thresholds for calling each integer copy number"""
    center: CenterMethod | None = None
    """Re-center the log2 ratio values using this estimator of the center or average value. ('median' if no argument given.)"""
    center_at: float | None = None
    """
    Subtract a constant number from all log2 ratios. For "manual" re-centering.

    When this parameter is set, the centering method should be left empty.
    """
    filter: FilterMethod | None = None
    """Merge segments flagged by the specified filter(s) with the adjacent segment(s)."""

    @model_validator(mode="after")
    def ensure_center_without_center_at(self) -> Self:
        if self.center_at is not None and self.center is not None:
            raise ValueError("'center' and 'center_at' parameters cannot be used together")
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
    """File listing the chromosomal ranges to display, as BED, interval list or 'chr:start-end' text (currently not implemented)"""
    chromosome: str | None = None
    """Name of the chromosome to display (whole genome if empty)"""
    gene: str | None = None
    """Name of gene or genes (comma-separated) to display (currently not implemented)"""
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


class ReportStats(enum.StrEnum):
    MEAN = "mean"
    MEDIAN = "median"
    MODE = "mode"
    T_TEST = "t-test"
    STDEV = "stdev"
    SEM = "sem"
    MAD = "mad"
    MSE = "mse"
    IQR = "iqr"
    BIVAR = "bivar"
    CI = "ci"
    PI = "pi"


class ReportSegmetrics(Report):
    alpha: float = 0.05
    """Level to estimate confidence and prediction intervals; use with --ci and --pi."""
    bootstrap: int = 100
    """Number of bootstrap iterations to estimate confidence interval; use with --ci."""
    smooth_bootstrap: bool = False
    """Apply Gaussian noise to bootstrap samples, a.k.a. smoothed bootstrap, to estimate confidence interval"""
    stats: list[ReportStats] = [
        ReportStats.MEAN,
        ReportStats.MEDIAN,
        ReportStats.MODE,
        ReportStats.T_TEST,
        ReportStats.STDEV,
        ReportStats.SEM,
        ReportStats.MAD,
        ReportStats.MSE,
        ReportStats.IQR,
        ReportStats.BIVAR,
        ReportStats.CI,
        ReportStats.PI,
    ]


class ReportGenemetrics(Report):
    alpha: float = 0.05
    """Level to estimate confidence and prediction intervals; use with --ci and --pi."""
    bootstrap: int = 100
    """Number of bootstrap iterations to estimate confidence interval; use with --ci."""
    threshold: float = 0.2
    """Copy number change threshold to report a gene gain/loss"""
    min_probes: int = 3
    """Minimum number of covered probes to report a gain/loss"""
    stats: list[ReportStats] = [
        ReportStats.MEAN,
        ReportStats.MEDIAN,
        ReportStats.MODE,
        ReportStats.T_TEST,
        ReportStats.STDEV,
        ReportStats.SEM,
        ReportStats.MAD,
        ReportStats.MSE,
        ReportStats.IQR,
        ReportStats.BIVAR,
        ReportStats.CI,
        ReportStats.PI,
    ]


class CnvkitToReference(SnappyModel):
    # Substep-secific parameters
    access: Access = Access()
    target: Target = Target()
    antitarget: Antitarget = Antitarget()

    coverage: Coverage = Coverage()

    metrics: Report = Report()
    segmetrics: ReportSegmetrics = ReportSegmetrics()
    genemetrics: ReportGenemetrics = ReportGenemetrics()

    # Generic parameters (used in different substeps & must agree)
    male_reference: bool = False
    """Create/use male reference (for shifting chrX & chrY)"""
    diploid_parx_genome: str | None = None
    """Considers the given human genome's PAR of chromosome X as autosomal. Example: 'grch38'"""

    cluster: bool = False
    """Calculate and store summary stats for clustered subsets of the normal samples with similar coverage profiles."""
    min_cluster_size: int = 4
    """Minimum cluster size to keep in reference profiles."""

    gc: bool = True
    """Skip GC correction."""
    edge: bool | None = None
    """Skip edge correction. Automatic selection when None (True for WGS & Panel, False for WES)"""
    rmask: bool = True
    """Skip RepeatMasker correction."""

    drop_low_coverage: bool = False
    """Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples."""


class Cnvkit(CnvkitToReference):
    fix: Fix = Fix()
    segment: Segment
    call: Call = Call()
    bintest: Bintest = Bintest()

    diagram: PlotDiagram = PlotDiagram()
    scatter: PlotScatter = PlotScatter()
