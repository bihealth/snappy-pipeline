import enum
import typing
from typing import Annotated

from pydantic import Field, model_validator  #  , validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, Parallel


class WgsCaller(enum.StrEnum):
    CNVKIT = "cnvkit"
    CONTROL_FREEC = "control_freec"


class WesCaller(enum.StrEnum):
    CNVKIT = "cnvkit"
    PURECN = "purecn"
    SEQUENZA = "sequenza"


class Tools(SnappyModel):
    wgs: Annotated[typing.List[WgsCaller], EnumField(WgsCaller, [])]
    """WGS calling tools"""

    wes: Annotated[typing.List[WesCaller], EnumField(WesCaller, [])]
    """WES calling tools"""


class Sex(enum.StrEnum):
    SAMPLESHEET = "samplesheet"
    """Obtain the sex from the samplesheet"""
    DIPLOID_ONLY = "diploid_only"
    """Compute CNV for diploid chromosomes only"""
    AUTO = "auto"
    """Automatic sex detection using X/Y coverage"""
    FEMALE = "female"
    """Assume all samples are female"""
    MALE = "male"
    """Assume all samples are male"""
    UNKNOWN = "unknown"
    """Sex is unknown"""


class SequencingMethod(enum.StrEnum):
    WES = "hybrid"
    PANEL = "amplicon"
    WGS = "wgs"


class LibraryKitDefinition(SnappyModel):
    """
    Mapping from enrichment kit to target region BED file, for either computing per--target
    region coverage or selecting targeted exons.

    The following will match both the stock IDT library kit and the ones
    with spike-ins seen fromr Yale genomics.  The path above would be
    mapped to the name "default".
      - name: IDT_xGen_V1_0
        pattern: "xGen Exome Research Panel V1\\.0*"
        path: "path/to/targets.bed"
    """

    name: Annotated[str, Field(examples=["IDT_xGen_V1_0"])]

    pattern: Annotated[str, Field(examples=["xGen Exome Research Panel V1\\.0*"])]

    path: Annotated[str, Field(examples=["path/to/targets.bed"])]


class PanelOfNormalsOrigin(enum.StrEnum):
    PREVIOUS_STEP = "previous_step"
    """Use (& possibly create) a panel of normals from the current cohort of normals in the panel_of_normals step"""
    STATIC = "static"
    """Use an panel of normals from another cohort or from public data"""


class PanelOfNormals(SnappyModel):
    enabled: bool = False
    origin: PanelOfNormalsOrigin = PanelOfNormalsOrigin.PREVIOUS_STEP
    path_panel_of_normals: str = "../panel_of_normals"
    """
    Path to panel of normals created in current project

    The panel of normals can be either a file (typically from another project),
    or from the current project's panel_of_normals step.

    In the latter case, the missing one(s) (in case there are more than one panel, or if there are WES & WGS)
    will be created when not present.
    The matching of genome release & exome baits is done on genome name & exome baits md5 checksum.
    These are computed in the panel of normals step, and saved with the panel itself.

    There is no such matching if a panel of normal file is provided. The panel of normals validity is left to the user.
    """


class Mutect2(Parallel):
    panel_of_normals: PanelOfNormals | None = None
    """
    Panel of normals created by the PanelOfNormals program.
    """

    germline_resource: str

    common_variants: str | None = ""
    """Common germline variants for contamination estimation"""

    arguments_for_purecn: bool = True
    """
    PureCN requires that Mutect2 be called with arguments:
    --genotype-germline-sites true --genotype-pon-sites true
    """

    extra_arguments: Annotated[
        typing.List[str],
        # AfterValidator(argument),
        Field(
            examples=[
                "--read-filter CigarContainsNoNOperator",
                "--annotation AssemblyComplexity BaseQuality",
            ]
        ),
    ] = []
    """
    List additional Mutect2 arguments.
    Each additional argument must be of the form:
    "--<argument name> <argument value>"
    For example, to filter reads prior to calling & to add annotations to the output vcf:
      - "--read-filter CigarContainsNoNOperator"
      - "--annotation AssemblyComplexity BaseQuality"
    """

    window_length: int = 300000000


class VariantTool(enum.StrEnum):
    MUTECT2 = "mutect2"


class Variant(SnappyModel):
    enabled: bool = False
    tool: VariantTool | None = None

    mutect2: Mutect2 | None = None


class Ascat(SnappyModel):
    pass
    """TODO: configure purity tools (except for PureCN)"""


class Sequenza(SnappyModel):
    pass


class ControlFreec(SnappyModel):
    pass


class PureCn(SnappyModel):
    panel_of_normals: PanelOfNormals
    """
    Panel of normals created by the NormalDB.R script.
    This is required even if the normal/tumor paired mode won't use it.
    """

    variants: VariantTool

    mappability: str = ""
    """
    GRCh38:
     /fast/work/groups/cubi/projects/biotools/static_data/app_support/PureCN/hg38/mappability.bw
    """

    reptiming: str = ""
    """Nothing for GRCh38"""

    seed: int = 1234567
    extra_commands: typing.Dict[str, typing.Any] = {
        "model": "betabin",
        "fun-segmentation": "PSCBS",
        "post-optimize": "",
    }
    """Recommended extra arguments for PureCN, extra_commands: {} to clear them all"""

    path_container: Annotated[
        str, Field(examples=["../panel_of_normals/work/containers/out/purecn.simg"])
    ]
    """Conda installation not working well, container is required"""

    path_intervals: Annotated[
        str,
        Field(
            examples=[
                "../panel_of_normals/output/purecn/out/<enrichement_kit_name>_<genome_name>.list"
            ]
        ),
    ]


class PurityTool(enum.StrEnum):
    ASCAT = "ascat"
    PURECN = "purecn"


class Purity(SnappyModel):
    enabled: bool = False

    ignore_samplesheet: bool = False
    """Discard purity values in samplesheet when they exist"""
    default_value: float | None = None
    """Purity value for all samples"""

    tool: PurityTool | None = None
    """Tool used for purity estimation, if not set, try samplesheet, otherwise default_value"""

    ascat: Ascat | None = None


class CnvkitSegmentationMethod(enum.StrEnum):
    CBS = "cbs"
    FLASSO = "flasso"
    HAAR = "haar"
    HMM = "hmm"
    HMM_TUMOR = "hmm-tumor"
    HMM_GERMLINE = "hmm-germline"
    NONE = "none"


class CnvkitCallingMethod(enum.StrEnum):
    THRESHOLD = "threshold"
    CLONAL = "clonal"
    NONE = "none"


class CnvkitCenterMethod(enum.StrEnum):
    MEAN = "mean"
    MEDIAN = "median"
    MODE = "mode"
    BIWEIGHT = "biweight"


class CnvkitFilterMethod(enum.StrEnum):
    AMPDEL = "ampdel"
    CN = "cn"
    CI = "ci"
    SEM = "sem"


class CnvkitAccess(SnappyModel):
    exclude: Annotated[
        str | None,
        Field(
            examples=[
                "/fast/work/groups/cubi/projects/biotools/static_data/app_support/cnvkit/access-5k-mappable.grch37.bed"
            ]
        ),
    ] = None
    """Regions accessible to mapping"""

    min_gap_size: int = 5000
    """Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together."""


class CnvkitTarget(SnappyModel):
    split: bool = False
    """Split large tiled intervals into smaller, consecutive targets."""
    avg_size: float = 800 / 3
    """Average size of split target bins (results are approximate)"""


class CnvkitAntitarget(SnappyModel):
    avg_size: float = 150000
    """Average size of split antitarget bins (results are approximate)"""
    min_size: float | None = None
    """Minimum size of antitarget bins (smaller regions are dropped). When missing, 1/16 avg size"""


class CnvkitCoverage(SnappyModel):
    count: bool = False
    """Get read depths by counting read midpoints within each bin."""
    min_mapq: int = 0
    """Minimum mapping quality score (phred scale 0-60) to count a read for coverage depth."""


class CnvkitReference(SnappyModel):
    cluster: bool = False
    """Calculate and store summary stats for clustered subsets of the normal samples with similar coverage profiles."""
    min_cluster_size: int = 4
    """Minimum cluster size to keep in reference profiles."""
    no_gc: bool = False
    """Skip GC correction."""
    no_edge: bool = None
    """Skip edge correction. Automatic selection when None (True for WGS & Panel, False for WES)"""
    no_rmask: bool = False
    """Skip RepeatMasker correction."""


class CnvkitFix(SnappyModel):
    cluster: bool = False
    """Compare and use cluster-specific values present in the reference profile."""
    no_gc: bool = False
    """Skip GC correction."""
    no_edge: bool = False
    """Skip edge correction."""
    no_rmask: bool = False
    """Skip RepeatMasker correction."""


class CnvkitSegment(SnappyModel):
    method: CnvkitSegmentationMethod = CnvkitSegmentationMethod.CBS
    """Segmentation method, or 'NONE' for chromosome arm-level averages as segments"""
    threshold: float = 0.0001
    """Significance threshold (p-value or FDR, depending on method) to accept breakpoints during segmentation. For HMM methods, this is the smoothing window size."""
    drop_low_coverage: bool = False
    """Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples."""
    drop_outliers: float = 10
    """Drop outlier bins more than this many multiples of the 95th quantile away from the average within a rolling window. Set to 0 for no outlier filtering."""
    smooth_cbs: bool = False

    min_variant_depth: float = 20
    """Minimum read depth for a SNV to be displayed in the b-allele frequency plot."""
    zygocity_freq: float = 0.25
    """Ignore VCF's genotypes (GT field) and instead infer zygosity from allele frequencies."""

    @model_validator(mode="after")
    def ensure_smooth_for_cbs_only(self) -> typing.Self:
        if self.smooth_cbs and self.method != CnvkitSegmentationMethod.CBS:
            raise ValueError("'smooth_cbs' option can be used only with 'CBS' segmentation method")
        return self


class CnvkitCall(SnappyModel):
    method: CnvkitCallingMethod = CnvkitCallingMethod.THRESHOLD
    """Calling method."""
    thresholds: str | None = None
    """Hard thresholds for calling each integer copy number, separated by commas"""
    center: CnvkitCenterMethod | None = CnvkitCenterMethod.MEDIAN
    """Re-center the log2 ratio values using this estimator of the center or average value. ('median' if no argument given.)"""
    center_at: float | None = None
    """Subtract a constant number from all log2 ratios. For "manual" re-centering."""
    filter: CnvkitFilterMethod | None = None
    """Merge segments flagged by the specified filter(s) with the adjacent segment(s)."""
    ploidy: float | None = 2
    """Ploidy of the sample cells."""
    drop_low_coverage: bool = False
    """Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples."""

    min_variant_depth: float = 20
    """Minimum read depth for a SNV to be displayed in the b-allele frequency calculation."""
    zygocity_freq: float = 0.25
    """Ignore VCF's genotypes (GT field) and instead infer zygosity from allele frequencies."""


class CnvkitBintest(SnappyModel):
    alpha: float = 0.005
    """Significance threhold."""
    target: bool = False
    """Test target bins only; ignore off-target bins."""


class CnvkitPlotDiagram(SnappyModel):
    threshold: float = 0.5
    """Copy number change threshold to label genes."""
    min_probes: int = 3
    """Minimum number of covered probes to label a gene."""
    no_shift_xy: bool = False


class CnvkitPlotScatter(SnappyModel):
    antitarget_marker: str | None = None
    """Plot antitargets using this symbol when plotting in a selected chromosomal region."""
    by_bin: bool = False
    """Plot data x-coordinates by bin indices instead of genomic coordinates."""
    segment_color: str | None = None
    """Plot segment lines in this color. Value can be any string accepted by matplotlib."""
    trend: bool = False
    """Draw a smoothed local trendline on the scatter plot."""
    y_max: float | None = None
    """y-axis upper limit."""
    y_min: float | None = None
    """y-axis lower limit."""
    fig_size: typing.Tuple[float, float] | None = None
    """Width and height of the plot in inches."""

    min_variant_depth: float = 20
    """Minimum read depth for a SNV to be displayed in the b-allele frequency calculation."""
    zygocity_freq: float = 0.25
    """Ignore VCF's genotypes (GT field) and instead infer zygosity from allele frequencies."""


class CnvkitPlot(SnappyModel):
    diagram: CnvkitPlotDiagram = CnvkitPlotDiagram()
    scatter: CnvkitPlotScatter = CnvkitPlotScatter()


class CnvkitReportMetrics(SnappyModel):
    drop_low_coverage: bool = False
    """Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples."""


class CnvkitReportSegmetrics(SnappyModel):
    drop_low_coverage: bool = False
    """Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples."""
    alpha: float = 0.05
    """Level to estimate confidence and prediction intervals; use with --ci and --pi."""
    bootstrap: int = 100
    """Number of bootstrap iterations to estimate confidence interval; use with --ci."""


class CnvkitReport(enum.StrEnum):
    METRICS = "metrics"
    SEGMETRICS = "segmetrics"


class Cnvkit(SnappyModel):
    panel_of_normals: PanelOfNormals | None = None

    variants: VariantTool | None = None

    purity: Purity
    """
    When present, purity estimates can be used for calling segments. The requested tool must be configured.
    Or the purity can be provided in the samplesheet, as an extra information attached to the library.

    Note that PureCN cannot be used to estimate purity for WGS samples (because PureCN is WES & Panel-only).
    TODO: This should be tested by a validation method, I don't know how to do (Till help!!)
    TODO: The exact name is not yet set.
    """

    access: CnvkitAccess = CnvkitAccess()
    target: CnvkitTarget = CnvkitTarget()
    antitarget: CnvkitAntitarget = CnvkitAntitarget()
    coverage: CnvkitCoverage = CnvkitCoverage()

    reference: CnvkitReference | None = None

    @model_validator(mode="after")
    def set_default_reference(self) -> typing.Self:
        if self.reference is None and not self.panel_of_normals.enabled:
            self.reference = CnvkitReference()
        return self

    fix: CnvkitFix = CnvkitFix()
    segment: CnvkitSegment = CnvkitSegment()
    call: CnvkitCall = CnvkitCall()
    bintest: CnvkitBintest = CnvkitBintest()

    use_male_reference: bool = False
    """Create/use a male reference. Must be identical to panel of normals creation, when using one"""

    plots: typing.List[CnvkitPlot] = []

    reports: typing.List[CnvkitReport] = []
    metrics: CnvkitReportMetrics | None = None

    # @validator("metrics")
    # def get_default_reference(cls, v, values) -> CnvkitReportMetrics | None:
    #     if v is None and "metrics" in values["reports"]:
    #         return CnvkitReportMetrics()
    #     return None

    segmetrics: CnvkitReportSegmetrics | None = None

    # @validator("segmetrics")
    # def get_default_reference(cls, v, values) -> CnvkitReportSegmetrics | None:
    #     if v is None and "segmetrics" in values["reports"]:
    #         return CnvkitReportSegmetrics()
    #     return None


class SomaticCnvCalling(SnappyStepModel):
    path_ngs_mapping: str
    """Path to bam files"""

    tools: Tools
    """Tools for WGS & WES data"""

    path_target_interval_list_mapping: typing.List[LibraryKitDefinition] | None = None

    sex: Sex = Sex.DIPLOID_ONLY

    cnvkit: Cnvkit
    purecn: PureCn | None = None
    sequenza: Sequenza | None = None
    control_freec: ControlFreec | None = None

    mutect2: Mutect2 | None = None

    default_ploidy: float | None = None

    # @model_validator(mode="after")
    # def ensure_single_pon_step(self) -> typing.Self:
    #     """
    #     I am not sure this is absolutely required.
    #     I am trying to avoid registering the panel_of_normals step when initializing SomaticCnvCalling
    #     """
    #     pon_steps = set()
    #     for tool in itertools.chain(self.tools.wgs, self.tools.wes):
    #         tool_config = getattr(self, tool)
    #         if (
    #             tool_config
    #             and getattr(tool_config, "use_panel_of_normals")
    #             and tool_config.use_panel_of_normals == PanelOfNormalsUse.PREVIOUS_STEP
    #         ):
    #             pon_steps.add(str(tool_config.panel_of_normals.panel_of_normals))
    #     if len(pon_steps) > 1:
    #         raise ValueError("Too many panel_of_normals steps")
    #     return self
