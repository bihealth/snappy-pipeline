import enum
from typing import Annotated

from pydantic import Field

from models import SnappyModel


class SegmentationMethod(enum.Enum):
    cbs = "cbs"
    flasso = "flasso"
    haar = "haar"
    hmm = "hmm"
    hmm_tumor = "hmm-tumor"
    hmm_germline = "hmm-germline"
    none = "none"


class CenterMode(enum.Enum):
    mean = "mean"
    median = "median"
    mode = "mode"
    biweight = "biweight"


class FilterMode(enum.Enum):
    ampdel = "ampdel"
    cn = "cn"
    ci = "ci"
    sem = "sem"


class CallingMethod(enum.Enum):
    threshold = "threshold"
    clonal = "clonal"
    none = ""


class Gender(enum.Enum):
    male = "male"
    female = "female"
    guess = ""


class Cnvkit(SnappyModel):
    path_target: Annotated[
        str, Field(examples=["../panel_of_normals/output/cnvkit.target/out/cnvkit.target.bed"])
    ]
    """Path to target regions"""

    path_antitarget: Annotated[
        str,
        Field(examples=["../panel_of_normals/output/cnvkit.antitarget/out/cnvkit.antitarget.bed"]),
    ]
    """Path to antitarget regions"""

    path_panel_of_normals: Annotated[
        str,
        Field(
            examples=[
                "../panel_of_normals/output/{mapper}.cnvkit.create_panel/out/{mapper}.cnvkit.panel_of_normals.cnn"
            ]
        ),
    ]
    """Path to panel of normals (reference)"""

    plot: bool = True
    """Generate plots (very slow)"""

    min_mapq: int = 0
    """[coverage] Mininum mapping quality score to count a read for coverage depth"""

    count: bool = False
    """[coverage] Alternative counting algorithm"""

    gc_correction: bool = True
    """[fix] Use GC correction"""

    edge_correction: bool = True
    """[fix] Use edge correction"""

    rmask_correction: bool = True
    """[fix] Use rmask correction"""
    # BCBIO uses
    # seg_method: haar
    # seg_threshold: 0.0001
    # -- OR
    # seg_method: cbs
    # seg_threshold: 0.000001
    segmentation_method: SegmentationMethod = SegmentationMethod.cbs
    """[segment] One of cbs, flasso, haar, hmm, hmm-tumor, hmm-germline, none"""

    segmentation_threshold: float = 0.000001
    """[segment] Significance threshold (hmm methods: smoothing window size)"""

    drop_low_coverage: bool = False
    """[segment, call, genemetrics] Drop very low coverage bins"""

    drop_outliers: 10
    """[segment] Drop outlier bins (0 for no outlier filtering)"""

    smooth_cbs: bool = True
    """[segment] Additional smoothing of CBS segmentation (WARNING- not the default value)"""

    center: CenterMode | float | None = None
    """[call] Either one of mean, median, mode, biweight, or a constant log2 ratio value."""

    filter: FilterMode | str = FilterMode.ampdel
    """
    [call] One of ampdel, cn, ci, sem (merging segments flagged with the specified filter),
    "" for no filtering
    """

    calling_method: CallingMethod = CallingMethod.threshold
    """[call] One of threshold, clonal, none"""

    call_thresholds: str = "-1.1,-0.25,0.2,0.7"
    """[call] Thresholds for calling integer copy number"""

    ploidy: int = 2
    """[call] Ploidy of sample cells"""
    purity: Annotated[float, Field(0, ge=0, le=1)]
    """[call] Estimated tumor cell fraction (0 for discarding tumor cell purity)"""

    gender: Gender = Gender.guess
    """
    [call, diagram] Specify the chromosomal sex of all given samples as male or female.
    Guess when missing
    """

    male_reference: bool = False
    """[call, diagram] Create male reference"""
    diagram_threshold: float = 0.5
    """[diagram] Copy number change threshold to label genes"""

    diagram_min_probes: int = 3
    """[diagram] Min number of covered probes to label genes"""

    shift_xy: bool = True
    """[diagram] Shift X & Y chromosomes according to sample sex"""

    breaks_min_probes: int = 1
    """[breaks] Min number of covered probes for a break inside the gene"""

    genemetrics_min_probes: int = 3
    """[genemetrics] Min number of covered probes to consider a gene"""

    genemetrics_threshold: float = 0.2
    """[genemetrics] Min abs log2 change to consider a gene"""

    genemetrics_alpha: float = 0.05
    """[genemetrics] Significance cutoff"""

    genemetrics_bootstrap: int = 100
    """[genemetrics] Number of bootstraps"""

    segmetrics_alpha: float = 0.05
    """[segmetrics] Significance cutoff"""

    segmetrics_bootstrap: int = 100
    """[segmetrics] Number of bootstraps"""

    smooth_bootstrap: bool = False
    """[segmetrics] Smooth bootstrap results"""
