import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel


class SegmentationMethod(enum.StrEnum):
    cbs = "cbs"
    flasso = "flasso"
    haar = "haar"
    hmm = "hmm"
    hmm_tumor = "hmm-tumor"
    hmm_germline = "hmm-germline"
    none = "none"


class CenterMode(enum.StrEnum):
    mean = "mean"
    median = "median"
    mode = "mode"
    biweight = "biweight"


class FilterMode(enum.StrEnum):
    ampdel = "ampdel"
    cn = "cn"
    ci = "ci"
    sem = "sem"


class CallingMethod(enum.StrEnum):
    threshold = "threshold"
    clonal = "clonal"
    none = ""


class Gender(enum.StrEnum):
    male = "male"
    female = "female"
    guess = ""


class Coverage(SnappyModel):
    min_mapq: int = 0
    """[coverage] Mininum mapping quality score to count a read for coverage depth"""

    count: bool = False
    """[coverage] Alternative counting algorithm"""


class Fix(SnappyModel):
    gc_correction: bool = True
    """[fix] Use GC correction"""

    edge_correction: bool = True
    """[fix] Use edge correction"""

    rmask_correction: bool = True
    """[fix] Use rmask correction"""


class Segment(SnappyModel):
    # BCBIO uses
    # seg_method: haar
    # seg_threshold: 0.0001
    # -- OR
    # seg_method: cbs
    # seg_threshold: 0.000001
    method: SegmentationMethod = SegmentationMethod.cbs
    """[segment] One of cbs, flasso, haar, hmm, hmm-tumor, hmm-germline, none"""

    threshold: float = 0.000001
    """[segment] Significance threshold (hmm methods: smoothing window size)"""

    drop_outliers: int = 10
    """[segment] Drop outlier bins (0 for no outlier filtering)"""

    smooth_cbs: bool = True
    """[segment] Additional smoothing of CBS segmentation (WARNING- not the default value)"""


class Call(SnappyModel):
    center: CenterMode | float = CenterMode.median
    """[call] Either one of mean, median, mode, biweight, or a constant log2 ratio value."""

    filter: FilterMode | str = FilterMode.ampdel
    """
    [call] One of ampdel, cn, ci, sem (merging segments flagged with the specified filter),
    "" for no filtering
    """

    method: CallingMethod = CallingMethod.threshold
    """[call] One of threshold, clonal, none"""

    thresholds: str = "-1.1,-0.25,0.2,0.7"
    """[call] Thresholds for calling integer copy number"""

    ploidy: int = 2
    """[call] Ploidy of sample cells"""

    purity: Annotated[float, Field(0, ge=0, le=1)]
    """[call] Estimated tumor cell fraction (0 for discarding tumor cell purity)"""


class Diagram(SnappyModel):
    threshold: float = 0.5
    """[diagram] Copy number change threshold to label genes"""

    min_probes: int = 3
    """[diagram] Min number of covered probes to label genes"""

    shift_xy: bool = True
    """[diagram] Shift X & Y chromosomes according to sample sex"""


class Genemetrics(SnappyModel):
    min_probes: int = 3
    """[genemetrics] Min number of covered probes to consider a gene"""

    threshold: float = 0.2
    """[genemetrics] Min abs log2 change to consider a gene"""

    alpha: float = 0.05
    """[genemetrics] Significance cutoff"""

    bootstrap: int = 100
    """[genemetrics] Number of bootstraps"""


class Segmetrics(SnappyModel):
    alpha: float = 0.05
    """[segmetrics] Significance cutoff"""

    bootstrap: int = 100
    """[segmetrics] Number of bootstraps"""

    smooth_bootstrap: bool = False
    """[segmetrics] Smooth bootstrap results"""


class Breaks(SnappyModel):
    min_probes: int = 1
    """[breaks] Min number of covered probes for a break inside the gene"""


class Target(SnappyModel):
    split: bool = False
    """Split large tiled intervals into smaller, consecutive targets."""
    avg_size: int | None = None
    """ Average size of split target bins (results are approximate). Default 266.6666667 is not allowed (ints only)."""


class Antitarget(SnappyModel):
    min_size: float | None = None
    """Minimum size of antitarget bins (smaller regions are dropped). Default: 1/16 avg size when not set"""
    avg_size: int = 150000
    """ Average size of split target bins (results are approximate). Default 150000."""


class Reference(Fix):
    min_cluster_size: int = 0
    """Minimum cluster size to keep in reference profiles. When 0, don't use cluster, default value = 4"""


class Autobin(SnappyModel):
    bp_per_bin: float = 100000.0
    """Desired average number of sequencing read bases mapped to each bin."""
    target_min_size: int = 20
    """Minimum size of target bins"""
    target_max_size: int = 20000
    """Maximum size of target bins"""
    antitarget_min_size: int = 500
    """Minimum size of antitarget bins"""
    antitarget_max_size: int = 500000
    """Maximum size of antitarget bins"""


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

    coverage: Coverage = Coverage()
    fix: Fix = Fix()
    segment: Segment = Segment()
    call: Call = Call()
    diagram: Diagram = Diagram()
    genemetrics: Genemetrics = Genemetrics()
    segmetrics: Segmetrics = Segmetrics()
    breaks: Breaks = Breaks()

    drop_low_coverage: bool = False
    """[segment, call, genemetrics] Drop very low coverage bins"""

    gender: Gender = Gender.guess
    """
    [call, diagram] Specify the chromosomal sex of all given samples as male or female.
    Guess when missing
    """

    male_reference: bool = False
    """[call, diagram] Create male reference"""


class PanelOfNormals(SnappyModel):
    path_target: Annotated[
        str, Field(examples=["../panel_of_normals/output/cnvkit.target/out/cnvkit.target.bed"])
    ]
    """Path to target regions"""

    path_access: Annotated[
        str | None,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/cnvkit/access-10kb.hg38.bed"
            ]
        ),
    ] = None
    """Path to accessible regions"""

    path_annotation: Annotated[
        str | None,
        Field(
            examples=[
                "/data/cephfs-1/work/projects/cubit/20.05/static_data/annotation/GENCODE/33/GRCh38/gencode.v33.annotation.gtf"
            ]
        ),
    ] = None
    """Path to accessible regions"""

    coverage: Coverage = Coverage()
    target: Target = Target()
    antitarget: Antitarget = Antitarget()
    reference: Reference = Reference()
    genemetrics: Genemetrics = Genemetrics()
    segmetrics: Segmetrics = Segmetrics()
    breaks: Breaks = Breaks()

    drop_low_coverage: bool = False
    """[genemetrics] Drop very low coverage bins"""

    gender: Gender = Gender.guess
    """
    [call, diagram] Specify the chromosomal sex of all given samples as male or female.
    Guess when missing
    """

    male_reference: bool = False
    """[call, diagram] Create male reference"""

    # Autobin is used is for targets in WGS mode
    bp_per_bin: float = 100000.0
    """Desired average number of sequencing read bases mapped to each bin."""
