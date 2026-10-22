import enum
from typing import Annotated, Literal

from pydantic import model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel, validators
from snappy_pipeline.models.cnvkit import Cnvkit, Diagram, Fix
from snappy_pipeline.workflows.abstract.protocol import (
    TOOL,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    cnvetti = "cnvetti"
    control_freec = "control_freec"
    cnvkit = "cnvkit"


class CnvettiPreset(SnappyModel):
    window_length: int
    count_kind: str
    segmentation: str
    normalization: str


class Cnvetti(ToolModel):
    window_length: int | None = None
    count_kind: str | None = None
    segmentation: str | None = None
    normalization: str | None = None
    presets: dict[str, CnvettiPreset] = {
        "deep_wgs": CnvettiPreset(
            **{
                "window_length": 200,
                "count_kind": "Coverage",
                "segmentation": "HaarSeg",
                "normalization": "MedianGcBinned",
            }
        )
    }
    preset: str = "deep_wgs"


class ControlFreecConvert(SnappyModel):
    org_obj: str = "org.Hs.eg.db::org.Hs.eg.db"

    tx_obj: str = "TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene"

    bs_obj: str = "BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5"

    amplification: int = 6
    """Smallest copy number that is called an amplification (GISTIC value 2); gains have 3 up to this number minus 1"""


class ExpectedControlFreecDatabase(SnappyModel):
    """Consumer-driven contract: the files of Control-FREEC."""

    chr_len: str
    """ChrLenFile"""

    mappability: str = ""
    """Mappability file, used when set"""


class ControlFreec(ToolModel):
    window_size: int = -1
    """set to a value >=0 you want a specific fixed window size"""

    ploidy: int = 2
    """Ploidy of the sample (``ploidy``)"""

    break_point_threshold: float = 0.8
    """Segmentation sensitivity, lower values give more breakpoints (``breakPointThreshold``)"""

    coefficient_of_variation: float = 0.05
    """Coefficient of variation to select the window size (``coefficientOfVariation``)"""

    contamination_adjustment: bool = True
    """Correct the copy numbers for normal cell contamination (``contaminationAdjustment``)"""

    contamination: float = 0.4
    """Fraction of normal cells in the tumor sample, used with ``contamination_adjustment`` (``contamination``)"""

    minimal_subclone_presence: float = 0.2
    """Minimal fraction of cells for a subclone to be reported (``minimalSubclonePresence``)"""

    unique_match: bool = True
    """Use only uniquely mapped reads when a mappability file is given (``uniqueMatch``)"""

    force_gc_content_normalization: int = 0
    """GC content normalization mode, 0 uses the polynomial fit only when it improves the fit (``forceGCcontentNormalization``)"""

    min_cna_length: int = 2
    """Minimal number of consecutive windows for a copy number alteration (``minCNAlength``)"""

    unknown_sex: Literal["XX", "XY"] = "XX"
    """Control-FREEC's ``sex`` for donors whose sex the sample sheet does not give; male donors get XY, female ones XX"""

    convert: ControlFreecConvert

    diagram: Diagram = Diagram()
    """Settings of the cnvkit diagram plot of the calls"""


# If defaults need to be overwritten, subclass the model and override the defaults
class CnvkitWgs(Cnvkit):
    fix: Fix = Fix(edge_correction=False)
    """Whole genomes have no baits, so there is no edge correction, as in ``cnvkit.py batch
    --method wgs`` and the whole genome panel of normals"""


class ExpectedCnvkitPanelOfNormals(SnappyModel):
    """Consumer-driven contract: the files of a ``cnvkit`` panel of normals task."""

    panel_of_normals: str
    """The reference (``.cnn``), from the normals or flat"""

    target: str
    """The target regions the panel was built with"""

    antitarget: str
    """The antitarget regions the panel was built with (empty for whole genomes)"""


class SomaticWgsCnvCallingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference = ""

    panel_of_normals: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"pon", "cnvkit"})),
        ExpectedPathSchema(ExpectedCnvkitPanelOfNormals),
    ] = ""
    """Needed for ``tool: cnvkit``: a ``panel_of_normals`` task with ``tool: cnvkit`` in whole
    genome mode (without ``depends_on.kits``), which also gives the target and antitarget
    regions"""

    database: Annotated[
        str,
        DataSignature(DataType.DATABASE, frozenset({TOOL})),
        ExpectedPathSchema(ExpectedControlFreecDatabase),
    ] = ""
    """Needed for ``tool: control_freec``: the chromosome lengths and, optionally, the
    mappability"""


TOOL_DEPENDENCIES = {
    Tool.cnvetti: ("reference",),
    Tool.cnvkit: ("panel_of_normals", "reference"),
    Tool.control_freec: ("database",),
}


class SomaticWgsCnvCalling(SnappyStepModel):
    depends_on: SomaticWgsCnvCallingDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    cnvetti: Cnvetti | None = None

    control_freec: ControlFreec | None = None

    cnvkit: CnvkitWgs | None = None

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
