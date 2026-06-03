"""Pydantic models for the unified ``variant_filtration`` step."""

from typing import Annotated, Literal, Self

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class ExpectedVariantVcf(BaseModel):
    vcf: str
    vcf_tbi: str


class Bcftools(SnappyModel):
    include: str = ""
    exclude: str = ""

    @model_validator(mode="after")
    def ensure_exactly_one(self) -> Self:
        if not self.include and not self.exclude:
            raise ValueError("Either 'include' or 'exclude' must be set for bcftools")
        if self.include and self.exclude:
            raise ValueError("Only one of 'include' or 'exclude' may be set for bcftools")
        return self


class Regions(SnappyModel):
    include: str = ""
    exclude: str = ""
    path_bed: Annotated[str, Field(deprecated="Use 'exclude' instead")] = ""

    @model_validator(mode="after")
    def ensure_exactly_one(self) -> Self:
        n = sum(bool(v) for v in (self.include, self.exclude, self.path_bed))
        if n == 0:
            raise ValueError("One of 'include', 'exclude', or 'path_bed' must be set for regions")
        if n > 1:
            raise ValueError(
                "Only one of 'include', 'exclude', or 'path_bed' may be set for regions"
            )
        return self


class Vembrane(SnappyModel):
    expressions: dict[str, str] = Field(
        default_factory=dict,
        examples=[
            {
                "silent": 'ANN["Consequence"] == ["synonymous_variant"]',
                "poor_support": (
                    '(FORMAT["DP"][SAMPLES[1]] < 50) or '
                    '(FORMAT["AD"][SAMPLES[1]][1] < 5) or '
                    '(FORMAT["AD"][SAMPLES[1]][1] / '
                    '(FORMAT["AD"][SAMPLES[1]][0] + FORMAT["AD"][SAMPLES[1]][1]) < 0.05)'
                ),
            }
        ],
    )
    mode: Literal["tag", "filter"] = "tag"
    """``tag`` uses ``vembrane tag``; ``filter`` uses ``vembrane filter``."""

    expression: str = ""
    """Single expression used in ``filter`` mode."""

    aux: dict[str, str] = {}
    """Mapping for ``--aux NAME=PATH`` files."""

    context: list[str] = []
    """Python statements for ``--context``."""

    context_files: list[str] = []
    """Paths for ``--context-file`` scripts."""

    ontology: str = ""
    """Optional ontology file passed via ``--ontology``."""

    extra_args: str = ""

    @model_validator(mode="after")
    def ensure_mode_specific_fields(self) -> Self:
        if self.mode == "tag":
            if not self.expressions:
                raise ValueError("tag mode requires 'expressions'")
            if self.expression:
                raise ValueError("tag mode does not use 'expression'")
        else:
            if not self.expression:
                raise ValueError("filter mode requires 'expression'")
            if self.expressions:
                raise ValueError("filter mode does not use 'expressions'")
        return self


class Dkfz(SnappyModel):
    """DKFZ bias filter – no configurable parameters."""


class Ebfilter(SnappyModel):
    ebfilter_threshold: float = 2.4
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: int = 20
    min_baseq: int = 15
    path_panel_of_normals_sample_list: str = ""


class RemoveTags(SnappyModel):
    tags: list[str]
    """FILTER tags to remove corresponding records from the VCF."""

    backend: Literal["bcftools", "vembrane"] = "bcftools"
    """Backend used for removing records tagged in ``FILTER``."""

    @model_validator(mode="after")
    def ensure_tags(self) -> Self:
        cleaned = [t for t in self.tags if t]
        if not cleaned:
            raise ValueError("tags must contain at least one non-empty tag")
        self.tags = cleaned
        return self


ToolLiteral = Literal[
    "bcftools",
    "vembrane",
    "regions",
    "dkfz",
    "ebfilter",
]

_BAM_TOOLS: frozenset[str] = frozenset({"dkfz", "ebfilter"})


class VariantFiltrationDependsOn(SnappyModel):
    variant: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ]

    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = ""


class VariantFiltration(SnappyStepModel):
    depends_on: VariantFiltrationDependsOn = Field(default_factory=VariantFiltrationDependsOn)

    tool: ToolLiteral
    """One task = one tool. Vembrane supports both ``tag`` and ``filter`` modes."""

    bcftools: Bcftools | None = None
    vembrane: Vembrane | None = None
    regions: Regions | None = None
    dkfz: Dkfz = Field(default_factory=Dkfz)
    ebfilter: Ebfilter | None = None

    @model_validator(mode="after")
    def validate_config(self) -> Self:
        if not self.depends_on.variant:
            raise ValueError("depends_on.variant must be set")

        if self.tool != "dkfz":
            tool_cfg = getattr(self, self.tool, None)
            if tool_cfg is None:
                raise ValueError(
                    f"Configuration block '{self.tool}:' must be provided when tool is '{self.tool}'"
                )

        if self.tool in _BAM_TOOLS and not self.depends_on.ngs_mapping:
            raise ValueError(f"depends_on.ngs_mapping is required when tool is '{self.tool}'")

        return self
