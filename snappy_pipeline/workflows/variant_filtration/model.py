"""Pydantic models for the unified ``variant_filtration`` step."""

import enum
from typing import Annotated, Literal, Self

from pydantic import Field, model_validator

from snappy_pipeline.models import (
    ConfigPath,
    EnumField,
    RequiredBy,
    SnappyModel,
    SnappyStepModel,
    ToolModel,
)
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    ExpectedVariants,
    Reference,
    Regions,
)


class Tool(enum.StrEnum):
    bcftools = "bcftools"
    vembrane = "vembrane"
    regions = "regions"
    dkfz = "dkfz"
    ebfilter = "ebfilter"
    varlociraptor = "varlociraptor"


class Bcftools(ToolModel):
    include: str = ""
    exclude: str = ""
    mode: Literal["tag", "filter"] = "tag"

    @model_validator(mode="after")
    def ensure_exactly_one(self) -> Self:
        if not self.include and not self.exclude:
            raise ValueError("Either 'include' or 'exclude' must be set for bcftools")
        if self.include and self.exclude:
            raise ValueError("Only one of 'include' or 'exclude' may be set for bcftools")
        return self


class RegionsFilter(ToolModel):
    action: Literal["include", "exclude"]
    """``exclude`` applies the filter to the variants inside the regions of ``depends_on.regions``,
    ``include`` to the variants outside them."""

    mode: Literal["tag", "filter"] = "tag"


class Vembrane(ToolModel):
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

    aux: dict[str, ConfigPath] = {}
    """Mapping for ``--aux NAME=PATH`` files."""

    context: list[str] = []
    """Python statements for ``--context``."""

    context_files: list[ConfigPath] = []
    """Paths for ``--context-file`` scripts."""

    ontology: ConfigPath = ""
    """Optional ontology file passed via ``--ontology``."""

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


class Dkfz(ToolModel):
    """DKFZ bias filter."""

    mode: Literal["tag", "filter"] = "tag"


class Ebfilter(ToolModel):
    ebfilter_threshold: float = 2.4
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: int = 20
    min_baseq: int = 15
    skip_filters: list[str] = ["germline", "weak_evidence", "OffExome"]
    """Variants with one of these ``FILTER`` IDs (that the VCF header declares) skip EBFilter"""

    skip_annotation: str = "stream_gene_variant"
    """With annotated variants, those with this string in their ``CSQ`` or ``ANN`` field skip EBFilter
    (``bcftools`` ``~`` match); empty for none"""

    path_panel_of_normals_sample_list: ConfigPath = ""
    mode: Literal["tag", "filter"] = "tag"


class VarlociraptorFdr(ToolModel):
    """FDR control of varlociraptor calls (``varlociraptor filter-calls control-fdr``); the
    ``extra_arguments`` go to that command (``--smart-retain-artifacts``, ``--minlen``, ...)"""

    events: Annotated[list[str], Field(min_length=1)]
    """The events of the scenario to control the FDR for, jointly:
    names such as ``somatic_tumor_high`` and ``somatic_tumor_medium`` of the somatic preset or ``present``,
    ``recessive`` and ``dominant`` of the germline one, which must exist in the calls,
    and patterns such as ``denovo_*`` for the germline preset's per-member events ``denovo_<member>``,
    which may match none (a pedigree without a child has no de novo event)"""

    fdr: Annotated[float, Field(gt=0, le=1)] = 0.05
    """The false discovery rate (``--fdr``)"""

    fdr_mode: Literal["local-smart", "local-strict", "global-smart", "global-strict"] = (
        "local-smart"
    )
    """``local`` controls the FDR per record, ``global`` over the whole call set; ``smart``
    controls the probability of presence and keeps records whose events reach half of it,
    ``strict`` the summed probability of the events (``--mode``)"""

    mode: Literal["tag", "filter"] = "filter"
    """``filter`` keeps only the records that pass the FDR control, ``tag`` keeps all records
    and marks the others with the task name in ``FILTER``"""

    variant_types: list[Literal["SNV", "MNV", "INS", "DEL", "BND", "INV", "DUP", "REP"]] = [
        "SNV",
        "MNV",
        "INS",
        "DEL",
        "BND",
        "INV",
        "DUP",
        "REP",
    ]
    """The variant types kept, each controlled on its own (``--var``)"""


class VariantFiltrationDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariants),
    ]

    alignments: Annotated[DnaAlignments, RequiredBy(Tool.dkfz, Tool.ebfilter)] = ""

    reference: Annotated[Reference, RequiredBy(Tool.dkfz, Tool.ebfilter)] = ""

    regions: Annotated[Regions, RequiredBy(Tool.regions)] = ""
    """Task with the BED file of the ``regions`` filter"""


class VariantFiltration(SnappyStepModel):
    depends_on: VariantFiltrationDependsOn

    tool: Annotated[Tool, EnumField(Tool)]
    """One task = one tool. Vembrane supports both ``tag`` and ``filter`` modes."""

    bcftools: Bcftools | None = None
    vembrane: Vembrane | None = None
    regions: RegionsFilter | None = None
    dkfz: Dkfz = Field(default_factory=Dkfz)
    ebfilter: Ebfilter | None = None
    varlociraptor: VarlociraptorFdr | None = None
