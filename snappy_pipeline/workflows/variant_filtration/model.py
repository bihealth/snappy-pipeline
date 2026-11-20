"""Pydantic models for the unified ``variant_filtration`` step."""

from typing import Annotated, Literal, Self

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import ConfigPath, RequiredBy, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
    Regions,
)
from snappy_pipeline.models.calls_report import CallsReport
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class ExpectedVariantVcf(BaseModel):
    vcf: str
    vcf_tbi: str


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
    """The events of the scenario to control the FDR for, jointly, such as ``somatic_tumor_high`` and ``somatic_tumor_medium`` of the somatic preset, or ``present``, ``denovo``, ``recessive`` and ``dominant`` of the germline one; a kind such as ``denovo`` stands for every per-member event ``denovo_<member>`` of the pedigree's scenario, and a kind the scenario lacks is skipped"""

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


ToolLiteral = Literal[
    "bcftools",
    "vembrane",
    "regions",
    "dkfz",
    "ebfilter",
    "varlociraptor",
]


class VariantFiltrationDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ]

    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
        RequiredBy("dkfz", "ebfilter"),
        RequiredBy(when="calls_report.igv"),
    ] = ""

    reference: Annotated[
        Reference, RequiredBy("dkfz", "ebfilter"), RequiredBy(when="calls_report.igv")
    ] = ""

    regions: Annotated[Regions, RequiredBy("regions")] = ""
    """Task with the BED file of the ``regions`` filter"""


class VariantFiltration(SnappyStepModel):
    depends_on: VariantFiltrationDependsOn

    calls_report: CallsReport | None = None
    """The calls of each VCF as a table and a datavzrd site for the report; unset: no report"""

    tool: ToolLiteral
    """One task = one tool. Vembrane supports both ``tag`` and ``filter`` modes."""

    bcftools: Bcftools | None = None
    vembrane: Vembrane | None = None
    regions: RegionsFilter | None = None
    dkfz: Dkfz = Field(default_factory=Dkfz)
    ebfilter: Ebfilter | None = None
    varlociraptor: VarlociraptorFdr | None = None

    @model_validator(mode="after")
    def validate_config(self) -> Self:
        if not self.depends_on.variants:
            raise ValueError("depends_on.variants must be set")

        if self.tool != "dkfz":
            tool_cfg = getattr(self, self.tool, None)
            if tool_cfg is None:
                raise ValueError(
                    f"Configuration block '{self.tool}:' must be provided when tool is '{self.tool}'"
                )

        return self
