from __future__ import annotations

import enum
from typing import Annotated, Any, TypedDict

from pydantic import BaseModel, ConfigDict, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToggleModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Features,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedGeneCounts
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants


class ExpectedCopyNumberCalls(BaseModel):
    """Consumer-driven contract for copy-number provider outputs used by cBioPortal export."""

    dnacopy_seg: str


class VariantAnnotationTool(enum.StrEnum):
    VEP = "vep"
    MEHARI = "mehari"


class CopyNumberTool(enum.StrEnum):
    CNVKIT = "cnvkit"
    SEQUENZA = "sequenza"
    CONTROL_FREEC = "Control_FREEC"
    """unsupported"""


class NcbiBuild(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class Vcf2Maf(SnappyModel):
    Center: str
    ncbi_build: NcbiBuild

    annotation_tool: VariantAnnotationTool = VariantAnnotationTool.VEP
    """Which annotation tool was used on the input VCF.

    The vcf2maf wrapper uses this to select the correct config for parsing
    VCF annotation fields (e.g. VEP vs Mehari format).
    """


class GenomeName(enum.StrEnum):
    hg19 = "hg19"  # GRCh37
    hg38 = "hg38"  # GRCh38
    mm10 = "mm10"  # GRCm38


class Expression(ToggleModel):
    """When missing, no expression data is uploaded to cBioPortal. The gene counts come from the
    STAR task in ``depends_on.alignments``."""


class CNA(ToggleModel):
    """When missing, no CNV data uploaded to portal. Access WES & WGS steps"""

    copy_number_tool: CopyNumberTool = CopyNumberTool.CNVKIT


class Study(SnappyModel):
    type_of_cancer: str
    """
    see http://oncotree.mskcc.org/#/home
    see also `curl https://oncotree.mskcc.org:443/api/tumorTypes | jq ".[].code"`
    """
    cancer_study_id: str
    """Usually: <type of cancer id>_<pi>_<year>"""
    study_description: str
    study_name: str
    study_name_short: str
    reference_genome: GenomeName = GenomeName.hg38


class ExtraInfos(TypedDict):
    name: str
    description: str
    datatype: str
    priority: str
    column: str


class CbioportalExportDependsOn(SnappyModel):
    #: RNA mapping task whose gene counts the expression export reads
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedGeneCounts),
    ] = ""
    copy_number: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedCopyNumberCalls),
    ] = ""
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", ("snv", "indel")})),
        ExpectedPathSchema(ExpectedSomaticVariants),
    ]

    features: Features = ""


class CbioportalExport(SnappyStepModel):
    depends_on: CbioportalExportDependsOn

    model_config = ConfigDict(
        extra="forbid",
    )

    """Annotation is mandatory, but filtration is optional, can happen before or after annotation"""

    path_gene_id_mappings: str
    """Mapping from pipeline gene ids to cBioPortal ids (HGNC symbols from GeneNexus)"""

    exclude_variant_with_flag: str = ""
    """Required to filter variants (typically by type)"""

    vcf2maf: Vcf2Maf

    expression: Expression = Expression()
    """Include mRNA expression data"""

    copy_number_alteration: CNA = CNA()
    """Include copy number alteration results"""

    @model_validator(mode="after")
    def validate_optional_dependencies(self):
        if self.copy_number_alteration.enabled and not self.depends_on.copy_number:
            raise ValueError("copy_number_alteration needs depends_on.copy_number")
        if self.expression.enabled and not self.depends_on.alignments:
            raise ValueError("expression needs depends_on.alignments")
        if (self.copy_number_alteration.enabled or self.expression.enabled) and not (
            self.depends_on.features
        ):
            raise ValueError("copy_number_alteration and expression need depends_on.features")
        return self

    study: Study

    patient_info: dict[str, Any] = {}
    """unimplemented"""

    sample_info: dict[str, Any] = {}
    """Implementation must be re-designed"""
