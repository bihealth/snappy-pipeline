from __future__ import annotations

import enum
from typing import Annotated, Any, TypedDict

from pydantic import BaseModel, ConfigDict

from snappy_pipeline.models import ConfigPath
from snappy_pipeline.models import RequiredBy, SnappyModel, SnappyStepModel, ToggleModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    ExpectedVariants,
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Features,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedGeneCounts


class ExpectedCopyNumberCalls(BaseModel):
    """Consumer-driven contract for copy-number provider outputs used by cBioPortal export."""

    dnacopy_seg: str


class NcbiBuild(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class Vcf2Maf(ToolModel):
    Center: str
    ncbi_build: NcbiBuild


class GenomeName(enum.StrEnum):
    hg19 = "hg19"  # GRCh37
    hg38 = "hg38"  # GRCh38
    mm10 = "mm10"  # GRCm38


class Expression(ToggleModel):
    """When missing, no expression data is uploaded to cBioPortal. The gene counts come from the
    STAR task in ``depends_on.alignments``."""


class CNA(ToggleModel):
    """When missing, no CNV data uploaded to portal. Access WES & WGS steps"""


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
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedGeneCounts),
        RequiredBy(when="expression"),
    ] = ""
    """RNA mapping task whose gene counts the expression export reads"""

    copy_number: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedCopyNumberCalls),
        RequiredBy(when="copy_number_alteration"),
    ] = ""
    """Somatic CNV calling task whose segments the copy number export reads"""

    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", ("snv", "indel"), "vep"})),
        ExpectedPathSchema(ExpectedVariants),
    ]
    """Somatic variants annotated by VEP, filtered or not; the conversion to MAF cannot read
    mehari's annotations yet"""

    features: Annotated[
        Features, RequiredBy(when="copy_number_alteration"), RequiredBy(when="expression")
    ] = ""


class CbioportalExport(SnappyStepModel):
    depends_on: CbioportalExportDependsOn

    model_config = ConfigDict(
        extra="forbid",
    )

    """Annotation is mandatory, but filtration is optional, can happen before or after annotation"""

    path_gene_id_mappings: ConfigPath
    """Mapping from pipeline gene ids to cBioPortal ids (HGNC symbols from GeneNexus)"""

    exclude_variant_with_flag: str = ""
    """Required to filter variants (typically by type)"""

    vcf2maf: Vcf2Maf

    expression: Expression = Expression()
    """Include mRNA expression data"""

    copy_number_alteration: CNA = CNA()
    """Include copy number alteration results"""

    study: Study

    sample_info: dict[str, Any] = {}
    """Implementation must be re-designed"""
