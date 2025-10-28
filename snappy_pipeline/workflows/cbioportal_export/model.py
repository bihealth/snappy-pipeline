from __future__ import annotations

import enum
from typing import Any, TypedDict

from pydantic import ConfigDict, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToggleModel


class MappingTool(enum.StrEnum):
    BWA = "bwa"
    BWA_MEM2 = "bwa_mem2"
    MBCS = "mbcs"
    MINIMAP2 = "minimap2"


class ExpressionTool(enum.StrEnum):
    STAR = "star"


class SomaticVariantCallingTool(enum.StrEnum):
    MUTECT2 = "mutect2"


class SomaticVariantAnnotationTool(enum.StrEnum):
    VEP = "vep"


class CopyNumberTool(enum.StrEnum):
    CNVKIT = "cnvkit"

    CONTROL_FREEC = "Control_FREEC"
    """unsupported"""


class NcbiBuild(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class Vcf2Maf(SnappyModel):
    Center: str
    ncbi_build: NcbiBuild
    # Remember to move path_gene_id_mappings option here when re-factoring the step


class GenomeName(enum.StrEnum):
    hg19 = "hg19"  # GRCh37
    hg38 = "hg38"  # GRCh38
    mm10 = "mm10"  # GRCm38


class Expression(ToggleModel):
    path_ngs_mapping: str = "../ngs_mapping"
    """When missing, no expression data is uploaded to cBioPortal"""

    expression_tool: ExpressionTool = ExpressionTool.STAR


class SomaticVariantStep(enum.StrEnum):
    ANNOTATION = "somatic_variant_annotation"
    FILTER = "somatic_variant_filtration"


class CNA(ToggleModel):
    path_copy_number: str | None = None
    """When missing, no CNV data uploaded to portal. Access WES & WGS steps"""

    copy_number_tool: CopyNumberTool = CopyNumberTool.CNVKIT

    @model_validator(mode="after")
    def ensure_path_set_when_enabled(self):
        if self.enabled and not self.path_copy_number:
            raise ValueError("Copy number path must be set when copy_number_alteration is enabled")
        return self


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


class CbioportalExport(SnappyStepModel):
    model_config = ConfigDict(
        extra="forbid",
    )

    path_somatic_variant: str
    """Annotation is mandatory, but filtration is optional, can happen before or after annotation"""

    mapping_tool: MappingTool = MappingTool.BWA

    somatic_variant_calling_tool: SomaticVariantCallingTool = SomaticVariantCallingTool.MUTECT2
    """mutect/scalpel combo unsupported"""

    somatic_variant_step: SomaticVariantStep = SomaticVariantStep.FILTER
    """Which pipeline step is used to compute signatures"""

    somatic_variant_annotation_tool: SomaticVariantAnnotationTool = SomaticVariantAnnotationTool.VEP

    is_filtered: bool = True
    """Is the vcf post-filtered"""

    path_gene_id_mappings: str
    """Mapping from pipeline gene ids to cBioPortal ids (HGNC symbols from GeneNexus)"""

    exclude_variant_with_flag: str = ""
    """Required to filter variants (typically by type)"""

    vcf2maf: Vcf2Maf

    expression: Expression = Expression()
    """Include mRNA expression data"""

    copy_number_alteration: CNA = CNA()
    """Include copy number alteration results"""

    study: Study

    patient_info: dict[str, Any] = {}
    """unimplemented"""

    sample_info: dict[str, Any] = {}
    """Implementation must be re-designed"""
    # sample_info: dict[str, ExtraInfos] = Field(
    #     {},
    #     examples=[
    #         {
    #             "tumor_mutational_burden": dict(
    #                 name="TMB",
    #                 description="Tumor mutational burden computed on CDS regions",
    #                 datatype="NUMBER",
    #                 priority="2",
    #                 column="TMB",
    #             )
    #         }
    #     ],
    # )
    # """Each additional sample column must have a name and a (possibly empty) config attached."""

    @model_validator(mode="after")
    def ensure_filtration_are_configured_correctly(self):
        if self.somatic_variant_step == SomaticVariantStep.FILTER:
            if not self.is_filtered:
                raise ValueError(
                    "When the input step is 'somatic_variant_filtration', the filtration status must be set to 'True'"
                )
        return self
