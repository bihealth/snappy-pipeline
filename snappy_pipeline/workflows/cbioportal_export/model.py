from __future__ import annotations

import enum
from typing import Annotated, TypedDict

from pydantic import ConfigDict, Field, model_validator

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


class FiltrationSchema(enum.StrEnum):
    unfiltered = "unfiltered"
    list = "list"
    sets = "sets"


class FilterSet(enum.StrEnum):
    NO_FILTER = "no_filter"
    DKFZ_ONLY = "dkfz_only"
    DKFZ_AND_EBFILTER = "dkfz_and_ebfilter"
    DKFZ_AND_EBFILTER_AND_OXOG = "dkfz_and_ebfilter_and_oxog"


class Filter(SnappyModel):
    pass


class DkfzFilter(Filter):
    pass


class EbFilter(Filter):
    ebfilter_threshold: float = 2.4
    shuffle_seed: int = 1
    panel_of_normals_size: int = 25
    min_mapq: int = 20
    min_baseq: int = 15


class Bcftools(Filter):
    include: str = ""
    exclude: str = ""


class Regions(Filter):
    path_bed: str = ""


class Protected(Filter):
    path_bed: str = ""


class CopyNumberTool(enum.StrEnum):
    CNVKIT = "cnvkit"

    CONTROL_FREEC = "Control_FREEC"
    """unsupported"""

    COPYWRITER = "CopywriteR"
    """unmaintained"""


class NcbiBuild(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class Vcf2Maf(SnappyModel):
    Center: str
    ncbi_build: NcbiBuild


class GenomeName(enum.StrEnum):
    hg19 = "hg19"  # GRCh37
    hg38 = "hg38"  # GRCh38
    mm10 = "mm10"  # GRCm38


class Expression(ToggleModel):
    path_ngs_mapping: str = "../ngs_mapping"
    """When missing, no expression data is uploaded to cBioPortal"""

    expression_tool: ExpressionTool = ExpressionTool.STAR


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
    reference_genome: GenomeName


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
    """REQUIRED (before or after filtration)"""

    mapping_tool: MappingTool = MappingTool.BWA

    somatic_variant_calling_tool: SomaticVariantCallingTool = SomaticVariantCallingTool.MUTECT2
    """mutect/scalpel combo unsupported"""

    somatic_variant_annotation_tool: SomaticVariantAnnotationTool = SomaticVariantAnnotationTool.VEP

    filtration_schema: FiltrationSchema = FiltrationSchema.list
    """Method of variant filtration (if any)"""

    filter_before_annotation: bool = False
    """Flag if filtration has been done before annotation"""

    # TODO: remove filtration by sets
    # When this is done, the deprecated options can be deleted, and the
    # signature computations will be done just on the output of whatever step
    # is used on input (somatic_variant_calling, somatic_variant_annotation or somatic_variant_filtration)

    filter_set: Annotated[FilterSet | None, Field(None, deprecated="use `filter_list` instead")]
    """
    DEPRECATED: use `filter_list` in `somatic_variant_filtration` instead.
    Must be set when the ``filtration_schema`` is ``sets``.
    """

    exon_list: Annotated[
        str | None,
        Field(
            "genome_wide",
            deprecated="Works together with filter_set, ignored when `filter_list` is selected",
        ),
    ]
    """
    DEPRECATED: use `filter_list` in `somatic_variant_filtration` instead.
    Must be set when the ``filtration_schema`` is ``sets``.
    """

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

    patient_info: None = None
    """unimplemented"""

    sample_info: None = None
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
    def ensure_tools_are_configured(self):
        if self.filtration_schema == FiltrationSchema.sets:
            if self.filter_set is None or self.exon_list is None:
                raise ValueError(
                    "'filter_set' & 'exon_list' must be set when the filtration schema is on sets"
                )
        return self
