from __future__ import annotations

from enum import Enum
from typing import Any, Annotated

from pydantic import ConfigDict, Field

from models import SnappyStepModel, SnappyModel


class MappingTool(Enum):
    BWA = "bwa"


class ExpressionTool(Enum):
    STAR = "star"


class SomaticVariantCallingTool(Enum):
    MUTECT2 = "mutect2"


class SomaticVariantAnnotationTool(Enum):
    VEP = "vep"


class FilterSet(Enum):
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


class CopyNumberTool(Enum):
    CNVKIT = "cnvkit"

    CONTROL_FREEC = "Control_FREEC"
    """unsupported"""

    COPYWRITER = "CopywriteR"
    """unmaintained"""


class NcbiBuild(Enum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class Vcf2Maf(SnappyModel):
    Center: str
    ncbi_build: NcbiBuild


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


class CbioportalExport(SnappyStepModel):
    model_config = ConfigDict(
        extra="forbid",
    )

    path_ngs_mapping: str | None = None
    """When missing, no expression data is uploaded to cBioPortal"""

    mapping_tool: MappingTool = MappingTool.BWA

    expression_tool: ExpressionTool = ExpressionTool.STAR

    path_somatic_variant: str
    """REQUIRED (before or after filtration)"""

    somatic_variant_calling_tool: SomaticVariantCallingTool = SomaticVariantCallingTool.MUTECT2
    """mutect/scalpel combo unsupported"""

    somatic_variant_annotation_tool: SomaticVariantAnnotationTool = SomaticVariantAnnotationTool.VEP

    filter_set: Annotated[FilterSet | None, Field(None, deprecated="use `filter_list` instead")]
    """
    DEPRECATED: use `filter_list instead`.
    Set it to an empty value when using annotated variants without filtration.
    """

    exon_list: Annotated[
        str | None,
        Field(
            "genome_wide",
            deprecated="Works together with filter_set, ignored when `filter_list` is selected",
        ),
    ]
    """
    DEPRECATED.
    Works together with filter_set, ignored when "filter_list" is selected
    """

    filter_list: list[Filter] = []

    exclude_variant_with_flag: str | None = None
    """Required for Copy Number Alterations"""

    path_copy_number: str | None = None
    """When missing, no CNV data uploaded to portal. Access WES & WGS steps"""

    copy_number_tool: CopyNumberTool = CopyNumberTool.CNVKIT

    path_gene_id_mappings: str
    """Mapping from pipeline gene ids to cBioPortal ids (HGNC symbols from GeneNexus)"""

    vcf2maf: Vcf2Maf

    study: Study

    patient_info: None
    """unimplemented"""

    sample_info: Any
    """Each additional sample column must have a name and a (possibly empty) config attached."""
