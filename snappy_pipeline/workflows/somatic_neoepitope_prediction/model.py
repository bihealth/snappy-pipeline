import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Format(enum.StrEnum):
    stringtie = "stringtie"
    snappy_custom = "snappy_custom"
    cufflinks = "cufflinks"
    kallisto = "kallisto"
    custom = "custom"


class Preparation(SnappyModel):
    format: Annotated[Format, EnumField(Format)] = Format.snappy_custom
    "The file format of the expression file to process. (stringtie, kallisto, cufflinks, snappy_custom, custom)"
    "Use `custom` to process file formats not explicitly supported."
    "The `custom` option requires the use of the --id-column and --expression-column arguments."
    path_features: str = ""
    """Gencode file path, required for star and snappy format"""
    mode: Annotated[str, Field(examples=["gene", "transcript"])] = "gene"
    """Determine whether the expression file contains gene or transcript TPM values."""
    full_vep_annotation: bool = True
    """The somatic_variant_annotation has been done in fully annotation mode"""
    id_column: str = ""
    """Gene/Transcript id column name. Need when using the `custom` format."""
    expression_column: str = ""
    """Expression column name. Need when using the `custom` format."""
    ignore_ensembl_id_version: bool = True
    """Ignore the ensemble id version"""
    max_depth: int = 4000


class SomaticNeoepitopePrediction(SnappyStepModel):
    path_somatic_variant_annotation: Annotated[
        str, Field(examples=["../somatic_variant_annotation"])
    ]
    path_rna_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]
    tools_somatic_variant_annotation: list[str] = ["vep"]
    tools_ngs_mapping: list[str] = []
    tools_somatic_variant_calling: list[str] = []
    tools_rna_mapping: list[str] = []
    tools_rna_mapping: list[str] = []
    """Deafult to those configured for ngs_mapping"""
    tools_ngs_mapping: list[str] = []
    """Deafult to those configured for ngs_mapping"""
    tools_somatic_variant_calling: list[str] = []
    """Deafult to those configured for somatic_variant_calling"""
    preparation: Preparation | None = None
