import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import ConfigPath
from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel, validators
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
    Regions,
)
from snappy_pipeline.workflows.common.sv_calling_model import GeneRegions
from snappy_pipeline.workflows.gcnv_model.model import ExpectedGcnvModel
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Tool(enum.StrEnum):
    gcnv = "gcnv"
    delly2 = "delly2"
    manta = "manta"
    melt = "melt"


class Gcnv(ToolModel):
    extra_arguments: list[str] = []
    """Additional arguments for ``GermlineCNVCaller``, the other GATK commands take none"""

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Delly2(ToolModel):
    map_qual: int = 1

    geno_qual: int = 5

    qual_tra: int = 20

    mad_cutoff: int = 9

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Manta(ToolModel):
    num_threads: int = 16

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Melt(ToolModel):
    me_refs_infix: str = Field(examples=["1KGP_Hg19"])

    me_types: list[str] = ["ALU", "LINE1", "SVA"]

    exclude_contigs: list[str] = ["hs37d5", "NC_007605"]
    """Contigs that IndivAnalysis ignores (MELT's ``-b``), here the GRCh37 decoy and EBV"""

    read_length: int = 150
    """Read length, for IndivAnalysis and GroupAnalysis (MELT's ``-r``)"""

    max_no_call_percent: int = 100
    """Percentage of samples with a no-call genotype that MakeVCF tolerates at a site (MELT's ``-j``)"""

    jar_file: ConfigPath = Field(examples=["MELT.jar"])

    me_refs_path: ConfigPath = Field(examples=["me_refs"])

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class SvCallingTargetedDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    gcnv_model: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"gcnv", "targeted"})),
        ExpectedPathSchema(ExpectedGcnvModel),
    ] = ""
    """Required for ``tool: gcnv``: a ``gcnv_model`` task with ``intervals: targeted``."""

    excluded_regions: Regions = ""
    """Task with a BED file of regions that ``delly2`` excludes"""

    gene_regions: GeneRegions = ""
    """Required for ``tool: melt``: task with the BED file of genes, such as ``add_bed_files/1KGP_Hg19/hg19.genes.bed`` of MELT"""


TOOL_DEPENDENCIES = {Tool.gcnv: ("gcnv_model",), Tool.melt: ("gene_regions",)}


class SvCallingTargeted(SnappyStepModel):
    depends_on: SvCallingTargetedDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    gcnv: Gcnv | None = None

    delly2: Delly2 | None = None

    manta: Manta | None = None

    melt: Melt | None = None

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
