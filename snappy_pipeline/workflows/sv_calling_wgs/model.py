import enum
from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import ConfigPath
from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators
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
    delly2 = "delly2"
    manta = "manta"
    popdel = "popdel"
    gcnv = "gcnv"
    melt = "melt"
    # sniffles2 = "sniffles2"  # exclude for now


class Gcnv(SnappyModel):
    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Delly2(SnappyModel):
    map_qual: int = 1

    geno_qual: int = 5

    qual_tra: int = 20

    mad_cutoff: int = 9

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Manta(SnappyModel):
    num_threads: int = 16

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Melt(SnappyModel):
    me_refs_infix: str = "1KGP_Hg19"
    me_types: list[str] = ["ALU", "LINE1", "SVA"]
    me_refs_path: ConfigPath = Field(examples=["/path/to/melt/support/files"])
    jar_file: ConfigPath = Field(examples=["/path/to/MELTv2.2.2/MELT.jar"])

    me_refs_path: ConfigPath = Field(examples=["me_refs"])

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Popdel(SnappyModel):
    window_size: int = 10000000

    max_sv_size: int = 20000
    """== padding"""

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class Sniffles2(SnappyModel):
    tandem_repeats: Annotated[
        str,
        Field(examples=["/path/to/sniffles2/trf/GRCh37/human_hs37d5.trf.bed"]),
    ]

    skip_libraries: list[str] = []
    """
    Skip processing of the following libraries.
    If the library is in family/pedigree then all of the family/pedigree will be skipped.
    """


class SvCallingWgsDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    gcnv_model: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"gcnv", "wgs"})),
        ExpectedPathSchema(ExpectedGcnvModel),
    ] = ""
    """Required for ``tool: gcnv``: a ``gcnv_model`` task with ``intervals: wgs``."""

    excluded_regions: Regions = ""
    """Task with a BED file of regions that ``delly2`` excludes"""

    gene_regions: GeneRegions = ""
    """Required for ``tool: melt``: task with the BED file of genes, such as ``add_bed_files/1KGP_Hg19/hg19.genes.bed`` of MELT"""


TOOL_DEPENDENCIES = {Tool.gcnv: ("gcnv_model",), Tool.melt: ("gene_regions",)}


class SvCallingWgs(SnappyStepModel):
    depends_on: SvCallingWgsDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    delly2: Delly2 | None = None

    manta: Manta | None = None

    popdel: Popdel | None = None

    gcnv: Gcnv | None = None

    melt: Melt | None = None

    sniffles2: Sniffles2 | None = None

    ignore_chroms: list[str] = [
        "^NC_007605$",
        "^hs37d5$",
        "^chrEBV$",
        "_decoy$",
        "^HLA-",
        "^chrUn_",
    ]
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""

    @model_validator(mode="after")
    def ensure_tools_are_configured(self):
        if not getattr(self, self.tool):
            raise ValueError(f"Tool {self.tool} not configured")
        return self

    @model_validator(mode="after")
    def validate_tool_dependencies(self):
        validators.require_tool_dependencies(self, TOOL_DEPENDENCIES)
        return self
