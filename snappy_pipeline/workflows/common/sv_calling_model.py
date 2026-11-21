"""Models shared by ``sv_calling_{targeted,wgs}``: the tool sections and field types"""

from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import ConfigPath, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedRegions,
)
from snappy_pipeline.workflows.gcnv_model.model import GcnvPostprocess

#: ``depends_on`` field type for the BED file of the genes of MELT
GeneRegions = Annotated[
    str, DataSignature(DataType.REGIONS, frozenset({"genes"})), ExpectedPathSchema(ExpectedRegions)
]

#: The ``skip_libraries`` field of every tool section; a field type rather than a base class, so
#: that it stays the last field of the section (pydantic puts inherited fields first)
SkipLibraries = Annotated[
    list[str],
    Field(
        description=(
            "Skip processing of the following libraries.\n"
            "If the library is in family/pedigree then all of the family/pedigree will be skipped."
        )
    ),
]


class Gcnv(ToolModel, GcnvPostprocess):
    extra_arguments: list[str] = []
    """Additional arguments for ``GermlineCNVCaller``, the other GATK commands take none"""

    skip_libraries: SkipLibraries = []


class Delly2(ToolModel):
    map_qual: int = 1

    geno_qual: int = 5

    qual_tra: int = 20

    mad_cutoff: int = 9

    skip_libraries: SkipLibraries = []


class Manta(ToolModel):
    num_threads: int = 16

    skip_libraries: SkipLibraries = []


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

    skip_libraries: SkipLibraries = []
