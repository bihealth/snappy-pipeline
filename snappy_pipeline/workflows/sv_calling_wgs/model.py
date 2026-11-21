import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import ReadBy, RequiredBy
from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    Reference,
    Regions,
)
from snappy_pipeline.workflows.common import sv_calling_model
from snappy_pipeline.workflows.common.sv_calling_model import (
    Delly2,
    Gcnv,
    GeneRegions,
    Manta,
    SkipLibraries,
)
from snappy_pipeline.workflows.gcnv_model.model import ExpectedGcnvModel


class Tool(enum.StrEnum):
    delly2 = "delly2"
    manta = "manta"
    popdel = "popdel"
    gcnv = "gcnv"
    melt = "melt"


class Melt(sv_calling_model.Melt):
    me_refs_infix: str = "1KGP_Hg19"


class Popdel(ToolModel):
    num_jobs: Annotated[int, Field(ge=1)] = 24
    """Number of jobs that share the windows in genome order, each at least one window and about equal bases.
    A reference with fewer windows than jobs is an error."""

    window_size: Annotated[int, Field(ge=1)] = 10000000
    """Length of the windows that the contigs are cut into from their start; popdel calls each window on its own."""

    max_sv_size: int = 20000
    """Bases added on both sides of each window, within its contig, so that every deletion up to this length lies whole in a padded window."""

    profile_regions: list[str] = [
        "chr1:35000000-36000000",
        "chr2:174000000-175000000",
        "chr3:36500000-37500000",
        "chr4:88000000-89000000",
        "chr5:38000000-39000000",
        "chr6:38000000-39000000",
        "chr7:38000000-39000000",
        "chr8:19000000-20000000",
        "chr9:19000000-20000000",
        "chr10:19000000-20000000",
        "chr11:19000000-20000000",
        "chr12:19000000-20000000",
        "chr13:25000000-26000000",
        "chr14:25000000-26000000",
        "chr15:25000000-26000000",
        "chr16:25000000-26000000",
        "chr17:31000000-32000000",
        "chr18:31000000-32000000",
        "chr19:31000000-32000000",
        "chr20:33000000-34000000",
        "chr21:21000000-22000000",
        "chr22:25000000-26000000",
    ]
    """
    Regions (``chrom:begin-end``) of ``popdel profile`` (``-i``), written with the ``chr`` prefix.
    The wrapper drops the prefix for a reference without it.
    """

    skip_libraries: SkipLibraries = []


class SvCallingWgsDependsOn(SnappyModel):
    alignments: DnaAlignments

    reference: Reference

    gcnv_model: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"gcnv", "wgs"})),
        ExpectedPathSchema(ExpectedGcnvModel),
        RequiredBy(Tool.gcnv),
    ] = ""
    """A ``gcnv_model`` task with ``intervals: wgs``."""

    excluded_regions: Annotated[Regions, ReadBy(Tool.delly2)] = ""
    """Task with a BED file of regions that ``delly2`` excludes"""

    gene_regions: Annotated[GeneRegions, RequiredBy(Tool.melt)] = ""
    """Task with the BED file of genes, such as ``add_bed_files/1KGP_Hg19/hg19.genes.bed`` of MELT"""


class SvCallingWgs(SnappyStepModel):
    depends_on: SvCallingWgsDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    delly2: Delly2 | None = None

    manta: Manta | None = None

    popdel: Popdel | None = None

    gcnv: Gcnv | None = None

    melt: Melt | None = None

    ignore_contigs: list[str] = [
        "^NC_007605$",
        "^hs37d5$",
        "^chrEBV$",
        "_decoy$",
        "^HLA-",
        "^chrUn_",
    ]
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""
