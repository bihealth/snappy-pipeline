import enum
from typing import Annotated

from snappy_pipeline.models import ReadBy, RequiredBy
from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    Reference,
    Regions,
)
from snappy_pipeline.workflows.common.sv_calling_model import (
    Delly2,
    Gcnv,
    GeneRegions,
    Manta,
    Melt,
)
from snappy_pipeline.workflows.gcnv_model.model import ExpectedGcnvModel


class Tool(enum.StrEnum):
    gcnv = "gcnv"
    delly2 = "delly2"
    manta = "manta"
    melt = "melt"


class SvCallingTargetedDependsOn(SnappyModel):
    alignments: DnaAlignments

    reference: Reference

    gcnv_model: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"gcnv", "targeted"})),
        ExpectedPathSchema(ExpectedGcnvModel),
        RequiredBy(Tool.gcnv),
    ] = ""
    """A ``gcnv_model`` task with ``intervals: targeted``."""

    excluded_regions: Annotated[Regions, ReadBy(Tool.delly2)] = ""
    """Task with a BED file of regions that ``delly2`` excludes"""

    gene_regions: Annotated[GeneRegions, RequiredBy(Tool.melt)] = ""
    """Task with the BED file of genes, such as ``add_bed_files/1KGP_Hg19/hg19.genes.bed`` of MELT"""


class SvCallingTargeted(SnappyStepModel):
    depends_on: SvCallingTargetedDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    gcnv: Gcnv | None = None

    delly2: Delly2 | None = None

    manta: Manta | None = None

    melt: Melt | None = None
