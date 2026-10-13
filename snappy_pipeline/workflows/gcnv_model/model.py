import enum
from typing import Annotated

from pydantic import BaseModel, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.gcnv import TargetIntervalEntry
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Intervals(enum.StrEnum):
    targeted = "targeted"
    wgs = "wgs"


class ExpectedGcnvModel(BaseModel):
    """Consumer-driven contract: the model files of a ``gcnv_model`` task for the case mode.

    ``calls`` is a template with the wildcard ``{shard}`` unless the consumer names the shard.
    """

    interval_list: str
    ploidy: str
    calls: str
    shards: str


class Gcnv(SnappyModel):
    path_par_intervals: str = ""
    """Path to interval block list with PAR region for contig calling."""

    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""

    path_target_interval_list_mapping: list[TargetIntervalEntry] = []
    """
    Required for ``intervals: targeted``, not allowed for ``intervals: wgs``.
    Defines one or more set of target intervals, one library kit each.
    Example:
     - name: "Agilent SureSelect Human All Exon V6"
       pattern: "Agilent SureSelect Human All Exon V6.*"
       path: "path/to/targets.bed"
    """


class GcnvModelDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference


class GcnvModel(SnappyStepModel):
    depends_on: GcnvModelDependsOn

    intervals: Annotated[Intervals, EnumField(Intervals)]
    """
    ``targeted`` bins the target regions of each library kit, ``wgs`` bins the whole genome
    (one kit, named ``wgs``).
    """

    gcnv: Gcnv

    @model_validator(mode="after")
    def check_target_interval_list_mapping(self):
        mapping = self.gcnv.path_target_interval_list_mapping
        if self.intervals == Intervals.targeted and not mapping:
            raise ValueError("intervals=targeted needs gcnv.path_target_interval_list_mapping")
        if self.intervals == Intervals.wgs and mapping:
            raise ValueError("intervals=wgs has no gcnv.path_target_interval_list_mapping")
        return self
