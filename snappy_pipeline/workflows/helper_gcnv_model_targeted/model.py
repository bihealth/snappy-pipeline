from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.models.gcnv import TargetIntervalEntry
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Gcnv(SnappyModel):
    path_par_intervals: str = ""
    """Path to interval block list with PAR region for contig calling."""

    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""

    path_target_interval_list_mapping: list[TargetIntervalEntry]
    """
    The following allows to define one or more set of target intervals.  This is only used by gcnv.
    Example:
     - name: "Agilent SureSelect Human All Exon V6"
       pattern: "Agilent SureSelect Human All Exon V6.*"
       path: "path/to/targets.bed"
    """


class HelperGcnvModelTargetedDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"


class HelperGcnvModelTargeted(SnappyStepModel):
    depends_on: HelperGcnvModelTargetedDependsOn = Field(
        default_factory=HelperGcnvModelTargetedDependsOn
    )

    gcnv: Gcnv
