from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Gcnv(SnappyModel):
    path_par_intervals: str = ""
    """Path to interval block list with PAR region for contig calling."""

    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""

    # NOTE: the wgs model do NOT need the path_target_interval_list_mapping


class HelperGcnvModelWgsDependsOn(SnappyModel):
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"


class HelperGcnvModelWgs(SnappyStepModel):
    depends_on: HelperGcnvModelWgsDependsOn = Field(default_factory=HelperGcnvModelWgsDependsOn)

    gcnv: Gcnv
