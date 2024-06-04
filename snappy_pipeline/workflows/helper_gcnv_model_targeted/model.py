from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.models.gcnv import TargetIntervalEntry


class Gcnv(SnappyModel):
    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""

    path_par_intervals: str = ""

    path_target_interval_list_mapping: list[TargetIntervalEntry] = []


class HelperGcnvModelTargeted(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    gcnv: Gcnv
