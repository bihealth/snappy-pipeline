from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.models.gcnv import TargetIntervalEntry


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


class HelperGcnvModelTargeted(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    gcnv: Gcnv
