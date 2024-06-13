from snappy_pipeline.models import SnappyModel, SnappyStepModel


class Gcnv(SnappyModel):
    path_par_intervals: str = ""
    """Path to interval block list with PAR region for contig calling."""

    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""

    #NOTE: the wgs model do NOT need the path_target_interval_list_mapping

class HelperGcnvModelWgs(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    gcnv: Gcnv
