from snappy_pipeline.models import SnappyModel, SnappyStepModel


class Gcnv(SnappyModel):
    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""


class HelperGcnvModelWgs(SnappyStepModel):
    path_ngs_mapping: str

    gcnv: Gcnv
