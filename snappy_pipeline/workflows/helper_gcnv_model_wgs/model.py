from pydantic import DirectoryPath

from snappy_pipeline.models import SnappyStepModel, SnappyModel


class Gcnv(SnappyModel):
    path_uniquely_mapable_bed: str
    """path to BED file with uniquely mappable regions."""


class HelperGcnvModelWgs(SnappyStepModel):
    path_ngs_mapping: DirectoryPath | str

    gcnv: Gcnv
