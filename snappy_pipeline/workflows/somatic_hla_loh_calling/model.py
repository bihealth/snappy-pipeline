from pydantic import DirectoryPath

from snappy_pipeline.models import SnappyStepModel


class SomaticHlaLohCalling(SnappyStepModel):
    path_ngs_mapping: DirectoryPath | str

    path_hla_typing: str

    path_somatic_purity_ploidy: str
