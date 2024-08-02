from snappy_pipeline.models import SnappyStepModel


class SomaticHlaLohCalling(SnappyStepModel):
    path_ngs_mapping: str = "../ngs_mapping"

    path_hla_typing: str

    path_somatic_purity_ploidy: str
