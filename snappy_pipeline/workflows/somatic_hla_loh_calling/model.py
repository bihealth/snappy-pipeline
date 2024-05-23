from models import SnappyStepModel


class SomaticHlaLohCalling(SnappyStepModel):
    path_ngs_mapping: str

    path_hla_typing: str

    path_somatic_purity_ploidy: str
