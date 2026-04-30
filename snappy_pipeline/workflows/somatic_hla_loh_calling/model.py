from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class SomaticHlaLohCallingDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"
    hla_typing: str = "hla_typing"


class SomaticHlaLohCalling(SnappyStepModel):
    depends_on: SomaticHlaLohCallingDependsOn = Field(default_factory=SomaticHlaLohCallingDependsOn)

    path_somatic_purity_ploidy: str
