from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants


class SomaticVariantSignaturesDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic"})),
        ExpectedPathSchema(ExpectedSomaticVariants),
    ]


class TabulateVcf(ToolModel):
    """Options of ``bcftools query``, which tabulates the somatic variants"""


class SomaticVariantSignatures(SnappyStepModel):
    depends_on: SomaticVariantSignaturesDependsOn

    tabulate_vcf: TabulateVcf = TabulateVcf()
