import enum

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


class SignatureSet(enum.StrEnum):
    COSMIC = "cosmic"
    NATURE2013 = "nature2013"


class TriCountsMethod(enum.StrEnum):
    DEFAULT = "default"
    EXOME = "exome"
    GENOME = "genome"
    EXOME2GENOME = "exome2genome"
    GENOME2EXOME = "genome2exome"


class DeconstructSigs(SnappyModel):
    """Options of ``deconstructSigs::whichSignatures``"""

    signatures: SignatureSet = SignatureSet.COSMIC
    """Reference signatures: the object ``signatures.<value>`` of the package (``signatures.ref``)"""

    tri_counts_method: TriCountsMethod = TriCountsMethod.DEFAULT
    """Normalization of the trinucleotide counts (``tri.counts.method``)"""


class SomaticVariantSignatures(SnappyStepModel):
    depends_on: SomaticVariantSignaturesDependsOn

    tabulate_vcf: TabulateVcf = TabulateVcf()

    deconstruct_sigs: DeconstructSigs = DeconstructSigs()
