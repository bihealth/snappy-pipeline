from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
    Regions,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class ExpectedSomaticCnvCalls(SnappyModel):
    """Consumer-driven contract: expected output keys from somatic CNV caller steps."""

    dnacopy_seg: str


class SomaticCnvCheckingDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]
    copy_number: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})),
        ExpectedPathSchema(ExpectedSomaticCnvCalls),
    ] = ""
    """Somatic CNV calling task whose segments (``dnacopy_seg``) are merged with the B-allele frequencies"""

    reference: Reference

    excluded_regions: Regions = ""
    """Task with a BED file of regions to be excluded"""


class SomaticCnvChecking(SnappyStepModel):
    depends_on: SomaticCnvCheckingDependsOn

    max_depth: int = 10000
    """Max depth for pileups"""

    min_depth: int = 20
    """Minimum depth for reference and alternative alleles to consider variant"""

    min_baf: Annotated[float, Field(0.4, ge=0, le=0.5)]
    """Maximum BAF to consider variant as heterozygous (between 0 & 1/2)"""
