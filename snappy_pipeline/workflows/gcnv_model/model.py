import enum
from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedIndexedRegions,
    ExpectedPathSchema,
    ExpectedRegions,
    Kits,
    Reference,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments


class Intervals(enum.StrEnum):
    targeted = "targeted"
    wgs = "wgs"


class IntervalMergingRule(enum.StrEnum):
    ALL = "ALL"
    OVERLAPPING_ONLY = "OVERLAPPING_ONLY"


class GcnvPostprocess(SnappyModel):
    """Settings of ``PostprocessGermlineCNVCalls``, for the model and the callers of a task"""

    autosomal_ref_copy_number: int = 2
    """Reference copy number of the autosomes (``--autosomal-ref-copy-number``)"""

    allosomal_contigs: list[str] = ["X", "Y"]
    """
    Allosomal contigs (``--allosomal-contig``), without the ``chr`` prefix.
    The wrapper adds the prefix if the contigs of the ploidy calls have it.
    """


class ExpectedGcnvModel(BaseModel):
    """Consumer-driven contract: the model files of a ``gcnv_model`` task for the case mode.

    ``calls`` is a template with the wildcard ``{shard}`` unless the consumer names the shard.
    """

    interval_list: str
    ploidy: str
    calls: str
    shards: str


class GcnvModelDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ]

    reference: Reference

    kits: Kits = ""
    """Target regions of the library kits, required for ``intervals: targeted``."""

    mappable_regions: Annotated[
        str,
        DataSignature(DataType.REGIONS, frozenset({"mappable"})),
        ExpectedPathSchema(ExpectedIndexedRegions),
    ]
    """Task with the BED file of uniquely mappable regions, with its ``.idx`` (plain BED) or
    ``.tbi`` (bgzipped BED) next to it"""

    par_regions: Annotated[
        str,
        DataSignature(DataType.REGIONS, frozenset({"par"})),
        ExpectedPathSchema(ExpectedRegions),
    ] = ""
    """Task with the BED file of the pseudo-autosomal regions, excluded from the contig ploidy"""


class GcnvModel(SnappyStepModel, ToolModel, GcnvPostprocess):
    depends_on: GcnvModelDependsOn

    intervals: Annotated[Intervals, EnumField(Intervals)]
    """
    ``targeted`` bins the target regions of each library kit, ``wgs`` bins the whole genome
    (one kit, named ``wgs``).
    """

    extra_arguments: list[str] = []
    """Additional arguments for ``GermlineCNVCaller``, the other GATK commands take none"""

    interval_merging_rule: IntervalMergingRule = IntervalMergingRule.OVERLAPPING_ONLY
    """
    Interval merging rule of all GATK commands that read intervals (``--interval-merging-rule``).
    The callers of the model use it for their coverage, too.
    """

    bin_length: int = 0
    """Only for ``intervals: targeted``: bin length of ``PreprocessIntervals`` (``--bin-length``), ``0`` keeps the target regions as they are"""

    padding: int = 0
    """Only for ``intervals: wgs``: padding of ``PreprocessIntervals`` (``--padding``)"""

    minimum_gc_content: float = 0.1
    """``FilterIntervals`` ``--minimum-gc-content``"""

    maximum_gc_content: float = 0.9
    """``FilterIntervals`` ``--maximum-gc-content``"""

    minimum_mappability: float = 0.9
    """``FilterIntervals`` ``--minimum-mappability``"""

    maximum_mappability: float = 1.0
    """``FilterIntervals`` ``--maximum-mappability``"""

    intervals_per_shard: int = 5000
    """Number of intervals per shard of ``GermlineCNVCaller`` (``IntervalListTools --SCATTER_CONTENT``)"""

    autosome_ploidy_prior: Annotated[list[float], Field(min_length=4, max_length=4)] = [
        0,
        0.01,
        0.98,
        0.01,
    ]
    """Prior probabilities of ploidy 0 to 3 for the autosomes (``DetermineGermlineContigPloidy`` priors table)"""

    x_ploidy_prior: Annotated[list[float], Field(min_length=4, max_length=4)] = [
        0.01,
        0.49,
        0.49,
        0.01,
    ]
    """Prior probabilities of ploidy 0 to 3 for the contig X"""

    y_ploidy_prior: Annotated[list[float], Field(min_length=4, max_length=4)] = [
        0.495,
        0.495,
        0.01,
        0,
    ]
    """Prior probabilities of ploidy 0 to 3 for the contig Y"""

    @model_validator(mode="after")
    def check_kits(self):
        if self.intervals == Intervals.targeted and not self.depends_on.kits:
            raise ValueError("intervals=targeted needs depends_on.kits")
        if self.intervals == Intervals.wgs and self.depends_on.kits:
            raise ValueError("intervals=wgs has no depends_on.kits")
        return self
