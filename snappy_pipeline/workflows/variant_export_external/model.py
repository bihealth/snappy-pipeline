from typing import Annotated

from snappy_pipeline.models import ConfigPath
from snappy_pipeline.models import RequiredBy, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    Exons,
    ExpectedAlignments,
    ExpectedPathSchema,
    ExpectedVcf,
    Kits,
    Reference,
)


class VariantExportExternalDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedVcf),
    ]
    """``external_data`` task with the (g)VCF of each library."""

    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
        RequiredBy(when="bam_available_flag"),
    ] = ""
    """``external_data`` task with the BAM file of each library"""

    reference: Reference

    kits: Annotated[Kits, RequiredBy(when="bam_available_flag")] = ""
    """Task with the ``targets`` of each kit"""

    exons: Exons = ""
    """Optional task with a BED file of exons; the export keeps only the small variants near them"""


class VariantExportExternal(SnappyStepModel):
    depends_on: VariantExportExternalDependsOn

    bam_available_flag: bool = False
    """BAM QC only possible if BAM files are present."""

    ignore_read_groups: bool = True
    """BAM QC: the target coverage of all read groups of a library together (``alfred qc --ignore``)"""

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str | None = None
    """How to merge VCF, used in `bcftools --merge` argument."""

    gvcf_option: bool = True
    """Flag to indicate if inputs are genomic VCFs."""

    release: str = "GRCh37"
    """genome release; default 'GRCh37'."""

    path_refseq_ser: ConfigPath
    """path to RefSeq .ser file."""

    path_ensembl_ser: ConfigPath
    """path to ENSEMBL .ser file."""

    path_db: ConfigPath
    """path to annotator DB file to use."""
