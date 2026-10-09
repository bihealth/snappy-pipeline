from typing import Annotated

from pydantic import BaseModel, FilePath, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    Reference,
)


class ExpectedExternalVcf(BaseModel):
    """Consumer-driven contract: the VCF of a library."""

    vcf: str


class ExpectedExternalBam(BaseModel):
    """Consumer-driven contract: the BAM file and index of a library."""

    bam: str
    bai: str


class VariantExportExternalDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedExternalVcf),
    ]
    """``external_data`` task with the (g)VCF of each library."""

    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedExternalBam),
    ] = ""
    """``external_data`` task with the BAM file of each library; needed for bam_available_flag."""

    reference: Reference


class TargetCoverageReport(SnappyModel):
    path_targets_bed: str | None = None
    """
    Mapping from enrichment kit to target region BED file, for either computing per target
    region coverage or selecting targeted exons. Only used if 'bam_available_flag' is True.
    It will not generate detailed reporting.
    """


class VariantExportExternal(SnappyStepModel):
    depends_on: VariantExportExternalDependsOn

    external_tool: str = "dragen"
    """external tool name."""

    bam_available_flag: bool = False
    """BAM QC only possible if BAM files are present."""

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str | None = None
    """How to merge VCF, used in `bcftools --merge` argument."""

    gvcf_option: bool = True
    """Flag to indicate if inputs are genomic VCFs."""

    release: str = "GRCh37"
    """genome release; default 'GRCh37'."""

    path_exon_bed: str = ""
    """Path to BED file with exons; used for reducing data to near-exon small variants."""

    path_refseq_ser: FilePath
    """path to RefSeq .ser file."""

    path_ensembl_ser: FilePath
    """path to ENSEMBL .ser file."""

    path_db: FilePath
    """path to annotator DB file to use."""

    target_coverage_report: TargetCoverageReport = TargetCoverageReport()

    @model_validator(mode="after")
    def bam_qc_needs_alignments(self):
        if self.bam_available_flag and not self.depends_on.alignments:
            raise ValueError("bam_available_flag needs depends_on.alignments")
        return self
