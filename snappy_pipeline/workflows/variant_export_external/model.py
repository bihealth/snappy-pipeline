from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyStepModel, SnappyModel


class TargetCoverageReport(SnappyModel):
    path_targets_bed: str | None = None
    """
    Mapping from enrichment kit to target region BED file, for either computing per target
    region coverage or selecting targeted exons. Only used if 'bam_available_flag' is True.
    It will not generate detailed reporting.
    """


class VariantExportExternal(SnappyStepModel):
    external_tool: str = "dragen"
    """external tool name."""

    bam_available_flag: bool
    """BAM QC only possible if BAM files are present."""

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str | None = None
    """How to merge VCF, used in `bcftools --merge` argument."""

    gvcf_option: bool = True
    """Flag to indicate if inputs are genomic VCFs."""

    search_paths: Annotated[list[str], Field(min_length=1)]
    """list of paths to VCF files."""

    search_patterns: Annotated[
        list[dict[str, str]],
        Field(examples=[{"vcf": "*.vcf.gz"}, {"bam": "*.bam"}, {"bai": "*.bam.bai"}], min_length=1),
    ]
    """list of search patterns"""

    release: str = "GRCh37"
    """genome release; default 'GRCh37'."""

    path_exon_bed: str  # FIXME old code: "path_exon_bed: null # REQUIRED: exon BED file to use"
    """Path to BED file with exons; used for reducing data to near-exon small variants."""

    path_refseq_ser: str
    """path to RefSeq .ser file."""

    path_ensembl_ser: str
    """path to ENSEMBL .ser file."""

    path_db: str
    """path to annotator DB file to use."""

    target_coverage_report: TargetCoverageReport = TargetCoverageReport()
