from typing import Annotated

from pydantic import DirectoryPath, Field, FilePath

from snappy_pipeline.models import SnappyStepModel


class WgsSvExportExternal(SnappyStepModel):
    tool_ngs_mapping: str | None = None
    """used to create output file prefix."""

    tool_sv_calling_wgs: str | None = None
    """used to create output file prefix."""

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str = "id"
    """How to merge VCF, used in `bcftools --merge` call."""

    search_paths: Annotated[list[DirectoryPath], Field(min_length=1)]
    """path to all VCF files."""

    search_patterns: Annotated[
        list[dict[str, str]], Field(examples=[{"vcf": "*/*.vcf.gz"}], min_length=1)
    ]
    """list of search pattern"""

    release: str = "GRCh37"

    path_refseq_ser: FilePath
    """path to RefSeq .ser file"""

    path_ensembl_ser: FilePath
    """path to ENSEMBL .ser file"""

    path_db: FilePath
    """path to annotator DB file to use"""

    varfish_server_compatibility: bool = False
    """
    build output compatible with
    varfish-server v1.2 (Anthenea) and early versions of the v2 (Bollonaster)
    """
