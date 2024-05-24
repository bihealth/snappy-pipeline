from typing import Annotated

from pydantic import Field

from models import SnappyStepModel


class WgsSvExportExternal(SnappyStepModel):
    tool_ngs_mapping: str | None = None
    """used to create output file prefix."""

    tool_sv_calling_wgs: str | None = None
    """used to create output file prefix."""

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str = "id"
    """How to merge VCF, used in `bcftools --merge` call."""

    search_paths: list[str]
    """path to all VCF files."""

    search_patterns: Annotated[list[dict[str, str]], Field(examples=[{"vcf": "*/*.vcf.gz"}])]
    """list of search pattern"""

    release: str = "GRCh37"

    path_refseq_ser: str
    """path to RefSeq .ser file"""

    path_ensembl_ser: str
    """path to ENSEMBL .ser file"""

    path_db: str
    """path to annotator DB file to use"""

    varfish_server_compatibility: bool = False
    """
    build output compatible with
    varfish-server v1.2 (Anthenea) and early versions of the v2 (Bollonaster)
    """
