from typing import Annotated

from pydantic import BaseModel, FilePath

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
)


class ExpectedExternalVcf(BaseModel):
    """Consumer-driven contract: the VCF of a library."""

    vcf: str


class WgsCnvExportExternalDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"cnv"})),
        ExpectedPathSchema(ExpectedExternalVcf),
    ]
    """``external_data`` task with the CNV VCF of each library."""


class WgsCnvExportExternal(SnappyStepModel):
    depends_on: WgsCnvExportExternalDependsOn

    tool_ngs_mapping: str | None = None
    """used to create output file prefix."""

    tool_wgs_cnv_calling: str | None = None
    """used to create output file prefix."""

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str = "id"
    """How to merge VCF, used in `bcftools --merge` call."""

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
