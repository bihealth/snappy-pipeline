from typing import Annotated

from pydantic import BaseModel

from snappy_pipeline.models import ConfigPath
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


class WgsSvExportExternalDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"sv"})),
        ExpectedPathSchema(ExpectedExternalVcf),
    ]
    """``external_data`` task with the SV VCF of each library."""

    reference: Reference


class WgsSvExportExternal(SnappyStepModel):
    depends_on: WgsSvExportExternalDependsOn

    merge_vcf_flag: bool = False
    """true if pedigree VCFs still need merging (not recommended)."""

    merge_option: str = "id"
    """How to merge VCF, used in `bcftools --merge` call."""

    release: str = "GRCh37"

    path_refseq_ser: ConfigPath
    """path to RefSeq .ser file"""

    path_ensembl_ser: ConfigPath
    """path to ENSEMBL .ser file"""

    path_db: ConfigPath
    """path to annotator DB file to use"""

    varfish_server_compatibility: bool = False
    """
    build output compatible with
    varfish-server v1.2 (Anthenea) and early versions of the v2 (Bollonaster)
    """
