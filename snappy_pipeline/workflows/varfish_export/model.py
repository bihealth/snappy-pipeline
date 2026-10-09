from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants


class VarfishExportDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = "ngs_mapping"
    """Used output of ngs_mapping is alignment quality control data"""

    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        ExpectedPathSchema(ExpectedGermlineVariants),
    ] = "variant_calling"
    """Used output of variant_calling is variant calls"""

    structural_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"})),
    ] = ""
    """Optional ``sv_calling_targeted`` or ``sv_calling_wgs`` task whose calls are exported"""


class VarfishExport(SnappyStepModel):
    """Configuration of the input path enables export from the corresponding pipeline step."""

    depends_on: VarfishExportDependsOn = Field(default_factory=VarfishExportDependsOn)

    # Optionally, you can override the exported mappers and variant callers by setting
    # the following variables.

    # The following configuration is used for parameterizing the output itself.
    release: str = "GRCh37"
    """The release of the genome reference that data has been aligned to."""

    path_exon_bed: str
    """Path to BED file with exons; used for reducing data to near-exon small variants."""

    path_mehari_db: str
    """Path to mehari database."""
