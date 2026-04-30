from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class VarfishExportDependsOn(SnappyModel):
    ngs_mapping: str = "ngs_mapping"
    variant_calling: str = "variant_calling"
    sv_calling_targeted: str = "sv_calling_targeted"
    sv_calling_wgs: str = "sv_calling_wgs"


class VarfishExport(SnappyStepModel):
    """Configuration of the input path enables export from the corresponding pipeline step."""

    """Used output of ngs_mapping is alignment quality control data"""

    """Used output of variant_calling is variant calls"""

    """Used output of targeted SV calling is variant calls"""

    """Used output of WGS SV calling is variant calls"""

    depends_on: VarfishExportDependsOn = Field(default_factory=VarfishExportDependsOn)

    # Optionally, you can override the exported mappers and variant callers by setting
    # the following variables.
    tools_ngs_mapping: list[str] = []
    """Can be used to override the exported mappers and variant callers"""

    tools_variant_calling: list[str] = []
    """Can be used to override the exported mappers and variant callers"""

    tools_sv_calling_targeted: list[str] = []
    """Can be used to override the exported mappers and variant callers"""

    tools_sv_calling_wgs: list[str] = []
    """Can be used to override the exported mappers and variant callers"""

    # The following configuration is used for parameterizing the output itself.
    release: str = "GRCh37"
    """The release of the genome reference that data has been aligned to."""

    path_exon_bed: str
    """Path to BED file with exons; used for reducing data to near-exon small variants."""

    path_mehari_db: str
    """Path to mehari database."""
