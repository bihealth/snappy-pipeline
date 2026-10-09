from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.reads import compile_search_pattern
from snappy_pipeline.workflows.abstract.protocol import DataType


class Produces(SnappyModel):
    """The DataSignature of the provided data."""

    type: DataType
    tags: list[str] = []


class ExternalData(SnappyStepModel):
    """Existing data from outside the project, per library or as project-wide files.

    Paths are absolute or relative to the project directory.
    """

    produces: Produces
    """What the data is, e.g. ``{type: variants, tags: [germline, snv, indel]}``."""

    files: dict[str, str] = {}
    """Project-wide files by output key, e.g. ``{fasta: /refs/GRCh38.fa, fai: ...}``."""

    search_paths: list[str] = []
    """Per-library data: directories that contain a folder named like each library."""

    search_patterns: list[dict[str, str]] = []
    """Regular expressions for the paths below a library folder, by output key, e.g.
    ``{vcf: '.+\\.vcf\\.gz', vcf_tbi: '.+\\.vcf\\.gz\\.tbi'}``. Reads (``type: raw``) need a
    ``left`` key and a ``(?P<readgroup>...)`` group, which pairs the mates."""

    @model_validator(mode="after")
    def check_mode(self):
        if bool(self.files) == bool(self.search_paths):
            raise ValueError(
                "Set either files (project-wide data) or search_paths with search_patterns "
                "(per-library data)"
            )
        if self.search_paths and not self.search_patterns:
            raise ValueError("search_paths need search_patterns")
        for pattern in self.search_patterns:
            compile_search_pattern(pattern, reads=self.produces.type == DataType.RAW)
        return self
