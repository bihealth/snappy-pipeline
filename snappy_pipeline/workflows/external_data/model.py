import re

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.reads import compile_search_pattern
from snappy_pipeline.workflows.abstract.protocol import DataType


class Produces(SnappyModel):
    """The DataSignature of the provided data."""

    type: DataType
    tags: list[str] = []


#: Kit names become path components (gCNV models are per kit)
KIT_NAME = re.compile(r"[A-Za-z0-9_.-]+")


class ExternalData(SnappyStepModel):
    """Existing data from outside the project: project-wide files, files per library or files
    per enrichment kit.

    Paths are absolute or relative to the project directory.
    """

    produces: Produces
    """What the data is, e.g. ``{type: variants, tags: [germline, snv, indel]}``."""

    files: dict[str, str | list[str]] = {}
    """Project-wide files by output key, e.g. ``{fasta: /refs/GRCh38.fa}``; a key may name a list
    of files."""

    search_paths: list[str] = []
    """Per-library data: directories that contain a folder named like each library."""

    search_patterns: list[dict[str, str]] = []
    """Regular expressions for the paths below a library folder, by output key, e.g.
    ``{vcf: '.+\\.vcf\\.gz', vcf_tbi: '.+\\.vcf\\.gz\\.tbi'}``. Reads (``type: raw``) need a
    ``left`` key and a ``(?P<readgroup>...)`` group, which pairs the mates."""

    kits: dict[str, dict[str, str]] = {}
    """Files per enrichment kit, by kit name and output key, e.g.
    ``{MedExome: {targets: /refs/MedExome.targets.bed}}``. A targeted library uses the files of
    the kit named exactly like its ``libraryKit`` in the sample sheet. Every kit has the same
    keys."""

    @model_validator(mode="after")
    def check_mode(self):
        if sum(map(bool, (self.files, self.search_paths, self.kits))) != 1:
            raise ValueError(
                "Set one of files (project-wide data), search_paths with search_patterns "
                "(per-library data) or kits (per-kit data)"
            )
        if bad := sorted(name for name in self.kits if not KIT_NAME.fullmatch(name)):
            raise ValueError(f"Kit names may only contain A-Z, a-z, 0-9, '_', '.' and '-': {bad}")
        if len({frozenset(files) for files in self.kits.values()}) > 1:
            raise ValueError("Every kit must have the same keys")
        if self.search_paths and not self.search_patterns:
            raise ValueError("search_paths need search_patterns")
        for pattern in self.search_patterns:
            compile_search_pattern(pattern, reads=self.produces.type == DataType.RAW)
        return self
