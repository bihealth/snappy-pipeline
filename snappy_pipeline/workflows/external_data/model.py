import re

from pydantic import model_validator

from snappy_pipeline.models import ConfigPath, LibraryFolders, SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataType


class Produces(SnappyModel):
    """The DataSignature of the provided data."""

    type: DataType
    tags: list[str] = []


#: Kit names become path components (gCNV models are per kit)
KIT_NAME = re.compile(r"[A-Za-z0-9_.-]+")


class ExternalData(SnappyStepModel, LibraryFolders):
    """Existing data from outside the project: project-wide files, files per library or files
    per enrichment kit.

    Paths are absolute or relative to the project directory. Per-library data is in a folder
    named like each library below ``search_paths``, found like the reads of a data set.
    """

    produces: Produces
    """What the data is, e.g. ``{type: variants, tags: [germline, snv, indel]}``."""

    files: dict[str, ConfigPath | list[ConfigPath]] = {}
    """Project-wide files by output key, e.g. ``{fasta: /refs/GRCh38.fa}``; a key may name a list
    of files."""

    kits: dict[str, dict[str, ConfigPath]] = {}
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
        return self

    def patterns_find_reads(self) -> bool:
        return self.produces.type == DataType.RAW
