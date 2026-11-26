import re

from pydantic import model_validator

from snappy_pipeline.models import ConfigPath, SnappyModel, SnappyStepModel
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

    Paths are absolute or relative to the project directory. Per-library data is listed in the
    files tables of ``library_files``, which ``snappy scan-files`` writes.
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

    library_files: list[ConfigPath] = []
    """Files tables with the files of each library: TSV files with one row per read group
    (``type: raw``: the columns ``library``, ``readgroup``, ``left`` and, for paired-end data,
    ``right``) or one row per library (``library`` and the output keys, such as ``bam`` and
    ``bai``). An empty cell is no value; relative paths are relative to the table's directory.
    ``snappy scan-files`` writes the first table."""

    mixed_se_pe: bool = False
    """For reads: accept read groups without a ``right`` file, so that single-end and paired-end
    data can be mixed"""

    @model_validator(mode="after")
    def check_mode(self):
        if sum(map(bool, (self.files, self.kits, self.library_files))) != 1:
            raise ValueError(
                "Set one of files (project-wide data), library_files (files tables per library) "
                "or kits (per-kit data)"
            )
        if self.mixed_se_pe and not self.patterns_find_reads():
            raise ValueError("mixed_se_pe is only for reads")
        if bad := sorted(name for name in self.kits if not KIT_NAME.fullmatch(name)):
            raise ValueError(f"Kit names may only contain A-Z, a-z, 0-9, '_', '.' and '-': {bad}")
        if len({frozenset(files) for files in self.kits.values()}) > 1:
            raise ValueError("Every kit must have the same keys")
        return self

    def patterns_find_reads(self) -> bool:
        """Whether the files are reads (rows per read group), or one file per key and library."""
        return self.produces.type == DataType.RAW
