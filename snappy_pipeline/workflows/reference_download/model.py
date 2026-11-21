import re
import warnings
from enum import StrEnum
from typing import Annotated

from pydantic import AnyUrl, BaseModel, BeforeValidator, TypeAdapter, model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class ExpectedReferenceDownloadFiles(BaseModel):
    """Consumer-driven contract for downloaded reference FASTA outputs."""

    fasta: str


class Source(StrEnum):
    ensembl = "ensembl"
    refseq = "refseq"
    ucsc = "ucsc"


class Molecule(StrEnum):
    dna = "dna"
    rna = "rna"
    protein = "protein"


class EnsemblDataType(StrEnum):
    dna = "dna"
    cdna = "cdna"
    cds = "cds"
    ncrna = "ncrna"
    pep = "pep"


#: An NCBI assembly accession, whose digits name the directories of its files
ASSEMBLY_ACCESSION = re.compile(r"GC[AF]_\d{9}\.\d+")


def validate_url_or_empty(v: str) -> str:
    """Validate that the string is either empty or a syntactically valid URL."""
    if not v:
        return ""

    ta = TypeAdapter(AnyUrl)
    return str(ta.validate_python(v))


class DownloadCommon(SnappyModel):
    url: Annotated[str, BeforeValidator(validate_url_or_empty)] = ""
    """Optional explicit URL override. When set, source-specific URL construction is skipped."""

    contigs: list[str] = []
    """Optional explicit list of contigs to keep."""

    contigs_regex: str = ""
    """Optional regular expression to keep contigs by name."""

    molecule: Molecule = Molecule.dna
    """Molecule of the downloaded sequences, a tag of the task's output."""

    @model_validator(mode="after")
    def validate_url_exclusivity(self) -> "DownloadCommon":
        """Warn the user if they provided source-specific parameters while url is set."""
        if self.url:
            # Downstream processing filters (molecule, contigs, etc.) are still allowed
            allowed_with_url = {"url", "molecule", "contigs", "contigs_regex"}
            conflicting_fields = self.model_fields_set - allowed_with_url

            if conflicting_fields:
                warnings.warn(
                    f"A custom 'url' is specified ('{self.url}'), which overrides standard construction. "
                    f"The following configuration parameters will be ignored: {', '.join(sorted(conflicting_fields))}",
                    UserWarning,
                    stacklevel=2,
                )
        return self


class Ensembl(DownloadCommon):
    species: str = "homo_sapiens"
    """Species identifier, usually lowercase with underscore for Ensembl."""

    build: str = "GRCh38"
    """Genome build identifier (e.g., GRCh38, GRCh37)."""

    assembly: str = "primary_assembly"
    """Assembly subset used when ``subset`` is empty."""

    version: str = "113"
    subset: str = "primary_assembly"
    file_name: str = ""
    datatype: EnsemblDataType = EnsemblDataType.dna
    chromosome: list[str] = []
    branch: str = ""

    @model_validator(mode="after")
    def validate_datatype_molecule(self) -> "Ensembl":
        # Skip subclass-specific validation rules if a custom URL override is used
        if self.url:
            return self

        if self.datatype == EnsemblDataType.dna and self.molecule != Molecule.dna:
            raise ValueError("Ensembl datatype 'dna' requires molecule='dna'")
        if (
            self.datatype
            in {
                EnsemblDataType.cdna,
                EnsemblDataType.cds,
                EnsemblDataType.ncrna,
            }
            and self.molecule != Molecule.rna
        ):
            raise ValueError("Ensembl datatype 'cdna/cds/ncrna' requires molecule='rna'")
        if self.datatype == EnsemblDataType.pep and self.molecule != Molecule.protein:
            raise ValueError("Ensembl datatype 'pep' requires molecule='protein'")
        if self.chromosome and self.datatype != EnsemblDataType.dna:
            raise ValueError("Ensembl chromosome selection is only supported for datatype='dna'")
        return self


class Refseq(DownloadCommon):
    assembly_accession: str = ""
    """Accession of the NCBI assembly, such as ``GCF_000001405.40``; required without ``url``."""

    assembly_name: str = ""
    """Name of the NCBI assembly, such as ``GRCh38.p14``; required without ``url``."""

    file_name: str = ""
    """File in the assembly's directory; empty: ``<assembly_accession>_<assembly_name>_genomic.fna.gz``."""

    @model_validator(mode="after")
    def validate_assembly(self) -> "Refseq":
        if self.url:
            return self
        if not self.assembly_accession or not self.assembly_name:
            raise ValueError("refseq needs assembly_accession and assembly_name, or a url")
        if not ASSEMBLY_ACCESSION.fullmatch(self.assembly_accession):
            raise ValueError(
                f"refseq assembly_accession {self.assembly_accession!r} is no NCBI assembly "
                "accession such as GCF_000001405.40"
            )
        return self


class Ucsc(DownloadCommon):
    db: str = "hg38"
    """UCSC genome database, such as ``hg38`` or ``mm39``."""

    file_name: str = ""
    """File in the database's ``bigZips`` directory; empty: ``<db>.fa.gz``."""


class ReferenceDownload(SnappyStepModel):
    source: Source
    """Source of the download; its section selects the file."""

    ensembl: Ensembl = Ensembl()
    """The Ensembl FASTA file, by species, build, release (``version``) and ``datatype``."""

    refseq: Refseq | None = None
    """The genomic FASTA file of an NCBI assembly; required for ``source: refseq``."""

    ucsc: Ucsc = Ucsc()
    """The FASTA file of a UCSC genome database."""

    @model_validator(mode="after")
    def validate_source_section(self) -> "ReferenceDownload":
        if getattr(self, self.source) is None:
            raise ValueError(
                f"source={self.source} requires explicit '{self.source}' config section"
            )
        return self
