import warnings
from enum import StrEnum
from typing import Annotated

from pydantic import AnyUrl, BaseModel, BeforeValidator, TypeAdapter, model_validator

from snappy_pipeline.models import ResolvablePath, SnappyModel, SnappyStepModel


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


def validate_url_or_empty(v: str) -> str:
    """Validate that the string is either empty or a syntactically valid URL."""
    if not v:
        return ""

    ta = TypeAdapter(AnyUrl)
    return str(ta.validate_python(v))


class DownloadCommon(SnappyModel):
    url: Annotated[str, BeforeValidator(validate_url_or_empty)] = ""
    """Optional explicit URL override. When set, source-specific URL construction is skipped."""

    species: str = "homo_sapiens"
    """Species identifier, usually lowercase with underscore for Ensembl."""

    build: str = "GRCh38"
    """Genome build identifier (e.g., GRCh38, hg38)."""

    assembly: str = "primary_assembly"
    """Assembly/subset indicator, source-specific semantics."""

    version: str = ""
    """Release/version string, source-specific semantics (e.g., Ensembl release)."""

    contigs: list[str] = []
    """Optional explicit list of contigs to keep."""

    contigs_regex: str = ""
    """Optional regular expression to keep contigs by name."""

    molecule: Molecule = Molecule.dna
    """Molecule class of the downloaded reference payload."""

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
    assembly_name: str = ""
    file_name: str = ""


class Ucsc(DownloadCommon):
    db: str = "hg38"
    file_name: str = ""


class ReferenceDownload(SnappyStepModel):
    source: Source = Source.ensembl

    path_output_fasta: ResolvablePath = ""
    """Optional output FASTA path override. When empty, standard task-local output path is used."""

    ensembl: Ensembl = Ensembl()
    refseq: Refseq = Refseq()
    ucsc: Ucsc = Ucsc()
