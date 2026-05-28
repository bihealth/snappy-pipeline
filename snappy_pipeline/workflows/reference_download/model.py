from enum import StrEnum

from pydantic import BaseModel, model_validator

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


class DownloadCommon(SnappyModel):
    url: str = ""
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


class Ensembl(DownloadCommon):
    version: str = "113"
    subset: str = "primary_assembly"
    file_name: str = ""
    datatype: EnsemblDataType = EnsemblDataType.dna
    chromosome: list[str] = []
    branch: str = ""

    @model_validator(mode="after")
    def validate_datatype_molecule(self):
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
