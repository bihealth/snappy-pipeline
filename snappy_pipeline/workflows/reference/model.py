from enum import StrEnum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Source(StrEnum):
    Ensembl = "Ensembl"
    NCBI = "NCBI"
    Custom = "Custom"


class DataType(StrEnum):
    dna = "dna"
    cds = "cds"
    cdna = "cdna"
    ncrna = "ncrna"
    pep = "pep"


class Region(SnappyModel):
    name: str
    start: int | None
    end: int | None


class Annotation(SnappyModel):
    reference: list[str] | None = None


class Reference(SnappyModel):
    description: str
    """Description of the reference."""

    source: Annotated[Source, EnumField(Source)]
    """Source of the reference."""

    custom_url: str | None = Field(
        None, examples=["file:///path/to/reference.fa", "http://example.com/reference.fa"]
    )
    """URL to custom reference. Only used when source is 'Custom'."""

    species: str = Field(examples=["Homo Sapiens"])
    """Species name."""

    taxon_id: str | int = Field(examples=[9606])
    """Taxon ID."""

    datatype: Annotated[DataType, EnumField(DataType)]
    """Data type of the reference."""

    release: str | int = Field(examples=[112])
    """Release of the reference."""

    build: str | None = Field(None, examples=["GRCh37", "GRCh38"])
    """Build of the reference."""

    branch: str | None = Field(None, examples=["grch37"])
    """Branch of the reference."""

    exclude_contigs: str | None = None
    """Regular expression to exclude contigs with"""

    regions: list[Region] | None = None
    """Regions of the reference."""

    additional_sequences: list[str] | None = None
    """List of local fasta files to add to the reference"""

    annotations: dict[str, Annotation] = {}


class ReferenceModel(SnappyStepModel):
    references: dict[str, Reference] = {
        "GRCh38-foo": Reference(
            source="Ensembl", species="Homo Sapiens", taxon_id=9606, datatype="dna", release=112
        )
    }
