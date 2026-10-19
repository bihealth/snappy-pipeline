from dataclasses import dataclass, field
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel

#: Reserved ``depends_on.reads`` value: search FASTQs in the data sets' search paths.
DATA_SETS = "data_sets"


class DataType(StrEnum):
    RAW = "raw"
    REGIONS = "regions"
    DATABASE = "database"
    INDEX = "index"
    ALIGNMENTS = "alignments"
    VARIANTS = "variants"
    EXPRESSION = "expression"
    TABULAR = "tabular"
    MODELS = "models"
    EXPORTS = "exports"
    QC = "qc"


@dataclass(frozen=True)
class DataSignature:
    type: DataType
    tags: frozenset = field(default_factory=frozenset)

    def satisfies(self, requirement: "DataSignature") -> bool:
        """
        Evaluates if the current signature fulfills the requirement.
        The base DataType must match.
        For tags:
        - string: must be present in self.tags
        - tuple: at least one element must be present in self.tags (OR constraint)
        - string starting with '-': element (without '-') must NOT be present in self.tags (NOT constraint)
        """
        if self.type != requirement.type:
            return False

        for req in requirement.tags:
            if isinstance(req, tuple):
                if not any(r in self.tags for r in req):
                    return False
            elif isinstance(req, str) and req.startswith("-"):
                if req[1:] in self.tags:
                    return False
            elif req not in self.tags:
                return False
        return True

    def with_tags(self, *tags: str) -> "DataSignature":
        """Return this signature with ``tags`` added."""
        return DataSignature(self.type, self.tags | frozenset(tags))

    def __str__(self) -> str:
        tags = sorted("|".join(t) if isinstance(t, tuple) else t for t in self.tags)
        return f"{self.type.value} [{', '.join(tags)}]" if tags else self.type.value


def select_signature(
    produced: "tuple[DataSignature, ...]", required: "DataSignature | None"
) -> "DataSignature | None":
    """Return the signature in ``produced`` that a consumer requiring ``required`` reads.

    A task can produce several signatures, e.g. mutect2 calling writes all calls and the PASS
    calls (tagged ``filtered``). When more than one satisfies ``required``, the filtered one wins,
    so consumers get PASS calls unless they require ``-filtered``. Returns ``None`` when nothing
    satisfies ``required``.
    """
    matching = [s for s in produced if required is None or s.satisfies(required)]
    if not matching:
        return None
    return max(matching, key=lambda s: "filtered" in s.tags)


@dataclass(frozen=True)
class ExpectedPathSchema:
    """Annotated metadata wrapper for expected upstream output-path schema."""

    schema: type[BaseModel]


class ExpectedReference(BaseModel):
    """The genome FASTA; tools find its ``.fai`` and ``.dict`` next to it."""

    fasta: str


class ExpectedFeatures(BaseModel):
    """The gene annotation."""

    gtf: str


class ExpectedVcf(BaseModel):
    """A VCF file; tools find its ``.tbi`` next to it."""

    vcf: str


class ExpectedKits(BaseModel):
    """The regions of one enrichment kit as BED files: its targets and its baits (probes)."""

    targets: str
    baits: str


class ExpectedRegions(BaseModel):
    """A BED file of regions."""

    bed: str


#: ``depends_on`` field types for the reference data that many steps read
Reference = Annotated[
    str,
    DataSignature(DataType.RAW, frozenset({"reference", "dna"})),
    ExpectedPathSchema(ExpectedReference),
]
Features = Annotated[
    str, DataSignature(DataType.RAW, frozenset({"features"})), ExpectedPathSchema(ExpectedFeatures)
]
Dbsnp = Annotated[
    str, DataSignature(DataType.VARIANTS, frozenset({"dbsnp"})), ExpectedPathSchema(ExpectedVcf)
]
#: Population allele frequencies, such as the af-only gnomAD VCF of the GATK resource bundle
GermlineResource = Annotated[
    str,
    DataSignature(DataType.VARIANTS, frozenset({"population"})),
    ExpectedPathSchema(ExpectedVcf),
]
#: Common germline variants, such as ``small_exac_common_3.vcf`` of the GATK resource bundle
CommonVariants = Annotated[
    str, DataSignature(DataType.VARIANTS, frozenset({"common"})), ExpectedPathSchema(ExpectedVcf)
]
#: ``depends_on`` field types for regions (``type: regions``), usually ``external_data`` tasks.
#: ``Regions`` accepts any BED file of regions, the others add a tag for what the regions are.
Regions = Annotated[str, DataSignature(DataType.REGIONS), ExpectedPathSchema(ExpectedRegions)]
Exons = Annotated[
    str, DataSignature(DataType.REGIONS, frozenset({"exons"})), ExpectedPathSchema(ExpectedRegions)
]
#: The files per enrichment kit: an ``external_data`` task with ``kits``, read with
#: ``BaseStep.library_kit_paths()``
Kits = Annotated[
    str, DataSignature(DataType.REGIONS, frozenset({"kits"})), ExpectedPathSchema(ExpectedKits)
]
