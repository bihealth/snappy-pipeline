from dataclasses import dataclass, field
from enum import StrEnum

from pydantic import BaseModel

#: Reserved ``depends_on.reads`` value: search FASTQs in the data sets' search paths.
DATA_SETS = "data_sets"


class DataType(StrEnum):
    RAW = "raw"
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
