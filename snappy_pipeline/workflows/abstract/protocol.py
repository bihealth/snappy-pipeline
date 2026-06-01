from dataclasses import dataclass, field
from enum import StrEnum

from pydantic import BaseModel


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


@dataclass(frozen=True)
class ExpectedPathSchema:
    """Annotated metadata wrapper for expected upstream output-path schema."""

    schema: type[BaseModel]
