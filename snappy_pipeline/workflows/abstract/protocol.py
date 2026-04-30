from dataclasses import dataclass, field
from enum import Enum


class DataType(Enum):
    RAW = "raw"
    ALIGNMENTS = "alignments"
    VARIANTS = "variants"


@dataclass(frozen=True)
class DataSignature:
    type: DataType
    tags: frozenset[str] = field(default_factory=frozenset)

    def satisfies(self, requirement: "DataSignature") -> bool:
        """
        Evaluates if the current signature fulfills the requirement.
        The base DataType must match, and the provider must possess
        all tags specified in the requirement.
        """
        if self.type != requirement.type:
            return False
        return requirement.tags.issubset(self.tags)
