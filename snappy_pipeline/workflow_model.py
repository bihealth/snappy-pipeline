import enum
from typing import Any, TypedDict

from pydantic import ConfigDict, Field, field_validator

from snappy_pipeline.models import ResolvablePath, SnappyModel, SnappyStepModel
from snappy_pipeline.reads import compile_search_pattern


class PathModel(SnappyModel):
    """
    A simple model holding a single resolvable path field.
    Uses ResolvablePath which ensures relative paths are resolved to absolute ones.
    """

    path: ResolvablePath = ""


class StaticDataConfig(SnappyModel):
    reference: PathModel
    cosmic: PathModel | None = None
    dbsnp: PathModel | None = None
    dbnsfp: PathModel | None = None
    features: PathModel | None = None


class SearchPattern(TypedDict):
    """Regular expressions for the paths below a library's folder, by mate.

    Each needs a ``(?P<readgroup>...)`` group, which pairs the mates of one read group (lane).
    """

    left: str
    right: str | None


class DataSetType(enum.StrEnum):
    MATCHED_CANCER = "matched_cancer"
    GERMLINE_VARIANTS = "germline_variants"
    GENERIC = "generic"


class NamingScheme(enum.StrEnum):
    ONLY_SECONDARY_ID = "only_secondary_id"
    SECONDARY_ID_PK = "secondary_id_pk"


class DataSet(SnappyModel):
    file: str = ""
    search_patterns: list[SearchPattern] = [
        SearchPattern(
            left=r"(?P<readgroup>.+)\.R1\.fastq\.gz", right=r"(?P<readgroup>.+)\.R2\.fastq\.gz"
        )
    ]
    search_paths: list[str] = ["../raw"]
    type: DataSetType = DataSetType.MATCHED_CANCER
    naming_scheme: NamingScheme = NamingScheme.SECONDARY_ID_PK
    is_background: bool = False
    mixed_se_pe: bool = False
    sodar_uuid: str | None = None
    sodar_title: str | None = None
    pedigree_field: str | None = None

    @field_validator("search_patterns")
    @classmethod
    def check_search_patterns(cls, patterns: list[SearchPattern]) -> list[SearchPattern]:
        for pattern in patterns:
            compile_search_pattern(pattern)
        return patterns


class TaskModel(SnappyModel):
    step: str
    """The workflow to use (e.g., 'ngs_mapping', 'somatic_variant_calling')"""

    name: str
    """The unique instance name for this pipeline execution (e.g., 'mapping_stringent')"""

    config: dict[str, Any] = Field(default_factory=dict)
    """The raw configuration dictionary for the step. Validated downstream."""


class ConfigModel(SnappyStepModel):
    model_config = ConfigDict(
        extra="allow",
        use_attribute_docstrings=True,
        use_enum_values=True,
    )

    static_data_config: StaticDataConfig
    tasks: list[TaskModel]
    data_sets: dict[str, DataSet]
