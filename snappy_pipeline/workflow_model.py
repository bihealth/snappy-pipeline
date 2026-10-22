import enum
from typing import Any

from pydantic import ConfigDict, Field

from snappy_pipeline.models import (
    ConfigPath,
    LibraryFolders,
    SnappyModel,
    SnappyStepModel,
)


class DataSetType(enum.StrEnum):
    MATCHED_CANCER = "matched_cancer"
    GERMLINE_VARIANTS = "germline_variants"
    GENERIC = "generic"


class NamingScheme(enum.StrEnum):
    ONLY_SECONDARY_ID = "only_secondary_id"
    SECONDARY_ID_PK = "secondary_id_pk"


class DataSet(LibraryFolders):
    """A sample sheet, and the folders with the reads of its libraries below ``search_paths``.

    The sheet's ``folderName`` names the folder of each library. A library's reads are searched
    only below the search paths of its own data set."""

    file: ConfigPath = ""
    """The sample sheet"""

    search_patterns: list[dict[str, str]] = [
        {"left": r"(?P<readgroup>.+)\.R1\.fastq\.gz", "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz"}
    ]
    """Regular expressions for the paths below a library's folder, by mate, with a
    ``(?P<readgroup>...)`` group that pairs the mates; the default finds ``*.R1.fastq.gz`` and
    ``*.R2.fastq.gz``"""

    type: DataSetType = DataSetType.MATCHED_CANCER
    naming_scheme: NamingScheme = NamingScheme.SECONDARY_ID_PK
    is_background: bool = False
    """Libraries that are no task's results, such as a reference cohort; mapping, panels of
    normals and the SV and gCNV callers use them, narrowed by a task's ``background_selection``"""
    sodar_uuid: str | None = None
    sodar_title: str | None = None
    pedigree_field: str | None = None


class TaskModel(SnappyModel):
    step: str
    """The workflow to use (e.g., 'ngs_mapping', 'variant_calling')"""

    name: str
    """The unique instance name for this pipeline execution (e.g., 'mapping_stringent')"""

    config: dict[str, Any] = Field(default_factory=dict)
    """The raw configuration dictionary for the step. Validated downstream."""


class Resources(SnappyModel):
    """Project-wide rules for the memory and runtime that jobs request."""

    retry_factor: float = 1.5
    """Memory and runtime grow by this factor with each retry of a failed job"""


class ConfigModel(SnappyStepModel):
    model_config = ConfigDict(
        extra="allow",
        use_attribute_docstrings=True,
        use_enum_values=True,
    )

    tasks: list[TaskModel]
    data_sets: dict[str, DataSet]
    pathvars: dict[str, str] = {}
    """Named roots for paths, used as ``<name>`` in any path of the configuration, with the
    defaults and semantics of Snakemake's ``pathvars``"""
    resources: Resources = Resources()
