from typing import Any

from pydantic import ConfigDict, Field

from snappy_pipeline.models import ConfigPath, SnappyModel


class TaskModel(SnappyModel):
    step: str
    """The workflow to use (e.g., 'ngs_mapping', 'variant_calling')"""

    name: str
    """The unique instance name for this pipeline execution (e.g., 'mapping_stringent')"""

    config: dict[str, Any] = Field(default_factory=dict)
    """The raw configuration dictionary for the step. Validated downstream."""


class Resources(SnappyModel):
    """Project-wide rules for the memory, runtime and temporary directory that jobs request."""

    retry_factor: float = 1.5
    """Memory and runtime grow by this factor with each retry of a failed job"""

    tmpdir: str = ""
    """The temporary directory of the jobs (their ``tmpdir`` resource), with ``$VARIABLES`` of the
    environment expanded; unset: ``$TMPDIR``, else the system's temporary directory"""


class ProjectInfo(SnappyModel):
    """The name of the project, and its identifier in an external data management system, such
    as the SODAR project UUID, which tracks where the raw data lives."""

    name: str
    external_id: str = ""
    description: str = ""


class ConfigModel(SnappyModel):
    """The project configuration, ``config.yaml``; ``snappy run`` adds ``task`` and ``frozen``."""

    model_config = ConfigDict(
        extra="allow",
        use_attribute_docstrings=True,
        use_enum_values=True,
    )

    project: ProjectInfo
    pathvars: dict[str, str] = {}
    """Named roots for paths, used as ``<name>`` in any path of the configuration, with the
    defaults and semantics of Snakemake's ``pathvars``"""
    sample_sheets: list[ConfigPath] = Field(default=["samplesheet.tsv"], validate_default=True)
    """The sample sheets (``snappy_pipeline.sample_sheet``), read together; ``[]`` for a project
    without libraries, such as one that only builds reference data"""
    resources: Resources = Resources()
    tasks: list[TaskModel]
