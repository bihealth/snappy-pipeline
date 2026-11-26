import enum
from typing import Any

from pydantic import ConfigDict, Field, model_validator

from snappy_pipeline.models import ConfigPath, SnappyModel


class DataSetType(enum.StrEnum):
    MATCHED_CANCER = "matched_cancer"
    GERMLINE_VARIANTS = "germline_variants"
    GENERIC = "generic"


class NamingScheme(enum.StrEnum):
    ONLY_SECONDARY_ID = "only_secondary_id"
    SECONDARY_ID_PK = "secondary_id_pk"


class DataSet(SnappyModel):
    """A sample sheet. The files of its libraries come from ``external_data`` tasks, whose files
    tables ``snappy scan-files`` writes from the sheet's ``folderName`` of each library."""

    file: ConfigPath
    """The sample sheet"""

    type: DataSetType = DataSetType.MATCHED_CANCER
    naming_scheme: NamingScheme = NamingScheme.SECONDARY_ID_PK
    is_background: bool = False
    """Libraries that are no task's results, such as a reference cohort; mapping, panels of
    normals and the SV and gCNV callers use them, narrowed by a task's ``background_selection``"""
    sodar_uuid: str | None = None
    sodar_title: str | None = None
    pedigree_field: str | None = None
    """The custom field of a germline or generic sheet that groups its donors into pedigrees and names them.
    Unset: the parent columns group them, and the ``familyId`` field, else the index donor (``FAM_`` on germline sheets), names them."""

    @model_validator(mode="after")
    def _pedigree_field_needs_pedigrees(self):
        if self.pedigree_field and self.type == DataSetType.MATCHED_CANCER:
            raise ValueError(
                "pedigree_field groups the donors of a germline or generic sheet; a matched_cancer sheet has no pedigrees"
            )
        return self


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
    data_sets: dict[str, DataSet] = {}
    """The sample sheets by name; a project without libraries, such as one that only builds
    reference data, has none"""
    resources: Resources = Resources()
    tasks: list[TaskModel]
