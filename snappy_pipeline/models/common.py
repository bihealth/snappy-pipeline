import enum

from typing import Annotated
from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel


class LibraryKitEntry(SnappyModel):
    """
    Mapping from enrichment kit to target region BED file, for either computing per--target
    region coverage or selecting targeted exons.

    The following will match both the stock IDT library kit and the ones
    with spike-ins seen fromr Yale genomics.  The path above would be
    mapped to the name "default".
      - name: IDT_xGen_V1_0
        pattern: "xGen Exome Research Panel V1\\.0*"
        path: "path/to/targets.bed"
    """

    name: Annotated[str, Field(examples=["IDT_xGen_V1_0"])]

    pattern: Annotated[str, Field(examples=["xGen Exome Research Panel V1\\.0*"])]

    path: Annotated[str, Field(examples=["path/to/targets.bed"])]


class LibraryKit(SnappyModel):
    path_target_interval_list_mapping: list[LibraryKitEntry] = []
    """Connects sample-based library kit in sample sheets with corresponding bed files"""


class SexValue(enum.StrEnum):
    MALE = "male"
    FEMALE = "female"


class SexOrigin(enum.StrEnum):
    AUTOMATIC = "auto"
    SAMPLESHEET = "samplesheet"
    CONFIG = "config"


class Sex(SnappyModel):
    source: SexOrigin = SexOrigin.AUTOMATIC
    """Where is the sex information taken from? auto (guessed from data), samplesheet or config (single value for the cohort)"""

    path_guess_sex: str = "../guess_sex"
    """
    Path to the ``guess_sex`` step, where the decision files can be found.
    It is generally used when the sex is ``automatic``, but some steps/tools may have their own way to determine sex
    """

    guess_sex_tool: str = "samtools"
    """
    Tool used to guess the sex during the ``guess_sex`` step.
    It is generally used when the sex is ``automatic``, but some steps/tools may have their own way to determine sex
    """

    cohort: SexValue | None = None
    """Sex of the cohort, used when the source is ``config``"""

    column_name: str | None = None
    """Column name of the sex information in the sample sheet, used when the source is ``samplesheet``"""

    @model_validator(mode="after")
    def ensure_valid_values(self):
        if self.source == SexOrigin.CONFIG and not self.cohort:
            raise ValueError("Undefined cohort sex value in configuration file")
        if self.source == SexOrigin.SAMPLESHEET and not self.column_name:
            raise ValueError("Undefined column name for sex information")
        # We leave open the possibility that some step/tool relies on different input when ``automatic`` is set
        # if self.source == SexOrigin.AUTOMATIC and (
        #     not self.path_guess_sex or not self.guess_sex_tool
        # ):
        #     raise ValueError("Path to or tool used by the 'guess_sex' step are missing")
        return self
