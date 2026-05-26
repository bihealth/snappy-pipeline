from pydantic import BaseModel

from snappy_pipeline.models import SnappyStepModel


class ExpectedLinkedRawFastq(BaseModel):
    """Consumer-driven contract for link_in raw FASTQ provider paths."""

    path: str


class LinkIn(SnappyStepModel):
    """Configuration for linking in external pre-processed FASTQ files."""

    path: str
    """Absolute path to an external directory containing pre-processed FASTQ files.

    When set, workflows that declare ``depends_on: {link_in: <this_task_name>}``
    will use this directory to locate FASTQ files instead of crawling ``data_sets``
    search paths.
    """
