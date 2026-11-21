"""Tests for ``reference_download``: a task produces the tags of its own source and molecule."""

import pytest
from pydantic import ValidationError

from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.reference_download import ReferenceDownloadWorkflow
from snappy_pipeline.workflows.reference_download.model import ReferenceDownload

DNA_REFERENCE = DataSignature(DataType.RAW, frozenset({"reference", "dna"}))


def test_source_is_required():
    with pytest.raises(ValidationError, match="source"):
        ReferenceDownload()


def test_the_task_produces_its_source_and_molecule():
    config = ReferenceDownload(source="ucsc", ucsc={"molecule": "rna"})

    produced = ReferenceDownloadWorkflow.task_produces(config, {})

    assert produced == (DataSignature(DataType.RAW, frozenset({"ucsc", "reference", "rna"})),)
    assert set(produced) <= set(ReferenceDownloadWorkflow.produces)


def test_a_protein_download_is_no_dna_reference():
    config = ReferenceDownload(source="ensembl", ensembl={"datatype": "pep", "molecule": "protein"})

    (produced,) = ReferenceDownloadWorkflow.task_produces(config, {})

    assert not produced.satisfies(DNA_REFERENCE)
