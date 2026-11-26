"""Tests for ``reference_download``: a task produces the tags of its own source and molecule."""

import pytest
from pydantic import ValidationError

from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.reference_download import ReferenceDownloadWorkflow
from snappy_pipeline.workflows.reference_download.model import ReferenceDownload

DNA_REFERENCE = DataSignature(DataType.SEQUENCES, frozenset({"reference", "dna"}))


def test_source_is_required():
    with pytest.raises(ValidationError, match="source"):
        ReferenceDownload()


def test_the_task_produces_its_source_and_molecule():
    config = ReferenceDownload(source="ucsc", ucsc={"molecule": "rna"})

    produced = ReferenceDownloadWorkflow.task_produces(config, {})

    assert produced == (DataSignature(DataType.SEQUENCES, frozenset({"ucsc", "reference", "rna"})),)
    assert set(produced) <= set(ReferenceDownloadWorkflow.produces)


def test_a_protein_download_is_no_dna_reference():
    config = ReferenceDownload(source="ensembl", ensembl={"datatype": "pep", "molecule": "protein"})

    (produced,) = ReferenceDownloadWorkflow.task_produces(config, {})

    assert not produced.satisfies(DNA_REFERENCE)


REFSEQ = {"assembly_accession": "GCF_000001405.40", "assembly_name": "GRCh38.p14"}


def test_refseq_needs_its_section():
    with pytest.raises(ValidationError, match="requires explicit 'refseq' config section"):
        ReferenceDownload(source="refseq")


@pytest.mark.parametrize("missing", ["assembly_accession", "assembly_name"])
def test_refseq_needs_the_assembly_without_url(missing):
    refseq = {key: value for key, value in REFSEQ.items() if key != missing}
    with pytest.raises(ValidationError, match="needs assembly_accession and assembly_name"):
        ReferenceDownload(source="refseq", refseq=refseq)


def test_refseq_rejects_a_malformed_accession():
    with pytest.raises(ValidationError, match="no NCBI assembly accession"):
        ReferenceDownload(source="refseq", refseq={**REFSEQ, "assembly_accession": "GRCh38"})


def test_a_url_replaces_the_refseq_assembly():
    config = ReferenceDownload(source="refseq", refseq={"url": "https://example.org/genome.fa.gz"})
    assert config.refseq.assembly_accession == ""


def test_ucsc_has_no_build():
    with pytest.raises(ValidationError, match="build"):
        ReferenceDownload(source="ucsc", ucsc={"build": "hg19"})


@pytest.mark.parametrize(
    "source, section",
    [("ensembl", {"species": "mus_musculus"}), ("refseq", REFSEQ), ("ucsc", {"db": "hg19"})],
)
def test_a_url_warns_about_the_fields_it_overrides(source, section):
    with pytest.warns(UserWarning, match=f"ignored: {', '.join(sorted(section))}"):
        ReferenceDownload(source=source, **{source: {"url": "https://example.org/x.fa", **section}})
