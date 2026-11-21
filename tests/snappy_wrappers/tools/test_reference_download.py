"""Tests for ``snappy_wrappers.tools.reference_download``"""

import gzip
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from snappy_wrappers.tools import reference_download

FASTA = ">chr1 first\nAC\nGT\n>chr2\nGG\n>chrM\nTT\n"


@pytest.mark.parametrize(
    "contigs, regex, kept",
    [([], "", ["chr1", "chr2", "chrM"]), (["chr2", "chrM"], "", ["chr2", "chrM"])]
    + [([], "^chr[0-9]$", ["chr1", "chr2"]), (["chr1", "chrM"], "M", ["chrM"])],
)
def test_filter_fasta_keeps_the_selected_contigs(tmp_path, contigs, regex, kept):
    (tmp_path / "in.fa").write_text(FASTA)
    reference_download.filter_fasta(tmp_path / "in.fa", tmp_path / "out.fa", contigs, regex)
    out = (tmp_path / "out.fa").read_text()
    assert [line[1:].split()[0] for line in out.splitlines() if line.startswith(">")] == kept
    assert out.startswith(">chr1 first\nAC\nGT\n") == ("chr1" in kept)


def test_fetch_gunzips_a_file_url(tmp_path):
    with gzip.open(tmp_path / "genome.fa.gz", "wt") as f:
        f.write(FASTA)
    fasta = tmp_path / "out" / "reference.fa"
    reference_download.fetch((tmp_path / "genome.fa.gz").as_uri(), str(fasta), ["chr2"], "")
    assert fasta.read_text() == ">chr2\nGG\n"


NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14"
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"


@pytest.mark.parametrize(
    "url, file_name, expected",
    [
        ("", "", f"{NCBI}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"),
        (
            "",
            "GCF_000001405.40_GRCh38.p14_rna.fna.gz",
            f"{NCBI}/GCF_000001405.40_GRCh38.p14_rna.fna.gz",
        ),
        ("https://example.org/genome.fa.gz", "", "https://example.org/genome.fa.gz"),
    ],
)
def test_refseq_url(url, file_name, expected):
    assert (
        reference_download.refseq_url(url, "GCF_000001405.40", "GRCh38.p14", file_name) == expected
    )


def test_refseq_url_needs_the_assembly():
    with pytest.raises(ValueError, match="assembly_accession and assembly_name"):
        reference_download.refseq_url("", "GCF_000001405.40", "", "")


@pytest.mark.parametrize(
    "url, db, file_name, expected",
    [
        ("", "hg19", "", f"{UCSC}/hg19/bigZips/hg19.fa.gz"),
        ("", "hg38", "hg38.analysisSet.fa.gz", f"{UCSC}/hg38/bigZips/hg38.analysisSet.fa.gz"),
        ("https://example.org/genome.fa.gz", "", "", "https://example.org/genome.fa.gz"),
    ],
)
def test_ucsc_url(url, db, file_name, expected):
    assert reference_download.ucsc_url(url, db, file_name) == expected


def test_ucsc_url_needs_the_db():
    with pytest.raises(ValueError, match="db"):
        reference_download.ucsc_url("", "", "")


@pytest.mark.parametrize(
    "version, build, branch, prefix",
    [
        ("113", "GRCh38", "", "pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38"),
        ("87", "GRCh37", "", "pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37"),
        ("75", "GRCh37", "/x/", "pub/x/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75"),
    ],
)
def test_ensembl_prefix(version, build, branch, prefix):
    url = reference_download.ensembl_prefix("Homo_Sapiens", version, build, "dna", branch)
    assert url == f"https://ftp.ensembl.org/{prefix}"


def test_ensembl_suffixes():
    suffixes = reference_download.ensembl_suffixes
    assert suffixes("dna", [], "", "primary_assembly", "") == [
        "dna.primary_assembly.fa.gz",
        "dna.toplevel.fa.gz",
    ]
    assert suffixes("dna", ["1", "X"], "", "", "") == [
        "dna.chromosome.1.fa.gz",
        "dna.chromosome.X.fa.gz",
    ]
    assert suffixes("cdna", [], "", "", "") == ["cdna.all.fa.gz"]
    assert suffixes("pep", [], "", "", "custom.fa.gz") == ["custom.fa.gz"]


WRAPPERS = Path(__file__).parents[3] / "snappy_wrappers/wrappers/reference_download"

#: Stands in for ``curl``: records the URL, the last argument, and writes ``$FASTA_GZ`` to
#: ``--output``
FAKE_CURL = """#!/usr/bin/env bash
echo "${@: -1}" >> "$RECORD"
while [[ $# -gt 0 ]]; do
    if [[ $1 == --output ]]; then cp "$FASTA_GZ" "$2"; fi
    shift
done
"""

ASSEMBLY = {"assembly_accession": "GCF_000001405.40", "assembly_name": "GRCh38.p14"}


@pytest.mark.parametrize(
    "source, args, url, kept",
    [
        (
            "refseq",
            {"url": "", **ASSEMBLY, "file_name": "", "contigs": ["chr2"], "contigs_regex": ""},
            f"{NCBI}/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
            ">chr2\nGG\n",
        ),
        (
            "ucsc",
            {"url": "", "db": "hg19", "file_name": "", "contigs": [], "contigs_regex": "M"},
            f"{UCSC}/hg19/bigZips/hg19.fa.gz",
            ">chrM\nTT\n",
        ),
        (
            "ucsc",
            {"url": "https://example.org/genome.fa.gz", "db": "hg38", "file_name": ""}
            | {"contigs": [], "contigs_regex": ""},
            "https://example.org/genome.fa.gz",
            FASTA,
        ),
    ],
)
def test_wrapper_downloads_the_url_of_its_params(tmp_path, monkeypatch, source, args, url, kept):
    with gzip.open(tmp_path / "genome.fa.gz", "wt") as f:
        f.write(FASTA)
    (tmp_path / "bin").mkdir()
    curl = tmp_path / "bin" / "curl"
    curl.write_text(FAKE_CURL)
    curl.chmod(0o755)
    record = tmp_path / "record.txt"
    monkeypatch.setenv("PATH", f"{tmp_path / 'bin'}:{os.environ['PATH']}")
    monkeypatch.setenv("RECORD", str(record))
    monkeypatch.setenv("FASTA_GZ", str(tmp_path / "genome.fa.gz"))
    fasta = tmp_path / "out" / "reference.fa"
    snakemake = SimpleNamespace(
        output=SimpleNamespace(fasta=str(fasta)), params=SimpleNamespace(args=args)
    )

    wrapper = WRAPPERS / source / "wrapper.py"
    namespace = {"__name__": "wrapper", "snakemake": snakemake}
    exec(compile(wrapper.read_text(), str(wrapper), "exec"), namespace)
    namespace["main"]()

    assert record.read_text().split() == [url]
    assert fasta.read_text() == kept
