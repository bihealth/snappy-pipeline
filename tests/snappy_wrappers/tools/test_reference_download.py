"""Tests for ``snappy_wrappers.tools.reference_download``"""

import gzip

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


def test_refseq_url():
    url = reference_download.refseq_url("", "GCF_000001405.40", "GRCh38.p14", "")
    assert url == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
        "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    )


def test_ucsc_url_falls_back_to_the_build():
    url = reference_download.ucsc_url("", "", "hg19", "")
    assert url == "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"


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
