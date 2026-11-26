# -*- coding: utf-8 -*-
"""Tests for ``snappy scan-files``, which writes the files tables of ``external_data`` tasks."""

from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main
from snappy_pipeline.reads import ReadDiscovery, read_files_tables

SHEET = """\
[Data]
patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tfolderName\thpoTerms
P001\t.\t.\tF\tY\tWGS\t%(folder1)s\t.
P002\t.\t.\tM\tN\tWGS\tP002\t.
"""

PATTERNS = {
    "reads": "{left: '(?P<readgroup>.+)_R1.fastq.gz', right: '(?P<readgroup>.+)_R2.fastq.gz'}",
    "bams": "{bam: '.+\\.bam', bai: '.+\\.bam\\.bai'}",
}

CONFIG = """\
project:
  name: scan
data_sets:
  first:
    file: sheet.tsv
    type: germline_variants
    naming_scheme: only_secondary_id
tasks:
  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]
      scan: {search_paths: [raw], search_patterns: [%(reads)s]}
  - name: bams
    step: external_data
    config:
      produces: {type: alignments, tags: [dna]}
      library_files: [tables/bams.tsv]
      scan: {search_paths: [aligned], search_patterns: [%(bams)s]}
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: refs/genome.fa}
"""


def _touch(root, *paths):
    for path in paths:
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        (root / path).touch()


def _project(tmp_path, folder1="P001"):
    (tmp_path / "config.yaml").write_text(CONFIG % PATTERNS)
    (tmp_path / "sheet.tsv").write_text(SHEET % {"folder1": folder1})
    _touch(
        tmp_path,
        "raw/P001/L2_R1.fastq.gz",
        "raw/P001/L2_R2.fastq.gz",
        "raw/P001/L1_R1.fastq.gz",
        "raw/P001/L1_R2.fastq.gz",
        "raw/batch2/P002/x_R1.fastq.gz",
        "raw/batch2/P002/x_R2.fastq.gz",
        "aligned/P001/P001.bam",
        "aligned/P001/P001.bam.bai",
        "aligned/P002/P002.bam",
        "aligned/P002/P002.bam.bai",
    )


def test_scan_files_writes_the_tables(tmp_path):
    _project(tmp_path)
    result = CliRunner().invoke(main, ["scan-files", "--directory", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert f"reads: 3 read groups of 2 libraries written to {tmp_path}/reads.tsv" in result.output
    assert f"bams: 2 rows of 2 libraries written to {tmp_path}/tables/bams.tsv" in result.output

    assert (tmp_path / "reads.tsv").read_text() == (
        "library\treadgroup\tleft\tright\n"
        "P001-N1-DNA1-WGS1\tL1\traw/P001/L1_R1.fastq.gz\traw/P001/L1_R2.fastq.gz\n"
        "P001-N1-DNA1-WGS1\tL2\traw/P001/L2_R1.fastq.gz\traw/P001/L2_R2.fastq.gz\n"
        "P002-N1-DNA1-WGS1\tx\traw/batch2/P002/x_R1.fastq.gz\traw/batch2/P002/x_R2.fastq.gz\n"
    )
    # Paths outside the table's directory are absolute
    assert (tmp_path / "tables/bams.tsv").read_text() == (
        "library\tbam\tbai\n"
        f"P001-N1-DNA1-WGS1\t{tmp_path}/aligned/P001/P001.bam\t{tmp_path}/aligned/P001/P001.bam.bai\n"
        f"P002-N1-DNA1-WGS1\t{tmp_path}/aligned/P002/P002.bam\t{tmp_path}/aligned/P002/P002.bam.bai\n"
    )
    table = read_files_tables([str(tmp_path / "reads.tsv")], reads=True)
    patterns = [
        {"left": r"(?P<readgroup>.+)_R1.fastq.gz", "right": r"(?P<readgroup>.+)_R2.fastq.gz"}
    ]
    for library, folder in (("P001-N1-DNA1-WGS1", "P001"), ("P002-N1-DNA1-WGS1", "P002")):
        found = ReadDiscovery().find([str(tmp_path / "raw")], folder, patterns)
        assert table.read_groups(library) == found


def test_scan_files_reports_libraries_without_files_and_writes_the_rest(tmp_path):
    _project(tmp_path)
    for path in (tmp_path / "raw/batch2/P002").iterdir():
        path.unlink()
    result = CliRunner().invoke(
        main, ["scan-files", "--directory", str(tmp_path), "--task", "reads"]
    )
    assert result.exit_code == 1
    assert (
        f"reads: P002-N1-DNA1-WGS1: Found no reads of 'P002' below {tmp_path}/raw" in result.output
    )
    assert (tmp_path / "reads.tsv").read_text().count("\n") == 3
    assert not (tmp_path / "tables").exists()


def test_scan_files_rejects_two_libraries_in_one_folder(tmp_path):
    _project(tmp_path, folder1="P002")
    result = CliRunner().invoke(
        main, ["scan-files", "--directory", str(tmp_path), "--task", "reads"]
    )
    assert result.exit_code == 1
    assert (
        "reads: the libraries P001-N1-DNA1-WGS1, P002-N1-DNA1-WGS1 share the folder 'P002'"
        in result.output
    )
    assert (tmp_path / "reads.tsv").read_text() == "library\treadgroup\tleft\tright\n"


def test_scan_files_names_tasks_without_a_scan(tmp_path):
    _project(tmp_path)
    result = CliRunner().invoke(
        main, ["scan-files", "--directory", str(tmp_path), "--task", "genome"]
    )
    assert result.exit_code == 1
    assert "Task 'genome' is no external_data task with a scan section" in result.output
