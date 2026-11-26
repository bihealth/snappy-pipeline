# -*- coding: utf-8 -*-
"""Tests for ``snappy scan-files``, which writes the files tables of ``external_data`` tasks."""

from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main
from snappy_pipeline.reads import ReadDiscovery, read_files_tables

SHEET = """\
donor_name\tsex\tlibrary_name\tlibrary_type\tfolder_name
P001\tfemale\tP001-N1-DNA1-WGS1\tWGS\t%(folder1)s
P002\tmale\tP002-N1-DNA1-WGS1\tWGS\t%(folder2)s
"""

CONFIG = """\
project:
  name: scan
sample_sheets: [sheet.tsv]
tasks:
  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]
  - name: bams
    step: external_data
    config:
      produces: {type: alignments, tags: [dna]}
      library_files: [tables/bams.tsv]
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: refs/genome.fa}
"""

BAM_ARGS = [
    "--task",
    "bams",
    "--search-path",
    "aligned",
    "--pattern",
    r"bam=.+\.bam",
    "--pattern",
    r"bai=.+\.bam\.bai",
]


def _touch(root, *paths):
    for path in paths:
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        (root / path).touch()


def _project(tmp_path, folder1="P001", folder2="P002"):
    (tmp_path / "config.yaml").write_text(CONFIG)
    (tmp_path / "sheet.tsv").write_text(SHEET % {"folder1": folder1, "folder2": folder2})
    _touch(
        tmp_path,
        "raw/P001/L2.R1.fastq.gz",
        "raw/P001/L2.R2.fastq.gz",
        "raw/P001/L1.R1.fastq.gz",
        "raw/P001/L1.R2.fastq.gz",
        "raw/batch2/P002/x.R1.fastq.gz",
        "raw/batch2/P002/x.R2.fastq.gz",
        "aligned/P001/P001.bam",
        "aligned/P001/P001.bam.bai",
        "aligned/P002/P002.bam",
        "aligned/P002/P002.bam.bai",
    )


def _scan(tmp_path, *args):
    return CliRunner().invoke(main, ["scan-files", "--directory", str(tmp_path), *args])


def test_scan_files_writes_the_reads_table_with_the_defaults(tmp_path):
    _project(tmp_path)
    result = _scan(tmp_path)
    assert result.exit_code == 0, result.output
    assert f"reads: 3 read groups of 2 libraries written to {tmp_path}/reads.tsv" in result.output
    assert "bams" not in result.output  # other data needs its patterns once

    assert (tmp_path / "reads.tsv").read_text() == (
        "# snappy scan-files --task reads --search-path raw\n"
        "library\treadgroup\tleft\tright\n"
        "P001-N1-DNA1-WGS1\tL1\traw/P001/L1.R1.fastq.gz\traw/P001/L1.R2.fastq.gz\n"
        "P001-N1-DNA1-WGS1\tL2\traw/P001/L2.R1.fastq.gz\traw/P001/L2.R2.fastq.gz\n"
        "P002-N1-DNA1-WGS1\tx\traw/batch2/P002/x.R1.fastq.gz\traw/batch2/P002/x.R2.fastq.gz\n"
    )
    table = read_files_tables([str(tmp_path / "reads.tsv")], reads=True)
    patterns = [
        {"left": r"(?P<readgroup>.+)\.R1\.fastq\.gz", "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz"}
    ]
    for library, folder in (("P001-N1-DNA1-WGS1", "P001"), ("P002-N1-DNA1-WGS1", "P002")):
        found = ReadDiscovery().find([str(tmp_path / "raw")], folder, patterns)
        assert table.read_groups(library) == found


def test_scan_files_takes_the_patterns_once_and_the_table_remembers_them(tmp_path):
    _project(tmp_path)
    result = _scan(tmp_path, *BAM_ARGS)
    assert result.exit_code == 0, result.output
    assert f"bams: 2 rows of 2 libraries written to {tmp_path}/tables/bams.tsv" in result.output
    # Paths outside the table's directory are absolute
    assert (tmp_path / "tables/bams.tsv").read_text() == (
        f"# snappy scan-files --task bams --search-path {tmp_path}/aligned "
        "'--pattern=bam=.+\\.bam' '--pattern=bai=.+\\.bam\\.bai'\n"
        "library\tbam\tbai\n"
        f"P001-N1-DNA1-WGS1\t{tmp_path}/aligned/P001/P001.bam\t{tmp_path}/aligned/P001/P001.bam.bai\n"
        f"P002-N1-DNA1-WGS1\t{tmp_path}/aligned/P002/P002.bam\t{tmp_path}/aligned/P002/P002.bam.bai\n"
    )
    before = (tmp_path / "tables/bams.tsv").read_text()
    _touch(tmp_path, "aligned/P003/P003.bam")
    result = _scan(tmp_path)  # rescans both tables the way they record
    assert result.exit_code == 0, result.output
    assert "reads: 3 read groups" in result.output and "bams: 2 rows" in result.output
    assert (tmp_path / "tables/bams.tsv").read_text() == before


def test_scan_files_needs_patterns_for_other_data_without_a_recorded_scan(tmp_path):
    _project(tmp_path)
    result = _scan(tmp_path, "--task", "bams")
    assert result.exit_code == 1
    assert (
        "Task 'bams': its table records no scan; give --search-path and --pattern" in result.output
    )
    result = _scan(tmp_path, "--pattern", "bam")
    assert result.exit_code == 1
    assert "--search-path and --pattern apply to one task; name it with --task" in result.output
    result = _scan(tmp_path, "--task", "bams", "--search-path", "aligned", "--pattern", "bam")
    assert result.exit_code == 1 and "--pattern takes key=regex, not 'bam'" in result.output


def test_scan_files_reports_libraries_without_files_and_writes_the_rest(tmp_path):
    _project(tmp_path)
    for path in (tmp_path / "raw/batch2/P002").iterdir():
        path.unlink()
    result = _scan(tmp_path, "--task", "reads")
    assert result.exit_code == 1
    assert (
        f"reads: P002-N1-DNA1-WGS1: Found no reads of 'P002' below {tmp_path}/raw" in result.output
    )
    assert (tmp_path / "reads.tsv").read_text().count("\n") == 4  # comment, header, two groups
    assert not (tmp_path / "tables").exists()


def test_scan_files_rejects_two_libraries_in_one_folder(tmp_path):
    # The first library's folder is the second library's name, whose folder defaults to it
    _project(tmp_path, folder1="P002-N1-DNA1-WGS1", folder2="")
    result = _scan(tmp_path, "--task", "reads")
    assert result.exit_code == 1
    assert (
        "reads: the libraries P001-N1-DNA1-WGS1, P002-N1-DNA1-WGS1 share the folder "
        "'P002-N1-DNA1-WGS1'" in result.output
    )
    assert (tmp_path / "reads.tsv").read_text().endswith("library\treadgroup\tleft\tright\n")


def test_scan_files_names_tasks_without_tables(tmp_path):
    _project(tmp_path)
    result = _scan(tmp_path, "--task", "genome")
    assert result.exit_code == 1
    assert "Task 'genome' is no external_data task with library_files" in result.output


def test_the_defaults_cover_flat_layouts_and_the_common_namings(tmp_path):
    (tmp_path / "config.yaml").write_text(CONFIG)
    (tmp_path / "sheet.tsv").write_text(SHEET % {"folder1": "P001", "folder2": "P002"})
    _touch(
        tmp_path,
        "raw/P001_S3_L001_R1_001.fastq.gz",  # bcl2fastq, flat below raw
        "raw/P001_S3_L001_R2_001.fastq.gz",
        "raw/P001_S3_L002_R1_001.fastq.gz",
        "raw/P001_S3_L002_R2_001.fastq.gz",
        "raw/batch2/P002_1.fq.gz",  # SRA style, in a batch folder
        "raw/batch2/P002_2.fq.gz",
    )
    result = _scan(tmp_path, "--task", "reads")
    assert result.exit_code == 0, result.output
    lines = (tmp_path / "reads.tsv").read_text().splitlines()
    assert (
        lines[0] == "# snappy scan-files --task reads --search-path raw"
    )  # the defaults are not recorded
    assert lines[2:] == [
        "P001-N1-DNA1-WGS1\tS3_L001\traw/P001_S3_L001_R1_001.fastq.gz\traw/P001_S3_L001_R2_001.fastq.gz",
        "P001-N1-DNA1-WGS1\tS3_L002\traw/P001_S3_L002_R1_001.fastq.gz\traw/P001_S3_L002_R2_001.fastq.gz",
        "P002-N1-DNA1-WGS1\tP002\traw/batch2/P002_1.fq.gz\traw/batch2/P002_2.fq.gz",
    ]


def test_a_file_that_two_libraries_claim_is_a_problem(tmp_path):
    # The second library's folder name starts like the first one's, so S1_A_... files fit both
    _project(tmp_path, folder1="S1", folder2="S1_A")
    _touch(tmp_path, "raw/S1_A_L1.R1.fastq.gz", "raw/S1_A_L1.R2.fastq.gz")
    result = _scan(tmp_path, "--task", "reads")
    assert result.exit_code == 1
    assert (
        f"reads: P002-N1-DNA1-WGS1: {tmp_path}/raw/S1_A_L1.R1.fastq.gz belongs to library "
        "'P001-N1-DNA1-WGS1' as well" in result.output
    )
