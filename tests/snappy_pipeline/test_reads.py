# -*- coding: utf-8 -*-
"""Tests for the central FASTQ discovery in ``snappy_pipeline.reads``."""

import pytest

from snappy_pipeline.reads import ReadDiscovery, compile_search_pattern, read_files_tables

PAIRED = {"left": r"(?P<readgroup>.+)\.R1\.fastq\.gz", "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz"}


def _touch(root, *paths):
    for path in paths:
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        (root / path).touch()


def test_finds_paired_read_groups_in_left_path_order(tmp_path):
    _touch(
        tmp_path,
        "batch/P001-N1/L2/x.R1.fastq.gz",
        "batch/P001-N1/L2/x.R2.fastq.gz",
        "batch/P001-N1/L1/x.R1.fastq.gz",
        "batch/P001-N1/L1/x.R2.fastq.gz",
        "batch/P001-N1/L1/x.R1.fastq.gz.md5",
        "batch/P002-N1/L1/x.R1.fastq.gz",
    )
    groups = ReadDiscovery().find([str(tmp_path)], "P001-N1", [PAIRED])

    assert [g.name for g in groups] == ["L1/x", "L2/x"]
    assert groups[0].left == str(tmp_path / "batch/P001-N1/L1/x.R1.fastq.gz")
    assert groups[0].right == str(tmp_path / "batch/P001-N1/L1/x.R2.fastq.gz")
    assert groups[0].relpaths == {"left": "L1/x.R1.fastq.gz", "right": "L1/x.R2.fastq.gz"}


def test_rebased_keeps_the_paths_below_the_library_folder(tmp_path):
    _touch(tmp_path, "P1/L1/x.R1.fastq.gz", "P1/L1/x.R2.fastq.gz")
    (group,) = ReadDiscovery().find([str(tmp_path)], "P1", [PAIRED])
    assert group.rebased("tasks/trim/output/P1/out").paths == {
        "left": "tasks/trim/output/P1/out/L1/x.R1.fastq.gz",
        "right": "tasks/trim/output/P1/out/L1/x.R2.fastq.gz",
    }


def test_walks_each_root_once(tmp_path, monkeypatch):
    _touch(tmp_path, "P1/x.R1.fastq.gz", "P1/x.R2.fastq.gz", "P2/x.R1.fastq.gz", "P2/x.R2.fastq.gz")
    walks = []
    original_walk = __import__("os").walk
    monkeypatch.setattr(
        "snappy_pipeline.reads.os.walk", lambda *a, **k: walks.append(a) or original_walk(*a, **k)
    )
    discovery = ReadDiscovery()
    discovery.find([str(tmp_path)], "P1", [PAIRED])
    discovery.find([str(tmp_path)], "P2", [PAIRED])
    assert len(walks) == 1


def test_single_end_needs_mixed_se_pe(tmp_path):
    _touch(tmp_path, "P1/x.R1.fastq.gz")
    with pytest.raises(ValueError, match="no right files for 'P1'; set mixed_se_pe"):
        ReadDiscovery().find([str(tmp_path)], "P1", [PAIRED])
    (group,) = ReadDiscovery().find([str(tmp_path)], "P1", [PAIRED], single_end=True)
    assert group.right is None


@pytest.mark.parametrize(
    "files, error",
    [
        (
            ["P1/a.R1.fastq.gz", "P1/a.R2.fastq.gz", "P1/b.R1.fastq.gz"],
            "mixed single-end and paired",
        ),
        (["P1/a.R2.fastq.gz"], "Read group 'a' of 'P1' has no left file"),
    ],
)
def test_errors(tmp_path, files, error):
    _touch(tmp_path, *files)
    with pytest.raises(ValueError, match=error):
        ReadDiscovery().find([str(tmp_path)], "P1", [PAIRED])


def test_duplicate_mates_across_roots(tmp_path):
    _touch(tmp_path, "a/P1/x.R1.fastq.gz", "b/P1/x.R1.fastq.gz")
    with pytest.raises(ValueError, match="Read group 'x' of 'P1' has two left files"):
        ReadDiscovery().find([str(tmp_path / "a"), str(tmp_path / "b")], "P1", [PAIRED])


def test_missing_library_suggests_similar_folders(tmp_path):
    _touch(tmp_path, "P001-N1-DNA1-WES1/x.R1.fastq.gz")
    discovery = ReadDiscovery()
    assert discovery.find([str(tmp_path)], "P001-N1-DNA1-WES2", [PAIRED]) == []
    assert discovery.missing([str(tmp_path)], "P001-N1-DNA1-WES2") == (
        f"Found no reads of 'P001-N1-DNA1-WES2' below {tmp_path}; "
        "similar folders: P001-N1-DNA1-WES1"
    )


def test_patterns_need_a_readgroup_group():
    with pytest.raises(ValueError, match=r"has no \(\?P<readgroup>...\) group"):
        compile_search_pattern({"left": r".*\.R1\.fastq\.gz"})


BAM = {"bam": r".+\.bam", "bai": r".+\.bam\.bai"}


def test_find_files_returns_one_file_per_key(tmp_path):
    _touch(tmp_path, "P1/x.bam", "P1/x.bam.bai", "P1/x.bam.md5", "P2/y.bam")
    assert ReadDiscovery().find_files([str(tmp_path)], "P1", [BAM]) == {
        "bam": str(tmp_path / "P1/x.bam"),
        "bai": str(tmp_path / "P1/x.bam.bai"),
    }


@pytest.mark.parametrize(
    "files, error",
    [
        (["P1/x.bam", "P1/y.bam"], "Two bam files for 'P1'"),
        (["P2/x.bam"], "Found no files of 'P1' below"),
    ],
)
def test_find_files_errors(tmp_path, files, error):
    _touch(tmp_path, *files)
    with pytest.raises(ValueError, match=error):
        ReadDiscovery().find_files([str(tmp_path)], "P1", [BAM])


# -- files tables -------------------------------------------------------------------------------


def _table(path, *lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_reads_table_gives_read_groups_in_left_path_order(tmp_path):
    table = _table(
        tmp_path / "reads.tsv",
        "library\treadgroup\tleft\tright",
        "P1\tL2\traw/P1/L2/x.R1.fastq.gz\traw/P1/L2/x.R2.fastq.gz",
        "P1\tL1\traw/P1/L1/x.R1.fastq.gz\traw/P1/L1/x.R2.fastq.gz",
        "",
        "P2\tL1\t/data/P2.R1.fastq.gz\t/data/P2.R2.fastq.gz",
    )
    files = read_files_tables([table], reads=True)

    assert files.keys == ("left", "right")
    assert files.libraries == ["P1", "P2"]
    groups = files.read_groups("P1")
    assert [g.name for g in groups] == ["L1", "L2"]
    assert groups[0].left == str(tmp_path / "raw/P1/L1/x.R1.fastq.gz")
    assert groups[0].right == str(tmp_path / "raw/P1/L1/x.R2.fastq.gz")
    assert groups[0].relpaths == {"left": "L1/x.R1.fastq.gz", "right": "L1/x.R2.fastq.gz"}
    assert files.read_groups("P2")[0].relpaths == {
        "left": "P2.R1.fastq.gz",
        "right": "P2.R2.fastq.gz",
    }
    assert files.read_groups("P3") == []
    assert len(files.paths()) == 6


def test_reads_table_skips_comment_lines(tmp_path):
    table = _table(
        tmp_path / "t.tsv",
        "# snappy scan-files --task reads",
        "library\treadgroup\tleft",
        "P1\tL1\tx",
    )
    assert read_files_tables([table], reads=True).libraries == ["P1"]
    with pytest.raises(ValueError, match=":2: the columns must start with"):
        read_files_tables([_table(tmp_path / "u.tsv", "# comment", "library\tleft")], reads=True)


def test_reads_table_expands_pathvars(tmp_path):
    from snakemake.pathvars import Pathvars

    table = _table(tmp_path / "t.tsv", "library\treadgroup\tleft", "P1\tL1\t<raw>/P1.fastq.gz")
    pathvars = Pathvars.with_defaults()
    pathvars.update(Pathvars.from_config({"pathvars": {"raw": "/data/raw"}}))
    files = read_files_tables([table], reads=True, pathvars=pathvars)
    assert files.read_groups("P1")[0].left == "/data/raw/P1.fastq.gz"


def test_reads_tables_merge_and_reject_a_path_listed_twice(tmp_path):
    a = _table(tmp_path / "a/reads.tsv", "library\treadgroup\tleft", "P1\tL1\tP1.fastq.gz")
    b = _table(tmp_path / "b/reads.tsv", "library\treadgroup\tleft", "P2\tL1\tP2.fastq.gz")
    assert read_files_tables([a, b], reads=True).libraries == ["P1", "P2"]
    c = _table(tmp_path / "a/more.tsv", "library\treadgroup\tleft", "P3\tL1\tP1.fastq.gz")
    with pytest.raises(
        ValueError, match=rf"more.tsv:2: {tmp_path}/a/P1.fastq.gz is also listed at .*reads.tsv:2"
    ):
        read_files_tables([a, c], reads=True)
    d = _table(tmp_path / "d/reads.tsv", "library\treadgroup\tleft\tright", "P4\tL1\tx.R1\tx.R2")
    with pytest.raises(ValueError, match="columns left, right differ from left"):
        read_files_tables([a, d], reads=True)


@pytest.mark.parametrize(
    "lines, error",
    [
        (["library\tleft"], ":1: the columns must start with library, readgroup, left"),
        (["library\treadgroup\tleft\tleft"], ":1: a column name is repeated"),
        (["library\treadgroup\tleft", "P1\tL1"], ":2: 2 cells, the header has 3"),
        (["library\treadgroup\tleft", "P1\t\tx"], ":2: readgroup is empty"),
        (["library\treadgroup\tleft", "\tL1\tx"], ":2: library is empty"),
        (["library\treadgroup\tleft", "P1\tL1\t"], ":2: left is empty"),
        (
            ["library\treadgroup\tleft", "P1\tL1\tx", "P1\tL1\ty"],
            ":3: read group 'L1' of 'P1' is also listed at .*:2",
        ),
        (
            ["library\treadgroup\tleft\tright", "P1\tL1\tx\ty", "P1\tL2\tz\t"],
            "Library 'P1' mixes single-end and paired-end",
        ),
        (
            ["library\treadgroup\tleft\tright", "P1\tL1\tx\t"],
            "Library 'P1' has no right files; set mixed_se_pe",
        ),
    ],
)
def test_reads_table_errors(tmp_path, lines, error):
    table = _table(tmp_path / "reads.tsv", *lines)
    with pytest.raises(ValueError, match=error):
        read_files_tables([table], reads=True)


def test_reads_table_single_end_with_mixed_se_pe(tmp_path):
    table = _table(
        tmp_path / "t.tsv", "library\treadgroup\tleft\tright", "P1\tL1\tx\t", "P2\tL1\ty\tz"
    )
    files = read_files_tables([table], reads=True, mixed_se_pe=True)
    assert files.read_groups("P1")[0].right is None
    assert files.read_groups("P2")[0].right == str(tmp_path / "z")


def test_reads_table_without_right_column_is_single_end(tmp_path):
    table = _table(tmp_path / "t.tsv", "library\treadgroup\tleft", "P1\tL1\tx")
    assert read_files_tables([table], reads=True).read_groups("P1")[0].right is None


def test_files_table_gives_one_row_of_files_per_library(tmp_path):
    table = _table(
        tmp_path / "bams.tsv", "library\tbam\tbai", "P1\tP1.bam\tP1.bam.bai", "P2\t/x/P2.bam\t"
    )
    files = read_files_tables([table], reads=False)
    assert files.keys == ("bam", "bai")
    assert files.files("P1") == {
        "bam": str(tmp_path / "P1.bam"),
        "bai": str(tmp_path / "P1.bam.bai"),
    }
    assert files.files("P2") == {"bam": "/x/P2.bam"}
    assert files.files("P3") == {}
    assert files.paths() == [str(tmp_path / "P1.bam"), str(tmp_path / "P1.bam.bai"), "/x/P2.bam"]


@pytest.mark.parametrize(
    "lines, error",
    [
        (["library"], ":1: the table needs a file column besides library"),
        (["bam\tlibrary"], ":1: the columns must start with library"),
        (["library\tbam", "P1\ta", "P1\tb"], ":3: library 'P1' is also listed at .*:2"),
        ([], ":1: the files table is empty"),
    ],
)
def test_files_table_errors(tmp_path, lines, error):
    table = _table(tmp_path / "bams.tsv", *lines) if lines else _table(tmp_path / "bams.tsv")
    with pytest.raises(ValueError, match=error):
        read_files_tables([table], reads=False)


def test_missing_files_table(tmp_path):
    with pytest.raises(ValueError, match="Files table .*/none.tsv: No such file"):
        read_files_tables([str(tmp_path / "none.tsv")], reads=True)


# -- files named after the library -----------------------------------------------------------


def test_files_named_after_the_library_are_found_without_a_folder(tmp_path):
    _touch(
        tmp_path,
        "P001-N1_S3_L001.R1.fastq.gz",
        "P001-N1_S3_L001.R2.fastq.gz",
        "batch2/P001-N1_S3_L002.R1.fastq.gz",
        "batch2/P001-N1_S3_L002.R2.fastq.gz",
        "P001-N1-DNA1.R1.fastq.gz",  # another library whose name starts like this one
        "P001-N1-DNA1.R2.fastq.gz",
        "P001-N10.R1.fastq.gz",
    )
    groups = ReadDiscovery().find([str(tmp_path)], "P001-N1", [PAIRED])
    assert [(g.name, g.relpaths["left"]) for g in groups] == [
        ("S3_L001", "P001-N1_S3_L001.R1.fastq.gz"),
        ("S3_L002", "P001-N1_S3_L002.R1.fastq.gz"),
    ]


def test_a_file_named_like_the_library_keeps_the_name_as_read_group(tmp_path):
    _touch(tmp_path, "P1.R1.fastq.gz", "P1.R2.fastq.gz")
    (group,) = ReadDiscovery().find([str(tmp_path)], "P1", [PAIRED])
    assert group.name == "P1"


def test_a_folder_named_after_the_library_takes_precedence_over_its_name(tmp_path):
    _touch(tmp_path, "P1/P1_L1.R1.fastq.gz", "P1/P1_L1.R2.fastq.gz")
    (group,) = ReadDiscovery().find([str(tmp_path)], "P1", [PAIRED])
    assert group.name == "L1" and group.relpaths["left"] == "P1_L1.R1.fastq.gz"


def test_find_files_matches_file_names_too(tmp_path):
    _touch(tmp_path, "P1.bam", "P1.bam.bai", "P10.bam")
    found = ReadDiscovery().find_files(
        [str(tmp_path)], "P1", [{"bam": r".+\.bam", "bai": r".+\.bai"}]
    )
    assert found == {"bam": str(tmp_path / "P1.bam"), "bai": str(tmp_path / "P1.bam.bai")}
