# -*- coding: utf-8 -*-
"""Tests for the central FASTQ discovery in ``snappy_pipeline.reads``."""

import pytest

from snappy_pipeline.reads import ReadDiscovery, compile_search_pattern

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
