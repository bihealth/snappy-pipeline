# -*- coding: utf-8 -*-
"""Tests for linking a task's result files into ``output/`` (``snappy_pipeline.results``)."""

import os

import pytest

from snappy_pipeline.results import link_results, work_paths


def _task(tmp_path):
    work = tmp_path / "work" / "L1" / "out"
    work.mkdir(parents=True)
    for name, text in (("L1.vcf.gz", "vcf"), ("L1.vcf.gz.md5", "md5"), ("L1.txt", "text")):
        (work / name).write_text(text)
    return tmp_path


def test_work_paths_mirror_the_result_files():
    assert work_paths(["output/L1/out/L1.vcf.gz"]) == ["work/L1/out/L1.vcf.gz"]
    with pytest.raises(ValueError, match="outside output/"):
        work_paths(["work/L1/out/L1.vcf.gz"])


def test_results_job_links_the_files_and_writes_the_manifest(tmp_path):
    task = _task(tmp_path)
    files = [str(task / "work/L1/out" / n) for n in ("L1.vcf.gz", "L1.vcf.gz.md5", "L1.txt")]

    link_results(str(task), files, str(task / "output/.manifest.tsv"))

    link = task / "output/L1/out/L1.vcf.gz"
    assert os.readlink(link) == "../../../work/L1/out/L1.vcf.gz"
    assert (task / "output/.manifest.tsv").read_text().splitlines() == [
        "path\tsize\tmd5",
        "L1/out/L1.txt\t4\t",
        "L1/out/L1.vcf.gz\t3\tL1/out/L1.vcf.gz.md5",
        "L1/out/L1.vcf.gz.md5\t3\t",
    ]


def test_results_job_removes_stale_links_but_no_other_files(tmp_path):
    task = _task(tmp_path)
    output = task / "output/L1/out"
    output.mkdir(parents=True)
    (output / "old.vcf.gz").symlink_to("../../../work/L1/out/old.vcf.gz")
    (output / "notes.txt").write_text("kept")

    link_results(str(task), [str(task / "work/L1/out/L1.txt")], str(task / "output/.manifest.tsv"))

    assert sorted(p.name for p in output.iterdir()) == ["L1.txt", "notes.txt"]


def test_a_lone_done_marker_links_its_directory(tmp_path):
    work = tmp_path / "work" / "L1" / "out"
    work.mkdir(parents=True)
    for name in (".done", "L1.R1.fastq.gz", "L1.R2.fastq.gz"):
        (work / name).write_text(name)
    (tmp_path / "output" / "L1" / "out").mkdir(parents=True)  # from linking file by file

    link_results(str(tmp_path), [str(work / ".done")], str(tmp_path / "output/.manifest.tsv"))

    assert os.readlink(tmp_path / "output/L1/out") == "../../work/L1/out"
    assert (tmp_path / "output/.manifest.tsv").read_text().splitlines()[1:] == ["L1/out\t33\t"]


def test_a_directory_result_is_linked_as_a_directory(tmp_path):
    work = tmp_path / "work" / "L1" / "out"
    work.mkdir(parents=True)
    (work / "L1.R1.fastq.gz").write_text("reads")
    (work / ".snakemake_timestamp").touch()

    link_results(str(tmp_path), [str(work)], str(tmp_path / "output/.manifest.tsv"))

    assert os.readlink(tmp_path / "output/L1/out") == "../../work/L1/out"
    assert (tmp_path / "output/.manifest.tsv").read_text().splitlines()[1:] == ["L1/out\t5\t"]
