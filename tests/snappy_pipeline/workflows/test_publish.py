# -*- coding: utf-8 -*-
"""Tests for linking published files with their md5 file and index (``PublishStepPart``)."""

import re
from types import SimpleNamespace

import pytest
from snakemake.io import regex_from_filepattern

from snappy_pipeline.workflows.abstract import PublishStepPart


def _part(prefix=""):
    part = object.__new__(PublishStepPart)
    part.parent = SimpleNamespace(task_path_prefix=lambda: prefix)
    return part


@pytest.mark.parametrize(
    "name, action",
    [
        ("L1.segments.txt", "run"),
        ("L1.segments.txt.md5", "run"),
        ("L1.tsv.gz", "run"),
        ("L1.vcf.gz", "indexed"),
        ("L1.vcf.gz.md5", "indexed"),
        ("L1.vcf.gz.tbi", "indexed"),
        ("L1.vcf.gz.tbi.md5", "indexed"),
        ("L1.bed.gz.tbi", "indexed"),
        (".done", "marker"),
        ("L1.done", "marker"),
    ],
)
def test_exactly_one_action_links_a_file(name, action):
    part = _part()
    matches = [
        a
        for a in PublishStepPart.actions
        for pattern in part.get_output_files(a)
        if re.match(regex_from_filepattern(pattern), f"output/L1/out/{name}")
    ]
    assert matches == [action]


def test_run_locally_links_the_file_with_its_companions(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    work = tmp_path / "tasks" / "calling" / "work" / "L1" / "out"
    work.mkdir(parents=True)
    for suffix in ("", ".md5", ".tbi", ".tbi.md5"):
        (work / f"L1.vcf.gz{suffix}").write_text(suffix)
    (tmp_path / "tasks" / "calling" / "output" / "L1" / "out").mkdir(parents=True)

    _part("tasks/calling/").run_locally("indexed", SimpleNamespace(path="L1/out", name="L1.vcf.gz"))

    output = tmp_path / "tasks" / "calling" / "output" / "L1" / "out"
    assert sorted(p.name for p in output.iterdir()) == [
        "L1.vcf.gz",
        "L1.vcf.gz.md5",
        "L1.vcf.gz.tbi",
        "L1.vcf.gz.tbi.md5",
    ]
    assert (output / "L1.vcf.gz.tbi").readlink().as_posix() == "../../../work/L1/out/L1.vcf.gz.tbi"
