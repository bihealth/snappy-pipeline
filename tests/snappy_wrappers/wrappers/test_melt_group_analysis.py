"""MELT GroupAnalysis must see the IndivAnalysis output of all pedigree members in one directory."""

import os
import stat
from pathlib import Path
from types import SimpleNamespace

from snappy_wrappers import snappy_wrapper

WRAPPER = Path(__file__).parents[3] / "snappy_wrappers/wrappers/melt/group_analysis/wrapper.py"

#: Stands in for ``java``: records the files below the ``-discoverydir`` argument
FAKE_JAVA = """#!/usr/bin/env bash
while [[ $# -gt 0 ]]; do
    if [[ $1 == -discoverydir ]]; then ls "$2" > "$RECORD"; fi
    shift
done
"""

LIBRARIES = ("index", "father", "mother")


def test_group_analysis_discoverydir_holds_all_members(tmp_path, monkeypatch):
    done = []
    for library in LIBRARIES:
        out = tmp_path / f"indiv_analysis.{library}.ALU" / "out"
        out.mkdir(parents=True)
        (out / f"{library}.ALU.tmp.bed").touch()
        (out / f"{library}.ALU.hum_breaks.sorted.bam").touch()
        (out / f".done.{library}").touch()
        done.append(str(out / f".done.{library}"))
    (tmp_path / "group_analysis" / "out").mkdir(parents=True)
    (tmp_path / "bin").mkdir()
    java = tmp_path / "bin" / "java"
    java.write_text(FAKE_JAVA)
    java.chmod(java.stat().st_mode | stat.S_IXUSR)
    record = tmp_path / "record.txt"
    monkeypatch.setenv("PATH", f"{tmp_path / 'bin'}:{os.environ['PATH']}")
    monkeypatch.setenv("RECORD", str(record))
    monkeypatch.setattr(snappy_wrapper, "write_declared_md5s", lambda snakemake: None)
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            reference="ref.fa",
            indiv_analysis=done,
            jar="MELT.jar",
            me_refs="refs",
            genes="genes.bed",
        ),
        output=SimpleNamespace(done=str(tmp_path / "group_analysis" / "out" / ".done")),
        log=SimpleNamespace(log=str(tmp_path / "log" / "group_analysis.log")),
        params=SimpleNamespace(args={"me_refs_infix": "1KGP_Hg19", "me_type": "ALU"}),
        resources=SimpleNamespace(mem_mb=4096),
        rule="melt_group_analysis",
        threads=1,
        wildcards=SimpleNamespace(),
    )

    exec(compile(WRAPPER.read_text(), str(WRAPPER), "exec"), {"snakemake": snakemake})

    assert record.read_text().split() == sorted(
        f"{library}.ALU.{ext}"
        for library in LIBRARIES
        for ext in ("tmp.bed", "hum_breaks.sorted.bam")
    )
