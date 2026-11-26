# -*- coding: utf-8 -*-
"""Tests for the ``umi_mapping`` step: its config model, read structures, contract and paths."""

from types import SimpleNamespace

import pydantic
import pytest

from snappy_pipeline.base import ConfigError
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.umi_mapping import (
    GROUP_STRATEGY,
    UmiMappingWorkflow,
    parse_read_structures,
)
from snappy_pipeline.workflows.umi_mapping.model import Consensus, UmiMapping

DEPENDS_ON = {"reads": "data_sets", "reference": "genome", "index": "bwa_mem3_index"}


def _config(consensus="duplex", depends_on=None, **fgumi):
    section = {consensus: {}} | fgumi
    return UmiMapping(
        tool="fgumi", consensus=consensus, fgumi=section, depends_on=DEPENDS_ON | (depends_on or {})
    )


def test_the_section_of_the_consensus_mode_is_explicit():
    with pytest.raises(pydantic.ValidationError, match="requires the section fgumi.duplex"):
        UmiMapping(tool="fgumi", consensus="duplex", fgumi={}, depends_on=DEPENDS_ON)
    with pytest.raises(pydantic.ValidationError, match="requires the section fgumi.simplex"):
        UmiMapping(tool="fgumi", consensus="simplex", fgumi={"duplex": {}}, depends_on=DEPENDS_ON)
    # The section of the other mode is ignored, so a generated config may carry both
    config = _config("simplex", duplex={})
    assert config.fgumi.simplex.filter_min_reads == 3
    assert config.fgumi.duplex.filter_min_reads == [3, 2, 2]


def test_correct_and_the_umis_task_come_together():
    with pytest.raises(pydantic.ValidationError, match="umis"):
        _config(correct={})
    with pytest.raises(pydantic.ValidationError, match="umis"):
        _config(depends_on={"umis": "kit_umis"})
    config = _config(correct={}, depends_on={"umis": "kit_umis"})
    assert config.fgumi.correct.max_mismatches == 2
    assert config.depends_on.umis == "kit_umis"


@pytest.mark.parametrize(
    "section",
    [
        {"duplex": {"min_reads": []}},
        {"duplex": {"min_reads": [1, 1, 1, 1]}},
        {"duplex": {"filter_min_reads": [3, 2, 2, 1]}},
        {"simplex": {"filter_min_reads": 0}},
    ],
)
def test_min_reads_take_one_to_three_values(section):
    mode = next(iter(section))
    with pytest.raises(pydantic.ValidationError):
        UmiMapping(tool="fgumi", consensus=mode, fgumi=section, depends_on=DEPENDS_ON)


def test_the_defaults_follow_the_consensus_mode():
    duplex, simplex = _config("duplex").fgumi.duplex, _config("simplex").fgumi.simplex
    assert (duplex.min_reads, duplex.filter_min_reads) == ([1], [3, 2, 2])
    assert duplex.require_single_strand_agreement is True
    assert (simplex.min_reads, simplex.filter_min_reads) == (1, 3)
    assert GROUP_STRATEGY == {Consensus.duplex: "paired", Consensus.simplex: "adjacency"}


@pytest.mark.parametrize(
    "value, structures",
    [
        ("8M+T 8M+T", ["8M+T", "8M+T"]),
        ("8M+T +T", ["8M+T", "+T"]),
        ("3M4S+T 3M4S+T", ["3M4S+T"] * 2),
    ],
)
def test_read_structures_are_split_per_file(value, structures):
    assert parse_read_structures("lib", value) == structures


@pytest.mark.parametrize("value", [None, float("nan"), "", ".", "8X+T", "M+T", "8M+T+T"])
def test_missing_or_malformed_read_structures_are_rejected(value):
    with pytest.raises(ConfigError, match="lib"):
        parse_read_structures("lib", value)


def test_the_alignments_carry_the_consensus_mode():
    for mode in ("duplex", "simplex"):
        (signature,) = UmiMappingWorkflow.task_produces(SimpleNamespace(consensus=mode), {})
        assert signature == DataSignature(
            DataType.ALIGNMENTS, frozenset({"dna", "consensus", mode})
        )
        assert signature in UmiMappingWorkflow.produces
    assert signature.satisfies(DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})))


def test_output_paths_are_those_of_a_mapping_task():
    paths = UmiMappingWorkflow.get_output_paths(None, library_name="L1")
    assert paths == {
        "bam": "work/L1/out/L1.bam",
        "bai": "work/L1/out/L1.bam.bai",
        "bamstats": "work/L1/report/bam_qc/L1.bam.bamstats.txt",
        "flagstats": "work/L1/report/bam_qc/L1.bam.flagstats.txt",
        "idxstats": "work/L1/report/bam_qc/L1.bam.idxstats.txt",
    }


def test_the_default_config_requires_the_consensus_mode():
    yaml = UmiMappingWorkflow.default_config_yaml()
    assert "consensus:" in yaml and "REQUIRED; Options: 'duplex', 'simplex'" in yaml
    assert "#umis: ''" in yaml and "REQUIRED for tool: fgumi with fgumi.correct" in yaml
