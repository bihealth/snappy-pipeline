# -*- coding: utf-8 -*-
"""Tests that the cnvkit wrappers put the options of their params on the command line."""

import runpy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

WRAPPERS = Path(__file__).resolve().parents[3] / "snappy_wrappers" / "wrappers"

CALL_ARGS = {"method": "threshold", "thresholds": "-1.1,-0.25,0.2,0.7", "filter": "ampdel"}


def _run(wrapper, args):
    """Run the wrapper script and return its globals and the command it runs."""
    snakemake = SimpleNamespace(params=SimpleNamespace(args=args))
    with patch("snappy_wrappers.snappy_wrapper.ShellWrapper") as shell:
        result = runpy.run_path(
            str(WRAPPERS / wrapper / "wrapper.py"), init_globals={"snakemake": snakemake}
        )
    return result, shell.return_value.run.call_args.args[0]


@pytest.mark.parametrize("count, flag", [(True, "--count"), (False, "")])
def test_coverage_counts_reads_on_request(count, flag):
    result, command = _run("cnvkit/coverage", {"min_mapq": 0, "count": count})

    assert result["count"] == flag
    assert "{count}" in command


@pytest.mark.parametrize("drop, flag", [(True, " --drop-low-coverage"), (False, "")])
def test_call_drops_low_coverage_bins_on_request(drop, flag):
    result, command = _run("cnvkit/call", {**CALL_ARGS, "drop_low_coverage": drop})

    assert result["drop_low_coverage"] == flag
    assert "{drop_low_coverage}" in command
