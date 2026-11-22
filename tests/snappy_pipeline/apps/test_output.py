# -*- coding: utf-8 -*-
"""Tests of the listing output: JSON, tab-separated lines and the table on a terminal."""

import json

import click
import pytest
from rich.console import Console

from snappy_pipeline.apps.impl import output
from snappy_pipeline.apps.impl.output import Column, print_listing, text

ROWS = [
    {
        "name": "mapping",
        "step": "ngs_mapping",
        "depends_on": {"reads": "data_sets", "kits": ["a", "b"]},
    },
    {"name": "genome", "step": "external_data", "depends_on": {}},
]
COLUMNS = [Column("name", "Task"), Column("step", "Step"), Column("depends_on", "Depends on")]


def test_text_renders_none_lists_and_mappings():
    assert text(None) == ""
    assert text(3) == "3"
    assert text(["a", "b"]) == "a, b"
    assert text({"reads": "data_sets", "kits": ["a", "b"]}) == "reads: data_sets, kits: a, b"


def test_json_keeps_every_key_of_the_rows(capsys):
    print_listing(ROWS, COLUMNS, "Tasks", as_json=True)

    assert json.loads(capsys.readouterr().out) == ROWS


def test_lines_on_a_pipe_and_with_plain(capsys, monkeypatch):
    print_listing(ROWS, COLUMNS, "Tasks")  # the captured output is no terminal

    lines = ["mapping\tngs_mapping\treads: data_sets, kits: a, b", "genome\texternal_data\t"]
    assert capsys.readouterr().out.splitlines() == lines
    monkeypatch.setattr(output, "Console", lambda: Console(force_terminal=True, width=80))
    print_listing(ROWS, COLUMNS, "Tasks", plain=True)
    assert capsys.readouterr().out.splitlines() == lines


def test_table_on_a_terminal(capsys, monkeypatch):
    monkeypatch.setattr(output, "Console", lambda: Console(force_terminal=True, width=80))

    print_listing(ROWS, COLUMNS, "Tasks")

    out = capsys.readouterr().out
    assert "Tasks" in out and "Depends on" in out and "ngs_mapping" in out
    assert "\t" not in out
    print_listing([], COLUMNS, "Tasks", empty="No tasks")
    assert "No tasks" in capsys.readouterr().out


def test_an_empty_listing_is_empty_on_a_pipe_and_as_json(capsys):
    print_listing([], COLUMNS, "Tasks", empty="No tasks")
    assert capsys.readouterr().out == ""
    print_listing([], COLUMNS, "Tasks", as_json=True)
    assert json.loads(capsys.readouterr().out) == []


def test_json_and_plain_exclude_each_other():
    with pytest.raises(click.UsageError, match="exclude each other"):
        print_listing(ROWS, COLUMNS, "Tasks", as_json=True, plain=True)
