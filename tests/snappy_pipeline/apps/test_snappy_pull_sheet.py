# -*- coding: utf-8 -*-
"""Tests for ``snappy pull-sheet``."""

import pytest
import requests
from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli, snappy_pull_sheet

#: An investigation without samples: one study and one assay with empty tables
INVESTIGATION = {
    "studies": {
        "study": {
            "study": {
                "top_header": [{"colspan": 0}, {"colspan": 0}],
                "field_header": [],
                "table_data": [],
            },
            "assays": {"assay": {"top_header": [], "field_header": [], "table_data": []}},
        }
    }
}

ARGS = ["pull-sheet", "--api-key", "SECRET", "--sodar-host", "sodar.example.org"]
ARGS += ["--project-uuid", "uuid", "--project-name", "project"]

#: The source columns of a study with one donor, and the donor's values
SOURCE_COLUMNS = ["Name", "Family", "Father", "Mother", "Sex", "Disease Status", "Batch"]
SOURCE_VALUES = ["P1", "FAM1", "0", "0", "male", "affected", "run-2024-A"]


def _investigation(columns, values):
    """An investigation with one study source (no assay) of the given source columns."""
    study = {
        "top_header": [{"colspan": len(columns)}, {"colspan": 0}],
        "field_header": [{"value": column} for column in columns],
        "table_data": [[{"value": value} for value in values]],
    }
    return {"studies": {"study": {**INVESTIGATION["studies"]["study"], "study": study}}}


class Response:
    def __init__(self, status_code=200, investigation=INVESTIGATION):
        self.status_code = status_code
        self.investigation = investigation

    def raise_for_status(self):
        if self.status_code != 200:
            raise requests.HTTPError(
                f"{self.status_code} Error for url: {snappy_pull_sheet.URL_TPL}", response=self
            )

    def json(self):
        return self.investigation


@pytest.mark.parametrize("status_code", [200, 403])
def test_pull_sheet_keeps_the_api_key_out_of_the_output(monkeypatch, status_code):
    calls = []

    def get(url, **kwargs):
        calls.append((url, kwargs))
        return Response(status_code)

    monkeypatch.setattr(snappy_pull_sheet.requests, "get", get)
    result = CliRunner().invoke(snappy_cli.main, ARGS)
    assert calls == [
        (
            "https://sodar.example.org/samplesheets/api/remote/get/uuid/SECRET",
            {"timeout": snappy_pull_sheet.TIMEOUT},
        )
    ]
    assert result.stderr.startswith(
        "Fetching the sample sheet of project uuid from sodar.example.org\n"
    )
    assert "SECRET" not in result.output
    if status_code == 200:
        assert result.exit_code == 0, result.output
        assert result.stdout.startswith(snappy_pull_sheet.HEADER + "\n")
    else:
        assert result.exit_code == 1
        assert result.exception.code == "Fetching the sample sheet failed: HTTPError 403"


@pytest.mark.parametrize("with_batch", [True, False])
def test_pull_sheet_writes_the_batch_of_sodar_verbatim_as_a_string(monkeypatch, with_batch):
    """A batch is an id such as ``run-2024-A``, not a number; without a batch column, ``.``."""
    end = None if with_batch else -1
    investigation = _investigation(SOURCE_COLUMNS[:end], SOURCE_VALUES[:end])
    monkeypatch.setattr(
        snappy_pull_sheet.requests, "get", lambda url, **kwargs: Response(200, investigation)
    )

    result = CliRunner().invoke(snappy_cli.main, ARGS)

    assert result.exit_code == 0, result.output
    assert result.stdout.startswith(snappy_pull_sheet.HEADER + "\n")
    batch = "run-2024-A" if with_batch else ""
    assert result.stdout.splitlines()[-1] == "\t".join(
        ["FAM1", "P1", "", "", "male", "true", "", "", "", "", batch, "", "", "project", "uuid"]
    )
