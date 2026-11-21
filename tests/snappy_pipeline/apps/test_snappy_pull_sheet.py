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


class Response:
    def __init__(self, status_code=200):
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code != 200:
            raise requests.HTTPError(
                f"{self.status_code} Error for url: {snappy_pull_sheet.URL_TPL}", response=self
            )

    def json(self):
        return INVESTIGATION


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
        assert result.stdout.startswith("[Metadata]\n") and "title\t.\n" in result.stdout
    else:
        assert result.exit_code == 1
        assert result.exception.code == "Fetching the sample sheet failed: HTTPError 403"
