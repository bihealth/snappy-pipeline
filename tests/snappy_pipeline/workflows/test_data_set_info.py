# -*- coding: utf-8 -*-
"""Tests for loading the sample sheets of the data sets (``DataSetInfo``)."""

import textwrap

import pytest

from snappy_pipeline.base import ConfigError
from snappy_pipeline.library_table import DataSetInfo

SHEET = textwrap.dedent(
    """
    [Data]
    patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tfolderName\thpoTerms
    P001\t.\t.\tF\tY\t{library_type}\tP001\t.
    """
).lstrip()


def _info(path):
    return DataSetInfo(
        name="first",
        sheet_path=str(path),
        sheet_type="germline_variants",
        is_background=False,
        naming_scheme="only_secondary_id",
    )


def test_sheet_errors_name_the_sheet_file(tmp_path):
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(SHEET.format(library_type="WESS"))
    with pytest.raises(ConfigError, match="^Sample sheet .*sheet.tsv: Invalid library type WESS"):
        _info(sheet)

    sheet.write_text(SHEET.format(library_type="WES"))
    assert [donor.name for donor in _info(sheet).sheet.bio_entities.values()] == ["P001"]


@pytest.mark.parametrize("library_type", ["", "."])
def test_sheet_without_library_types_is_a_config_error(tmp_path, library_type):
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(SHEET.format(library_type=library_type))
    with pytest.raises(ConfigError, match="^Sample sheet .*sheet.tsv: no row has a libraryType"):
        _info(sheet)


def test_missing_sheet_is_a_config_error(tmp_path):
    with pytest.raises(ConfigError, match="Could not find sample sheet file"):
        _info(tmp_path / "missing.tsv")
