# -*- coding: utf-8 -*-
"""Tests for ``snappy sheet check``."""

from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main


def test_sheet_check_reports_the_sheets_or_their_problems(tmp_path):
    runner = CliRunner()
    path = tmp_path / "s.tsv"
    path.write_text(
        "donor_name\tsex\tsample_name\tis_tumor\tlibrary_name\tlibrary_type\n"
        "P1\tfemale\tP1-N1\tfalse\tP1-N1-DNA1-WES1\tWES\n"
        "P1\tfemale\tP1-T1\ttrue\tP1-T1-DNA1-WES1\tWES\n"
    )
    result = runner.invoke(main, ["sheet", "check", str(path)])
    assert result.exit_code == 0, result.output
    assert f"2 libraries of 1 donors in {path}" in result.output

    (tmp_path / "bad.tsv").write_text("donor_name\tlibrary_name\tlibrary_type\nP1\tL1\tRNA\n")
    result = runner.invoke(main, ["sheet", "check", str(tmp_path / "bad.tsv")])
    assert result.exit_code == 1
    assert "bad.tsv:2: library_type: 'RNA' is no known value" in result.output


def test_sheet_check_lists_synonyms_and_fixes_them(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text("donor_name\tsex\tlibrary_name\tlibrary_type\nP1\tF\tL1\texome\n")
    result = CliRunner().invoke(main, ["sheet", "check", str(path)])
    assert result.exit_code == 0, result.output
    assert f"{path}:2: sex: 'F' is read as 'female'" in result.output
    assert (
        "2 cells are synonyms of their value; snappy sheet check --fix rewrites them"
        in result.output
    )
    result = CliRunner().invoke(main, ["sheet", "check", "--fix", str(path)])
    assert result.exit_code == 0, result.output
    assert f"Rewrote 2 cells of {path}" in result.output
    assert path.read_text() == "donor_name\tsex\tlibrary_name\tlibrary_type\nP1\tfemale\tL1\tWES\n"
