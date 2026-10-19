# -*- coding: utf-8 -*-
"""Tests for ``snappy check-files``."""

from click.testing import CliRunner

from snappy_pipeline.apps.snappy_cli import main

CONFIG = """\
pathvars:
  refs: refs
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: <refs>/genome.fa}
data_sets: {}
"""


def test_check_files_lists_missing_files(tmp_path):
    (tmp_path / "config.yaml").write_text(CONFIG)
    result = CliRunner().invoke(main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.output == f"genome\t{tmp_path}/refs/genome.fa\n"


def test_check_files_succeeds_when_all_files_exist(tmp_path):
    (tmp_path / "config.yaml").write_text(CONFIG)
    (tmp_path / "refs").mkdir()
    (tmp_path / "refs" / "genome.fa").touch()
    result = CliRunner().invoke(main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 0
