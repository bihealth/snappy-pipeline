# -*- coding: utf-8 -*-
"""Tests for ``snappy`` CLI subcommands, completion, and task listing."""

import textwrap

from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_steps_match_workflow_registry():
    """STEPS used by the CLI must match the WORKFLOW_REGISTRY keys."""
    assert set(snappy_cli.STEPS) == set(WORKFLOW_REGISTRY)


def test_snappy_run_help(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check whether the call to ``snappy run --help`` works."""
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_cli", fake_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", fake_fs, mocker)
    m = mocker.MagicMock()
    mocker.patch("snappy_pipeline.apps.snappy_cli.snakemake_main", m)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["run", "--help"])
    assert result.exit_code == 0
    m.assert_not_called()


def test_snappy_completion():
    """Check whether ``snappy completion`` generates autocompletion code."""
    runner = CliRunner()
    result_zsh = runner.invoke(snappy_cli.main, ["completion", "zsh"])
    assert result_zsh.exit_code == 0
    assert "_SNAPPY_COMPLETE" in result_zsh.output

    result_bash = runner.invoke(snappy_cli.main, ["completion", "bash"])
    assert result_bash.exit_code == 0
    assert "_SNAPPY_COMPLETE" in result_bash.output


def test_snappy_task_steps():
    """Check whether ``snappy task steps`` lists available workflow steps."""
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "steps"])
    assert result.exit_code == 0
    assert "Available Workflow Steps in SNAPPY" in result.output
    assert "ngs_mapping" in result.output


def test_snappy_task_list_available():
    """Check whether ``snappy task list --available`` lists available steps."""
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "list", "--available"])
    assert result.exit_code == 0
    assert "Available Workflow Steps in SNAPPY" in result.output
    assert "ngs_mapping" in result.output


def test_task_param_type_shell_complete():
    """Check shell completion helper for TaskParamType."""
    param_type = snappy_cli.TaskParamType()
    completions = param_type.shell_complete(None, None, "ngs_")
    names = [c.value for c in completions]
    assert "ngs_mapping" in names

    completions_prefix = param_type.shell_complete(None, None, "my_custom_task=ngs_")
    names_prefix = [c.value for c in completions_prefix]
    assert "my_custom_task=ngs_mapping" in names_prefix


def test_complete_task_names(germline_sheet_fake_project_ngs_mapping_fs, mocker):
    """Check task name completion from project config.yaml."""
    fake_fs = germline_sheet_fake_project_ngs_mapping_fs
    patch_module_fs("snappy_pipeline.apps.snappy_cli", fake_fs, mocker)

    ctx = mocker.MagicMock()
    ctx.params = {"directory": "/project-dir/ngs_mapping"}

    completions = snappy_cli._complete_task_names(ctx, None, "")
    names = [c.value for c in completions]
    assert "ngs_mapping" in names


#: A mapping task reading the FASTQ files of one germline library
CHECK_FILES_CONFIG = """
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: refs/genome.fa}
  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: refs/genome}
  - name: mapping
    step: ngs_mapping
    config:
      depends_on: {reads: data_sets, reference: genome, index: bwa_index}
      tool: bwa
      bwa: {}
data_sets:
  first:
    file: sheet.tsv
    search_paths: [raw]
    search_patterns:
      - {left: '(?P<readgroup>.+)_R1.fastq.gz', right: '(?P<readgroup>.+)_R2.fastq.gz'}
    type: germline_variants
    naming_scheme: only_secondary_id
"""

CHECK_FILES_SHEET = """
[Data]
patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tfolderName\thpoTerms
P001\t.\t.\tF\tY\tWGS\tP001\t.
"""


def test_check_files_lists_missing_files_and_libraries_without_reads(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG)
    (tmp_path / "sheet.tsv").write_text(textwrap.dedent(CHECK_FILES_SHEET).lstrip())
    (tmp_path / "raw").mkdir()
    for name in ("genome.fa", "genome.fa.fai", "genome.dict"):
        (tmp_path / "refs").mkdir(exist_ok=True)
        (tmp_path / "refs" / name).touch()

    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.output.splitlines() == [
        f"mapping\tP001-N1-DNA1-WGS1\tFound no reads of 'P001' below {tmp_path}/raw",
    ]

    (tmp_path / "raw" / "P001").mkdir()
    for mate in ("R1", "R2"):
        (tmp_path / "raw" / "P001" / f"P001_{mate}.fastq.gz").touch()
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 0, result.output

    (tmp_path / "refs" / "genome.dict").unlink()
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.output.splitlines() == [f"genome\t{tmp_path}/refs/genome.dict"]


def test_check_files_prints_configuration_errors_without_traceback(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG.replace("tool: bwa", "tool: bwa2"))
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.exception is None or isinstance(result.exception, SystemExit)
