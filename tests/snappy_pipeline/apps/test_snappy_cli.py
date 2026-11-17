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


def test_task_add_keeps_the_style_of_the_config_and_comments_the_new_task(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG.replace("tool: bwa", 'tool: "bwa"'))
    result = CliRunner().invoke(
        snappy_cli.main, ["task", "add", "--directory", str(tmp_path), "qc=ngs_data_qc"]
    )
    assert result.exit_code == 0, result.output
    text = (tmp_path / "config.yaml").read_text()
    # The user's quotes and sequence style stay, the backup holds the previous file
    assert 'tool: "bwa"' in text
    assert (tmp_path / "config.yaml.bak").read_text() == CHECK_FILES_CONFIG.replace(
        "tool: bwa", 'tool: "bwa"'
    )
    lines = text.splitlines()
    assert "  - name: qc" in lines and "    step: ngs_data_qc" in lines
    # The comments of the task block sit at the indentation of their keys
    assert "      #library_selection:" in lines
    assert "      #depends_on:" in lines
    assert "        # Alignments of DNA libraries" in lines
    assert any(line.startswith("        #alignments: ''") for line in lines)
    assert "      tool:" in {line.split("#")[0].rstrip() for line in lines}
    # The file still loads, with both tasks
    assert [
        task["name"] for task in snappy_cli.load_configfile(str(tmp_path / "config.yaml"))["tasks"]
    ] == [
        "genome",
        "bwa_index",
        "mapping",
        "qc",
    ]


def test_task_add_minimal_writes_the_required_fields_only(tmp_path):
    (tmp_path / "config.yaml").write_text("tasks: []\n")
    result = CliRunner().invoke(
        snappy_cli.main,
        [
            "task",
            "add",
            "--directory",
            str(tmp_path),
            "--config-mode",
            "minimal",
            "mapping=ngs_mapping",
        ],
    )
    assert result.exit_code == 0, result.output
    text = (tmp_path / "config.yaml").read_text()
    assert "bwa:" not in text and "library_selection" not in text
    assert "      tool:" in {line.split("#")[0].rstrip() for line in text.splitlines()}
    # No descriptions, only the REQUIRED tags
    assert not any(line.lstrip().startswith("#") for line in text.splitlines())
    assert "# REQUIRED for tool: star" in text and "Options:" not in text
    config = snappy_cli.load_configfile(str(tmp_path / "config.yaml"))
    assert config["tasks"][0]["config"]["depends_on"] == {
        "reads": None,
        "index": None,
        "reference": None,
        "features": "",
        "exons": "",
    }


def test_init_writes_the_static_sections_before_the_tasks(tmp_path):
    project = tmp_path / "project"
    result = CliRunner().invoke(
        snappy_cli.main,
        ["init", str(project), "--task", "genome=external_data", "--task", "mapping=ngs_mapping"],
    )
    assert result.exit_code == 0, result.output
    lines = (project / "config.yaml").read_text().splitlines()
    # The tasks end the file, so that an added task never lands before a section comment
    data_sets, tasks = lines.index("data_sets: {}"), lines.index("tasks:")
    assert data_sets < tasks < lines.index("  - name: genome") < lines.index("  - name: mapping")
    assert not any(line.startswith("#") for line in lines[tasks:])
    config = snappy_cli.load_configfile(str(project / "config.yaml"))
    assert [task["name"] for task in config["tasks"]] == ["genome", "mapping"]
    assert config["data_sets"] == {}


def test_init_minimal_points_to_the_docs_instead_of_explaining(tmp_path):
    project = tmp_path / "project"
    result = CliRunner().invoke(
        snappy_cli.main,
        ["init", str(project), "--config-mode", "minimal", "--task", "mapping=ngs_mapping"],
    )
    assert result.exit_code == 0, result.output
    text = (project / "config.yaml").read_text()
    lines = text.splitlines()
    # The comments are the header, which names the docs; the task has only its REQUIRED tags
    comments = [line for line in lines if line.lstrip().startswith("#")]
    assert comments == lines[: len(comments)] and "snappy-pipeline.readthedocs.io" in text
    assert "data_sets: {}" in lines and "  - name: mapping" in lines
    assert "# REQUIRED" in text and "Options:" not in text
