# -*- coding: utf-8 -*-
"""Tests for ``snappy`` CLI subcommands, completion, and task listing."""

import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli


def test_snappy_run_help(mocker):
    """Check whether the call to ``snappy run --help`` works."""
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


def test_snappy_task_list_available():
    """``snappy task list --available`` lists the steps: lines on a pipe, objects as JSON."""
    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["task", "list", "--available"])
    assert result.exit_code == 0
    assert "ngs_mapping\tRead mapping" in result.output
    result = runner.invoke(snappy_cli.main, ["task", "list", "--available", "--json"])
    steps = {row["step"]: row["description"] for row in json.loads(result.output)}
    assert steps["ngs_mapping"].startswith("Read mapping")


def test_task_list_shows_name_step_tool_and_dependencies(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG)
    runner = CliRunner()

    result = runner.invoke(snappy_cli.main, ["task", "list", "--directory", str(tmp_path)])

    assert result.exit_code == 0
    assert result.output.splitlines() == [
        "genome\texternal_data\t\t",
        "bwa_index\texternal_data\t\t",
        "reads\texternal_data\t\t",
        "mapping\tngs_mapping\tbwa\treads: reads, reference: genome, index: bwa_index",
    ]
    result = runner.invoke(
        snappy_cli.main, ["task", "list", "--json", "--directory", str(tmp_path)]
    )
    assert json.loads(result.output)[3] == {
        "name": "mapping",
        "step": "ngs_mapping",
        "tool": "bwa",
        "depends_on": {"reads": "reads", "reference": "genome", "index": "bwa_index"},
    }
    (tmp_path / "config.yaml").write_text("tasks: []\n")
    result = runner.invoke(snappy_cli.main, ["task", "list", "--directory", str(tmp_path)])
    assert result.exit_code == 0 and result.output == ""
    result = runner.invoke(
        snappy_cli.main, ["task", "list", "--json", "--plain", "--directory", str(tmp_path)]
    )
    assert result.exit_code == 2
    assert "excludeeachother" in re.sub(r"\x1b\[[0-9;]*m|[\s│╭╮╰╯─]", "", result.output)


def test_task_param_type_shell_complete():
    """Check shell completion helper for TaskParamType."""
    param_type = snappy_cli.TaskParamType()
    completions = param_type.shell_complete(None, None, "ngs_")
    names = [c.value for c in completions]
    assert "ngs_mapping" in names

    completions_prefix = param_type.shell_complete(None, None, "my_custom_task=ngs_")
    names_prefix = [c.value for c in completions_prefix]
    assert "my_custom_task=ngs_mapping" in names_prefix


def test_complete_task_names(tmp_path, mocker):
    """The task names of config.yaml complete ``--task``; without a config, the step names."""
    ctx = mocker.MagicMock()
    ctx.params = {"directory": str(tmp_path)}
    assert [c.value for c in snappy_cli._complete_task_names(ctx, None, "")] == list(
        snappy_cli.STEPS
    )

    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG)
    assert [c.value for c in snappy_cli._complete_task_names(ctx, None, "")] == [
        "bwa_index",
        "genome",
        "mapping",
        "reads",
    ]
    assert [c.value for c in snappy_cli._complete_task_names(ctx, None, "ma")] == ["mapping"]


#: A mapping task reading the FASTQ files of one germline library from a files table
CHECK_FILES_CONFIG = """
project:
  name: check
sample_sheets: [sheet.tsv]
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
  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]
  - name: mapping
    step: ngs_mapping
    config:
      depends_on: {reads: reads, reference: genome, index: bwa_index}
      tool: bwa
      bwa: {}
"""

CHECK_FILES_SHEET = """
donor_name\tsex\tlibrary_name\tlibrary_type\tfolder_name
P001\tfemale\tP001-N1-DNA1-WGS1\tWGS\tP001
"""


def test_check_files_lists_missing_files_and_libraries_without_reads(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG)
    (tmp_path / "sheet.tsv").write_text(textwrap.dedent(CHECK_FILES_SHEET).lstrip())
    (tmp_path / "raw").mkdir()
    for name in ("genome.fa", "genome.fa.fai", "genome.dict"):
        (tmp_path / "refs").mkdir(exist_ok=True)
        (tmp_path / "refs" / name).touch()

    # Before the scan, the files table is missing
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.output.splitlines() == [
        f"file\treads\t{tmp_path}/reads.tsv\tdoes not exist",
        f"reads\treads\t{tmp_path}/reads.tsv\tdoes not exist; run `snappy scan-files`",
    ]

    (tmp_path / "raw" / "P001").mkdir()
    for mate in ("R1", "R2"):
        (tmp_path / "raw" / "P001" / f"P001_{mate}.fastq.gz").touch()
    patterns = [
        "--pattern",
        "left=(?P<readgroup>.+)_R1.fastq.gz",
        "--pattern",
        "right=(?P<readgroup>.+)_R2.fastq.gz",
    ]
    result = CliRunner().invoke(
        snappy_cli.main, ["scan-files", "--directory", str(tmp_path)] + patterns
    )
    assert result.exit_code == 0, result.output
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 0, result.output

    # A listed FASTQ file that goes missing, and a sidecar of the genome
    (tmp_path / "raw" / "P001" / "P001_R2.fastq.gz").unlink()
    (tmp_path / "refs" / "genome.dict").unlink()
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.output.splitlines() == [
        f"file\tgenome\t{tmp_path}/refs/genome.dict\tdoes not exist",
        f"file\treads\t{tmp_path}/raw/P001/P001_R2.fastq.gz\tdoes not exist",
    ]


def test_rule_rows_name_the_task_step_part_and_action():
    part = lambda name, *actions: SimpleNamespace(name=name, actions=actions)  # noqa: E731
    instances = {
        "mapping": SimpleNamespace(
            name="ngs_mapping",
            sub_steps={"bwa": part("bwa", "run"), "qc_summary": part("qc_summary", "run")},
        ),
        "mapping_qc": SimpleNamespace(
            name="ngs_data_qc", sub_steps={"bam_collect_doc": part("bam_collect_doc", "run")}
        ),
    }
    names = [
        "all",
        "mapping_ngs_mapping_bwa_run",
        "mapping_ngs_mapping_qc_summary_run",
        "mapping_qc_ngs_data_qc_bam_collect_doc_run",
        "mapping_qc_results",
        "mapping_results",
    ]

    rows = snappy_cli._rule_rows(names, instances)

    assert [tuple(row.values()) for row in rows] == [
        ("all", "", "", "", ""),
        ("mapping_ngs_mapping_bwa_run", "mapping", "ngs_mapping", "bwa", "run"),
        ("mapping_ngs_mapping_qc_summary_run", "mapping", "ngs_mapping", "qc_summary", "run"),
        (
            "mapping_qc_ngs_data_qc_bam_collect_doc_run",
            "mapping_qc",
            "ngs_data_qc",
            "bam_collect_doc",
            "run",
        ),
        ("mapping_qc_results", "mapping_qc", "ngs_data_qc", "", ""),
        ("mapping_results", "mapping", "ngs_mapping", "", ""),
    ]


def test_rules_lists_the_rules_of_the_project(tmp_path):
    """The names come from Snakemake, so this runs in a subprocess like the dry runs."""
    _write_project(tmp_path / "project", CHECK_FILES_CONFIG)

    result = _snappy(tmp_path / "project", "rules", "--plain")

    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert "mapping_ngs_mapping_bwa_run\tmapping\tngs_mapping\tbwa\trun" in lines
    assert "mapping_results\tmapping\tngs_mapping\t\t" in lines
    assert "all\t\t\t\t" in lines
    assert "[snappy]" not in result.stdout  # the library summary of the Snakefile is left out
    result = _snappy(tmp_path / "project", "rules", "--part", "qc_summary", "--names")
    assert result.stdout.splitlines() == ["mapping_ngs_mapping_qc_summary_run"]
    result = _snappy(tmp_path / "project", "rules", "--task", "genome", "--json")
    assert json.loads(result.stdout) == []


#: Two tasks whose contract check fails in ``load_project()`` (``task_produces``)
FDR_OF_BCFTOOLS_CALLS = """
  - name: calling
    step: variant_calling
    config:
      depends_on: {alignments: mapping, reference: genome}
      tool: bcftools_call
      bcftools_call: {assembly: GRCh38}
  - name: fdr
    step: variant_filtration
    config:
      depends_on: {variants: calling}
      tool: varlociraptor
      varlociraptor: {events: [present]}
"""


def _write_project(directory: Path, config: str) -> None:
    """Write ``config`` with the sheet, the reads and the reference files it names."""
    directory.mkdir()
    (directory / "config.yaml").write_text(config)
    (directory / "sheet.tsv").write_text(textwrap.dedent(CHECK_FILES_SHEET).lstrip())
    (directory / "raw" / "P001").mkdir(parents=True)
    (directory / "refs").mkdir()
    for name in ("raw/P001/P001_R1.fastq.gz", "raw/P001/P001_R2.fastq.gz"):
        (directory / name).touch()
    (directory / "reads.tsv").write_text(
        "library\treadgroup\tleft\tright\n"
        "P001-N1-DNA1-WGS1\tP001\traw/P001/P001_R1.fastq.gz\traw/P001/P001_R2.fastq.gz\n"
    )
    for name in ("genome.fa", "genome.fa.fai", "genome.dict"):
        (directory / "refs" / name).touch()


def _snappy(directory: Path, *args: str) -> subprocess.CompletedProcess[str]:
    """Run a ``snappy`` command in a subprocess, as Snakemake keeps global state."""
    env = os.environ.copy()
    root = str(Path(__file__).resolve().parents[3])
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [root, env.get("PYTHONPATH")]))
    # The project directory goes to snappy, before anything that goes to Snakemake after --
    cmd = [sys.executable, "-m", "snappy_pipeline.apps.snappy_cli", args[0]]
    cmd += ["--directory", str(directory), *args[1:]]
    return subprocess.run(cmd, text=True, capture_output=True, check=False, env=env)


def _snappy_run_dry(directory: Path) -> subprocess.CompletedProcess[str]:
    """Run ``snappy run -n`` in a subprocess."""
    return _snappy(directory, "run", "-n", "--", "--cores", "1", "--quiet")


@pytest.mark.parametrize(
    "config, message",
    [
        pytest.param(
            CHECK_FILES_CONFIG.replace(
                "      tool: bwa\n",
                "      library_selection: \"library_type = 'WGS'\"\n      tool: bwa\n",
            ),
            "[snappy] Error: Task 'mapping' (step ngs_mapping): library_selection expression "
            "\"library_type = 'WGS'\" is invalid: ",
            id="constructor",
        ),
        pytest.param(
            CHECK_FILES_CONFIG + FDR_OF_BCFTOOLS_CALLS,
            "[snappy] Error: Task 'fdr' (step variant_filtration): tool: varlociraptor controls "
            "the FDR of the calls of a variant_calling task with tool: varlociraptor, but "
            "depends_on.variants produces variants tagged germline, indel, snv",
            id="task_produces",
        ),
    ],
)
def test_run_prints_setup_errors_without_snakemake_framing(tmp_path, config, message):
    """A ``ConfigError`` of a workflow constructor or of ``task_produces`` is printed by the main
    Snakefile as one line; Snakemake adds neither its location nor traceback lines."""
    _write_project(tmp_path / "project", config)
    result = _snappy_run_dry(tmp_path / "project")
    output = result.stdout + result.stderr
    assert result.returncode == 1, output
    errors = [line for line in output.splitlines() if line.startswith("[snappy] Error:")]
    assert len(errors) == 1 and errors[0].startswith(message), output
    assert 'in file "' not in output and '  File "' not in output, output


def test_check_files_prints_configuration_errors_without_traceback(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG.replace("tool: bwa", "tool: bwa2"))
    result = CliRunner().invoke(snappy_cli.main, ["check-files", "--directory", str(tmp_path)])
    assert result.exit_code == 1
    assert result.exception is None or isinstance(result.exception, SystemExit)


def test_task_add_keeps_the_style_of_the_config_and_comments_the_new_task(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG.replace("tool: bwa", 'tool: "bwa"'))
    result = CliRunner().invoke(
        snappy_cli.main,
        [
            "task",
            "add",
            "--directory",
            str(tmp_path),
            "--config-mode",
            "commented",
            "qc=ngs_data_qc",
        ],
    )
    assert result.exit_code == 0, result.output
    text = (tmp_path / "config.yaml").read_text()
    # One line per task, no diff; the reminder once
    assert "Added task qc (step ngs_data_qc)" in result.output
    assert "Applying the following change" not in result.output
    assert result.output.count("REQUIRED fields") == 1
    # The user's quotes and sequence style stay, and there is no backup copy
    assert 'tool: "bwa"' in text
    assert not (tmp_path / "config.yaml.bak").exists()
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
        "reads",
        "mapping",
        "qc",
    ]


def test_task_add_full_writes_every_field_active(tmp_path):
    (tmp_path / "config.yaml").write_text("tasks: []\n")
    result = CliRunner().invoke(
        snappy_cli.main,
        ["task", "add", "--directory", str(tmp_path), "--config-mode", "full", "ngs_mapping"],
    )
    assert result.exit_code == 0, result.output
    lines = (tmp_path / "config.yaml").read_text().splitlines()
    assert "  - name: ngs_mapping" in lines and "    step: ngs_mapping" in lines
    # The optional fields are set, not commented out
    assert "      library_selection:" in lines and "      bwa:" in lines
    assert not any(line.lstrip().startswith("#library_selection") for line in lines)


def test_task_add_asks_for_the_tasks_on_a_terminal(tmp_path, answers):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG)
    answers(
        "variant_filtration",  # step
        "mapping",  # a taken name, asked again
        "filtration",
        "vembrane",  # tool
        snappy_cli.ADD_PROVIDER,  # variants: no task of the project calls variants yet
        "variant_calling",  # the step of the new task
        "calling",  # its name
        "mutect2",  # its tool
        "mapping",  # its alignments
        "genome",  # its reference
        False,  # no further task
    )

    result = CliRunner().invoke(
        snappy_cli.main, ["task", "add", "--interactive", "--directory", str(tmp_path)]
    )

    assert result.exit_code == 0, result.output
    assert answers.pending == []
    text = (tmp_path / "config.yaml").read_text()
    assert text.index("- name: calling\n") < text.index("- name: filtration\n")
    assert "tool: mutect2" in text and "alignments: mapping" in text and "reference: genome" in text
    assert "tool: vembrane" in text and "variants: calling" in text
    # The task table closes the questions
    assert "filtration\tvariant_filtration\tvembrane\tvariants: calling" in result.output


def test_task_add_without_tasks_needs_a_terminal(tmp_path):
    (tmp_path / "config.yaml").write_text(CHECK_FILES_CONFIG)

    result = CliRunner().invoke(snappy_cli.main, ["task", "add", "--directory", str(tmp_path)])

    assert result.exit_code == 2
    assert "runonaterminal" in re.sub(r"\x1b\[[0-9;]*m|[\s│╭╮╰╯─]", "", result.output)


def test_task_add_without_a_config_fails_with_a_hint(tmp_path):
    result = CliRunner().invoke(
        snappy_cli.main, ["task", "add", "--directory", str(tmp_path), "ngs_mapping"]
    )
    assert result.exit_code == 2
    # rich-click colors the output on a CI runner and wraps it in a box
    plain = re.sub(r"\x1b\[[0-9;]*m", "", result.output).replace("│", "")
    assert "createaprojectwithsnappyinit" in "".join(plain.split())


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
    sheets, tasks = lines.index("sample_sheets:"), lines.index("tasks:")
    assert sheets < tasks < lines.index("  - name: genome") < lines.index("  - name: mapping")
    assert not any(line.startswith("#") for line in lines[tasks:])
    config = snappy_cli.load_configfile(str(project / "config.yaml"))
    # The reads task of the template comes first
    assert [task["name"] for task in config["tasks"]] == ["reads", "genome", "mapping"]
    assert config["sample_sheets"] == ["samplesheet.tsv"]


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
    assert "sample_sheets:" in lines and "  - name: mapping" in lines
    assert "# REQUIRED" in text and "Options:" not in text
