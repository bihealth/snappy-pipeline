# -*- coding: utf-8 -*-
"""Tests for ``snappy init`` subcommand."""

import re
import shutil

import pytest
from biomedsheets import io_tsv
from click.testing import CliRunner
from snakemake.common.configfile import load_configfile

from snappy_pipeline.apps import snappy_cli
from snappy_pipeline.orchestration import load_project


def test_init_default_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(snappy_cli.main, ["init", "--no-git"])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "README.md").exists() and (tmp_path / "config.yaml").exists()


def test_init_puts_a_named_project_into_a_directory_of_its_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(snappy_cli.main, ["init", "--no-git", "--project-name", "trio"])

    assert result.exit_code == 0, result.output
    assert (tmp_path / "trio" / "config.yaml").exists()
    assert not (tmp_path / "config.yaml").exists()
    assert 'name: "trio"' in (tmp_path / "trio" / "config.yaml").read_text()
    # The current directory, if it is named like the project; a given directory always
    monkeypatch.chdir(tmp_path / "trio")
    (tmp_path / "trio" / "config.yaml").unlink()
    result = CliRunner().invoke(snappy_cli.main, ["init", "--no-git", "--project-name", "trio"])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "trio" / "config.yaml").exists() and not (
        tmp_path / "trio" / "trio"
    ).exists()
    result = CliRunner().invoke(
        snappy_cli.main, ["init", str(tmp_path / "other"), "--no-git", "--project-name", "trio"]
    )
    assert result.exit_code == 0, result.output
    assert 'name: "trio"' in (tmp_path / "other" / "config.yaml").read_text()


def test_init_fails_on_a_duplicate_task_but_keeps_the_others(tmp_path):
    project = tmp_path / "project"
    result = CliRunner().invoke(
        snappy_cli.main,
        ["init", str(project), "--no-git", "--task", "ngs_mapping", "--task", "ngs_mapping"],
    )
    assert result.exit_code == 1
    config = load_configfile(str(project / "config.yaml"))
    assert [task["name"] for task in config["tasks"]] == ["reads", "ngs_mapping"]


needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


@needs_git
def test_init_tracks_only_the_configuration_in_git(tmp_path, git):
    project = tmp_path / "project"
    result = CliRunner().invoke(
        snappy_cli.main, ["init", str(project), "--task", "mapping=ngs_mapping"]
    )
    assert result.exit_code == 0, result.output
    assert git(project, "ls-files").split() == [".gitignore", "README.md", "config.yaml"]
    assert git(project, "log", "--format=%s").splitlines() == ["Initialize the snappy project"]
    # The task added by init is part of the commit
    assert git(project, "status", "--porcelain") == ""
    assert "mapping" in (project / "config.yaml").read_text()
    assert git(project, "remote") == ""
    # Results and sample sheets stay out
    (project / "samplesheet.tsv").write_text("")
    (project / "output").mkdir()
    (project / "output" / "result.txt").write_text("")
    ignored = git(project, "status", "--porcelain", "--ignored").splitlines()
    assert ignored == ["!! output/", "!! samplesheet.tsv"]


@needs_git
def test_init_adds_the_remote_without_pushing(tmp_path, git):
    project = tmp_path / "project"
    url = "git@example.org:lab/project.git"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--remote", url])
    assert result.exit_code == 0, result.output
    assert git(project, "remote", "get-url", "origin").strip() == url


@needs_git
def test_init_inside_a_repository_stages_there(tmp_path, git):
    git(tmp_path, "init", "-q")
    project = tmp_path / "project"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project)])
    assert result.exit_code == 0, result.output
    assert not (project / ".git").exists()
    staged = git(tmp_path, "diff", "--cached", "--name-only").split()
    assert staged == ["project/.gitignore", "project/README.md", "project/config.yaml"]


def test_init_without_git(tmp_path):
    project = tmp_path / "project"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--no-git"])
    assert result.exit_code == 0, result.output
    assert not (project / ".git").exists() and not (project / ".gitignore").exists()


SHEET_READERS = {
    "generic": io_tsv.read_generic_tsv_sheet,
    "germline": io_tsv.read_germline_tsv_sheet,
    "cancer": io_tsv.read_cancer_tsv_sheet,
}


@pytest.mark.parametrize("sheet_type", ["generic", "germline", "cancer"])
def test_init_creates_a_sample_sheet_that_the_config_names(tmp_path, sheet_type):
    project = tmp_path / "project"
    result = CliRunner().invoke(
        snappy_cli.main,
        ["init", str(project), "--no-git", "--sheet-type", sheet_type, "--external-id", "42"],
    )
    assert result.exit_code == 0, result.output
    # An empty sheet of the format, with the project as title
    with open(project / "samplesheet.tsv") as f:
        sheet = SHEET_READERS[sheet_type](f, "samplesheet.tsv")
    assert sheet.bio_entities == {} and sheet.json_data["title"] == "project"
    # The batch is a string field of the library, with a column in the data header
    assert sheet.json_data["extraInfoDefs"]["ngsLibrary"]["batchNo"]["type"] == "string"
    header = (project / "samplesheet.tsv").read_text().partition("[Data]\n")[2].splitlines()[0]
    assert "\tfolderName\tbatchNo\t" in header
    # The configuration loads and names the sheet with its type
    config = load_configfile(str(project / "config.yaml"))
    loaded = load_project(config, str(project))
    assert (loaded.model.project.name, loaded.model.project.external_id) == ("project", "42")
    data_set = loaded.model.data_sets["first_batch"]
    assert data_set.type == snappy_cli.SHEET_TYPES[sheet_type]
    assert data_set.file == str(project / "samplesheet.tsv")
    assert "- external id: 42" in (project / "README.md").read_text()


def test_init_keeps_an_existing_sample_sheet(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    (project / "samplesheet.tsv").write_text("mine\n")
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--no-git"])
    assert result.exit_code == 0, result.output
    assert (project / "samplesheet.tsv").read_text() == "mine\n"


@needs_git
def test_init_tracks_the_sample_sheet_only_on_request(tmp_path, git):
    project = tmp_path / "project"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--track-sheet"])
    assert result.exit_code == 0, result.output
    assert git(project, "ls-files").split() == [
        ".gitignore",
        "README.md",
        "config.yaml",
        "samplesheet.tsv",
    ]
    assert "!/samplesheet.tsv\n" in (project / ".gitignore").read_text()


@needs_git
def test_init_asks_for_what_no_option_gives(tmp_path, git, answers):
    project = tmp_path / "project"
    answers(
        "My project",  # name
        "sodar-uuid",  # external id
        "A trio",  # description
        "germline",  # sheet type
        "minimal",  # configuration style
        True,  # track the sheet
        "git@example.org:lab/project.git",  # remote
        False,  # the remote is not approved for the sheet
        True,  # add tasks now
        "ngs_mapping",  # step
        "mapping",  # name
        "bwa_mem3",  # tool
        "reads",  # reads: the external_data task of the template
        snappy_cli.LEAVE_EMPTY,  # index: nothing provides it yet, filled in later
        snappy_cli.ADD_PROVIDER,  # reference: add a task that provides it now
        "external_data",  # its step
        "genome",  # its name
        False,  # no further task
    )

    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--interactive"])

    assert result.exit_code == 0, result.output
    assert answers.pending == []
    assert answers.asked[3:5] == ["Sample sheet type", "Configuration style"]
    assert answers.asked[8:12] == ["Add tasks now?", "Step", "Task name", "Tool"]
    assert [q.split(":")[0] for q in answers.asked[12:15]] == ["reads", "index", "reference"]
    assert answers.asked[15:17] == ["Step of the task", "Task name"]
    text = (project / "config.yaml").read_text()
    assert 'name: "My project"' in text and 'external_id: "sodar-uuid"' in text
    assert 'description: "A trio"' in text and "type: germline_variants" in text
    # The genome task was added while the mapping task was asked for, so it comes first
    assert text.index("- name: genome\n") < text.index("- name: mapping\n")
    assert "type: raw\n" in text and "- dna\n" in text and "- reference\n" in text
    assert re.search(r"\n      files:\n        fasta:\s+# REQUIRED\n", text)  # to fill in
    assert "tool: bwa_mem3\n" in text and "reads: reads\n" in text
    assert "reference: genome\n" in text
    assert re.search(r"\n\s+index:\s+# REQUIRED", text)  # left to be filled
    # The block is that of bwa_mem3: its section, no keys that only STAR reads
    assert "\n      bwa_mem3: {}\n" in text and "features" not in text and "exons" not in text
    # The remote is not approved for the sheet, so it stays untracked
    assert git(project, "ls-files").split() == [".gitignore", "README.md", "config.yaml"]
    assert git(project, "remote", "get-url", "origin").strip() == "git@example.org:lab/project.git"


def test_init_without_a_terminal_asks_nothing(tmp_path):
    project = tmp_path / "project"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--no-git"], input="x\n")
    assert result.exit_code == 0, result.output
    assert "Project name" not in result.output
    assert 'name: "project"' in (project / "config.yaml").read_text()
