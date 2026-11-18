# -*- coding: utf-8 -*-
"""Tests for ``snappy init`` subcommand."""

import shutil

import pytest
from biomedsheets import io_tsv
from click.testing import CliRunner
from snakemake.common.configfile import load_configfile

from snappy_pipeline.apps import snappy_cli
from snappy_pipeline.orchestration import load_project
from tests.snappy_pipeline.workflows.conftest import patch_module_fs


def test_init_with_positional_directory(germline_sheet_fake_noproject_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_noproject_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["init", "project"])
    assert result.exit_code == 0

    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/README.md")
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/project/config.yaml")


def test_init_default_cwd(germline_sheet_fake_noproject_fs, mocker):
    patch_module_fs("snappy_pipeline.apps.snappy_cli", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("snappy_pipeline.apps.impl.fsmanip", germline_sheet_fake_noproject_fs, mocker)
    patch_module_fs("shutil", germline_sheet_fake_noproject_fs, mocker)

    runner = CliRunner()
    result = runner.invoke(snappy_cli.main, ["init"])
    assert result.exit_code == 0

    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/README.md")
    assert germline_sheet_fake_noproject_fs.os.path.exists("/projects/config.yaml")


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
        ["init", str(project), "--no-git", "--sheet-type", sheet_type, "--project-id", "42"],
    )
    assert result.exit_code == 0, result.output
    # An empty sheet of the format, with the project as title
    with open(project / "samplesheet.tsv") as f:
        sheet = SHEET_READERS[sheet_type](f, "samplesheet.tsv")
    assert sheet.bio_entities == {} and sheet.json_data["title"] == "project"
    # The configuration loads and names the sheet with its type
    config = load_configfile(str(project / "config.yaml"))
    loaded = load_project(config, str(project))
    assert (loaded.model.project.name, loaded.model.project.id) == ("project", "42")
    data_set = loaded.model.data_sets["first_batch"]
    assert data_set.type == snappy_cli.SHEET_TYPES[sheet_type]
    assert data_set.file == str(project / "samplesheet.tsv")
    assert "- project id: 42" in (project / "README.md").read_text()


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
def test_init_asks_for_what_no_option_gives(tmp_path, git):
    project = tmp_path / "project"
    answers = [
        "My project",  # name
        "sodar-uuid",  # id
        "A trio",  # description
        "germline",  # sheet type
        "y",  # add tasks
        "mapping=ngs_mapping",
        "",  # no more tasks
        "y",  # track the sheet
        "git@example.org:lab/project.git",  # remote
        "n",  # the remote is not approved for the sheet
    ]
    result = CliRunner().invoke(
        snappy_cli.main, ["init", str(project), "--interactive"], input="\n".join(answers) + "\n"
    )
    assert result.exit_code == 0, result.output
    text = (project / "config.yaml").read_text()
    assert 'name: "My project"' in text and 'id: "sodar-uuid"' in text
    assert 'description: "A trio"' in text and "type: germline_variants" in text
    assert "\n  - name: mapping\n" in text
    # The remote is not approved for the sheet, so it stays untracked
    assert git(project, "ls-files").split() == [".gitignore", "README.md", "config.yaml"]
    assert git(project, "remote", "get-url", "origin").strip() == "git@example.org:lab/project.git"


def test_init_without_a_terminal_asks_nothing(tmp_path):
    project = tmp_path / "project"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--no-git"], input="x\n")
    assert result.exit_code == 0, result.output
    assert "Project name" not in result.output
    assert 'name: "project"' in (project / "config.yaml").read_text()
