# -*- coding: utf-8 -*-
"""Tests for ``snappy init`` subcommand."""

import os
import shutil
import subprocess

import pytest
from click.testing import CliRunner

from snappy_pipeline.apps import snappy_cli
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


GIT_ENV = {
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_AUTHOR_NAME": "snappy",
    "GIT_AUTHOR_EMAIL": "snappy@example.org",
    "GIT_COMMITTER_NAME": "snappy",
    "GIT_COMMITTER_EMAIL": "snappy@example.org",
}

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


@pytest.fixture
def git_env(monkeypatch):
    for key, value in GIT_ENV.items():
        monkeypatch.setenv(key, value)


def git(directory, *args) -> str:
    return subprocess.run(
        ["git", "-C", str(directory), *args], capture_output=True, text=True, check=True
    ).stdout


@needs_git
def test_init_tracks_only_the_configuration_in_git(tmp_path, git_env):
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
def test_init_adds_the_remote_without_pushing(tmp_path, git_env):
    project = tmp_path / "project"
    url = "git@example.org:lab/project.git"
    result = CliRunner().invoke(snappy_cli.main, ["init", str(project), "--remote", url])
    assert result.exit_code == 0, result.output
    assert git(project, "remote", "get-url", "origin").strip() == url


@needs_git
def test_init_inside_a_repository_stages_there(tmp_path, git_env):
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
