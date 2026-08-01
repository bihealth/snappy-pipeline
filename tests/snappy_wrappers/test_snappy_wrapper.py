# -*- coding: utf-8 -*-
"""Tests for the base wrapper classes in snappy_wrappers.snappy_wrapper."""

import os

import pytest

from snappy_wrappers.snappy_wrapper import PythonWrapper, SnappyWrapper


class _Log:
    def __init__(self, path):
        self.log = path
        self.conda_list = path + ".conda_list.txt"
        self.conda_info = path + ".conda_info.txt"


class _FakeSnakemake:
    def __init__(self, log_path):
        self.log = _Log(log_path)
        self.input = None
        self.output = None


@pytest.mark.parametrize(
    "log_path",
    ["/project-dir/logs/rule/sample/rule.log", "relative/logs/rule.log"],
)
def test_header_always_tees_into_snakemake_log(log_path):
    """The wrapper header must always tee into the Snakemake log file.

    Regression test: the header previously only piped output into the log file
    when a TTY was present, so under the SLURM executor the Snakemake-declared
    log stayed empty while the output went to the SLURM job log only.
    """
    snakemake = _FakeSnakemake(log_path)
    header = SnappyWrapper.header.format(snakemake=snakemake)
    assert f'rm -f "{log_path}"' in header
    assert f'tee -a "{log_path}"' in header
    assert "tty" not in header.lower()


def test_python_wrapper_logging_context_redirects_output(tmp_path):
    log_path = str(tmp_path / "work" / "logs" / "rule.log")
    snakemake = _FakeSnakemake(log_path)
    wrapper = PythonWrapper(snakemake)

    with wrapper.logging_context():
        print("hello stdout")
        print("hello stderr", file=__import__("sys").stderr)

    assert os.path.exists(log_path)
    with open(log_path, "rt") as f:
        content = f.read()
    assert "hello stdout" in content
    assert "hello stderr" in content


def test_python_wrapper_logging_context_appends(tmp_path):
    log_path = str(tmp_path / "rule.log")
    snakemake = _FakeSnakemake(log_path)
    wrapper = PythonWrapper(snakemake)

    with wrapper.logging_context():
        print("first run")

    with wrapper.logging_context():
        print("second run")

    with open(log_path, "rt") as f:
        content = f.read()
    assert content.count("first run") == 1
    assert content.count("second run") == 1


def test_python_wrapper_logging_context_requires_log(tmp_path):
    wrapper = PythonWrapper(_FakeSnakemake(str(tmp_path / "rule.log")))
    wrapper.snakemake.log = None
    with pytest.raises(AttributeError):
        with wrapper.logging_context():
            pass


def test_python_wrapper_compute_log_md5(mocker, tmp_path):
    log_path = str(tmp_path / "rule.log")
    wrapper = PythonWrapper(_FakeSnakemake(log_path))
    shell_mock = mocker.patch("snappy_wrappers.snappy_wrapper.shell")

    wrapper.compute_log_md5()

    assert shell_mock.call_count == 1
    cmd = shell_mock.call_args.args[0]
    assert os.path.realpath(log_path) in cmd
    assert ".md5" in cmd


def test_python_wrapper_write_conda_info(mocker, tmp_path):
    log_path = str(tmp_path / "rule.log")
    wrapper = PythonWrapper(_FakeSnakemake(log_path))
    shell_mock = mocker.patch("snappy_wrappers.snappy_wrapper.shell")

    wrapper.write_conda_info()

    cmds = [call.args[0] for call in shell_mock.call_args_list]
    assert len(cmds) == 4  # conda list, list.md5, conda info, info.md5
    assert cmds[0] == "conda list > {}.conda_list.txt".format(log_path)
    assert cmds[2] == "conda info > {}.conda_info.txt".format(log_path)
    assert "md5sum $f > $f.md5" in cmds[1]
    assert "md5sum $f > $f.md5" in cmds[3]
