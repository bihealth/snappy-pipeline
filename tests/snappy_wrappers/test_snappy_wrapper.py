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


def test_python_wrapper_run_logs_output_and_writes_md5_files(mocker, tmp_path):
    log_path = str(tmp_path / "rule.log")
    (tmp_path / "rule.log").write_text("previous job\n")
    output = tmp_path / "out.txt"
    snakemake = _FakeSnakemake(log_path)
    snakemake.output = [str(output), str(tmp_path / "out.txt.md5")]
    shell_mock = mocker.patch("snappy_wrappers.snappy_wrapper.shell")

    def main():
        print("hello from main")
        output.write_text("result\n")

    PythonWrapper(snakemake).run(main)

    assert (tmp_path / "rule.log").read_text() == "hello from main\n"
    cmds = [call.args[0] for call in shell_mock.call_args_list]
    assert cmds[0].startswith("conda list") and cmds[2].startswith("conda info")
    md5_of = [cmd for cmd in cmds if "md5sum $f > $f.md5" in cmd]
    assert len(md5_of) == 4  # conda_list, conda_info, out.txt, log
    assert str(output.name) in md5_of[2] and os.path.realpath(log_path) in md5_of[3]


def test_python_wrapper_logging_context_captures_subprocess_output(tmp_path):
    import subprocess

    log_path = str(tmp_path / "rule.log")
    with PythonWrapper(_FakeSnakemake(log_path)).logging_context():
        subprocess.run(["echo", "from a subprocess"], check=True)

    assert "from a subprocess" in (tmp_path / "rule.log").read_text()


def test_base_classes_run_on_the_oldest_python_that_snakemake_uses_for_wrappers():
    """Snakemake runs a wrapper with its environment's Python from MIN_PY_VERSION on."""
    import ast
    import pathlib

    from snakemake.script import MIN_PY_VERSION

    import snappy_wrappers.snappy_wrapper as module

    ast.parse(pathlib.Path(module.__file__).read_text(), feature_version=MIN_PY_VERSION)
