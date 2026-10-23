# -*- coding: utf-8 -*-
"""Tests for the base wrapper classes in snappy_wrappers.snappy_wrapper."""

import hashlib
import os

import pytest

from snappy_wrappers.snappy_wrapper import (
    PythonWrapper,
    ShellWrapper,
    SnappyWrapper,
    write_declared_md5s,
)


class _Log:
    def __init__(self, path):
        self.log = path


class _FakeSnakemake:
    def __init__(self, log_path):
        self.log = _Log(log_path)
        self.input = None
        self.output = None
        self.params = None
        self.threads = 4


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
    header = SnappyWrapper.header.format(snakemake=snakemake, tmpdir="/tmp/job")
    assert f'rm -f "{log_path}"' in header
    assert f'tee -a "{log_path}"' in header
    assert "tty" not in header.lower()


def test_header_gives_openmp_and_blas_the_threads_of_the_job():
    # Snakemake's shell() sets them to None when a wrapper calls it
    header = SnappyWrapper.header.format(snakemake=_FakeSnakemake("rule.log"), tmpdir="/tmp/job")
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        assert f"{variable}=4" in header


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


def _md5sum(path):
    return f"{hashlib.md5(path.read_bytes()).hexdigest()}  {path.name}\n"


def test_write_declared_md5s_hashes_only_outputs_with_a_declared_md5(tmp_path):
    result, intermediate = tmp_path / "out.vcf.gz", tmp_path / "out.tsv"
    result.write_text("result\n")
    intermediate.write_text("intermediate\n")
    snakemake = _FakeSnakemake(str(tmp_path / "rule.log"))
    snakemake.output = [str(result), str(result) + ".md5", str(intermediate)]

    write_declared_md5s(snakemake)

    assert (tmp_path / "out.vcf.gz.md5").read_text() == _md5sum(result)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["out.tsv", "out.vcf.gz", "out.vcf.gz.md5"]


def test_write_declared_md5s_replaces_a_link(tmp_path):
    upstream, output = tmp_path / "in.vcf.gz.md5", tmp_path / "out.vcf.gz"
    upstream.write_text("upstream\n")
    output.write_text("filtered\n")
    (tmp_path / "out.vcf.gz.md5").symlink_to(upstream)
    snakemake = _FakeSnakemake(str(tmp_path / "rule.log"))
    snakemake.output = [str(output), str(output) + ".md5"]

    write_declared_md5s(snakemake)

    assert upstream.read_text() == "upstream\n"
    assert not (tmp_path / "out.vcf.gz.md5").is_symlink()
    assert (tmp_path / "out.vcf.gz.md5").read_text() == _md5sum(output)


def test_python_wrapper_run_logs_output_and_writes_md5_files(tmp_path):
    log_path = str(tmp_path / "rule.log")
    (tmp_path / "rule.log").write_text("previous job\n")
    output = tmp_path / "out.txt"
    snakemake = _FakeSnakemake(log_path)
    snakemake.output = [str(output), str(tmp_path / "out.txt.md5")]

    def main():
        print("hello from main")
        output.write_text("result\n")

    PythonWrapper(snakemake).run(main)

    assert (tmp_path / "rule.log").read_text() == "hello from main\n"
    assert (tmp_path / "out.txt.md5").read_text() == _md5sum(output)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["out.txt", "out.txt.md5", "rule.log"]


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


def test_shell_wrapper_removes_its_exported_temp_directory(tmp_path):
    """The command sees an exported ``TMPDIR``, which is gone after the run, also when the
    command sets its own ``EXIT`` trap."""
    snakemake = _FakeSnakemake(str(tmp_path / "rule.log"))
    snakemake.output = []
    seen = tmp_path / "seen.txt"
    ShellWrapper(snakemake).run(
        f'trap "true" EXIT\nsh -c \'echo "$TMPDIR"\' > {seen}\ntouch "$TMPDIR/scratch"'
    )
    tmpdir = seen.read_text().strip()
    assert tmpdir
    assert not os.path.exists(tmpdir)
