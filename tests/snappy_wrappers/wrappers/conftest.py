"""Support code for the snappy pipeline wrapper tests."""

from collections import namedtuple
from itertools import chain
import os
import shutil
import subprocess
from unittest.mock import MagicMock

from pyfakefs import fake_filesystem
import pytest
import yaml

FORCE_RUN = os.environ.get("FORCE_RUN", "false") == "true"
DIFF_MASTER = os.environ.get("DIFF_MASTER", "false") == "true"
DIFF_LAST_COMMIT = os.environ.get("DIFF_LAST_COMMIT", "false") == "true"

if FORCE_RUN or DIFF_MASTER or DIFF_LAST_COMMIT:
    compare = "HEAD^" if DIFF_LAST_COMMIT else "origin/master"

    # check if wrapper is modified compared to master
    DIFF_FILES = set(
        subprocess.check_output(["git", "diff", compare, "--name-only"]).decode().split("\n")
    )

#: Allow running containerized (local only).
CONTAINERIZED = os.environ.get("CONTAINERIZED", "false") == "true"


def run_workflow(wrapper, test_dir, cmd, tmpdir, check_log=None):
    d = str(tmpdir)
    origdir = os.getcwd()
    dst = os.path.join(d, "master")
    os.makedirs(dst, exist_ok=True)

    def copy(pth, src):
        shutil.copy(os.path.join(pth, src), os.path.join(dst, pth))

    used_wrappers = []
    wrapper_file = "used_wrappers.yaml"
    if os.path.exists(os.path.join(wrapper, wrapper_file)):
        # is meta wrapper
        with open(os.path.join(wrapper, wrapper_file), "r") as wf:
            wf = yaml.safe_load(wf)
            used_wrappers = wf["wrappers"]
    else:
        used_wrappers.append(wrapper)

    for w in used_wrappers:
        success = False
        for ext in ("py", "R", "Rmd"):
            script = "wrapper." + ext
            if os.path.exists(os.path.join(w, script)):
                os.makedirs(os.path.join(dst, w), exist_ok=True)
                copy(w, script)
                success = True
                break
        assert success, "No wrapper script found for {}".format(w)
        copy(w, "environment.yaml")

    if (
        not FORCE_RUN
        and (DIFF_MASTER or DIFF_LAST_COMMIT)
        and not any(
            any(f.startswith(w) for f in DIFF_FILES) for w in chain(used_wrappers, [wrapper])
        )
    ):
        raise Skipped("wrappers not modified")

    testdir = os.path.join(d, "test")
    # pkgdir = os.path.join(d, "pkgs")
    shutil.copytree(os.path.join(os.path.dirname(__file__), test_dir), testdir)
    # prepare conda package dir
    # os.makedirs(pkgdir)
    # switch to test directory
    os.chdir(testdir)
    if os.path.exists(".snakemake"):
        shutil.rmtree(".snakemake")
    cmd = cmd + ["--wrapper-prefix", "file://{}/master/".format(d), "--conda-cleanup-pkgs"]

    if CONTAINERIZED:
        # run snakemake in container
        cmd = [
            "sudo",
            "docker",
            "run",
            "-it",
            "-v",
            "{}:{}".format(os.getcwd(), "/workdir/test"),
            "-w",
            "/workdir",
            "snakemake/snakemake",
            " ".join(cmd),
        ]

    # env = dict(os.environ)
    # env["CONDA_PKGS_DIRS"] = pkgdir
    try:
        subprocess.check_call(cmd)
    except Exception as e:
        # go back to original directory
        os.chdir(origdir)
        logfiles = [
            os.path.join(d, f)
            for d, _, files in os.walk(os.path.join(testdir, "logs"))
            for f in files
        ]
        for path in logfiles:
            with open(path) as f:
                msg = "###### Logfile: " + path + " ######"
                print(msg, "\n")
                print(f.read())
                print("#" * len(msg))
        if check_log is not None:
            for f in logfiles:
                check_log(open(f).read())
        else:
            raise e
    finally:
        # cleanup environments to save disk space
        subprocess.check_call(
            r"for env in `conda env list | grep -P '\.snakemake/conda' | "
            r"cut -f1 | tr -d ' '`; do conda env remove --prefix $env; done",
            shell=True,  # nosec
        )
        # go back to original directory
        os.chdir(origdir)


@pytest.fixture
def fai_file_content():
    """Returns FAI file content based on hs37d5 (chromosome 1 only)."""
    return "1\t249250621\t52\t60\t61"


@pytest.fixture
def fake_fs():
    """Return ``namedtuple`` with fake file system objects."""
    klass = namedtuple("FakeFsBundle", "fs os open inter_process_lock")
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    fake_lock = MagicMock()
    return klass(fs=fake_fs, os=fake_os, open=fake_open, inter_process_lock=fake_lock)


@pytest.fixture
def variant_caller_fake_fs(fake_fs, fai_file_content):
    """Return fake file system setup with files for the somatic variant calling workflow."""
    # Create work directory
    fake_fs.fs.makedirs("/work", exist_ok=True)
    # Create static files
    fake_fs.fs.create_file("/path/to/ref.fa", create_missing_dirs=True)
    fake_fs.fs.create_file("/path/to/ref.fa.fai", contents=fai_file_content)
    return fake_fs


def patch_module_fs(module_name, fake_fs, mocker):
    """Helper function to mock out the file-system related things in the module with the given
    name using the given fake_fs and pytest-mock mocker
    """
    mocker.patch("{}.os".format(module_name), fake_fs.os)
    mocker.patch("{}.open".format(module_name), fake_fs.open, create=True)
    mocker.patch("{}.os".format(module_name), fake_fs.os)
