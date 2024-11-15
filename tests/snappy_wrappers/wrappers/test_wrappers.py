import pytest

from .conftest import run_workflow

skip = pytest.mark.xfail(run=False)


@skip
def test_bwa_mem(tmpdir):
    run_workflow(
        "snappy_wrappers/wrappers/bwa",
        "bwa_mem_pe",
        [
            "snakemake",
            "--cores",
            "1",
            "--software-deployment-method",
            "conda",
            "--conda-frontend",
            "mamba",
        ],
        tmpdir=tmpdir,
    )
