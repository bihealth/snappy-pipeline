# from .conftest import run_workflow, skip_if_not_modified
#
#
# @skip_if_not_modified
# def test_bwa_mem(tmpdir):
#     run_workflow(
#         "snappy_wrappers/wrappers/bwa",
#         "bwa_mem_pe",
#         ["snakemake", "--cores", "1", "--use-conda", "--conda-frontend", "mamba"],
#         tmpdir=tmpdir,
#     )
