# -*- coding: utf-8 -*-
"""Code for testing mutect_par/run wrapper"""
from pathlib import Path

from snappy_wrappers.wrappers.mutect2_par.run.parallel_mutect2_wrapper import ParallelMutect2Wrapper

from .conftest import patch_module_fs


def test_mutect2_wrapper_run_construct_merge_rule(snakemake_obj, somatic_variant_fake_fs, mocker):
    """Tests ParallelMutect2Wrapper.construct_merge_rule()"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", somatic_variant_fake_fs, mocker)
    wrapper_par = ParallelMutect2Wrapper(snakemake=snakemake_obj)
    # Define expected
    data_path = (Path(__file__).parent / "data/mutect2_par.snakemake").resolve()
    with open(data_path, "r") as f:
        expected = f.read()
    # Get actual and assert
    actual = wrapper_par.construct_merge_rule()
    assert actual == expected
