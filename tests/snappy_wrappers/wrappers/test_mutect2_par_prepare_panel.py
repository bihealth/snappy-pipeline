# -*- coding: utf-8 -*-
"""Code for testing mutect2_par/run wrapper"""
import importlib.machinery
import os
import re
from pathlib import Path
import tempfile
import textwrap
import types

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import InputFiles, Log, OutputFiles, Params, Resources, Wildcards
from snakemake.script import Snakemake

from snappy_wrappers.wrappers.mutect2_par.prepare_panel.parallel_prepare_panel import (
    ParallelMutect2Wrapper,
)

from .conftest import mock_settings_env_vars, patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          cosmic:
            path: /path/to/cosmic.vcf.gz
          dbsnp:
            path: /path/to/dbsnp.vcf.gz
        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          panel_of_normals:
            tools: ['mutect2']
            mutect2:
              germline_resource: REQUIRED # Germline variants resource (same as panel of normals)
              # Parallelization configuration
              num_cores: 2              # number of cores to use locally
              window_length: 50000000   # split input into windows of this size, each triggers a job
              num_jobs: 500             # number of windows to process in parallel
              use_profile: true         # use Snakemake profile for parallel processing
              restart_times: 5          # number of times to re-launch jobs in case of failure
              max_jobs_per_second: 2    # throttling of job creation
              max_status_checks_per_second: 10   # throttling of status checks
              debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
              keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
              job_mult_memory: 2        # memory multiplier
              job_mult_time: 3          # running time multiplier
              merge_mult_memory: 4      # memory multiplier for merging
              merge_mult_time: 5        # running time multiplier for merging
              ignore_chroms:            # patterns of chromosome names to ignore
              - NC_007605    # herpes virus
              - hs37d5       # GRCh37 decoy
              - chrEBV       # Eppstein-Barr Virus
              - '*_decoy'    # decoy contig
              - 'HLA-*'      # HLA genes
              - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def snakemake_output_dict():
    """Returns dictionary that defined snakemake.output"""
    return {
        "vcf": "work/bwa.mutect2/out/bwa.mutect2.P001-N1-DNA1-WGS1.vcf.gz",
        "vcf_md5": "work/bwa.mutect2/out/bwa.mutect2.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "vcf_tbi": "work/bwa.mutect2/out/bwa.mutect2.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "vcf_tbi_md5": "work/bwa.mutect2/out/bwa.mutect2.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
    }


@pytest.fixture
def snakemake_obj(minimal_config, snakemake_output_dict):
    """Returns Snakemake object."""
    # Define helper variables
    rule_name = "panel_of_normals_mutect2_prepare_panel"
    threads = 2
    bench_iteration = 2
    script_dir = "/work"
    input_dict = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }

    log_base_name = "work/bwa.mutect2/log/bwa.mutect2.P001-N1-DNA1-WGS1"
    log_dict = {
        "conda_info": log_base_name + ".conda_info.txt",
        "conda_info_md5": log_base_name + ".conda_info.txt.md5",
        "conda_list": log_base_name + ".conda_list.txt",
        "conda_list_md5": log_base_name + ".conda_list.txt.md5",
        "log": log_base_name + ".log",
        "log_md5": log_base_name + ".log.md5",
    }
    wildcards_dict = {
        "mapper": "bwa",
        "normal_library": "P001-N1-DNA1-WGS1",
    }
    params_dict = {}

    # Define Snakemake class input
    input_ = InputFiles(fromdict=input_dict)
    output_ = OutputFiles(fromdict=snakemake_output_dict)
    params_ = Params(fromdict=params_dict)
    log_ = Log(fromdict=log_dict)
    wildcards_ = Wildcards(fromdict=wildcards_dict)
    resources_ = Resources(fromdict={})

    return Snakemake(
        rulename=rule_name,
        threads=threads,
        bench_iteration=bench_iteration,
        input_=input_,
        output=output_,
        log=log_,
        params=params_,
        wildcards=wildcards_,
        config=minimal_config,
        scriptdir=script_dir,
        resources=resources_,
    )


@pytest.fixture
def construct_preamble_module(
    snakemake_obj, variant_caller_fake_fs, mocker, mock_settings_env_vars
):
    """Returns ParallelMutect2Wrapper.construct_preamble() as a module"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", variant_caller_fake_fs, mocker)

    # Get methods as string
    wrapper_par = ParallelMutect2Wrapper(snakemake=snakemake_obj)
    preamble_str = wrapper_par.construct_preamble()

    # Remove Snakemake syntax
    clear_preamble_str = ""
    for line in preamble_str.split("\n"):
        if "rule all:" in line:
            continue
        elif "input: **" in line:
            continue
        elif "shell." in line:
            continue
        clear_preamble_str += line + "\n"

    # Push content to temp file
    tmp = tempfile.NamedTemporaryFile(suffix=".py")
    with open(tmp.name, "w") as f:
        f.write(clear_preamble_str)

    # Load as module
    module_name = os.path.basename(tmp.name).replace(".py", "")
    loader = importlib.machinery.SourceFileLoader(module_name, tmp.name)
    module = types.ModuleType(loader.name)
    loader.exec_module(module)
    return module


def test_mutect2_wrapper_prepare_panel_joint_chunks(snakemake_obj, variant_caller_fake_fs, mocker):
    """Tests ParallelMutect2Wrapper.joint_chunks()"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", variant_caller_fake_fs, mocker)
    wrapper_par = ParallelMutect2Wrapper(snakemake=snakemake_obj)
    # Force level_one
    wrapper_par.merge_block_size = 3
    # Define expected
    data_path = (Path(__file__).parent / "data/mutect2_par_prepare_panel.snakemake").resolve()
    with open(data_path, "r", encoding="utf8") as f:
        expected = f.read()
    # Remove wrapper lines that contain path to script
    remove_wrapper = re.compile(" +wrapper: [^\n]+\n")
    # Get actual and assert
    actual = remove_wrapper.sub("", wrapper_par.joint_chunks())
    assert actual == expected


def test_mutect2_wrapper_prepare_panel_preamble_multiply_time(construct_preamble_module):
    """Tests Parallel Preamble multiply_time()"""
    # Constant factor
    factor = 10
    # Define (input, expected) pair
    input_expected_pairs = (
        ("00:01:00", "0-00:10:00"),
        ("01:00:00", "0-10:00:00"),
        ("12:00:00", "5-00:00:00"),
        ("01-00:00:00", "10-00:00:00"),
    )

    # Test all pairs
    for pair in input_expected_pairs:
        _input = pair[0]
        _expected = pair[1]
        actual = construct_preamble_module.multiply_time(day_time_str=_input, factor=factor)
        msg = (
            f"For input '{_input}' * {factor} expected output '{_expected}'. "
            f"Received instead: {actual}"
        )
        assert actual == _expected, msg

    # Test invalid time
    with pytest.raises(ValueError):
        construct_preamble_module.multiply_time(day_time_str="_not_a_valid_time", factor=factor)


def test_mutect2_wrapper_prepare_panel_preamble_multiply_memory(construct_preamble_module):
    """Tests Parallel Preamble multiply_memory()"""
    # Constant factor
    factor = 10
    # Define (input, expected) pair
    input_expected_pairs = (
        ("1024k", 10),
        ("1M", 10),
        ("16G", 160000),
        ("1T", 10000000),
        ("12", 120),
    )
    # Test all pairs
    for pair in input_expected_pairs:
        _input = pair[0]
        _expected = pair[1]
        actual = construct_preamble_module.multiply_memory(memory_str=_input, factor=factor)
        msg = (
            f"For input '{_input}' * {factor} expected output '{_expected}'. "
            f"Received instead: {actual}"
        )
        assert actual == _expected, msg


def test_mutect2_wrapper_prepare_panel_preamble_resource_chunk_threads(construct_preamble_module):
    """Tests Parallel Preamble resource_chunk_threads() - Chunks always get a single thread"""
    expected = 1
    actual = construct_preamble_module.resource_chunk_threads(wildcards=None)
    assert actual == expected


def test_mutect2_wrapper_prepare_panel_preamble_resource_chunk_memory(construct_preamble_module):
    """
    Tests Parallel Preamble resource_chunk_memory() - Baseline memory defined in Snakemake object
    """
    # Define (input, expected) pair
    input_expected_pairs = (
        (1, 28000),
        (2, 56000),
        (3, 84000),
        (4, 112000),
        (5, 140000),
    )
    # Test all pairs
    for pair in input_expected_pairs:
        _input = pair[0]
        _expected = pair[1]
        actual = construct_preamble_module.resource_chunk_memory(wildcards=None, attempt=_input)
        msg = f"For input '{_input}' expected output '{_expected}'. " f"Received instead: {actual}"
        assert actual == _expected, msg


def test_mutect2_wrapper_prepare_panel_preamble_resource_chunk_time(construct_preamble_module):
    """
    Tests Parallel Preamble resource_chunk_time() - Baseline time defined in Snakemake object
    """
    # Define (input, expected) pair
    input_expected_pairs = (
        (1, "0-12:00:00"),
        (2, "1-00:00:00"),
        (3, "1-12:00:00"),
        (4, "2-00:00:00"),
        (5, "2-12:00:00"),
    )
    # Test all pairs
    for pair in input_expected_pairs:
        _input = pair[0]
        _expected = pair[1]
        actual = construct_preamble_module.resource_chunk_time(wildcards=None, attempt=_input)
        msg = f"For input '{_input}' expected output '{_expected}'. " f"Received instead: {actual}"
        assert actual == _expected, msg


def test_mutect2_wrapper_prepare_panel_preamble_resource_chunk_partition(construct_preamble_module):
    """Tests Parallel Preamble resource_chunk_partition()"""
    expected = "medium"
    actual = construct_preamble_module.resource_chunk_partition(wildcards=None)
    assert actual == expected


def test_mutect2_wrapper_prepare_panel_preamble_resource_merge_threads(construct_preamble_module):
    """Tests Parallel Preamble resource_merge_threads() - Merge always get a single thread"""
    expected = 1
    actual = construct_preamble_module.resource_merge_threads(wildcards=None)
    assert actual == expected


def test_mutect2_wrapper_prepare_panel_preamble_resource_merge_memory(construct_preamble_module):
    """Tests Parallel Preamble resource_merge_memory() - Merge always get the same value"""
    expected = "8G"
    actual = construct_preamble_module.resource_merge_memory(wildcards=None)
    assert actual == expected


def test_mutect2_wrapper_prepare_panel_preamble_resource_merge_time(construct_preamble_module):
    """Tests Parallel Preamble resource_merge_time() - Merge always get the same value"""
    expected = "20:00:00"
    actual = construct_preamble_module.resource_merge_time(wildcards=None)
    assert actual == expected


def test_mutect2_wrapper_prepare_panel_preamble_resource_merge_partition(construct_preamble_module):
    """Tests Parallel Preamble resource_merge_partition()"""
    expected = "medium"
    actual = construct_preamble_module.resource_merge_partition(wildcards=None)
    assert actual == expected
