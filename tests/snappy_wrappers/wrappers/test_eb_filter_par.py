# -*- coding: utf-8 -*-
"""Code for testing mutect2_par/run wrapper"""
from pathlib import Path
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import InputFiles, Log, OutputFiles, Params, Resources, Wildcards
from snakemake.script import Snakemake

from snappy_wrappers.wrappers.eb_filter_par.parallel_eb_filter import ParallelEasyBayesFilterWrapper

from .conftest import patch_module_fs


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
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_variant_filtration:
            tools_ngs_mapping: ['bwa']
            tools_somatic_variant_calling: ['mutect2']
            eb_filter:
              shuffle_seed: 1
              panel_of_normals_size: 25
              min_mapq: 20
              min_baseq: 15
              # Parallelization configuration
              window_length: 10000000   # split input into windows of this size, each triggers a job
              num_jobs: 500             # number of windows to process in parallel
              use_profile: true         # use Snakemake profile for parallel processing
              restart_times: 5          # number of times to re-launch jobs in case of failure
              max_jobs_per_second: 2    # throttling of job creation
              max_status_checks_per_second: 10   # throttling of status checks
              debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
              keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
              job_mult_memory: 1        # memory multiplier
              job_mult_time: 1          # running time multiplier
              merge_mult_memory: 1      # memory multiplier for merging
              merge_mult_time: 1        # running time multiplier for merging
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
    output_base_name = (
        "bwa.mutect2.jannovar_annotate_somatic_vcf.dkfz_bias_filter."
        "eb_filter.P001-T1-DNA1-WGS1/out/bwa.mutect2.jannovar_annotate_somatic_vcf."
        "dkfz_bias_filter.eb_filter.P001-T1-DNA1-WGS1"
    )
    return {
        "vcf": output_base_name + "vcf.gz",
        "vcf_md5": output_base_name + "vcf.gz.md5",
        "vcf_tbi": output_base_name + "vcf.gz.tbi",
        "vcf_tbi_md5": output_base_name + "vcf.gz.tbi.md5",
    }


@pytest.fixture
def snakemake_obj(minimal_config, snakemake_output_dict):
    """Returns Snakemake object."""
    # Define helper variables
    rule_name = "somatic_variant_calling_mutect2_run"
    threads = 2
    bench_iteration = 2
    script_dir = "/work"
    base_out = (
        "work/bwa.mutect2.jannovar_annotate_somatic_vcf.dkfz_bias_filter.P001-T1-DNA1-WGS1/out/"
        "bwa.mutect2.jannovar_annotate_somatic_vcf.dkfz_bias_filter.P001-T1-DNA1-WGS1"
    )
    input_dict = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "txt": "work/bwa.eb_filter.panel_of_normals/out/bwa.eb_filter.panel_of_normals.txt",
    }

    log_base_name = "work/bwa.mutect.P001-T1-DNA1-WGS1/log/bwa.mutect.P001-T1-DNA1-WGS1"
    log_dict = {
        "conda_info": log_base_name + ".conda_info.txt",
        "conda_info_md5": log_base_name + ".conda_info.txt.md5",
        "conda_list": log_base_name + ".conda_list.txt",
        "conda_list_md5": log_base_name + ".conda_list.txt.md5",
        "log": log_base_name + ".log",
        "log_md5": log_base_name + ".log.md5",
    }
    wildcards_dict = {"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"}
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


def test_eb_filter_wrapper_construct_parallel_rules(snakemake_obj, variant_caller_fake_fs, mocker):
    """Tests ParallelEasyBayesFilterWrapper.construct_parallel_rules()"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", variant_caller_fake_fs, mocker)
    wrapper_par = ParallelEasyBayesFilterWrapper(snakemake=snakemake_obj)
    # Define expected
    data_path = (Path(__file__).parent / "data/eb_filter_par_parallel.snakemake").resolve()
    with open(data_path, "r", encoding="utf8") as f:
        expected = f.read()
    # Get actual and assert if `rule chunk_0` is the same
    # Note: It is not feasible to test all chunks as the `wrapper` will be set to a local file
    _tmp_actual = list(wrapper_par.construct_parallel_rules())[0]
    actual = (
        "\n".join(
            [line for line in _tmp_actual.split("\n") if not ("wrapper" in line or line == "")]
        )
        + "\n"
    )
    assert actual == expected


def test_eb_filter_wrapper_construct_merge(snakemake_obj, variant_caller_fake_fs, mocker):
    """Tests ParallelEasyBayesFilterWrapper.construct_merge_rule()"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", variant_caller_fake_fs, mocker)
    wrapper_par = ParallelEasyBayesFilterWrapper(snakemake=snakemake_obj)
    # Define expected
    data_path = (Path(__file__).parent / "data/eb_filter_par_merge.snakemake").resolve()
    with open(data_path, "r", encoding="utf8") as f:
        expected = f.read()
    # Get actual and assert
    actual = wrapper_par.construct_merge_rule()
    print(actual)
    assert actual == expected
