# -*- coding: utf-8 -*-
"""Code for testing gatk_hc_par/run wrapper"""
from pathlib import Path
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import InputFiles, Log, OutputFiles, Params, Resources, Wildcards
from snakemake.script import Snakemake

from snappy_wrappers.wrappers.gatk_hc_par.parallel_gatk_hc import ParallelGatkHcWrapper

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

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          variant_calling:
            tools:
            - gatk
            gatk_hc:
              # Parallelization configuration
              num_cores: 2               # number of cores to use locally
              window_length: 3500000     # split input into windows, each triggers a job
              num_jobs: 500              # number of windows to process in parallel
              restart_times: 5           # number of times to re-launch jobs in case of failure
              max_jobs_per_second: 2     # throttling of job creation
              max_status_checks_per_second: 10   # throttling of status checks
              debug_trunc_tokens: 0      # truncation to first N tokens (0 for none)
              keep_tmpdir: never         # keep temporary directory, {always, never, onerror}
              job_mult_memory: 1         # memory multiplier
              job_mult_time: 1           # running time multiplier
              merge_mult_memory: 1       # memory multiplier for merging
              merge_mult_time: 1         # running time multiplier for merging
              ignore_chroms:             # patterns of chromosome names to ignore
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
            type: germline_variants
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def snakemake_obj(minimal_config):
    """Returns Snakemake object."""
    # Define helper variables
    rule_name = "variant_calling_gatk_hc_run"
    threads = 2
    bench_iteration = 2
    scriptdir = "/work"
    input_dict = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
    }
    output_base_name = (
        "work/{mapper}.gatk_hc.{index_library_name}/out/{mapper}.gatk_hc.{index_library_name}"
    )
    output_dict = {
        "vcf_tbi": output_base_name + ".vcf.gz.tbi",
        "vcf_tbi_md5": output_base_name + ".vcf.gz.tbi.md5",
        "vcf": output_base_name + ".vcf.gz",
        "vcf_md5": output_base_name + ".vcf.gz.md5",
    }
    log_base_name = "work/bwa.gatk_hc.P001-T1-DNA1-WGS1/log/bwa.gatk_hc.P001-T1-DNA1-WGS1"
    log_dict = {
        "conda_info": log_base_name + ".conda_info.txt",
        "conda_info_md5": log_base_name + ".conda_info.txt.md5",
        "conda_list": log_base_name + ".conda_list.txt",
        "conda_list_md5": log_base_name + ".conda_list.txt.md5",
        "log": log_base_name + ".log",
        "log_md5": log_base_name + ".log.md5",
    }
    wildcards_dict = {"mapper": "bwa", "index_library_name": "P001-N1-DNA1-WGS1"}
    params_dict = {"step_key": "variant_calling", "caller_key": "gatk_hc"}

    # Define Snakemake class input
    input_ = InputFiles(fromdict=input_dict)
    output_ = OutputFiles(fromdict=output_dict)
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
        scriptdir=scriptdir,
        resources=resources_,
    )


def test_gatk_hc_wrapper_run_construct_merge_rule(snakemake_obj, variant_caller_fake_fs, mocker):
    """Tests ParallelGatkHcWrapper.construct_merge_rule()"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", variant_caller_fake_fs, mocker)
    wrapper_par = ParallelGatkHcWrapper(snakemake=snakemake_obj)
    # Define expected
    data_path = (Path(__file__).parent / "data/gatk_hc_par.snakemake").resolve()
    with open(data_path, "r", encoding="utf8") as f:
        expected = f.read()
    # Get actual and assert
    actual = wrapper_par.construct_merge_rule()
    assert actual == expected
