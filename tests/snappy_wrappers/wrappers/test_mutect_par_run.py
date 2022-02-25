# -*- coding: utf-8 -*-
"""Code for testing mutect_par/run wrapper"""
from pathlib import Path
import textwrap

import pytest
from ruamel import yaml as yaml_ruamel
from snakemake.io import InputFiles, Log, OutputFiles, Params, Resources, Wildcards
from snakemake.script import Snakemake

from snappy_wrappers.wrappers.mutect_par.parallel_mutect_wrapper import ParallelMutectWrapper

from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config_mutect():
    """Return YAML parsing result for configuration"""
    return yaml_ruamel.round_trip_load(
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
          somatic_variant_calling:
            tools:
            - mutect
            mutect:
              # Parallelization configuration
              drmaa_snippet: ''          # value to pass in as additional DRMAA arguments
              num_cores: 2               # number of cores to use locally
              window_length: 3500000     # split input into windows of this size, each triggers a job
              num_jobs: 500              # number of windows to process in parallel
              use_drmaa: true            # use drmaa for parallel processing
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
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def snakemake_obj(minimal_config_mutect):
    """Returns Snakemake object."""
    # Define helper variables
    rule_name = "somatic_variant_calling_mutect_run"
    threads = 2
    bench_iteration = 2
    scriptdir = "/work"
    input_dict = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    output_base_name = "work/{mapper}.mutect.{tumor_library}/out/{mapper}.mutect.{tumor_library}"
    output_dict = {
        "tbi": output_base_name + ".vcf.gz.tbi",
        "tbi_md5": output_base_name + ".vcf.gz.tbi.md5",
        "vcf": output_base_name + ".vcf.gz",
        "vcf_md5": output_base_name + ".vcf.gz.md5",
        "full_tbi": output_base_name + ".full.vcf.gz.tbi",
        "full_tbi_md5": output_base_name + ".full.vcf.gz.tbi.md5",
        "full": output_base_name + ".full.vcf.gz",
        "full_md5": output_base_name + ".full.vcf.gz.md5",
        "txt": output_base_name + ".txt",
        "txt_md5": output_base_name + ".txt.md5",
        "wig": output_base_name + ".wig",
        "wig_md5": output_base_name + ".wig.md5",
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
    params_dict = {"normal_lib_name": "P001-N1-DNA1-WGS1"}

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
        config=minimal_config_mutect,
        scriptdir=scriptdir,
        resources=resources_,
    )


def test_mutect_wrapper_run_construct_merge_rule(snakemake_obj, somatic_variant_fake_fs, mocker):
    """Tests ParallelMutectWrapper.construct_merge_rule()"""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", somatic_variant_fake_fs, mocker)
    wrapper_par = ParallelMutectWrapper(snakemake=snakemake_obj)
    # Define expected
    data_path = (Path(__file__).parent / "data/mutect_par.snakemake").resolve()
    with open(data_path, "r") as f:
        expected = f.read()
    # Get actual and assert
    actual = wrapper_par.construct_merge_rule()
    print(actual)
    assert actual == expected
