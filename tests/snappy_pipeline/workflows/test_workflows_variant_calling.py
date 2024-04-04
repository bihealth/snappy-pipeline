# -*- coding: utf-8 -*-
"""Tests for the variant_calling workflow module code"""

import textwrap
from copy import deepcopy

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            bwa:
              path_index: /path/to/bwa/index.fa
            target_cov_report:
              path_target_interval_list_mapping:
              - name: "Agilent SureSelect Human All Exon V6"
                pattern: "Agilent SureSelect Human All Exon V6*"
                path: "path/to/targets.bed"

          variant_calling:
            baf_file_generation:
              enabled: true
            tools:
            - bcftools_call
            - gatk3_hc
            - gatk3_ug

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
def variant_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return VariantCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_calling", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return VariantCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for BcftoolsCallStepPart ------------------------------------------------------------------


def test_bcftools_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("bcftools_call", "run")(wildcards)
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        ],
    }
    assert actual == expected


def test_bcftools_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.bcftools_call.{library_name}/out/{mapper}.bcftools_call.{library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    expected["output_links"] = [path.replace("work/", "output/") for path in expected.values()]
    for ext in ("log", "conda_info.txt", "conda_list.txt", "wrapper.py", "environment.yaml"):
        for full_ext in (ext, f"{ext}.md5"):
            base = "output/{mapper}.bcftools_call.{library_name}/log/{mapper}.bcftools_call.{library_name}.bcftools_call_run"
            expected["output_links"].append(f"{base}.{full_ext}")
    # Get actual
    actual = variant_calling_workflow.get_output_files("bcftools_call", "run")
    assert actual == expected


def test_bcftools_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.bcftools_call.{library_name}/log/{mapper}.bcftools_call.{library_name}.bcftools_call_run"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=True)
    # Get actual
    actual = variant_calling_workflow.get_log_file("bcftools_call", "run")
    assert actual == expected


def test_bcftools_call_step_part_get_resource(variant_calling_workflow):
    """Tests BcftoolsCallStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "61440M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_calling_workflow.get_resource("bcftools_call", "run", resource)
        assert actual == expected, msg_error


# Tests for GatkHaplotypeCallerStepPart -----------------------------------------------------------


def test_gatk3_hc_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("gatk3_hc", "run")(wildcards)
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        ],
    }
    assert actual == expected


def test_gatk3_hc_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.gatk3_hc.{library_name}/out/{mapper}.gatk3_hc.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    expected["output_links"] = [
        "output/{mapper}.gatk3_hc.{library_name}/out/{mapper}.gatk3_hc.{library_name}.vcf.gz",
        "output/{mapper}.gatk3_hc.{library_name}/out/{mapper}.gatk3_hc.{library_name}.vcf.gz.md5",
        "output/{mapper}.gatk3_hc.{library_name}/out/{mapper}.gatk3_hc.{library_name}.vcf.gz.tbi",
        "output/{mapper}.gatk3_hc.{library_name}/out/{mapper}.gatk3_hc.{library_name}.vcf.gz.tbi.md5",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.log",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.log.md5",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.conda_info.txt",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.conda_info.txt.md5",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.conda_list.txt",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.conda_list.txt.md5",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.wrapper.py",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.wrapper.py.md5",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.environment.yaml",
        "output/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run.environment.yaml.md5",
    ]
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk3_hc", "run")
    assert actual == expected


def test_gatk3_hc_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk3_hc.{library_name}/log/{mapper}.gatk3_hc.{library_name}.gatk3_hc_run"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=True)
    # Get actual
    actual = variant_calling_workflow.get_log_file("gatk3_hc", "run")
    assert actual == expected


def test_gatk3_hc_step_part_get_resource(variant_calling_workflow):
    """Tests GatkHaplotypeCallerStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "88G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_calling_workflow.get_resource("gatk3_hc", "run", resource)
        assert actual == expected, msg_error


# Tests for GatkUnifiedGenotyperStepPart ----------------------------------------------------------


def test_gatk3_ug_step_part_get_input_files(variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = variant_calling_workflow.get_input_files("gatk3_ug", "run")(wildcards)
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "bam": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        ],
    }
    assert actual == expected


def test_gatk3_ug_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.gatk3_ug.{library_name}/out/{mapper}.gatk3_ug.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    expected["output_links"] = [
        "output/{mapper}.gatk3_ug.{library_name}/out/{mapper}.gatk3_ug.{library_name}.vcf.gz",
        "output/{mapper}.gatk3_ug.{library_name}/out/{mapper}.gatk3_ug.{library_name}.vcf.gz.md5",
        "output/{mapper}.gatk3_ug.{library_name}/out/{mapper}.gatk3_ug.{library_name}.vcf.gz.tbi",
        "output/{mapper}.gatk3_ug.{library_name}/out/{mapper}.gatk3_ug.{library_name}.vcf.gz.tbi.md5",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.log",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.log.md5",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.conda_info.txt",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.conda_info.txt.md5",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.conda_list.txt",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.conda_list.txt.md5",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.wrapper.py",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.wrapper.py.md5",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.environment.yaml",
        "output/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run.environment.yaml.md5",
    ]
    # Get actual
    actual = variant_calling_workflow.get_output_files("gatk3_ug", "run")
    assert actual == expected


def test_gatk3_ug_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.gatk3_ug.{library_name}/log/{mapper}.gatk3_ug.{library_name}.gatk3_ug_run"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=True)
    # Get actual
    actual = variant_calling_workflow.get_log_file("gatk3_ug", "run")
    assert actual == expected


def test_gatk3_ug_step_part_get_resource(variant_calling_workflow):
    """Tests GatkHaplotypeCallerStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "88G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_calling_workflow.get_resource("gatk3_ug", "run", resource)
        assert actual == expected, msg_error


# Tests for BcftoolsStatsStepPart ------------------------------------------------------------------


def test_bcftools_stats_step_part_get_input_files(variant_calling_workflow):
    # Define expected
    vcf_file = (
        "work/{mapper}.{var_caller}.{index_library_name}/out/"
        "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
    )
    expected = {"vcf": vcf_file}
    # Get actual
    actual = variant_calling_workflow.get_input_files("bcftools_stats", "run")
    assert actual == expected


def test_bcftools_stats_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.{index_library_name}/report/bcftools_stats/"
        "{mapper}.{var_caller}.{index_library_name}.{donor_library_name}"
    )
    expected = {
        "txt": base_name_out + ".txt",
        "txt_md5": base_name_out + ".txt.md5",
    }
    expected["output_links"] = [value.replace("work/", "output/") for value in expected.values()]
    for ext in ("log", "conda_info.txt", "conda_list.txt", "wrapper.py", "environment.yaml"):
        token = "{mapper}.{var_caller}.{index_library_name}"
        line = f"output/{token}/log/{token}.{{donor_library_name}}.bcftools_stats_run.{ext}"
        expected["output_links"].append(line)
        expected["output_links"].append(f"{line}.md5")
    # Get actual
    actual = variant_calling_workflow.get_output_files("bcftools_stats", "run")
    assert actual == expected


def test_bcftools_stats_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    expected = {
        "conda_info": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.conda_list.txt.md5",
        "env_yaml": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.environment.yaml",
        "env_yaml_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.environment.yaml.md5",
        "log": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.log",
        "log_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.log.md5",
        "wrapper": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.wrapper.py",
        "wrapper_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.bcftools_stats_run.wrapper.py.md5",
    }
    # Get actual
    actual = variant_calling_workflow.get_log_file("bcftools_stats", "run")
    assert actual == expected


def test_bcftools_stats_step_part_get_resource(variant_calling_workflow):
    """Tests BcftoolsStatsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "1024M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_calling_workflow.get_resource("bcftools_stats", "run", resource)
        assert actual == expected, msg_error


# Tests for JannovarStatisticsStepPart -------------------------------------------------------------


def test_jannovar_stats_step_part_get_input_files(variant_calling_workflow):
    # Define expected
    vcf_file = (
        "work/{mapper}.{var_caller}.{index_library_name}/out/"
        "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
    )
    expected = {"vcf": vcf_file}
    # Get actual
    actual = variant_calling_workflow.get_input_files("jannovar_stats", "run")
    assert actual == expected


def test_jannovar_stats_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.{index_library_name}/report/jannovar_stats/"
        "{mapper}.{var_caller}.{index_library_name}"
    )
    expected = {
        "report": base_name_out + ".txt",
        "report_md5": base_name_out + ".txt.md5",
    }
    expected["output_links"] = [value.replace("work/", "output/") for value in expected.values()]
    for ext in ("log", "conda_info.txt", "conda_list.txt", "wrapper.py", "environment.yaml"):
        token = "{mapper}.{var_caller}.{index_library_name}"
        line = f"output/{token}/log/{token}.jannovar_stats_run.{ext}"
        expected["output_links"].append(line)
        expected["output_links"].append(f"{line}.md5")
    # Get actual
    actual = variant_calling_workflow.get_output_files("jannovar_stats", "run")
    assert actual == expected


def test_jannovar_stats_step_part_get_log_file(variant_calling_workflow):
    expected = {
        "conda_info": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.conda_list.txt.md5",
        "env_yaml": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.environment.yaml",
        "env_yaml_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.environment.yaml.md5",
        "log": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.log",
        "log_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.log.md5",
        "wrapper": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.wrapper.py",
        "wrapper_md5": "work/{mapper}.{var_caller}.{index_library_name}/log/{mapper}.{var_caller}.{index_library_name}.jannovar_stats_run.wrapper.py.md5",
    }
    actual = variant_calling_workflow.get_log_file("jannovar_stats", "run")
    assert actual == expected


def test_jannovar_stats_stats_step_part_get_resource(variant_calling_workflow):
    """Tests JannovarStatisticsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "04:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_calling_workflow.get_resource("jannovar_stats", "run", resource)
        assert actual == expected, msg_error


# Tests for BafFileGeneration ---------------------------------------------------------------------


def test_baf_file_generation_step_part_get_input_files(variant_calling_workflow):
    # Define expected
    vcf_file = (
        "work/{mapper}.{var_caller}.{index_library_name}/out/"
        "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
    )
    expected = {"vcf": vcf_file}
    # Get actual
    actual = variant_calling_workflow.get_input_files("baf_file_generation", "run")
    assert actual == expected


def test_baf_file_generation_step_part_get_output_files(variant_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.{var_caller}.{index_library_name}/report/baf/"
        "{mapper}.{var_caller}.{index_library_name}.{donor_library_name,[^\\.]+}.baf"
    )
    expected = {
        "bw": base_name_out + ".bw",
        "bw_md5": base_name_out + ".bw.md5",
    }
    expected["output_links"] = [value.replace("work/", "output/") for value in expected.values()]
    for ext in ("log", "conda_info.txt", "conda_list.txt", "wrapper.py", "environment.yaml"):
        token = "{mapper}.{var_caller}.{index_library_name}"
        line = f"output/{token}/log/{token}.{{donor_library_name}}.baf_file_generation_run.{ext}"
        expected["output_links"].append(line)
        expected["output_links"].append(f"{line}.md5")
    # Get actual
    actual = variant_calling_workflow.get_output_files("baf_file_generation", "run")
    assert actual == expected


def test_baf_file_generation_step_part_get_log_file(variant_calling_workflow):
    # Define expected
    tpl = (
        "work/{mapper}.{var_caller}.{index_library_name}/log/"
        "{mapper}.{var_caller}.{index_library_name}.{donor_library_name}.baf_file_generation_run"
    )
    key_ext = {
        "log": ".log",
        "log_md5": ".log.md5",
        "conda_info": ".conda_info.txt",
        "conda_info_md5": ".conda_info.txt.md5",
        "conda_list": ".conda_list.txt",
        "conda_list_md5": ".conda_list.txt.md5",
        "env_yaml": ".environment.yaml",
        "env_yaml_md5": ".environment.yaml.md5",
        "wrapper": ".wrapper.py",
        "wrapper_md5": ".wrapper.py.md5",
    }
    expected = {key: tpl + ext for key, ext in key_ext.items()}
    # Get actual
    actual = variant_calling_workflow.get_log_file("baf_file_generation", "run")
    assert actual == expected


def test_baf_file_generation_step_part_get_resource(variant_calling_workflow):
    """Tests BafFileGeneration.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "1024M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_calling_workflow.get_resource("baf_file_generation", "run", resource)
        assert actual == expected, msg_error


# Tests for VariantCallingWorkflow ----------------------------------------------------------------


def test_variant_calling_workflow(variant_calling_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = [
        "baf_file_generation",
        "bcftools_call",
        "bcftools_roh",
        "bcftools_stats",
        "gatk3_hc",
        "gatk3_ug",
        "gatk4_hc_gvcf",
        "gatk4_hc_joint",
        "jannovar_stats",
        "write_pedigree",
    ]
    assert sorted(variant_calling_workflow.sub_steps.keys()) == sorted(expected)
    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
    ]
    base_out = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{step}.{ext}"
    )
    expected += [
        base_out.format(i=i, ext=ext, mapper=mapper, var_caller=var_caller, step=step)
        for i in (1, 4)  # only for indices
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
            "environment.yaml",
            "environment.yaml.md5",
            "wrapper.py",
            "wrapper.py.md5",
        )
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
        for step in (
            f"{var_caller}_run",
            "jannovar_stats_run",
            "bcftools_roh_run",
        )
    ]
    base_out = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.P00{k}-N1-DNA1-WGS1.{step}.{ext}"
    )
    expected += [
        base_out.format(i=i, k=i + j, ext=ext, mapper=mapper, var_caller=var_caller, step=step)
        for i in (1, 4)  # only for indices
        for j in (0, 1, 2)  # all donors
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
            "environment.yaml",
            "environment.yaml.md5",
            "wrapper.py",
            "wrapper.py.md5",
        )
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
        for step in ("bcftools_stats_run", "baf_file_generation_run")
    ]
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/report/"
        "bcftools_stats/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.P00{t}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (1, 2), (1, 3), (4, 4), (4, 5), (4, 6))
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
        for ext in ("txt", "txt.md5")
    ]
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/report/"
        "jannovar_stats/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
        for ext in ("txt", "txt.md5")
    ]
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/report/"
        "baf/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.P00{t}-N1-DNA1-WGS1.baf.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (1, 2), (1, 3), (4, 4), (4, 5), (4, 6))
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
        for ext in ("bw", "bw.md5")
    ]
    tpl = (
        "output/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1/report/"
        "roh/{mapper}.{var_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in (1, 4)
        for mapper in ("bwa",)
        for var_caller in (
            "bcftools_call",
            "gatk3_hc",
            "gatk3_ug",
        )
        for ext in ("txt", "txt.md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(variant_calling_workflow.get_result_files()))
    assert actual == expected


def test_variant_calling_custom_pedigree_field(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_trio_plus_sheet_fake_fs,
    mocker,
):
    """Tests VariantCallingWorkflow object pre-configured with germline trio plus sheet
    and custom pedigree field"""
    # Initialise variables
    index_standard_pedigree_list = ["P001-N1-DNA1-WGS1", "P004-N1-DNA1-WGS1", "P007-N1-DNA1-WGS1"]
    index_custom_field_pedigree_list = ["P001-N1-DNA1-WGS1", "P004-N1-DNA1-WGS1"]

    # Create alternative configuration file
    local_minimal_config = deepcopy(minimal_config)
    local_minimal_config["data_sets"]["first_batch"]["file"] = "sheet_trio_plus.tsv"

    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_trio_plus_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_trio_plus_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_calling", germline_trio_plus_sheet_fake_fs, mocker
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}

    # Construct the workflow object - should work, standard pedigree join
    vcw = VariantCallingWorkflow(
        dummy_workflow,
        local_minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )
    # Expects a single shortcut sheet
    assert len(vcw.shortcut_sheets) == 1
    sheet = vcw.shortcut_sheets[0]

    # P007 will be considered as a different index as the not defining pedigree based on `familyId`
    assert len(sheet.index_ngs_library_to_pedigree) == 3
    # Index name as expected
    assert all(
        [
            index in index_standard_pedigree_list
            for index in sheet.index_ngs_library_to_pedigree.keys()
        ]
    )

    # Construct the workflow object - should work, custom pedigree join
    local_minimal_config["data_sets"]["first_batch"]["pedigree_field"] = "familyId"
    vcw = VariantCallingWorkflow(
        dummy_workflow,
        local_minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )
    # Expects a single shortcut sheet
    assert len(vcw.shortcut_sheets) == 1
    sheet = vcw.shortcut_sheets[0]

    # P007 will be included in the same pedigree as P00{4,5,6}, i.e., all labelled 'family2'
    assert len(sheet.index_ngs_library_to_pedigree) == 2
    # Index name as expected
    assert all(
        [
            index in index_custom_field_pedigree_list
            for index in sheet.index_ngs_library_to_pedigree.keys()
        ]
    )

    # Construct the workflow object - should fail, custom pedigree field not defined
    with pytest.raises(Exception):
        local_minimal_config["data_sets"]["first_batch"]["pedigree_field"] = "_field_not_defined_"
        VariantCallingWorkflow(
            dummy_workflow,
            local_minimal_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )
