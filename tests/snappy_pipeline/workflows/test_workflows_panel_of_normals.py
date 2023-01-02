# -*- coding: utf-8 -*-
"""Tests for the panel_of_normals workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.panel_of_normals import PanelOfNormalsWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (cancer) configuration"""
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
              tools: ['mutect2', 'cnvkit']
              mutect2:
                  germline_resource: /path/to/germline_resource.vcf
                  path_normals_list: ""
              cnvkit:
                  path_excluded_regions: ""
                  path_target_regions: /path/to/target_regions.bed
                  path_normals_list: ""

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
def panel_of_normals_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return PanelOfNormalsWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return PanelOfNormalsWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for Mutect2StepPart ------------------------------------------------------------------------


def test_mutect2_step_part_get_input_files_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_input_files_prepare_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "normal_library": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "prepare_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_input_files_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "normals": [
            "work/bwa.mutect2.prepare_panel/out/P001-N1-DNA1-WGS1.vcf.gz",
            "work/bwa.mutect2.prepare_panel/out/P002-N1-DNA1-WGS1.vcf.gz",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "create_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_output_files_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_prepare_panel()"""
    expected = {
        "vcf": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.md5",
        "vcf_tbi": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.tbi",
        "vcf_tbi_md5": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.tbi.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2.create_panel/out/{mapper}.mutect2.panel_of_normals"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_output_files("mutect2", "create_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_prepare_panel()"""
    base_name_out = "work/{mapper}.mutect2.prepare_panel/log/{normal_library}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2.create_panel/log/{mapper}.mutect2.panel_of_normals"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "create_panel")
    assert actual == expected


def test_mutect2_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests Mutect2StepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    create_panel_expected_dict = {
        "threads": 2,
        "time": "08:00:00",
        "memory": "30G",
        "partition": "medium",
    }
    prepare_panel_expected_dict = {
        "threads": 2,
        "time": "3-00:00:00",
        "memory": "8G",
        "partition": "medium",
    }

    # Evaluate action `create_panel`
    for resource, expected in create_panel_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'create_panel'."
        actual = panel_of_normals_workflow.get_resource("mutect2", "create_panel", resource)
        assert actual == expected, msg_error

    # Evaluate action `prepare_panel`
    for resource, expected in prepare_panel_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'prepare_panel'."
        actual = panel_of_normals_workflow.get_resource("mutect2", "prepare_panel", resource)
        assert actual == expected, msg_error


# Tests for CnvkitStepPart ------------------------------------------------------------------------


# def test_cnvkit_step_part_get_input_files_autobin(panel_of_normals_workflow):
#     """Tests CnvkitStepPart._get_input_files_autobin()"""
#     wildcards = Wildcards(
#         fromdict={
#             "mapper": "bwa",
#         }
#     )
#     expected = {
#         "access": "work/cnvkit.access/out/cnvkit.access.bed",
#         "bams": [
#             "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
#             "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
#         ],
#         "bais": [
#             "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
#             "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
#         ],
#     }
#     actual = panel_of_normals_workflow.get_input_files("cnvkit", "autobin")(wildcards)
#     assert actual == expected


def test_cnvkit_step_part_get_input_files_antitarget(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_input_files_antitarget()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "normal_library": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "access": "work/cnvkit.access/out/cnvkit.access.bed",
        "target": "work/cnvkit.target/out/cnvkit.target.bed",
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "antitarget")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_coverage(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "normal_library": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "target": "work/cnvkit.target/out/cnvkit.target.bed",
        "antitarget": "work/cnvkit.antitarget/out/cnvkit.antitarget.bed",
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "coverage")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_create_panel(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "target": [
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P001-N1-DNA1-WGS1.target.cnn",
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P002-N1-DNA1-WGS1.target.cnn",
        ],
        "antitarget": [
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P001-N1-DNA1-WGS1.antitarget.cnn",
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P002-N1-DNA1-WGS1.antitarget.cnn",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "create_panel")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_sex(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_input_files_sex()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "target": [
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P001-N1-DNA1-WGS1.target.cnn",
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P002-N1-DNA1-WGS1.target.cnn",
        ],
        "antitarget": [
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P001-N1-DNA1-WGS1.antitarget.cnn",
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P002-N1-DNA1-WGS1.antitarget.cnn",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "sex")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_metrics(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_input_files_metrics()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "target": [
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P001-N1-DNA1-WGS1.target.cnn",
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P002-N1-DNA1-WGS1.target.cnn",
        ],
        "antitarget": [
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P001-N1-DNA1-WGS1.antitarget.cnn",
            "work/bwa.cnvkit.coverage/out/bwa.cnvkit.coverage.P002-N1-DNA1-WGS1.antitarget.cnn",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "metrics")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_output_files_access(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_access()"""
    expected = {
        "access": "work/cnvkit.access/out/cnvkit.access.bed",
        "access_md5": "work/cnvkit.access/out/cnvkit.access.bed.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "access")
    assert actual == expected


# def test_cnvkit_step_part_get_output_files_autobin(panel_of_normals_workflow):
#     """Tests CvnkitStepPart._get_output_files_autobin()"""
#     expected = {
#         "target": "work/{mapper}.cnvkit.autobin/out/{mapper}.cnvkit.autobin.target.bed",
#         "target_md5": "work/{mapper}.cnvkit.autobin/out/{mapper}.cnvkit.autobin.target.bed.md5",
#         "antitarget": "work/{mapper}.cnvkit.autobin/out/{mapper}.cnvkit.autobin.antitarget.bed",
#         "antitarget_md5": "work/{mapper}.cnvkit.autobin/out/{mapper}.cnvkit.autobin.antitarget.bed.md5",
#     }
#     actual = panel_of_normals_workflow.get_output_files("cnvkit", "autobin")
#     assert actual == expected


def test_cnvkit_step_part_get_output_files_target(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_target()"""
    expected = {
        "target": "work/cnvkit.target/out/cnvkit.target.bed",
        "target_md5": "work/cnvkit.target/out/cnvkit.target.bed.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "target")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_antitarget(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_antitarget()"""
    expected = {
        "antitarget": "work/cnvkit.antitarget/out/cnvkit.antitarget.bed",
        "antitarget_md5": "work/cnvkit.antitarget/out/cnvkit.antitarget.bed.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "antitarget")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_coverage(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_coverage()"""
    expected = {
        "target": "work/{mapper}.cnvkit.coverage/out/{mapper}.cnvkit.coverage.{normal_library}.target.cnn",
        "target_md5": "work/{mapper}.cnvkit.coverage/out/{mapper}.cnvkit.coverage.{normal_library}.target.cnn.md5",
        "antitarget": "work/{mapper}.cnvkit.coverage/out/{mapper}.cnvkit.coverage.{normal_library}.antitarget.cnn",
        "antitarget_md5": "work/{mapper}.cnvkit.coverage/out/{mapper}.cnvkit.coverage.{normal_library}.antitarget.cnn.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_create_panel()"""
    expected = {
        "panel": "work/{mapper}.cnvkit.create_panel/out/{mapper}.cnvkit.panel_of_normals.cnn",
        "panel_md5": "work/{mapper}.cnvkit.create_panel/out/{mapper}.cnvkit.panel_of_normals.cnn.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "create_panel")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_sex(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_sex()"""
    expected = {
        "table": "work/{mapper}.cnvkit.create_panel/report/{mapper}.cnvkit.sex.tsv",
        "table_md5": "work/{mapper}.cnvkit.create_panel/report/{mapper}.cnvkit.sex.tsv.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "sex")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_metrics(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_metrics()"""
    expected = {
        "table": "work/{mapper}.cnvkit.create_panel/report/{mapper}.cnvkit.metrics.tsv",
        "table_md5": "work/{mapper}.cnvkit.create_panel/report/{mapper}.cnvkit.metrics.tsv.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "metrics")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_access(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_access()"""
    base_name_out = "work/cnvkit.access/log/cnvkit.access"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "access")
    assert actual == expected


# def test_cnvkit_step_part_get_log_file_autobin(panel_of_normals_workflow):
#     """Tests CvnkitStepPart._get_log_files_autobin()"""
#     base_name_out = "work/{mapper}.cnvkit.autobin/log/{mapper}.cnvkit.autobin"
#     expected = get_expected_log_files_dict(base_out=base_name_out)
#     actual = panel_of_normals_workflow.get_log_file("cnvkit", "autobin")
#     assert actual == expected


def test_cnvkit_step_part_get_log_file_target(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_target()"""
    base_name_out = "work/cnvkit.target/log/cnvkit.target"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "target")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_antitarget(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_antitarget()"""
    base_name_out = "work/cnvkit.antitarget/log/cnvkit.antitarget"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "antitarget")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_coverage(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_coverage()"""
    base_name_out = "work/{mapper}.cnvkit.coverage/log/{mapper}.cnvkit.coverage.{normal_library}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_create_panel()"""
    base_name_out = "work/{mapper}.cnvkit.create_panel/log/{mapper}.cnvkit.panel_of_normals"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "create_panel")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_sex(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_sex()"""
    base_name_out = "work/{mapper}.cnvkit.create_panel/log/{mapper}.cnvkit.sex"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "sex")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_metrics(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_metrics()"""
    base_name_out = "work/{mapper}.cnvkit.create_panel/log/{mapper}.cnvkit.metrics"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "metrics")
    assert actual == expected


def test_cnvkit_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests CvnkitStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    access_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "8G",
        "partition": "medium",
    }
    # autobin_expected_dict = {
    #     "threads": 2,
    #     "time": "02:00:00",
    #     "memory": "16G",
    #     "partition": "medium",
    # }
    target_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "8G",
        "partition": "medium",
    }
    antitarget_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "8G",
        "partition": "medium",
    }
    coverage_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "16G",
        "partition": "medium",
    }
    reference_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "16G",
        "partition": "medium",
    }
    sex_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "16G",
        "partition": "medium",
    }
    metrics_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "16G",
        "partition": "medium",
    }

    # Evaluate action `access`
    for resource, expected in access_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'access'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "access", resource)
        assert actual == expected, msg_error

    # Evaluate action `autobin`
    # for resource, expected in autobin_expected_dict.items():
    #     msg_error = f"Assertion error for resource '{resource}' for action 'autobin'."
    #     actual = panel_of_normals_workflow.get_resource("cnvkit", "autobin", resource)
    #     assert actual == expected, msg_error

    # Evaluate action `target`
    for resource, expected in target_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'target'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "target", resource)
        assert actual == expected, msg_error

    # Evaluate action `antitarget`
    for resource, expected in antitarget_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'antitarget'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "antitarget", resource)
        assert actual == expected, msg_error

    # Evaluate action `coverage`
    for resource, expected in coverage_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'coverage'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "coverage", resource)
        assert actual == expected, msg_error

    # Evaluate action `create_panel`
    for resource, expected in reference_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'create_panel'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "create_panel", resource)
        assert actual == expected, msg_error

    # Evaluate action `sex`
    for resource, expected in reference_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'sex'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "sex", resource)
        assert actual == expected, msg_error

    # Evaluate action `metrics`
    for resource, expected in reference_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'metrics'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "metrics", resource)
        assert actual == expected, msg_error


# PanelOfNormalsWorkflow  --------------------------------------------------------------------------


def test_panel_of_normals_workflow(panel_of_normals_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnvkit", "link_out", "mutect2"]
    actual = list(sorted(panel_of_normals_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = []
    tpl = "output/{mapper}.mutect2.create_panel/out/{mapper}.mutect2.panel_of_normals.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
    ]
    # add log files
    tpl = "output/{mapper}.mutect2.create_panel/log/{mapper}.mutect2.panel_of_normals"
    for mapper in ("bwa",):
        expected += get_expected_log_files_dict(tpl.format(mapper=mapper)).values()
    # Now for basic cnvkit files (panel of normal only)
    tpl = "output/{mapper}.cnvkit.create_panel/out/{mapper}.cnvkit.panel_of_normals.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext) for ext in ("cnn", "cnn.md5") for mapper in ("bwa",)
    ]
    tpl = "output/cnvkit.access/out/cnvkit.access.{ext}"
    expected += [tpl.format(ext=ext) for ext in ("bed", "bed.md5")]
    # tpl = "output/{mapper}.cnvkit.autobin/out/{mapper}.cnvkit.autobin.{ext}"
    # expected += [
    #     tpl.format(mapper=mapper, ext=ext)
    #     for ext in ("target.bed", "target.bed.md5", "antitarget.bed", "antitarget.bed.md5")
    #     for mapper in ("bwa",)
    # ]
    tpl = "output/cnvkit.target/out/cnvkit.target.{ext}"
    expected += [tpl.format(ext=ext) for ext in ("bed", "bed.md5")]
    tpl = "output/cnvkit.antitarget/out/cnvkit.antitarget.{ext}"
    expected += [tpl.format(ext=ext) for ext in ("bed", "bed.md5")]
    tpl = "output/{mapper}.cnvkit.create_panel/report/{mapper}.cnvkit.sex.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext) for ext in ("tsv", "tsv.md5") for mapper in ("bwa",)
    ]
    tpl = "output/{mapper}.cnvkit.create_panel/report/{mapper}.cnvkit.metrics.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext) for ext in ("tsv", "tsv.md5") for mapper in ("bwa",)
    ]
    # add log files
    tpl = "output/{mapper}.cnvkit.create_panel/log/{mapper}.cnvkit.panel_of_normals"
    for mapper in ("bwa",):
        expected += get_expected_log_files_dict(tpl.format(mapper=mapper)).values()
    expected += get_expected_log_files_dict("output/cnvkit.access/log/cnvkit.access").values()
    # tpl = "output/{mapper}.cnvkit.autobin/log/{mapper}.cnvkit.autobin"
    # for mapper in ("bwa",):
    #     expected += get_expected_log_files_dict(tpl.format(mapper=mapper)).values()
    expected += get_expected_log_files_dict("output/cnvkit.target/log/cnvkit.target").values()
    expected += get_expected_log_files_dict(
        "output/cnvkit.antitarget/log/cnvkit.antitarget"
    ).values()
    tpl = "output/{mapper}.cnvkit.create_panel/log/{mapper}.cnvkit.sex"
    for mapper in ("bwa",):
        expected += get_expected_log_files_dict(tpl.format(mapper=mapper)).values()
    tpl = "output/{mapper}.cnvkit.create_panel/log/{mapper}.cnvkit.metrics"
    for mapper in ("bwa",):
        expected += get_expected_log_files_dict(tpl.format(mapper=mapper)).values()
    expected = list(sorted(expected))
    actual = list(sorted(panel_of_normals_workflow.get_result_files()))
    assert actual == expected
