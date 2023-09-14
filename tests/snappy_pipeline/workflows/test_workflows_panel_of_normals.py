# -*- coding: utf-8 -*-
"""Tests for the panel_of_normals workflow module code"""

from collections import OrderedDict
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
              tools: ['mutect2', 'cnvkit', 'access', 'purecn']
              path_ngs_mapping: ../ngs_mapping
              mutect2:
                  germline_resource: /path/to/germline_resource.vcf
                  path_normals_list: ""
              cnvkit:
                  path_excluded_regions: ""
                  path_target_regions: /path/to/regions.bed  # WES mode
                  path_normals_list: ""
              purecn:
                  path_normals_list: ""
                  path_bait_regions: /path/to/baits/regions.bed

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
            "work/bwa.mutect2/out/bwa.mutect2.P001-N1-DNA1-WGS1.prepare.vcf.gz",
            "work/bwa.mutect2/out/bwa.mutect2.P002-N1-DNA1-WGS1.prepare.vcf.gz",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "create_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_output_files_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_prepare_panel()"""
    expected = {
        "vcf": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz",
        "vcf_md5": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz.md5",
        "vcf_tbi": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz.tbi",
        "vcf_tbi_md5": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz.tbi.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_output_files("mutect2", "create_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_prepare_panel()"""
    base_name_out = "work/{mapper}.mutect2/log/{mapper}.mutect2.{normal_library}.prepare"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "create_panel")
    assert actual == expected


def test_mutect2_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests Mutect2StepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    create_panel_expected_dict = {
        "threads": 2,
        "time": "48:00:00",
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


def test_cnvkit_step_part_get_input_files_target(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_input_files_target()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "target")(wildcards)
    assert actual == {}


def test_cnvkit_step_part_get_input_files_antitarget(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_input_files_antitarget()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "target": "work/bwa.cnvkit/out/bwa.cnvkit.target.bed",
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
        "target": "work/bwa.cnvkit/out/bwa.cnvkit.target.bed",
        "antitarget": "work/bwa.cnvkit/out/bwa.cnvkit.antitarget.bed",
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
            "work/bwa.cnvkit/out/bwa.cnvkit.P001-N1-DNA1-WGS1.targetcoverage.cnn",
            "work/bwa.cnvkit/out/bwa.cnvkit.P002-N1-DNA1-WGS1.targetcoverage.cnn",
        ],
        "antitarget": [
            "work/bwa.cnvkit/out/bwa.cnvkit.P001-N1-DNA1-WGS1.antitargetcoverage.cnn",
            "work/bwa.cnvkit/out/bwa.cnvkit.P002-N1-DNA1-WGS1.antitargetcoverage.cnn",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "create_panel")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_report(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_input_files_report()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "target": [
            "work/bwa.cnvkit/out/bwa.cnvkit.P001-N1-DNA1-WGS1.targetcoverage.cnn",
            "work/bwa.cnvkit/out/bwa.cnvkit.P002-N1-DNA1-WGS1.targetcoverage.cnn",
        ],
        "antitarget": [
            "work/bwa.cnvkit/out/bwa.cnvkit.P001-N1-DNA1-WGS1.antitargetcoverage.cnn",
            "work/bwa.cnvkit/out/bwa.cnvkit.P002-N1-DNA1-WGS1.antitargetcoverage.cnn",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "report")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_output_files_target(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_target()"""
    expected = {
        "target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed",
        "target_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "target")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_antitarget(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_antitarget()"""
    expected = {
        "antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed",
        "antitarget_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "antitarget")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_coverage(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_coverage()"""
    expected = {
        "target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.targetcoverage.cnn",
        "target_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.targetcoverage.cnn.md5",
        "antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.antitargetcoverage.cnn",
        "antitarget_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.antitargetcoverage.cnn.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_create_panel()"""
    expected = {
        "panel": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.cnn",
        "panel_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.cnn.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "create_panel")
    assert actual == expected


def test_cnvkit_step_part_get_output_files_report(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_output_files_report()"""
    expected = {
        "sex": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.sex.tsv",
        "sex_md5": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.sex.tsv.md5",
        "metrics": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.metrics.tsv",
        "metrics_md5": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.metrics.tsv.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("cnvkit", "report")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_target(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_target()"""
    base_name_out = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.target"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "target")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_antitarget(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_antitarget()"""
    base_name_out = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.antitarget"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "antitarget")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_coverage(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_coverage()"""
    base_name_out = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.{normal_library}.coverage"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_create_panel()"""
    base_name_out = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.panel_of_normals"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "create_panel")
    assert actual == expected


def test_cnvkit_step_part_get_log_file_report(panel_of_normals_workflow):
    """Tests CvnkitStepPart._get_log_files_report()"""
    base_name_out = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.report"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "report")
    assert actual == expected


def test_cnvkit_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests CvnkitStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
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
        "threads": 8,
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
    report_expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "16G",
        "partition": "medium",
    }

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

    # Evaluate action `report`
    for resource, expected in report_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'report'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "report", resource)
        assert actual == expected, msg_error


# Tests for AccessStepPart -------------------------------------------------------------------------


def test_access_step_part_get_input_files_run(panel_of_normals_workflow):
    """Tests AccessStepPart._get_input_files_run()"""
    assert panel_of_normals_workflow.get_input_files("access", "run") == None


def test_access_step_part_get_output_files_run(panel_of_normals_workflow):
    """Tests AccessStepPart._get_output_files_run()"""
    expected = {
        "access": "work/cnvkit.access/out/cnvkit.access.bed",
        "access_md5": "work/cnvkit.access/out/cnvkit.access.bed.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("access", "run")
    assert actual == expected


def test_access_step_part_get_log_file_run(panel_of_normals_workflow):
    """Tests AccessStepPart._get_log_file_run()"""
    expected = get_expected_log_files_dict(base_out="work/cnvkit.access/log/cnvkit.access")
    actual = panel_of_normals_workflow.get_log_file("access", "run")
    assert actual == expected


def test_access_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests AccessStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {
        "threads": 2,
        "time": "02:00:00",
        "memory": "8G",
        "partition": "medium",
    }
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'run'."
        actual = panel_of_normals_workflow.get_resource("access", "run", resource)
        assert actual == expected, msg_error


# Tests for PureCnStepPart -------------------------------------------------------------------------


def test_purecn_step_part_get_params_install(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_params_install()"""
    expected = {
        "packages": [
            OrderedDict([("name", "TxDb.Hsapiens.UCSC.hg19.knownGene"), ("repo", "bioconductor")]),
            OrderedDict([("name", "lima1/PSCBS"), ("repo", "github")]),
        ]
    }
    actual = panel_of_normals_workflow.get_params("purecn", "install")
    assert actual == expected


def test_purecn_step_part_get_output_files_install(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_install()"""
    expected = {"done": "work/R_packages/out/PureCN.done"}
    actual = panel_of_normals_workflow.get_output_files("purecn", "install")
    assert actual == expected


def test_purecn_step_part_get_log_file_install(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_install()"""
    expected = get_expected_log_files_dict(base_out="work/R_packages/log/PureCN")
    actual = panel_of_normals_workflow.get_log_file("purecn", "install")
    assert actual == expected


def test_purecn_step_part_get_input_files_prepare(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_input_files_prepare()"""
    expected = {"packages": "work/R_packages/out/PureCN.done"}
    actual = panel_of_normals_workflow.get_input_files("purecn", "prepare")
    assert actual == expected


def test_purecn_step_part_get_output_files_prepare(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_prepare()"""
    expected = {
        "intervals": "work/PureCN/out/unknown_unknown.list",
        "optimized": "work/PureCN/out/unknown_unknown.bed.gz",
        "tbi": "work/PureCN/out/unknown_unknown.bed.gz.tbi",
        "intervals_md5": "work/PureCN/out/unknown_unknown.list.md5",
        "optimized_md5": "work/PureCN/out/unknown_unknown.bed.gz.md5",
        "tbi_md5": "work/PureCN/out/unknown_unknown.bed.gz.tbi.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("purecn", "prepare")
    assert actual == expected


def test_purecn_step_part_get_log_file_prepare(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_prepare()"""
    expected = get_expected_log_files_dict(base_out="work/PureCN/log/unknown_unknown")
    actual = panel_of_normals_workflow.get_log_file("purecn", "prepare")
    assert actual == expected


def test_purecn_step_part_get_input_files_coverage(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "normal_library": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "packages": "work/R_packages/out/PureCN.done",
        "intervals": "work/PureCN/out/unknown_unknown.list",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    actual = panel_of_normals_workflow.get_input_files("purecn", "coverage")(wildcards)
    assert actual == expected


def test_purecn_step_part_get_output_files_coverage(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_coverage()"""
    expected = {
        "coverage": "work/{mapper}.PureCN.{normal_library}/out/{mapper}.{normal_library}_coverage_loess.txt.gz"
    }
    actual = panel_of_normals_workflow.get_output_files("purecn", "coverage")
    assert actual == expected


def test_purecn_step_part_get_log_file_coverage(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_coverage()"""
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.PureCN.{normal_library}/log/{mapper}.{normal_library}"
    )
    actual = panel_of_normals_workflow.get_log_file("purecn", "coverage")
    assert actual == expected


def test_purecn_step_part_get_input_files_create_panel(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "normals": [
            "work/bwa.PureCN.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1_coverage_loess.txt.gz",
            "work/bwa.PureCN.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1_coverage_loess.txt.gz",
        ],
        "packages": "work/R_packages/out/PureCN.done",
    }
    actual = panel_of_normals_workflow.get_input_files("purecn", "create_panel")(wildcards)
    assert actual == expected


def test_purecn_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_create_panel()"""
    expected = {
        "rds": "work/{mapper}.PureCN/out/{mapper}.PureCN.panel_of_normals.rds",
        "rds_md5": "work/{mapper}.PureCN/out/{mapper}.PureCN.panel_of_normals.rds.md5",
        "bed": "work/{mapper}.PureCN/out/{mapper}.PureCN.panel_of_normals.bed",
        "png": "work/{mapper}.PureCN/out/{mapper}.PureCN.panel_of_normals.png",
    }
    actual = panel_of_normals_workflow.get_output_files("purecn", "create_panel")
    assert actual == expected


def test_purecn_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_create_panel()"""
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.PureCN/log/{mapper}.PureCN.panel_of_normals"
    )
    actual = panel_of_normals_workflow.get_log_file("purecn", "create_panel")
    assert actual == expected


def test_purecn_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests PureCnStepPart.get_resource_usage() for all actions"""
    expected = {
        "coverage": {"threads": 1, "memory": "24G", "time": "04:00:00"},
        "prepare": {"threads": 1, "memory": "24G", "time": "04:00:00"},
        "create_panel": {"threads": 1, "memory": "24G", "time": "04:00:00"},
    }
    for action, resources in expected.items():
        for resource, value in resources.items():
            actual = panel_of_normals_workflow.get_resource("purecn", action, resource)
            assert actual == value


# PanelOfNormalsWorkflow  --------------------------------------------------------------------------


def test_panel_of_normals_workflow(panel_of_normals_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["access", "cnvkit", "link_out", "mutect2", "purecn"]
    actual = list(sorted(panel_of_normals_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = []
    tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
    ]
    # add log files
    tpl = "output/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals"
    for mapper in ("bwa",):
        expected += get_expected_log_files_dict(base_out=tpl.format(mapper=mapper)).values()

    # Now for basic cnvkit files (panel of normal only)
    tpl = "output/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext) for ext in ("cnn", "cnn.md5") for mapper in ("bwa",)
    ]
    tpl = "output/{mapper}.cnvkit/out/{mapper}.cnvkit.{substep}.{ext}"
    for substep in ("target", "antitarget"):
        expected += [
            tpl.format(substep=substep, mapper=mapper, ext=ext)
            for ext in ("bed", "bed.md5")
            for mapper in ("bwa",)
        ]
    tpl = "output/{mapper}.cnvkit/report/{mapper}.cnvkit.{substep}.{ext}"
    for substep in ("sex", "metrics"):
        expected += [
            tpl.format(substep=substep, mapper=mapper, ext=ext)
            for ext in ("tsv", "tsv.md5")
            for mapper in ("bwa",)
        ]
    # add log files
    tpl = "output/{mapper}.cnvkit/log/{mapper}.cnvkit.{substep}"
    for substep in ("target", "antitarget", "panel_of_normals", "report"):
        for mapper in ("bwa",):
            expected += get_expected_log_files_dict(
                base_out=tpl.format(mapper=mapper, substep=substep)
            ).values()

    # Access
    tpl = "output/cnvkit.access/out/cnvkit.access.{ext}"
    expected += [tpl.format(ext=ext) for ext in ("bed", "bed.md5")]
    expected += get_expected_log_files_dict(
        base_out="output/cnvkit.access/log/cnvkit.access"
    ).values()

    # PureCN
    tpl = "output/{mapper}.PureCN/out/{mapper}.PureCN.panel_of_normals.rds{chksum}"
    expected += [
        tpl.format(mapper=mapper, chksum=chksum) for mapper in ("bwa",) for chksum in ("", ".md5")
    ]
    expected += get_expected_log_files_dict(
        base_out="output/{mapper}.PureCN/log/{mapper}.PureCN.panel_of_normals".format(mapper="bwa")
    ).values()
    tpl = "output/PureCN/out/unknown_unknown.{ext}{chksum}"
    expected += [
        tpl.format(ext=ext, chksum=chksum)
        for ext in ("list", "bed.gz", "bed.gz.tbi")
        for chksum in ("", ".md5")
    ]
    expected += get_expected_log_files_dict(
        base_out="output/PureCN/log/unknown_unknown".format(mapper="bwa")
    ).values()

    expected = list(sorted(expected))
    actual = list(sorted(panel_of_normals_workflow.get_result_files()))
    assert actual == expected
