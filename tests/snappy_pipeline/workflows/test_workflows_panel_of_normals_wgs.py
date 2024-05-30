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
            bwa:
              path_index: /path/to/bwa/index.fa

          panel_of_normals:
            path_ngs_mapping: NGS_MAPPING/
            tools: ['cnvkit']
            cnvkit:
              path_target_regions: ""  # WGS mode
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
    aligner_indices_fake_fs,
    mocker,
):
    """Return PanelOfNormalsWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
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


# Tests for CnvkitStepPart ------------------------------------------------------------------------


def test_cnvkit_step_part_get_input_files_target(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_input_files_target()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "bams": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        ],
        "bais": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "target")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_antitarget(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_input_files_antitarget()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    actual = panel_of_normals_workflow.get_input_files("cnvkit", "antitarget")(wildcards)
    assert actual == {}


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
        "logs": [
            "work/bwa.cnvkit/log/bwa.cnvkit.P001-N1-DNA1-WGS1.coverage.log",
            "work/bwa.cnvkit/log/bwa.cnvkit.P001-N1-DNA1-WGS1.coverage.conda_list.txt",
            "work/bwa.cnvkit/log/bwa.cnvkit.P001-N1-DNA1-WGS1.coverage.conda_info.txt",
            "work/bwa.cnvkit/log/bwa.cnvkit.P002-N1-DNA1-WGS1.coverage.log",
            "work/bwa.cnvkit/log/bwa.cnvkit.P002-N1-DNA1-WGS1.coverage.conda_list.txt",
            "work/bwa.cnvkit/log/bwa.cnvkit.P002-N1-DNA1-WGS1.coverage.conda_info.txt",
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
        "log": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.merged.tar.gz",
        "log_md5": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.merged.tar.gz.md5",
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
        actual = panel_of_normals_workflow.get_resource("cnvkit", "target", resource)()
        assert actual == expected, msg_error

    # Evaluate action `antitarget`
    for resource, expected in antitarget_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'antitarget'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "antitarget", resource)()
        assert actual == expected, msg_error

    # Evaluate action `coverage`
    for resource, expected in coverage_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'coverage'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "coverage", resource)()
        assert actual == expected, msg_error

    # Evaluate action `create_panel`
    for resource, expected in reference_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'create_panel'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "create_panel", resource)()
        assert actual == expected, msg_error

    # Evaluate action `report`
    for resource, expected in report_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'report'."
        actual = panel_of_normals_workflow.get_resource("cnvkit", "report", resource)()
        assert actual == expected, msg_error


# PanelOfNormalsWorkflow  --------------------------------------------------------------------------


def test_panel_of_normals_workflow(panel_of_normals_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["access", "cnvkit", "link_out", "mutect2", "purecn"]
    actual = list(sorted(panel_of_normals_workflow.sub_steps.keys()))
    assert actual == expected
    expected = []

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
    # add merged log
    tpl = "output/{mapper}.cnvkit/log/{mapper}.cnvkit.merged.tar.gz{chksum}"
    for mapper in ("bwa",):
        for chksum in ("", ".md5"):
            expected += [tpl.format(mapper=mapper, chksum=chksum)]

    expected = list(sorted(expected))
    actual = list(sorted(panel_of_normals_workflow.get_result_files()))
    assert actual == expected
