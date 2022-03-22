# -*- coding: utf-8 -*-
"""Tests for the panel_of_normals workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.panel_of_normals import PanelOfNormalsWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (cancer) configuration"""
    return yaml.round_trip_load(
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
            "tumor_library": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "normal_bam": ["NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam"],
        "normal_bai": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai"
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "prepare_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_input_files_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "tumor_library": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "txt": "work/{mapper}.mutect2.select_panel.txt",
        "vcf": [
            "work/bwa.mutect2.prepare_panel/out/P001-N1-DNA1-WGS1.vcf.gz",
            "work/bwa.mutect2.prepare_panel/out/P002-N1-DNA1-WGS1.vcf.gz",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "create_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_output_files_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_prepare_panel()"""
    # TODO: Potential extension error in output files, `vcf.tbi.gz` instead of `vcf.gz.tbi`.
    # TODO: If confirmed, create expected using `get_expected_output_vcf_files_dict()`
    expected = {
        "vcf": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.md5",
        "tbi": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.tbi.gz",
        "tbi_md5": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.tbi.md5",
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
    base_name_out = "work/{mapper}.mutect2.prepare_panel/log/{mapper}.mutect2.{normal_library}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2.create_panel/log/{mapper}.mutect2.create_panel"
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
        "partition": None,
    }
    prepare_panel_expected_dict = {
        "threads": 2,
        "time": "3-00:00:00",
        "memory": "3.7G",
        "partition": None,
    }
    run_expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": None}

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

    # Evaluate action `run`
    for resource, expected in run_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'run'."
        actual = panel_of_normals_workflow.get_resource("mutect2", "run", resource)
        assert actual == expected, msg_error


# PanelOfNormalsWorkflow  --------------------------------------------------------------------------


def test_panel_of_normals_workflow(panel_of_normals_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "mutect2"]
    actual = list(sorted(panel_of_normals_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = "output/{mapper}.mutect2.create_panel/out/{mapper}.mutect2.panel_of_normals.{ext}"
    expected = [
        tpl.format(mapper=mapper, ext=ext)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
    ]
    # add log files
    tpl = "output/{mapper}.mutect2.create_panel/log/{mapper}.mutect2.create_panel.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext)
        for ext in (
            "conda_info.txt",
            "conda_list.txt",
            "log",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
            "log.md5",
        )
        for mapper in ("bwa",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(panel_of_normals_workflow.get_result_files()))
    assert actual == expected
