# -*- coding: utf-8 -*-
"""Tests for the repeat_expansion workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.repeat_expansion import RepeatExpansionWorkflow

from .conftest import patch_module_fs


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

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fasta
          repeat_expansion:
            repeat_catalog: DUMMY
            repeat_annotation: DUMMY

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
def repeat_expansion_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return RepeatExpansionWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return RepeatExpansionWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_repeat_expansion_workflow_files(repeat_expansion_workflow):
    """Tests RepeatExpansionWorkflow.get_result_files()

    Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results for ExpansionHunter.
    """
    # Define expected
    pattern_json_out = (
        "output/bwa.expansionhunter_annotated.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.expansionhunter_annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    pattern_vcf_out = (
        "output/bwa.expansionhunter.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.expansionhunter.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        pattern_json_out.format(i=i, ext=ext)
        for i in range(1, 7)  # all donors: P001 - P006
        for ext in (
            "json",
            "json.md5",
        )
    ]
    expected += [
        pattern_vcf_out.format(i=i, ext=ext)
        for i in range(1, 7)  # all donors: P001 - P006
        for ext in (
            "vcf",
            "vcf.md5",
        )
    ]
    # Get actual
    actual = repeat_expansion_workflow.get_result_files()
    assert sorted(actual) == sorted(expected)


# Tests for ExpansionHunterStepPart (run) ----------------------------------------------------------


def test_expansionhunter_run_step_part_get_input_files(repeat_expansion_workflow):
    """Tests ExpansionHunterStepPart._get_input_files_run()"""
    # Define expected
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "reference": "/path/to/ref.fa",
        "repeat_catalog": "DUMMY",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = repeat_expansion_workflow.get_input_files("expansionhunter", "run")(wildcards)
    assert actual == expected


def test_expansionhunter_run_step_part_get_output_files(repeat_expansion_workflow):
    """Tests ExpansionHunterStepPart._get_output_files_run()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.expansionhunter.{library_name}/out/"
        "{mapper}.expansionhunter.{library_name}"
    )
    expected = {
        "json": pattern_out + ".json",
        "vcf": pattern_out + ".vcf",
        "vcf_md5": pattern_out + ".vcf.md5",
    }
    # Get actual
    actual = repeat_expansion_workflow.get_output_files("expansionhunter", "run")
    assert actual == expected


def test_expansionhunter_run_step_part_get_log_file(repeat_expansion_workflow):
    """Tests RepeatExpansionWorkflow._get_log_files_run()"""
    # Define expected
    expected = (
        "work/{mapper}.expansionhunter.{library_name}/log/"
        "{mapper}.expansionhunter.{library_name}.log"
    )
    # Get actual
    actual = repeat_expansion_workflow.get_log_file("expansionhunter", "run")
    assert actual == expected


def test_expansionhunter_run_step_part_get_args(repeat_expansion_workflow):
    """Tests RepeatExpansionWorkflow.get_args()"""
    # P001
    expected = "female"
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    actual = repeat_expansion_workflow.get_args("expansionhunter", "run")(wildcards)
    assert actual.get("sex") == expected, "P001-N1-DNA1-WGS1 associated with a female patient."
    # P002
    expected = "male"
    wildcards = Wildcards(fromdict={"library_name": "P002-N1-DNA1-WGS1"})
    actual = repeat_expansion_workflow.get_args("expansionhunter", "run")(wildcards)
    assert actual.get("sex") == expected, "P001-N1-DNA1-WGS1 associated with a male patient."


def test_expansionhunter_step_part_get_resource_usage(repeat_expansion_workflow):
    """Tests RepeatExpansionWorkflow.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = repeat_expansion_workflow.get_resource("expansionhunter", "run", resource)()
        assert actual == expected, msg_error


# Tests for ExpansionHunterStepPart (annotate) -----------------------------------------------------


def test_expansionhunter_annotate_step_part_get_input_files(repeat_expansion_workflow):
    """Tests ExpansionHunterStepPart._get_input_files_annotate()"""
    # Define expected
    json_out = (
        "work/{mapper}.expansionhunter.{library_name}/out/"
        "{mapper}.expansionhunter.{library_name}.json"
    )
    expected = [json_out]

    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = repeat_expansion_workflow.get_input_files("expansionhunter", "annotate")(wildcards)
    assert actual == expected


def test_expansionhunter_annotate_step_part_get_output_files(repeat_expansion_workflow):
    """Tests ExpansionHunterStepPart._get_output_files_annotate()"""
    # Define expected
    json_out = (
        "work/{mapper}.expansionhunter_annotated.{library_name}/out/"
        "{mapper}.expansionhunter_annotated.{library_name}"
    )
    expected = {"json": json_out + ".json", "json_md5": json_out + ".json.md5"}
    # Get actual
    actual = repeat_expansion_workflow.get_output_files("expansionhunter", "annotate")
    assert actual == expected


def test_expansionhunter_annotate_step_part_get_resource_usage(repeat_expansion_workflow):
    """Tests ExpansionHunterStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = repeat_expansion_workflow.get_resource("expansionhunter", "annotate", resource)()
        assert actual == expected, msg_error
