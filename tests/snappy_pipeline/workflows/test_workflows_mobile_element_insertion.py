# -*- coding: utf-8 -*-
"""Tests for the mobile_element_insertion workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.mobile_element_insertion import MEIWorkflow

from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    return yaml.round_trip_load(
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
def mei_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return MEIWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return MEIWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_mei_workflow_files(mei_workflow):
    """Tests MEIWorkflow::get_result_files()
    Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results for scramble.
    """
    # Define expected
    pattern_out = (
        "output/bwa.scramble_annotated.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.scramble_annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        pattern_out.format(i=i, ext=ext)
        for i in range(1, 7)  # all donors: P001 - P006
        for ext in (
            "json",
            "json.md5",
        )
    ]
    # Get actual
    actual = mei_workflow.get_result_files()
    assert sorted(actual) == sorted(expected)


# Tests for ScrambleStepPart (cluster) -------------------------------------------------------------


def test_scramble_cluster_step_part_get_input_files(mei_workflow):
    """Tests ScrambleStepPart::_get_input_files_cluster()"""
    # Define expected
    expected = ["NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam"]
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = mei_workflow.get_input_files("scramble", "cluster")(wildcards)
    assert actual == expected


def test_scramble_cluster_step_part_get_output_files(mei_workflow):
    """Tests ScrambleStepPart::_get_output_files_cluster()"""
    # Define expected
    pattern_out = "work/{mapper}.scramble.{library_name}/out/{mapper}.scramble.{library_name}"
    expected = {"txt": pattern_out + "_cluster.txt"}
    # Get actual
    actual = mei_workflow.get_output_files("scramble", "cluster")
    assert actual == expected


def test_scramble_cluster_step_part_get_log_file(mei_workflow):
    """Tests ScrambleStepPart::_get_log_files_cluster()"""
    # Define expected
    expected = (
        "work/{mapper}.scramble.{library_name}/log/{mapper}.scramble.{library_name}_cluster.log"
    )
    # Get actual
    actual = mei_workflow.get_log_file("scramble", "cluster")
    assert actual == expected


# Tests for ScrambleStepPart (analysis) ------------------------------------------------------------


def test_scramble_analysis_step_part_get_input_files(mei_workflow):
    """Tests ScrambleStepPart::_get_input_files_analysis()"""
    # Define expected
    expected = ["work/bwa.scramble.P001-N1-DNA1-WGS1/out/bwa.scramble.P001-N1-DNA1-WGS1.txt"]
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = mei_workflow.get_input_files("scramble", "analysis")(wildcards)
    assert actual == expected


def test_scramble_analysis_step_part_get_output_files(mei_workflow):
    """Tests ScrambleStepPart::_get_output_files_analysis()"""
    # Define expected
    pattern_out = "work/{mapper}.scramble.{library_name}/out/{mapper}.scramble.{library_name}"
    expected = {"txt": pattern_out + "_MEIs.txt"}
    # Get actual
    actual = mei_workflow.get_output_files("scramble", "analysis")
    assert actual == expected


def test_scramble_analysis_step_part_get_log_file(mei_workflow):
    """Tests ScrambleStepPart::_get_log_files_analysis()"""
    # Define expected
    expected = (
        "work/{mapper}.scramble.{library_name}/log/{mapper}.scramble.{library_name}_analysis.log"
    )
    # Get actual
    actual = mei_workflow.get_log_file("scramble", "analysis")
    assert actual == expected


def test_scramble_analysis_step_part_get_parameters(mei_workflow):
    """Tests ScrambleStepPart::_get_analysis_parameters()"""
    expected = {"rscript": "REQUIRED/SCRAMble.R", "mei_refs": "resources/MEI_consensus_seqs.fa"}
    # Get actual
    actual = mei_workflow.get_params("scramble", "analysis")(None)
    assert actual == expected


# Tests for ScrambleStepPart (annotate) ------------------------------------------------------------


def test_scramble_annotate_step_part_get_input_files(mei_workflow):
    """Tests ScrambleStepPart::_get_input_files_annotate()"""
    # Define expected
    expected = ["work/bwa.scramble.P001-N1-DNA1-WGS1/out/bwa.scramble.P001-N1-DNA1-WGS1_MEIs.txt"]
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = mei_workflow.get_input_files("scramble", "annotate")(wildcards)
    assert actual == expected


def test_scramble_annotate_step_part_get_output_files(mei_workflow):
    """Tests ScrambleStepPart::_get_output_files_annotate()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.scramble_annotated.{library_name}/out/"
        "{mapper}.scramble_annotated.{library_name}"
    )
    expected = {"json": pattern_out + ".json", "json_md5": pattern_out + ".json.md5"}
    # Get actual
    actual = mei_workflow.get_output_files("scramble", "annotate")
    assert actual == expected


def test_scramble_annotate_step_part_get_log_file(mei_workflow):
    """Tests ScrambleStepPart::_get_log_files_annotate()"""
    # Define expected
    expected = (
        "work/{mapper}.scramble_annotated.{library_name}/log/"
        "{mapper}.scramble_annotated.{library_name}.log"
    )
    # Get actual
    actual = mei_workflow.get_log_file("scramble", "annotate")
    assert actual == expected
