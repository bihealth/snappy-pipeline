# -*- coding: utf-8 -*-
"""Tests for the hla_typing workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow

from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    return yaml.round_trip_load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa

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
def hla_typing_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return HlaTypingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Construct the workflow object
    return HlaTypingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for OptiTypeStepPart ----------------------------------------------------------------------


def test_optitype_step_part_get_input_files(hla_typing_workflow):
    """Tests OptiTypeStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = hla_typing_workflow.get_input_files("optitype", "run")
    assert actual == expected


def test_optitype_step_part_get_output_files(hla_typing_workflow):
    """Tests OptiTypeStepPart.get_output_files()"""
    expected = {
        "cov_pdf": "work/optitype.{library_name}/out/optitype.{library_name}_coverage_plot.pdf",
        "tsv": "work/optitype.{library_name}/out/optitype.{library_name}_result.tsv",
        "txt": "work/optitype.{library_name}/out/optitype.{library_name}.txt",
        "txt_md5": "work/optitype.{library_name}/out/optitype.{library_name}.txt.md5",
    }
    actual = hla_typing_workflow.get_output_files("optitype", "run")
    assert actual == expected


def test_optitype_step_part_get_seq_type_rna(hla_typing_workflow):
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-RNA1-mRNA_seq1"})
    sinput = hla_typing_workflow.substep_dispatch("optitype", "get_args", "run")(wildcards)
    assert sinput["seq_type"] == "rna"


def test_optitype_step_part_get_seq_type_dna(hla_typing_workflow):
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1"})
    sinput = hla_typing_workflow.substep_dispatch("optitype", "get_args", "run")(wildcards)
    assert sinput["seq_type"] == "dna"


def test_optitype_step_part_get_args_input(hla_typing_workflow):
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1"})
    sinput = hla_typing_workflow.substep_dispatch("optitype", "get_args", "run")(wildcards)
    assert sinput["input"]["reads_left"] == [
        "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001-T1-DNA1-WGS1_R1.fastq.gz"
    ]
    assert sinput["input"]["reads_right"] == [
        "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001-T1-DNA1-WGS1_R2.fastq.gz"
    ]


def test_optitype_step_part_get_log_file(hla_typing_workflow):
    """Tests OptiTypeStepPart.get_log_file()"""
    expected = "work/optitype.{library_name}/log/snakemake.hla_typing.log"
    actual = hla_typing_workflow.get_log_file("optitype", "run")
    assert actual == expected


def test_optitype_step_part_get_resource_usage(hla_typing_workflow):
    """Tests OptiTypeStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 6, "time": "40:00:00", "memory": "45000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = hla_typing_workflow.get_resource("optitype", "run", resource)
        assert actual == expected, msg_error


# Tests for ArcasHlaStepPart ----------------------------------------------------------------------


def test_arcashla_step_part_get_input_files(hla_typing_workflow):
    """Tests ArcasHlaStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected_keys = ("ref_done", "bam")
    expected_ref = "work/arcashla.prepare_reference/out/.done"
    actual = hla_typing_workflow.get_input_files("arcashla", "run")(wildcards)
    assert all([key in expected_keys for key in actual])
    assert actual.get("ref_done") == expected_ref


def test_arcashla_step_part_get_output_files(hla_typing_workflow):
    """Tests ArcasHlaStepPart.get_output_files()"""
    expected = {
        "txt": "work/star.arcashla.{library_name}/out/star.arcashla.{library_name}.txt",
        "txt_md5": "work/star.arcashla.{library_name}/out/star.arcashla.{library_name}.txt.md5",
    }
    actual = hla_typing_workflow.get_output_files("arcashla", "run")
    assert actual == expected


def test_arcashla_step_part_get_log_file(hla_typing_workflow):
    """Tests ArcasHlaStepPart.get_log_file()"""
    expected = "work/arcashla.{library_name}/log/snakemake.hla_typing.log"
    actual = hla_typing_workflow.get_log_file("arcashla", "run")
    assert actual == expected


def test_arcashla_step_part_get_resource_usage(hla_typing_workflow):
    """Tests ArcasHlaStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 4, "time": "60:00:00", "memory": "15000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = hla_typing_workflow.get_resource("arcashla", "run", resource)
        assert actual == expected, msg_error


# Tests for HlaTypingWorkflow ---------------------------------------------------------------------


def test_hla_typing_workflow(hla_typing_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = ["arcashla", "link_in", "link_out", "optitype"]
    actual = list(sorted(hla_typing_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = [
        "output/optitype.P001-N1-DNA1-WGS1/out/optitype.P001-N1-DNA1-WGS1.txt",
        "output/optitype.P001-N1-DNA1-WGS1/out/optitype.P001-N1-DNA1-WGS1.txt.md5",
        "output/optitype.P001-T1-DNA1-WGS1/out/optitype.P001-T1-DNA1-WGS1.txt",
        "output/optitype.P001-T1-DNA1-WGS1/out/optitype.P001-T1-DNA1-WGS1.txt.md5",
        "output/optitype.P001-T1-RNA1-mRNA_seq1/out/optitype.P001-T1-RNA1-mRNA_seq1.txt",
        "output/optitype.P001-T1-RNA1-mRNA_seq1/out/optitype.P001-T1-RNA1-mRNA_seq1.txt.md5",
        "output/optitype.P002-N1-DNA1-WGS1/out/optitype.P002-N1-DNA1-WGS1.txt",
        "output/optitype.P002-N1-DNA1-WGS1/out/optitype.P002-N1-DNA1-WGS1.txt.md5",
        "output/optitype.P002-T1-DNA1-WGS1/out/optitype.P002-T1-DNA1-WGS1.txt",
        "output/optitype.P002-T1-DNA1-WGS1/out/optitype.P002-T1-DNA1-WGS1.txt.md5",
        "output/optitype.P002-T1-DNA1-WGS2/out/optitype.P002-T1-DNA1-WGS2.txt",
        "output/optitype.P002-T1-DNA1-WGS2/out/optitype.P002-T1-DNA1-WGS2.txt.md5",
        "output/optitype.P002-T2-DNA1-WGS1/out/optitype.P002-T2-DNA1-WGS1.txt",
        "output/optitype.P002-T2-DNA1-WGS1/out/optitype.P002-T2-DNA1-WGS1.txt.md5",
        "output/optitype.P002-T2-RNA1-mRNA_seq1/out/optitype.P002-T2-RNA1-mRNA_seq1.txt",
        "output/optitype.P002-T2-RNA1-mRNA_seq1/out/optitype.P002-T2-RNA1-mRNA_seq1.txt.md5",
    ]
    actual = hla_typing_workflow.get_result_files()
    assert actual == expected
