# -*- coding: utf-8 -*-
"""Tests for the homologous_recombination_deficiency module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.homologous_recombination_deficiency import (
    HomologousRecombinationDeficiencyWorkflow,
)

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs

__author__ = "Eric Blanc"


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
            bwa:
              path_index: /path/to/bwa/index.fasta
          homologous_recombination_deficiency:
            tools: ['scarHRD']
            path_ngs_mapping: ../ngs_mapping  # REQUIRED

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
def homologous_recombination_deficiency_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return HomologousRecombinationDeficiencyWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return HomologousRecombinationDeficiencyWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for ScarHRDStepPart ------------------------------------------------------------------


def test_scarHRD_step_part_get_input_files_run(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_input_files() - run"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "lib_path": "work/R_packages/out/.done",
        "gc": "work/static_data/out/grch37_50.wig.gz",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = homologous_recombination_deficiency_workflow.get_input_files("scarHRD", "run")(
        wildcards
    )
    assert actual == expected


def test_scarHRD_step_part_get_output_files_run(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_output_files() - run"""
    # Define expected
    base_name_out = "work/{mapper}.scarHRD.{library_name}/out/{mapper}.scarHRD.{library_name}"
    expected = {
        "sequenza": base_name_out + ".seqz.gz",
        "scarHRD": base_name_out + ".json",
    }
    # Get actual
    actual = homologous_recombination_deficiency_workflow.get_output_files("scarHRD", "run")
    assert actual == expected


def test_scarHRD_step_part_get_log_file_run(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_log_file() - run"""
    base_name = "work/{mapper}.scarHRD.{library_name}/log/{mapper}.scarHRD.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = homologous_recombination_deficiency_workflow.get_log_file("scarHRD", "run")
    assert actual == expected


def test_scarHRD_step_part_get_resource_usage_run(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_resource() - run"""
    # Define expected
    expected_dict = {"threads": 2, "time": "48:00:00", "memory": "32G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = homologous_recombination_deficiency_workflow.get_resource(
            "scarHRD", "run", resource
        )
        assert actual == expected, msg_error


def test_scarHRD_step_part_get_output_files_install(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_output_files() - install"""
    # Define expected
    expected = {"lib_path": "work/R_packages/out/.done"}
    # Get actual
    actual = homologous_recombination_deficiency_workflow.get_output_files("scarHRD", "install")
    assert actual == expected


def test_scarHRD_step_part_get_log_file_install(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_log_file() - install"""
    base_name = "work/R_packages/log/R_packages"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = homologous_recombination_deficiency_workflow.get_log_file("scarHRD", "install")
    assert actual == expected


def test_scarHRD_step_part_get_resource_usage_install(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_resource() - install"""
    # Define expected
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "4096M", "partition": "short"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = homologous_recombination_deficiency_workflow.get_resource(
            "scarHRD", "install", resource
        )
        assert actual == expected, msg_error


def test_scarHRD_step_part_get_output_files_gcreference(
    homologous_recombination_deficiency_workflow,
):
    """Tests ScarHRDStepPart.get_output_files() - gcreference"""
    # Define expected
    expected = {"gc": "work/static_data/out/grch37_50.wig.gz"}
    # Get actual
    actual = homologous_recombination_deficiency_workflow.get_output_files("scarHRD", "gcreference")
    assert actual == expected


def test_scarHRD_step_part_get_log_file_gcreference(homologous_recombination_deficiency_workflow):
    """Tests ScarHRDStepPart.get_log_file() - gcreference"""
    base_name = "work/static_data/log/grch37_50"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = homologous_recombination_deficiency_workflow.get_log_file("scarHRD", "gcreference")
    assert actual == expected


# Tests for SomaticMsiCallingWorkflow --------------------------------------------------------------


def test_homologous_recombination_deficiency_workflow(homologous_recombination_deficiency_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "scarHRD"]
    assert list(sorted(homologous_recombination_deficiency_workflow.sub_steps.keys())) == expected
    # Check result file construction
    expected = [
        "output/bwa.scarHRD.P001-T1-DNA1-WGS1/out/bwa.scarHRD.P001-T1-DNA1-WGS1.json",
        "output/bwa.scarHRD.P001-T1-DNA1-WGS1/out/bwa.scarHRD.P001-T1-DNA1-WGS1.seqz.gz",
        "output/bwa.scarHRD.P002-T1-DNA1-WGS1/out/bwa.scarHRD.P002-T1-DNA1-WGS1.json",
        "output/bwa.scarHRD.P002-T1-DNA1-WGS1/out/bwa.scarHRD.P002-T1-DNA1-WGS1.seqz.gz",
        "output/bwa.scarHRD.P002-T2-DNA1-WGS1/out/bwa.scarHRD.P002-T2-DNA1-WGS1.json",
        "output/bwa.scarHRD.P002-T2-DNA1-WGS1/out/bwa.scarHRD.P002-T2-DNA1-WGS1.seqz.gz",
    ]
    actual = set(homologous_recombination_deficiency_workflow.get_result_files())
    expected = set(expected)
    assert actual == expected
