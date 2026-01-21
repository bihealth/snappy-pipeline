# -*- coding: utf-8 -*-
"""Tests for the targeted_seq_mei_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.targeted_seq_mei_calling import MeiWorkflow

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

          targeted_seq_mei_calling:
            path_ngs_mapping: ../ngs_mapping
            scramble:
              blast_ref: /path/to/blast_ref.fa

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
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return MeiWorkflow object pre-configured with germline sheet"""
    # Create reference genome file
    germline_sheet_fake_fs.fs.create_file(
        "/path/to/blast_ref.fa",
        contents=">1\nGATACA",
        create_missing_dirs=True,
    )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.targeted_seq_mei_calling", germline_sheet_fake_fs, mocker
    )

    # Construct the workflow object
    return MeiWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for ScrambleStepPart (cluster) -------------------------------------------------------------


def test_scramble_cluster_step_part_get_input_files(mei_workflow):
    """Tests ScrambleStepPart._get_input_files_cluster()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    expected = ["../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam"]
    actual = mei_workflow.get_input_files("scramble", "cluster")(wildcards)
    assert actual == expected


def test_scramble_cluster_step_part_get_output_files(mei_workflow):
    """Tests ScrambleStepPart._get_output_files_cluster()"""
    pattern_out = "work/{mapper}.scramble.{library_name}/out/{mapper}.scramble.{library_name}"
    expected = {"txt": pattern_out + "_cluster.txt"}
    actual = mei_workflow.get_output_files("scramble", "cluster")
    assert actual == expected


def test_scramble_cluster_step_part_get_log_file(mei_workflow):
    """Tests ScrambleStepPart._get_log_files_cluster()"""
    expected = (
        "work/{mapper}.scramble.{library_name}/log/{mapper}.scramble.{library_name}_cluster.log"
    )
    actual = mei_workflow.get_log_file("scramble", "cluster")
    assert actual == expected


# Tests for ScrambleStepPart (analysis) ------------------------------------------------------------


def test_scramble_analysis_step_part_get_input_files(mei_workflow):
    """Tests ScrambleStepPart._get_input_files_analysis()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    expected = [
        "work/bwa.scramble.P001-N1-DNA1-WGS1/out/bwa.scramble.P001-N1-DNA1-WGS1_cluster.txt"
    ]
    actual = mei_workflow.get_input_files("scramble", "analysis")(wildcards)
    assert actual == expected


def test_scramble_analysis_step_part_get_output_files(mei_workflow):
    """Tests ScrambleStepPart._get_output_files_analysis()"""
    pattern_out = "work/{mapper}.scramble.{library_name}/out/{mapper}.scramble.{library_name}"
    expected = {
        "txt": pattern_out + "_MEIs.txt",
        "txt_md5": pattern_out + "_MEIs.txt.md5",
        "vcf": pattern_out + ".vcf",
        "vcf_gz": pattern_out + ".vcf.gz",
        "vcf_gz_md5": pattern_out + ".vcf.gz.md5",
        "vcf_tbi": pattern_out + ".vcf.gz.tbi",
        "vcf_tbi_md5": pattern_out + ".vcf.gz.tbi.md5",
    }
    actual = mei_workflow.get_output_files("scramble", "analysis")
    assert actual == expected


def test_scramble_analysis_step_part_get_log_file(mei_workflow):
    """Tests ScrambleStepPart._get_log_files_analysis()"""
    expected = (
        "work/{mapper}.scramble.{library_name}/log/{mapper}.scramble.{library_name}_analysis.log"
    )
    actual = mei_workflow.get_log_file("scramble", "analysis")
    assert actual == expected


def test_scramble_analysis_step_part_get_args(mei_workflow):
    """Tests ScrambleStepPart._get_analysis_args()"""
    expected = {
        "reference_genome": "/path/to/blast_ref.fa",
        "mei_refs": None,
        "n_cluster": 5,
        "mei_score": 50,
        "indel_score": 80,
        "mei_polya_frac": 0.75,
    }
    actual = mei_workflow.get_args("scramble", "analysis")(None)
    assert actual == expected


def test_scramble_analysis_step_part_get_resource_usage(mei_workflow):
    """Tests ScrambleStepPart.get_resource_usage()"""
    expected_dict = {"threads": 1, "time": "06:00:00", "memory": "8192M", "partition": "medium"}
    for action in ("analysis", "cluster"):
        for resource, expected in expected_dict.items():
            msg_error = (
                f"Assertion error for resource '{resource}' associated with action '{action}'."
            )
            actual = mei_workflow.get_resource("scramble", action, resource)()
            assert actual == expected, msg_error


# Tests for MeiWorkflow      -----------------------------------------------------------------------


def test_mei_workflow_files(mei_workflow):
    """Tests MeiWorkflow.get_result_files()
    Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results for scramble.
    """
    pattern_out = (
        "output/bwa.scramble.P00{i}-N1-DNA1-WGS1/out/bwa.scramble.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        pattern_out.format(i=i, ext=ext)
        for i in range(1, 7)  # all donors: P001 - P006
        for ext in (
            "vcf.gz",
            "vcf.gz.md5",
            "vcf.gz.tbi",
            "vcf.gz.tbi.md5",
        )
    ]
    actual = mei_workflow.get_result_files()
    assert sorted(actual) == sorted(expected)
