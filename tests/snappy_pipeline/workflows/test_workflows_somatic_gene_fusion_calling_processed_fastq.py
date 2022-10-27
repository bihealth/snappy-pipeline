# -*- coding: utf-8 -*-
"""Tests for the somatic_gene_fusion_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_gene_fusion_calling import SomaticGeneFusionCallingWorkflow

from .conftest import patch_module_fs

# TODO: implement tests for all steps with `get_args()`.


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
          somatic_gene_fusion_calling:
              path_link_in: /preprocess
              tools: ['fusioncatcher', 'jaffa', 'arriba']
              fusioncatcher:
                data_dir: REQUIRED   # REQUIRED
              pizzly:
                kallisto_index: REQUIRED    # REQUIRED
                transcripts_fasta: REQUIRED # REQUIRED
                annotations_gtf: REQUIRED       # REQUIRED
              hera:
                path_index: REQUIRED   # REQUIRED
                path_genome: REQUIRED  # REQUIRED
              star_fusion:
                path_ctat_resource_lib: REQUIRED
              defuse:
                path_dataset_directory: REQUIRED
              arriba:
                path_index: /path/to/star/index
                features: /path/to/features.gtf

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
def somatic_gene_fusion_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs_path_link_in,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticGeneFusionCallingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs_path_link_in, mocker)
    # Patch out files for aligner indices
    patch_module_fs(
        "snappy_pipeline.workflows.somatic_gene_fusion_calling", aligner_indices_fake_fs, mocker
    )
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SomaticGeneFusionCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FusioncatcherStepPart ------------------------------------------------------------------


def test_fusioncatcher_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests FusioncatcherStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("fusioncatcher", "run")
    assert actual == expected


def test_fusioncatcher_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests FusioncatcherStepPart.get_output_files()"""
    expected = {"done": "work/fusioncatcher.{library_name}/out/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_output_files("fusioncatcher", "run")
    assert actual == expected


def test_fusioncatcher_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests FusioncatcherStepPart.get_log_file()"""
    expected = "work/fusioncatcher.{library_name}/log/snakemake.gene_fusion_calling.log"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("fusioncatcher", "run")
    assert actual == expected


def test_fusioncatcher_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests FusioncatcherStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 4, "time": "5-00:00:00", "memory": "30000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("fusioncatcher", "run", resource)
        assert actual == expected, msg_error


# Tests for JaffaStepPart --------------------------------------------------------------------------


def test_jaffa_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests JaffaStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("jaffa", "run")
    assert actual == expected


def test_jaffa_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests JaffaStepPart.get_output_files()"""
    expected = {"done": "work/jaffa.{library_name}/out/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_output_files("jaffa", "run")
    assert actual == expected


def test_jaffa_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests JaffaStepPart.get_log_file()"""
    expected = "work/jaffa.{library_name}/log/snakemake.gene_fusion_calling.log"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("jaffa", "run")
    assert actual == expected


def test_jaffa_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests JaffaStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 4, "time": "5-00:00:00", "memory": "163840M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("jaffa", "run", resource)
        assert actual == expected, msg_error


# Tests for PizzlyStepPart -------------------------------------------------------------------------


def test_pizzly_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests PizzlyStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("pizzly", "run")
    assert actual == expected


def test_pizzly_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests PizzlyStepPart.get_output_files()"""
    expected = {"done": "work/pizzly.{library_name}/out/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_output_files("pizzly", "run")
    assert actual == expected


def test_pizzly_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests PizzlyStepPart.get_log_file()"""
    expected = "work/pizzly.{library_name}/log/snakemake.gene_fusion_calling.log"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("pizzly", "run")
    assert actual == expected


def test_pizzly_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests PizzlyStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 4, "time": "5-00:00:00", "memory": "81920M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("pizzly", "run", resource)
        assert actual == expected, msg_error


# Tests for StarFusionStepPart ---------------------------------------------------------------------


def test_star_fusion_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests StarFusionStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("star_fusion", "run")
    assert actual == expected


def test_star_fusion_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests StarFusionStepPart.get_output_files()"""
    expected = {"done": "work/star_fusion.{library_name}/out/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_output_files("star_fusion", "run")
    assert actual == expected


def test_star_fusion_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests StarFusionStepPart.get_log_file()"""
    expected = "work/star_fusion.{library_name}/log/snakemake.gene_fusion_calling.log"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("star_fusion", "run")
    assert actual == expected


def test_star_fusion_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests StarFusionStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 4, "time": "5-00:00:00", "memory": "122880M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("star_fusion", "run", resource)
        assert actual == expected, msg_error


# DefuseStepPart -----------------------------------------------------------------------------------


def test_defuse_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests DefuseStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("defuse", "run")
    assert actual == expected


def test_defuse_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests DefuseStepPart.get_output_files()"""
    expected = {"done": "work/defuse.{library_name}/out/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_output_files("defuse", "run")
    assert actual == expected


def test_defuse_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests DefuseStepPart.get_log_file()"""
    expected = "work/defuse.{library_name}/log/snakemake.gene_fusion_calling.log"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("defuse", "run")
    assert actual == expected


def test_defuse_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests DefuseStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 8, "time": "5-00:00:00", "memory": "81920M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("defuse", "run", resource)
        assert actual == expected, msg_error


# HeraStepPart -------------------------------------------------------------------------------------


def test_hera_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests HeraStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("hera", "run")
    assert actual == expected


def test_hera_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests HeraStepPart.get_output_files()"""
    expected = {"done": "work/hera.{library_name}/out/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_output_files("hera", "run")
    assert actual == expected


def test_hera_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests HeraStepPart.get_log_file()"""
    expected = "work/hera.{library_name}/log/snakemake.gene_fusion_calling.log"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("hera", "run")
    assert actual == expected


def test_hera_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests HeraStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 8, "time": "5-00:00:00", "memory": "163840M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("hera", "run", resource)
        assert actual == expected, msg_error


# ArribaStepPart -----------------------------------------------------------------------------------


def test_arriba_step_part_get_input_files(somatic_gene_fusion_calling_workflow):
    """Tests ArribaStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = somatic_gene_fusion_calling_workflow.get_input_files("arriba", "run")
    assert actual == expected


def test_arriba_step_part_get_args(somatic_gene_fusion_calling_workflow):
    """Tests ArribaStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-RNA1-mRNA_seq1"})
    expected = {
        "input": {
            "reads_left": [
                "work/input_links/P001-T1-RNA1-mRNA_seq1/FCXXXXXX/L001/out/P001_R1.fastq.gz"
            ],
            "reads_right": [
                "work/input_links/P001-T1-RNA1-mRNA_seq1/FCXXXXXX/L001/out/P001_R2.fastq.gz"
            ],
        }
    }
    actual = somatic_gene_fusion_calling_workflow.get_args("arriba", "run")(wildcards)
    assert actual == expected


def test_arriba_step_part_get_output_files(somatic_gene_fusion_calling_workflow):
    """Tests ArribaStepPart.get_output_files()"""
    key_ext = (
        ("fusions", "fusions.tsv"),
        ("discarded", "discarded_fusions.tsv.gz"),
    )
    expected = {"done": "work/arriba.{library_name}/out/.done"}
    prefix = "work/arriba.{library_name}/out/arriba.{library_name}."
    for key, ext in key_ext:
        expected[key] = prefix + ext
        expected[key + "_md5"] = prefix + ext + ".md5"
    actual = somatic_gene_fusion_calling_workflow.get_output_files("arriba", "run")
    assert actual == expected


def test_arriba_step_part_get_log_file(somatic_gene_fusion_calling_workflow):
    """Tests ArribaStepPart.get_log_file()"""
    expected = {}
    key_ext = (
        ("log", "log"),
        ("conda_info", "conda_info.txt"),
        ("conda_list", "conda_list.txt"),
    )
    prefix = "work/arriba.{library_name}/log/arriba.{library_name}."
    for key, ext in key_ext:
        expected[key] = prefix + ext
        expected[key + "_md5"] = prefix + ext + ".md5"
    key_ext = (
        ("out", "Log.out"),
        ("final", "Log.final.out"),
        ("std", "Log.std.out"),
        ("SJ", "SJ.out.tab"),
    )
    prefix = "work/arriba.{library_name}/log/"
    for key, ext in key_ext:
        expected[key] = prefix + ext
        expected[key + "_md5"] = prefix + ext + ".md5"
    actual = somatic_gene_fusion_calling_workflow.get_log_file("arriba", "run")
    assert actual == expected


def test_arriba_step_part_get_resource_usage(somatic_gene_fusion_calling_workflow):
    """Tests ArribaStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 8, "time": "24:00:00", "memory": "65536M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_gene_fusion_calling_workflow.get_resource("arriba", "run", resource)
        assert actual == expected, msg_error


# Tests for SomaticGeneFusionCallingWorkflow -------------------------------------------------------


# def test_somatic_gene_fusion_calling_workflow(somatic_gene_fusion_calling_workflow):
#     """Test simple functionality of the workflow"""
#     # Check created sub steps
#     expected = ["link_out", "mantis"]
#     assert list(sorted(somatic_gene_fusion_calling_workflow.sub_steps.keys())) == expected
#     # Check result file construction
#     expected = [
#         "output/mantis.bwa.P001-T1-DNA1-WGS1/out/mantis.bwa.P001-T1-DNA1-WGS1_results.txt",
#         "output/mantis.bwa.P001-T1-DNA1-WGS1/out/mantis.bwa.P001-T1-DNA1-WGS1_results.txt.status",
#         "output/mantis.bwa.P002-T1-DNA1-WGS1/out/mantis.bwa.P002-T1-DNA1-WGS1_results.txt",
#         "output/mantis.bwa.P002-T1-DNA1-WGS1/out/mantis.bwa.P002-T1-DNA1-WGS1_results.txt.status",
#         "output/mantis.bwa.P002-T2-DNA1-WGS1/out/mantis.bwa.P002-T2-DNA1-WGS1_results.txt",
#         "output/mantis.bwa.P002-T2-DNA1-WGS1/out/mantis.bwa.P002-T2-DNA1-WGS1_results.txt.status",
#     ]
#     actual = set(somatic_gene_fusion_calling_workflow.get_result_files())
#     expected = set(expected)
#     assert actual == expected
