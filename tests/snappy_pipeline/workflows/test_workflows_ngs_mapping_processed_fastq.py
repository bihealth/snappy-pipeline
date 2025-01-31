# -*- coding: utf-8 -*-
"""Tests for the ngs_mapping workflow module code using preprocessed FASTQ"""

import io
import textwrap
from collections import OrderedDict
from copy import deepcopy

import pytest
import ruamel.yaml as ruamel_yaml
from biomedsheets.io_tsv import read_cancer_tsv_sheet
from snakemake.io import Wildcards

from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .common import get_expected_log_files_dict
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
          features:
            path: /path/to/features.gtf

        step_config:
          ngs_mapping:
            path_link_in: "/preprocess"
            tools:
              dna: ['bwa']
              rna: ['star']
            target_coverage_report:
              path_target_interval_list_mapping:
              - pattern: "Agilent SureSelect Human All Exon V6.*"
                name: Agilent_SureSelect_Human_All_Exon_V6
                path: path/to/SureSelect_Human_All_Exon_V6_r2.bed
            bwa:
              path_index: /path/to/bwa/index.fasta.amb
            star:
              path_index: /path/to/star/index
              transcriptome: true
              out_filter_intron_motifs: ""
              out_sam_strand_field: ""
            minimap2:
              mapping_threads: 16
              path_index: /path/to/minimap2/index
            bam_collect_doc:
              enabled: true

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {"left": "*/*/*_R1.fastq.gz", "right": "*/*/*_R2.fastq.gz"}
            search_paths: ['/path']
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def ngs_mapping_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs_path_link_in,
    aligner_indices_fake_fs,
    mocker,
):
    """Return NgsMappingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs_path_link_in, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Construct the workflow object
    return NgsMappingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def get_expected_output_files_dict(bam_base_out: str, report_base_out: str, log_base_out: str):
    """Helper function.

    :param bam_base_out: Expected name pattern of BAM associated files without extension.
    For example if the full path would be '/path/to/step_part.bam', argument should be
    '/path/to/step_part'.

    :param report_base_out: Expected name pattern of report associated files without extension.
    For example if the full path would be '/path/to/step_report.bam.bamstats.txt', argument should
    be '/path/to/step_report'.

    :param log_base_out: Expected name pattern for log file.

    :return: Returns dictionary with expected path for BAM and report associated files based on the
    provided input.
    """

    def work_to_out(s: str) -> str:
        return s.replace("work/", "output/")

    expected = {
        "bam": f"{bam_base_out}.bam",
        "bam_md5": f"{bam_base_out}.bam.md5",
        "bam_bai": f"{bam_base_out}.bam.bai",
        "bam_bai_md5": f"{bam_base_out}.bam.bai.md5",
        "output_links": list(
            map(
                work_to_out,
                [
                    f"{bam_base_out}.bam",
                    f"{bam_base_out}.bam.bai",
                    f"{bam_base_out}.bam.md5",
                    f"{bam_base_out}.bam.bai.md5",
                    f"{report_base_out}.bam.bamstats.txt",
                    f"{report_base_out}.bam.flagstats.txt",
                    f"{report_base_out}.bam.idxstats.txt",
                    f"{report_base_out}.bam.bamstats.txt.md5",
                    f"{report_base_out}.bam.flagstats.txt.md5",
                    f"{report_base_out}.bam.idxstats.txt.md5",
                    f"{log_base_out}.mapping.log",
                    f"{log_base_out}.mapping.log.md5",
                    f"{log_base_out}.mapping.conda_info.txt",
                    f"{log_base_out}.mapping.conda_info.txt.md5",
                    f"{log_base_out}.mapping.conda_list.txt",
                    f"{log_base_out}.mapping.conda_list.txt.md5",
                    f"{log_base_out}.mapping.wrapper.py",
                    f"{log_base_out}.mapping.wrapper.py.md5",
                    f"{log_base_out}.mapping.environment.yaml",
                    f"{log_base_out}.mapping.environment.yaml.md5",
                ],
            )
        ),
        "report_bamstats_txt": f"{report_base_out}.bam.bamstats.txt",
        "report_bamstats_txt_md5": f"{report_base_out}.bam.bamstats.txt.md5",
        "report_flagstats_txt": f"{report_base_out}.bam.flagstats.txt",
        "report_flagstats_txt_md5": f"{report_base_out}.bam.flagstats.txt.md5",
        "report_idxstats_txt": f"{report_base_out}.bam.idxstats.txt",
        "report_idxstats_txt_md5": f"{report_base_out}.bam.idxstats.txt.md5",
    }

    # Return
    return expected


# Test for project validation ----------------------------------------------------------------------


def test_extraction_type_check(
    ngs_mapping_workflow,
    cancer_sheet_tsv,
):
    """Tests extraction type check method."""
    # Create mix data sample sheet
    cancer_sheet_io = io.StringIO(cancer_sheet_tsv)
    cancer_sheet = read_cancer_tsv_sheet(cancer_sheet_io)

    # Evaluate if both DNA and RNA are True
    dna_bool, rna_bool = ngs_mapping_workflow.extraction_type_check(sample_sheet=cancer_sheet)
    assert dna_bool, "Sample sheet contains both DNA and RNA."
    assert rna_bool, "Sample sheet contains both DNA and RNA."


def test_project_validation_cancer(ngs_mapping_workflow, cancer_sheet_tsv, minimal_config):
    """Tests project validation method in ngs mapping workflow"""
    # Convert yaml to dict
    minimal_config_dict = deepcopy(minimal_config)
    minimal_config_dict = dict(minimal_config_dict)
    minimal_config_dict = minimal_config_dict["step_config"].get("ngs_mapping", OrderedDict())
    config = ngs_mapping_workflow.config_model_class(**minimal_config_dict)

    # Create germline sample sheet
    cancer_sheet_io = io.StringIO(cancer_sheet_tsv)
    cancer_sheet = read_cancer_tsv_sheet(cancer_sheet_io)

    # Method returns None without exception, cause DNA sample sheet and DNA tool defined in config
    out = ngs_mapping_workflow.validate_project(config=config, sample_sheets_list=[cancer_sheet])
    assert out is None, "No exception expected: DNA sample sheet and DNA tool defined in config."


# Tests for BwaStepPart ----------------------------------------------------------------------------


def test_bwa_step_part_get_args(ngs_mapping_workflow):
    """Tests BaseStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    base_name = "work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/"
    expected = {
        "input": {
            "reads_left": [base_name + "P001_R1.fastq.gz"],
            "reads_right": [base_name + "P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-N1-DNA1-WGS1",
        "path_index": "/path/to/bwa/index.fasta",
        "num_threads_align": 16,
        "num_threads_trimming": 8,
        "num_threads_bam_view": 4,
        "num_threads_bam_sort": 4,
        "memory_bam_sort": "4G",
        "trim_adapters": False,
        "mask_duplicates": True,
        "split_as_secondary": False,
        "extra_args": [],
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_args("bwa", "run")(wildcards)
    assert actual == expected


def test_bwa_step_part_get_input_files(ngs_mapping_workflow):
    """Tests BaseStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = "work/input_links/P001-N1-DNA1-WGS1/.done"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("bwa", "run")(wildcards)
    assert actual == expected


def test_bwa_step_part_get_output_files(ngs_mapping_workflow):
    """Tests BaseStepPart.get_output_files()"""
    # Define expected
    bam_base_out = "work/bwa.{library_name}/out/bwa.{library_name}"
    report_base_out = "work/bwa.{library_name}/report/bam_qc/bwa.{library_name}"
    log_base_out = "work/bwa.{library_name}/log/bwa.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out, log_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("bwa", "run")
    assert actual == expected


def test_bwa_step_part_get_log_file(ngs_mapping_workflow):
    """Tests BaseStepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/bwa.{library_name}/log/bwa.{library_name}",
        infix="mapping",
        extended=True,
    )
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("bwa", "run")
    assert actual == expected


def test_bwa_step_part_get_resource(ngs_mapping_workflow):
    """Tests BaseStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "3-00:00:00", "memory": "73728M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("bwa", "run", resource)()
        assert actual == expected, msg_error


# Tests for StarStepPart --------------------------------------------------------------------------


def test_star_step_part_get_args(ngs_mapping_workflow):
    """Tests StarStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-RNA1-mRNA_seq1"})
    expected = {
        "input": {
            "reads_left": ["work/input_links/P001-T1-RNA1-mRNA_seq1/FCXXXXXX/L001/out/P001_R1.fastq.gz"],
            "reads_right": ["work/input_links/P001-T1-RNA1-mRNA_seq1/FCXXXXXX/L001/out/P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-T1-RNA1-mRNA_seq1",
        "path_index": "/path/to/star/index",
        "features": "/path/to/features.gtf",
        "num_threads_align": 16,
        "num_threads_trimming": 8,
        "num_threads_bam_view": 4,
        "num_threads_bam_sort": 4,
        "memory_bam_sort": "4G",
        "genome_load": "NoSharedMemory",
        "raw_star_options": "",
        "align_intron_max": 1000000,
        "align_intron_min": 20 ,
        "align_mates_gap_max": 1000000,
        "align_sjdb_overhang_min": 1,
        "align_sj_overhang_min": 8,
        "out_filter_mismatch_n_max": 999,
        "out_filter_mismatch_n_over_l_max": 0.04,
        "out_filter_multimap_n_max": 20,
        "out_filter_type": "BySJout",
        "out_filter_intron_motifs": "",
        "out_sam_strand_field": "",
        "transcriptome": True,
        "trim_adapters": False,
        "mask_duplicates": False,
        "include_unmapped": True,
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_args("star", "run")(wildcards)
    assert actual == expected


def test_star_step_part_get_input_files(ngs_mapping_workflow):
    """Tests StarStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-RNA1-mRNA_seq1"})
    expected = "work/input_links/P001-T1-RNA1-mRNA_seq1/.done"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("star", "run")(wildcards)
    assert actual == expected


def test_star_step_part_get_output_files(ngs_mapping_workflow):
    """Tests StarStepPart.get_output_files()"""
    # Define expected
    bam_base_out = "work/star.{library_name}/out/star.{library_name}"
    report_base_out = "work/star.{library_name}/report/bam_qc/star.{library_name}"
    log_base_out = "work/star.{library_name}/log/star.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out, log_base_out)
    expected["gene_counts"] = "work/star.{library_name}/out/star.{library_name}.GeneCounts.tab"
    expected["gene_counts_md5"] = (
        "work/star.{library_name}/out/star.{library_name}.GeneCounts.tab.md5"
    )
    expected["junctions"] = "work/star.{library_name}/out/star.{library_name}.Junctions.tab"
    expected["junctions_md5"] = "work/star.{library_name}/out/star.{library_name}.Junctions.tab.md5"
    expected["transcriptome"] = (
        "work/star.{library_name}/out/star.{library_name}.toTranscriptome.bam"
    )
    expected["transcriptome_md5"] = (
        "work/star.{library_name}/out/star.{library_name}.toTranscriptome.bam.md5"
    )
    expected["output_links"].extend(
        [
            "output/star.{library_name}/out/star.{library_name}.Junctions.tab",
            "output/star.{library_name}/out/star.{library_name}.Junctions.tab.md5",
            "output/star.{library_name}/out/star.{library_name}.toTranscriptome.bam",
            "output/star.{library_name}/out/star.{library_name}.toTranscriptome.bam.md5",
        ]
    )
    expected["output_links"].sort()
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("star", "run")
    actual["output_links"].sort()
    assert actual == expected


def test_star_step_part_get_log_file(ngs_mapping_workflow):
    """Tests StarStepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/star.{library_name}/log/star.{library_name}",
        infix="mapping",
        extended=True,
    )
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("star", "run")
    assert actual == expected


def test_star_step_part_get_resource(ngs_mapping_workflow):
    """Tests StarStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "56G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("star", "run", resource)()
        assert actual == expected, msg_error


# Tests for StrandednessStepPart -------------------------------------------------------------------


def test_strandedness_step_part_infer_get_input_files(ngs_mapping_workflow):
    """Tests StrandednessStepPart.get_input_files()"""
    # Define expected
    expected = {
        "bam": "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_input_files("strandedness", "infer")
    assert actual == expected


def test_strandedness_step_part_counts_get_input_files(ngs_mapping_workflow):
    """Tests StrandednessStepPart.get_input_files()"""
    # Define expected
    expected = {
        "decision": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.decision.json",
        "counts": "work/{mapper}.{library_name}/out/{mapper}.{library_name}.GeneCounts.tab",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_input_files("strandedness", "counts")
    assert actual == expected


def test_strandedness_step_part_infer_get_output_files(ngs_mapping_workflow):
    """Tests StrandednessStepPart.get_output_files()"""
    # Define expected
    expected = {
        "tsv": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.infer.txt",
        "tsv_md5": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.infer.txt.md5",
        "decision": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.decision.json",
        "decision_md5": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.decision.json.md5",
        "output": "output/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.decision.json",
        "output_md5": "output/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.decision.json.md5",
        "log": "output/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.log",
        "log_md5": "output/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.log.md5",
        "conda_list": "output/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_list.txt",
        "conda_list_md5": "output/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_list.txt.md5",
        "conda_info": "output/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_info.txt",
        "conda_info_md5": "output/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_info.txt.md5",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("strandedness", "infer")
    assert actual == expected


def test_strandedness_step_part_counts_get_output_files(ngs_mapping_workflow):
    """Tests StrandednessStepPart.get_output_files()"""
    # Define expected
    expected = {
        "counts": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.GeneCounts.tab",
        "counts_md5": "work/{mapper}.{library_name}/strandedness/{mapper}.{library_name}.GeneCounts.tab.md5",
        "output": "output/{mapper}.{library_name}/out/{mapper}.{library_name}.GeneCounts.tab",
        "output_md5": "output/{mapper}.{library_name}/out/{mapper}.{library_name}.GeneCounts.tab.md5",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("strandedness", "counts")
    assert actual == expected


def test_strandedness_step_part_infer_get_log_file(ngs_mapping_workflow):
    """Tests StrandednessStepPart.get_log_file()"""
    # Define expected
    expected = {
        "log": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.log",
        "log_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.log.md5",
        "conda_info": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.strandedness.conda_list.txt.md5",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("strandedness", "infer")
    assert actual == expected


def test_strandedness_step_part_infer_get_resource(ngs_mapping_workflow):
    """Tests StrandednessStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("strandedness", "infer", resource)()
        assert actual == expected, msg_error


# Tests for Minimap2StepPart -----------------------------------------------------------------------


def test_minimap2_step_part_get_args(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": {
            "reads_left": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R1.fastq.gz"],
            "reads_right": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-N1-DNA1-WGS1",
        "mapping_threads": 16,
        "path_index": "/path/to/minimap2/index",
        "extra_infos": {
            "libraryType": "WGS",
            "libraryKit": "Agilent SureSelect Human All Exon V6",
            "folderName": "P001-N1-DNA1-WGS1",
            "extractionType": "DNA",
            "seqPlatform": "Illumina",
        },
        "library_name": "P001-N1-DNA1-WGS1",
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_args("minimap2", "run")(wildcards)
    assert actual == expected


def test_minimap2_step_part_get_input_files(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = "work/input_links/P001-N1-DNA1-WGS1/.done"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("minimap2", "run")(wildcards)
    assert actual == expected


def test_minimap2_step_part_get_output_files(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_output_files()"""
    # Define expected
    bam_base_out = "work/minimap2.{library_name}/out/minimap2.{library_name}"
    report_base_out = "work/minimap2.{library_name}/report/bam_qc/minimap2.{library_name}"
    log_base_out = "work/minimap2.{library_name}/log/minimap2.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out, log_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("minimap2", "run")
    assert actual == expected


def test_minimap2_step_part_get_log_file(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/minimap2.{library_name}/log/minimap2.{library_name}",
        infix="mapping",
        extended=True,
    )
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("minimap2", "run")
    assert actual == expected


def test_minimap2_step_part_get_resource(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "56G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("minimap2", "run", resource)()
        assert actual == expected, msg_error


# Tests for ExternalStepPart -----------------------------------------------------------------------


# TODO(holtgrewe): the fake file system should be setup with fake BAM files.
def test_external_step_part_get_args(ngs_mapping_workflow):
    """Tests ExternalStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {"input": [], "platform": "EXTERNAL", "sample_name": "P001-N1-DNA1-WGS1"}
    # Get actual and assert
    actual = ngs_mapping_workflow.get_args("external", "run")(wildcards)
    assert actual == expected


def test_external_step_part_get_input_files(ngs_mapping_workflow):
    """Tests ExternalStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = "work/input_links/P001-N1-DNA1-WGS1/.done"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("external", "run")(wildcards)
    assert actual == expected


def test_external_step_part_get_output_files(ngs_mapping_workflow):
    """Tests ExternalStepPart.get_output_files()"""
    # Define expected
    bam_base_out = "work/external.{library_name}/out/external.{library_name}"
    report_base_out = "work/external.{library_name}/report/bam_qc/external.{library_name}"
    log_base_out = "work/external.{library_name}/log/external.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out, log_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("external", "run")
    assert actual == expected


def test_external_step_part_get_resource(ngs_mapping_workflow):
    """Tests ExternalStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "00:10:00", "memory": "1G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("external", "run", resource)()
        assert actual == expected, msg_error


# Tests for TargetCoverageReportStepPart ----------------------------------------------------------


def test_target_coverage_report_step_part_run_get_input_files(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_input_files() - run"""
    # Define expected
    expected = {
        "bam": "work/bwa.library/out/bwa.library.bam",
        "bai": "work/bwa.library/out/bwa.library.bam.bai",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "library"})
    actual = ngs_mapping_workflow.get_input_files("target_coverage_report", "run")(wildcards)
    assert actual == expected


def test_target_coverage_report_step_part_get_output_files(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_output_files() - run"""
    expected = {
        "json": "work/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz",
        "json_md5": "work/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz.md5",
        "output_links": [
            "output/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz",
            "output/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.log",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.log.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_info.txt",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_info.txt.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_list.txt",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_list.txt.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.wrapper.py",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.wrapper.py.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.environment.yaml",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.environment.yaml.md5",
        ],
    }
    assert ngs_mapping_workflow.get_output_files("target_coverage_report", "run") == expected


def test_target_coverage_report_step_part_run_get_log_file(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_log_file() - run"""
    expected = "work/{mapper}.{library_name}/log/snakemake.target_coverage.log"
    expected = {
        "log": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.log",
        "log_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.log.md5",
        "conda_info": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.conda_list.txt.md5",
        "wrapper": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.wrapper.py",
        "wrapper_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.wrapper.py.md5",
        "env_yaml": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.environment.yaml",
        "env_yaml_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report.environment.yaml.md5",
    }
    assert ngs_mapping_workflow.get_log_file("target_coverage_report", "run") == expected


def test_target_coverage_report_step_part_get_resource(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_resource()"""
    # Define expected
    actions = ("run",)
    expected_dict = {"threads": 2, "time": "04:00:00", "memory": "20G", "partition": "medium"}
    # Evaluate
    for action in actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = ngs_mapping_workflow.get_resource("target_coverage_report", action, resource)()
            assert actual == expected, msg_error


# Tests for NgsMappingWorkflow --------------------------------------------------------------------


def test_ngs_mapping_workflow_steps(ngs_mapping_workflow):
    """Tests simple functionality of the workflow: checks if sub steps are created,
    i.e., the tools associated with gene expression quantification."""
    # Check created sub steps
    expected = [
        "bam_collect_doc",
        "bwa",
        "bwa_mem2",
        "external",
        "link_in",
        "mbcs",
        "minimap2",
        "ngs_chew",
        "star",
        "strandedness",
        "target_coverage_report",
    ]
    actual = ngs_mapping_workflow.sub_steps.keys()
    assert sorted(actual) == sorted(expected)


def test_ngs_mapping_workflow_files(ngs_mapping_workflow):
    """Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools, namely: bwa, external, link_in, minimap2, star,
    target_coverage_report.
    """
    # Expected library names
    dna = (
        "P001-N1-DNA1-WGS1",
        "P001-T1-DNA1-WGS1",
        "P002-N1-DNA1-WGS1",
        "P002-T1-DNA1-WGS1",
        "P002-T2-DNA1-WGS1",
    )
    rna = ("P001-T1-RNA1-mRNA_seq1", "P002-T2-RNA1-mRNA_seq1")

    # Check result file construction (dna/bwa)
    expected = [
        "output/bwa.{library_name}/out/bwa.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in ("bam", "bam.bai", "bam.md5", "bam.bai.md5")
        for library_name in dna
    ]
    for infix in ("bam_collect_doc", "mapping", "target_cov_report", "ngs_chew_fingerprint"):
        expected += [
            "output/bwa.{library_name}/log/bwa.{library_name}.{ext}".format(
                library_name=library_name, ext=ext
            )
            for ext in (
                f"{infix}.log",
                f"{infix}.log.md5",
                f"{infix}.conda_info.txt",
                f"{infix}.conda_info.txt.md5",
                f"{infix}.conda_list.txt",
                f"{infix}.conda_list.txt.md5",
                f"{infix}.environment.yaml",
                f"{infix}.environment.yaml.md5",
                f"{infix}.wrapper.py",
                f"{infix}.wrapper.py.md5",
            )
            for library_name in dna
        ]
    bam_stats_text_out = (
        "output/bwa.{library_name}/report/bam_qc/bwa.{library_name}.bam.{stats}.{ext}"
    )
    expected += [
        bam_stats_text_out.format(library_name=library_name, stats=stats, ext=ext)
        for ext in ("txt", "txt.md5")
        for library_name in dna
        for stats in ("bamstats", "flagstats", "idxstats")
    ]
    expected += [
        "output/bwa.{library_name}/report/cov/bwa.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in (
            "cov.bw",
            "cov.bw.md5",
            "cov.vcf.gz",
            "cov.vcf.gz.md5",
            "cov.vcf.gz.tbi",
            "cov.vcf.gz.tbi.md5",
            "mq.bw",
            "mq.bw.md5",
        )
        for library_name in dna
    ]
    expected += [
        "output/bwa.{library_name}/report/fingerprint/bwa.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in ("npz", "npz.md5")
        for library_name in dna
    ]
    expected += [
        "output/bwa.{library_name}/report/alfred_qc/bwa.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in ("alfred.json.gz", "alfred.json.gz.md5")
        for library_name in dna
    ]

    # Check result file construction (rna/star)
    expected += [
        "output/star.{library_name}/out/star.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in (
            "bam",
            "bam.bai",
            "bam.md5",
            "bam.bai.md5",
            "GeneCounts.tab",
            "GeneCounts.tab.md5",
            "Junctions.tab",
            "Junctions.tab.md5",
            "toTranscriptome.bam",
            "toTranscriptome.bam.md5",
        )
        for library_name in rna
    ]
    for infix in ("mapping", "ngs_chew_fingerprint"):
        expected += [
            "output/star.{library_name}/log/star.{library_name}.{ext}".format(
                library_name=library_name, ext=ext
            )
            for ext in (
                f"{infix}.log",
                f"{infix}.log.md5",
                f"{infix}.conda_info.txt",
                f"{infix}.conda_info.txt.md5",
                f"{infix}.conda_list.txt",
                f"{infix}.conda_list.txt.md5",
                f"{infix}.environment.yaml",
                f"{infix}.environment.yaml.md5",
                f"{infix}.wrapper.py",
                f"{infix}.wrapper.py.md5",
            )
            for library_name in rna
        ]
    bam_stats_text_out = (
        "output/star.{library_name}/report/bam_qc/star.{library_name}.bam.{stats}.{ext}"
    )
    expected += [
        bam_stats_text_out.format(library_name=library_name, stats=stats, ext=ext)
        for ext in ("txt", "txt.md5")
        for library_name in rna
        for stats in ("bamstats", "flagstats", "idxstats")
    ]
    expected += [
        "output/star.{library_name}/report/fingerprint/star.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in ("npz", "npz.md5")
        for library_name in rna
    ]
    expected += [
        "output/star.{library_name}/strandedness/star.{library_name}.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in ("decision.json", "decision.json.md5")
        for library_name in rna
    ]
    expected += [
        "output/star.{library_name}/log/star.{library_name}.strandedness.{ext}".format(
            library_name=library_name, ext=ext
        )
        for ext in (
            "log",
            "log.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
        )
        for library_name in rna
    ]

    assert sorted(ngs_mapping_workflow.get_result_files()) == sorted(expected)
