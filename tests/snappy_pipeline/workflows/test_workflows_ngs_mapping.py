# -*- coding: utf-8 -*-
"""Tests for the ngs_mapping workflow module code"""

from collections import OrderedDict
from copy import deepcopy
import io
import textwrap

from biomedsheets.io_tsv import read_generic_tsv_sheet, read_germline_tsv_sheet
import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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
            target_coverage_report:
              path_target_interval_list_mapping:
              - pattern: "Agilent SureSelect Human All Exon V6.*"
                name: Agilent_SureSelect_Human_All_Exon_V6
                path: path/to/SureSelect_Human_All_Exon_V6_r2.bed
            compute_coverage_bed: true
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
            pedigree_field: pedigree_field
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
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return NgsMappingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
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


def get_expected_output_files_dict(bam_base_out, report_base_out):
    """Helper method.

    :param bam_base_out: Expected name pattern of BAM associated files without extension.
    For example if the full path would be '/path/to/step_part.bam', argument should be
    '/path/to/step_part'.
    :type bam_base_out: str

    :param report_base_out: Expected name pattern of report associated files without extension.
    For example if the full path would be '/path/to/step_report.bam.bamstats.html', argument should
    be '/path/to/step_report'.

    :return: Returns dictionary with expected path for BAM and report associated files based on the
    provided input.
    """
    # Define expected
    expected = {
        "bam": bam_base_out + ".bam",
        "bam_bai": bam_base_out + ".bam.bai",
        "bam_bai_md5": bam_base_out + ".bam.bai.md5",
        "bam_md5": bam_base_out + ".bam.md5",
        "report_bamstats_html": report_base_out + ".bam.bamstats.html",
        "report_bamstats_html_md5": report_base_out + ".bam.bamstats.html.md5",
        "report_bamstats_txt": report_base_out + ".bam.bamstats.txt",
        "report_bamstats_txt_md5": report_base_out + ".bam.bamstats.txt.md5",
        "report_flagstats_txt": report_base_out + ".bam.flagstats.txt",
        "report_flagstats_txt_md5": report_base_out + ".bam.flagstats.txt.md5",
        "report_idxstats_txt": report_base_out + ".bam.idxstats.txt",
        "report_idxstats_txt_md5": report_base_out + ".bam.idxstats.txt.md5",
    }
    # Return
    return expected


# Test for project validation ----------------------------------------------------------------------


def test_extraction_type_check(
    ngs_mapping_workflow,
    germline_sheet_tsv,
    generic_rna_sheet_tsv,
    generic_mix_extraction_sheet_tsv,
):
    """Tests extraction type check method."""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(germline_sheet_tsv)
    dna_sheet = read_germline_tsv_sheet(germline_sheet_io)

    # Create rna sample sheet
    rna_sheet_io = io.StringIO(generic_rna_sheet_tsv)
    rna_sheet = read_generic_tsv_sheet(rna_sheet_io)

    # Create mix data sample sheet
    mix_sheet_io = io.StringIO(generic_mix_extraction_sheet_tsv)
    mix_sheet = read_generic_tsv_sheet(mix_sheet_io)

    # Evaluate if only DNA is True
    dna_bool, rna_bool = ngs_mapping_workflow.extraction_type_check(sample_sheet=dna_sheet)
    assert dna_bool, "Germline extraction type are set to DNA by default."
    assert not rna_bool, "No RNA sample was included in the sample sheet."

    # Evaluate if only RNA is True
    dna_bool, rna_bool = ngs_mapping_workflow.extraction_type_check(sample_sheet=rna_sheet)
    assert not dna_bool, "No DNA sample was included in the sample sheet."
    assert rna_bool, "Only RNA samples were included in the sample sheet."

    # Evaluate if both DNA and RNA are True
    dna_bool, rna_bool = ngs_mapping_workflow.extraction_type_check(sample_sheet=mix_sheet)
    assert dna_bool, "Sample sheet contains both DNA and RNA."
    assert rna_bool, "Sample sheet contains both DNA and RNA."


def test_project_validation_germline(
    ngs_mapping_workflow, germline_sheet_tsv, generic_rna_sheet_tsv, minimal_config
):
    """Tests project validation method in ngs mapping workflow"""
    # Convert yaml to dict
    minimal_config_dict = deepcopy(minimal_config)
    minimal_config_dict = dict(minimal_config_dict)
    minimal_config_dict = minimal_config_dict["step_config"].get("ngs_mapping", OrderedDict())

    # Create germline sample sheet
    germline_sheet_io = io.StringIO(germline_sheet_tsv)
    germline_sheet = read_germline_tsv_sheet(germline_sheet_io)

    # Create rna sample sheet
    rna_sheet_io = io.StringIO(generic_rna_sheet_tsv)
    rna_sheet = read_generic_tsv_sheet(rna_sheet_io)

    # Method returns None without exception, cause DNA sample sheet and DNA tool defined in config
    out = ngs_mapping_workflow.validate_project(
        config_dict=minimal_config_dict, sample_sheets_list=[germline_sheet]
    )
    assert out is None, "No exception expected: DNA sample sheet and DNA tool defined in config."

    # Exception raised cause no RNA mapper defined in config
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(
            config_dict=minimal_config_dict, sample_sheets_list=[rna_sheet]
        )
    error_msg = "RNA sample provided, but config only contains DNA mapper."
    assert exec_info.value.args[0] is not None, error_msg

    # Exception raised cause only DNA mapper defined in config
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(
            config_dict=minimal_config_dict, sample_sheets_list=[germline_sheet, rna_sheet]
        )
    error_msg = "DNA and RNA sample provided, but config only contains DNA mapper."
    assert exec_info.value.args[0] is not None, error_msg

    # Update config and remove RNA exception
    minimal_config_dict["tools"]["rna"] = ["rna_mapper"]
    out = ngs_mapping_workflow.validate_project(
        config_dict=minimal_config_dict, sample_sheets_list=[germline_sheet, rna_sheet]
    )
    error_msg = (
        "No exception expected: DNA, RNA sample sheet and respective tools defined in config."
    )
    assert out is None, error_msg

    # Update config and introduce DNA exception
    minimal_config_dict["tools"]["dna"] = []
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(
            config_dict=minimal_config_dict, sample_sheets_list=[germline_sheet, rna_sheet]
        )
    error_msg = "DNA and RNA sample provided, but config only contains RNA mapper."
    assert exec_info.value.args[0] is not None, error_msg

    # Exception raised cause no DNA mapper defined in config
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(
            config_dict=minimal_config_dict, sample_sheets_list=[germline_sheet]
        )
    error_msg = "DNA and RNA sample provided, but config only contains RNA mapper."
    assert exec_info.value.args[0] is not None, error_msg


# Tests for BwaStepPart ----------------------------------------------------------------------------


def test_bwa_step_part_get_args(ngs_mapping_workflow):
    """Tests BaseStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": {
            "reads_left": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz"],
            "reads_right": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-N1-DNA1-WGS1",
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
    expected = get_expected_output_files_dict(bam_base_out, report_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("bwa", "run")
    assert actual == expected


def test_bwa_step_part_get_log_file(ngs_mapping_workflow):
    """Tests BaseStepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/bwa.{library_name}/log/bwa.{library_name}"
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
        actual = ngs_mapping_workflow.get_resource("bwa", "run", resource)
        assert actual == expected, msg_error


# Tests for StarStepPart --------------------------------------------------------------------------


def test_star_step_part_get_args(ngs_mapping_workflow):
    """Tests StarStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": {
            "reads_left": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz"],
            "reads_right": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-N1-DNA1-WGS1",
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_args("star", "run")(wildcards)
    assert actual == expected


def test_star_step_part_get_input_files(ngs_mapping_workflow):
    """Tests StarStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = "work/input_links/P001-N1-DNA1-WGS1/.done"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("star", "run")(wildcards)
    assert actual == expected


def test_star_step_part_get_output_files(ngs_mapping_workflow):
    """Tests StarStepPart.get_output_files()"""
    # Define expected
    bam_base_out = "work/star.{library_name}/out/star.{library_name}"
    report_base_out = "work/star.{library_name}/report/bam_qc/star.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("star", "run")
    assert actual == expected


def test_star_step_part_get_log_file(ngs_mapping_workflow):
    """Tests StarStepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/star.{library_name}/log/star.{library_name}"
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
        actual = ngs_mapping_workflow.get_resource("star", "run", resource)
        assert actual == expected, msg_error


# Tests for Minimap2StepPart -----------------------------------------------------------------------


def test_minimap2_step_part_get_args(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": {
            "reads_left": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz"],
            "reads_right": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-N1-DNA1-WGS1",
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
    expected = get_expected_output_files_dict(bam_base_out, report_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("minimap2", "run")
    assert actual == expected


def test_minimap2_step_part_get_log_file(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/minimap2.{library_name}/log/minimap2.{library_name}"
    )
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("minimap2", "run")
    assert actual == expected


def test_minimap2_step_part_get_resource(ngs_mapping_workflow):
    """Tests Minimap2StepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "4-00:00:00", "memory": "56G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("minimap2", "run", resource)
        assert actual == expected, msg_error


# Tests for NgmlrStepPart -----------------------------------------------------------------------


def test_ngmlr_step_part_get_args(ngs_mapping_workflow):
    """Tests NgmlrStepPart.get_args()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": {
            "reads_left": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz"],
            "reads_right": ["work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz"],
        },
        "platform": "ILLUMINA",
        "sample_name": "P001-N1-DNA1-WGS1",
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_args("ngmlr", "run")(wildcards)
    assert actual == expected


def test_ngmlr_step_part_get_input_files(ngs_mapping_workflow):
    """Tests NgmlrStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = "work/input_links/P001-N1-DNA1-WGS1/.done"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("ngmlr", "run")(wildcards)
    assert actual == expected


def test_ngmlr_step_part_get_output_files(ngs_mapping_workflow):
    """Tests NgmlrStepPart.get_output_files()"""
    # Define expected
    bam_base_out = "work/ngmlr.{library_name}/out/ngmlr.{library_name}"
    report_base_out = "work/ngmlr.{library_name}/report/bam_qc/ngmlr.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out)
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("ngmlr", "run")
    assert actual == expected


def test_ngmlr_step_part_get_log_file(ngs_mapping_workflow):
    """Tests NgmlrStepPart.get_log_file()"""
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/ngmlr.{library_name}/log/ngmlr.{library_name}"
    )
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("ngmlr", "run")
    assert actual == expected


def test_ngmlr_step_part_get_resource(ngs_mapping_workflow):
    """Tests NgmlrStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "4-00:00:00", "memory": "50G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("ngmlr", "run", resource)
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
    expected = get_expected_output_files_dict(bam_base_out, report_base_out)
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
        actual = ngs_mapping_workflow.get_resource("external", "run", resource)
        assert actual == expected, msg_error


# Tests for GatkPostBamStepPart -------------------------------------------------------------------


def test_gatk_post_bam_step_part_get_input_files(ngs_mapping_workflow):
    """Tests GatkPostBamStepPart.get_input_files()"""
    expected = "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
    actual = ngs_mapping_workflow.get_input_files("gatk_post_bam", "run")
    assert actual == expected


def test_gatk_post_bam_step_part_get_output_files(ngs_mapping_workflow):
    """Tests GatkPostBamStepPart.get_output_files()"""
    # Define expected
    base_name_out = "work/{mapper}.{library_name}/out/{mapper}.{library_name}"
    expected = {
        "bai_md5_realigned": base_name_out + ".realigned.bam.bai.md5",
        "bai_md5_recalibrated": base_name_out + ".realigned.recalibrated.bam.bai.md5",
        "bai_realigned": base_name_out + ".realigned.bam.bai",
        "bai_recalibrated": base_name_out + ".realigned.recalibrated.bam.bai",
        "bam_md5_realigned": base_name_out + ".realigned.bam.md5",
        "bam_md5_recalibrated": base_name_out + ".realigned.recalibrated.bam.md5",
        "bam_realigned": base_name_out + ".realigned.bam",
        "bam_recalibrated": base_name_out + ".realigned.recalibrated.bam",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("gatk_post_bam", "run")
    assert actual == expected


def test_gatk_post_bam_step_part_get_log_file(ngs_mapping_workflow):
    """Tests GatkPostBamStepPart.get_log_file()"""
    expected = "work/{mapper}.{library_name}/log/snakemake.gatk_post_bam.log"
    actual = ngs_mapping_workflow.get_log_file("gatk_post_bam", "run")
    assert actual == expected


def test_gatk_post_bam_step_part_get_resource(ngs_mapping_workflow):
    """Tests GatkPostBamStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "2-00:00:00", "memory": "16G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("gatk_post_bam", "run", resource)
        assert actual == expected, msg_error


# Tests for LinkOutBamStepPart --------------------------------------------------------------------


def test_link_out_bam_step_part_get_input_files(ngs_mapping_workflow):
    """Tests LinkOutBamStepPart.get_input_files()"""
    # Define expected
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "library"})
    expected = [
        "work/bwa.library/out/bwa.library.bam",
        "work/bwa.library/out/bwa.library.bam.bai",
        "work/bwa.library/out/bwa.library.bam.md5",
        "work/bwa.library/out/bwa.library.bam.bai.md5",
    ]
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("link_out_bam", "run")(wildcards)
    assert actual == expected


def test_link_out_bam_step_part_get_output_files(ngs_mapping_workflow):
    """Tests LinkOutBamStepPart.get_output_files()"""
    # Define expected
    expected = [
        "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam",
        "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam.bai",
        "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam.md5",
        "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam.bai.md5",
    ]
    # Get actual and assert
    actual = ngs_mapping_workflow.get_output_files("link_out_bam", "run")
    assert actual == expected


def test_link_out_bam_step_part_get_shell_cmd(ngs_mapping_workflow):
    """Tests LinkOutBamStepPart.get_shell_cmd()"""
    # Define expected
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "library"})
    b_name = "bwa.library/out/bwa.library.bam"
    expected = textwrap.dedent(
        rf"""
       test -L output/{b_name} || ln -sr work/{b_name} output/{b_name}
       test -L output/{b_name}.bai || ln -sr work/{b_name}.bai output/{b_name}.bai
       test -L output/{b_name}.md5 || ln -sr work/{b_name}.md5 output/{b_name}.md5
       test -L output/{b_name}.bai.md5 || ln -sr work/{b_name}.bai.md5 output/{b_name}.bai.md5
       """
    ).strip()
    # Get actual and assert
    actual = ngs_mapping_workflow.get_shell_cmd("link_out_bam", "run", wildcards)
    assert actual == expected


# Tests for PicardHsMetricsStepPart ----------------------------------------------------------------


def test_picard_hs_metrics_step_part_get_input_files(ngs_mapping_workflow):
    """Tests PicardHsMetricsStepPart.get_input_files()"""
    # Define expected
    expected = {
        "bam": "work/{mapper_lib}/out/{mapper_lib}.bam",
        "bai": "work/{mapper_lib}/out/{mapper_lib}.bam.bai",
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("picard_hs_metrics", "run")
    assert actual == expected


def test_picard_hs_metrics_step_part_get_output_files(ngs_mapping_workflow):
    """Tests PicardHsMetricsStepPart.get_output_files()"""
    # Define expected
    expected = {
        "txt": "work/{mapper_lib}/report/picard_hs_metrics/{mapper_lib}.txt",
        "txt_md5": "work/{mapper_lib}/report/picard_hs_metrics/{mapper_lib}.txt.md5",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_output_files("picard_hs_metrics", "run")
    assert actual == expected


def test_picard_hs_metrics_step_part_get_log_file(ngs_mapping_workflow):
    """Tests PicardHsMetricsStepPart.get_log_file()"""
    # Define expected
    expected = "work/{mapper_lib}/log/snakemake.picard_hs_metrics.log"
    # Get actual
    actual = ngs_mapping_workflow.get_log_file("picard_hs_metrics", "run")
    assert actual == expected


def test_picard_hs_metrics_step_part_get_resource(ngs_mapping_workflow):
    """Tests PicardHsMetricsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "04:00:00", "memory": "20G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("picard_hs_metrics", "run", resource)
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
    wildcards = Wildcards(fromdict={"mapper_lib": "bwa.library"})
    actual = ngs_mapping_workflow.get_input_files("target_coverage_report", "run")(wildcards)
    assert actual == expected


def test_target_coverage_report_step_part_collect_get_input_files(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_input_files() - collect"""
    # Define expected
    expected = [
        "work/bwa.P001-N1-DNA1-WGS1/report/cov_qc/bwa.P001-N1-DNA1-WGS1.txt",
        "work/bwa.P002-N1-DNA1-WGS1/report/cov_qc/bwa.P002-N1-DNA1-WGS1.txt",
        "work/bwa.P003-N1-DNA1-WGS1/report/cov_qc/bwa.P003-N1-DNA1-WGS1.txt",
        "work/bwa.P004-N1-DNA1-WGS1/report/cov_qc/bwa.P004-N1-DNA1-WGS1.txt",
        "work/bwa.P005-N1-DNA1-WGS1/report/cov_qc/bwa.P005-N1-DNA1-WGS1.txt",
        "work/bwa.P006-N1-DNA1-WGS1/report/cov_qc/bwa.P006-N1-DNA1-WGS1.txt",
    ]
    # Get actual
    wildcards = Wildcards(fromdict={"mapper_lib": "bwa"})
    actual = ngs_mapping_workflow.get_input_files("target_coverage_report", "collect")(wildcards)
    assert actual == expected


def test_target_coverage_report_step_part_get_output_files(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_output_files() - run"""
    expected = {
        "txt": "work/{mapper_lib}/report/cov_qc/{mapper_lib}.txt",
        "txt_md5": "work/{mapper_lib}/report/cov_qc/{mapper_lib}.txt.md5",
    }
    assert ngs_mapping_workflow.get_output_files("target_coverage_report", "run") == expected


def test_target_coverage_report_step_part_run_get_log_file(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_log_file() - run"""
    expected = "work/{mapper_lib}/log/snakemake.target_coverage.log"
    assert ngs_mapping_workflow.get_log_file("target_coverage_report", "run") == expected


def test_target_coverage_report_step_part_collect_get_log_file(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_log_file() - collect"""
    expected = "work/target_cov_report/log/snakemake.target_coverage.log"
    assert ngs_mapping_workflow.get_log_file("target_coverage_report", "collect") == expected


def test_target_coverage_report_step_part_get_resource(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_resource()"""
    # Define expected
    actions = ("run", "collect")
    expected_dict = {"threads": 2, "time": "04:00:00", "memory": "20G", "partition": "medium"}
    # Evaluate
    for action in actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = ngs_mapping_workflow.get_resource("target_coverage_report", action, resource)
            assert actual == expected, msg_error


# Tests for GenomeCoverageReportStepPart ----------------------------------------------------------


def test_genome_coverage_report_step_part_get_input_files(ngs_mapping_workflow):
    """Tests GenomeCoverageReportStepPart.get_input_files()"""
    # Define expected
    expected = {
        "bam": "work/{mapper_lib}/out/{mapper_lib}.bam",
        "bai": "work/{mapper_lib}/out/{mapper_lib}.bam.bai",
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_input_files("genome_coverage_report", "run")
    assert actual == expected


def test_genome_coverage_report_step_part_get_output_files(ngs_mapping_workflow):
    """Tests GenomeCoverageReportStepPart.get_output_files()"""
    # Define expected
    expected = {
        "bed": "work/{mapper_lib}/report/coverage/{mapper_lib}.bed.gz",
        "tbi": "work/{mapper_lib}/report/coverage/{mapper_lib}.bed.gz.tbi",
    }
    # Get actual and assert
    actual = ngs_mapping_workflow.get_output_files("genome_coverage_report", "run")
    assert actual == expected


def test_genome_coverage_report_step_part_get_log_file(ngs_mapping_workflow):
    """Tests GenomeCoverageReportStepPart.get_log_file()"""
    # Define expected
    expected = "work/{mapper_lib}/log/snakemake.genome_coverage.log"
    # Get actual and assert
    actual = ngs_mapping_workflow.get_log_file("genome_coverage_report", "run")
    assert actual == expected


def test_genome_coverage_report_step_part_get_shell_cmd(ngs_mapping_workflow):
    """Tests GenomeCoverageReportStepPart.get_shell_cmd()"""
    cmd = ngs_mapping_workflow.get_shell_cmd("genome_coverage_report", "run", None)
    # The shell command returns a static blob and thus we only check for some important keywords.
    # The actual functionality has to be tested in an integration/system test.
    assert "samtools depth" in cmd
    assert "bgzip -c" in cmd
    assert "tabix {output.bed}" in cmd


def test_genome_coverage_report_step_part_get_resource(ngs_mapping_workflow):
    """Tests GenomeCoverageReportStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "04:00:00", "memory": "4G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = ngs_mapping_workflow.get_resource("genome_coverage_report", "run", resource)
        assert actual == expected, msg_error


# Tests for NgsMappingWorkflow --------------------------------------------------------------------


def test_ngs_mapping_workflow_steps(ngs_mapping_workflow):
    """Tests simple functionality of the workflow: checks if sub steps are created,
    i.e., the tools associated with gene expression quantification."""
    # Check created sub steps
    expected = [
        "bwa",
        "external",
        "gatk_post_bam",
        "genome_coverage_report",
        "link_in",
        "link_out",
        "link_out_bam",
        "minimap2",
        "ngmlr",
        "picard_hs_metrics",
        "star",
        "target_coverage_report",
    ]
    actual = list(sorted(ngs_mapping_workflow.sub_steps.keys()))
    assert actual == expected


def test_ngs_mapping_workflow_files(ngs_mapping_workflow):
    """Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools, namely: bwa, external, gatk_post_bam,
    genome_coverage_report, link_in, link_out, link_out_bam, minimap2, ngmlr, picard_hs_metrics,
    star, target_coverage_report."""
    # Check result file construction
    expected = [
        "output/bwa.P00{i}-N1-DNA1-WGS1/out/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(i=i, ext=ext)
        for i in range(1, 7)
        for ext in ("bam", "bam.bai", "bam.md5", "bam.bai.md5")
    ]
    expected += [
        "output/bwa.P00{i}-N1-DNA1-WGS1/log/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(i=i, ext=ext)
        for i in range(1, 7)
        for ext in (
            "log",
            "conda_info.txt",
            "conda_list.txt",
            "log.md5",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
        )
    ]
    bam_stats_text_out = (
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/bam_qc/bwa.P00{i}-N1-DNA1-WGS1.bam.{stats}.{ext}"
    )
    expected += [
        bam_stats_text_out.format(i=i, stats=stats, ext=ext)
        for ext in ("txt", "txt.md5")
        for i in range(1, 7)
        for stats in ("bamstats", "flagstats", "idxstats")
    ]
    bam_stats_html_out = (
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/bam_qc/bwa.P00{i}-N1-DNA1-WGS1.bam.bamstats.{ext}"
    )
    expected += [
        bam_stats_html_out.format(i=i, ext=ext) for ext in ("html", "html.md5") for i in range(1, 7)
    ]
    expected += [
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/cov_qc/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(
            i=i, ext=ext
        )
        for ext in ("txt", "txt.md5")
        for i in range(1, 7)
    ]
    expected += [
        "output/target_cov_report/out/target_cov_report.txt",
        "output/target_cov_report/out/target_cov_report.txt.md5",
    ]
    expected += [
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/coverage/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(
            i=i, ext=ext
        )
        for i in range(1, 7)
        for ext in ("bed.gz", "bed.gz.tbi")
    ]
    assert sorted(ngs_mapping_workflow.get_result_files()) == sorted(expected)
