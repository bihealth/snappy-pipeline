# -*- coding: utf-8 -*-
"""Tests for the ngs_mapping workflow module code"""

import io
import textwrap
from collections import OrderedDict
from copy import deepcopy

import pytest
import ruamel.yaml as ruamel_yaml
from biomedsheets.io_tsv import read_generic_tsv_sheet, read_germline_tsv_sheet
from snakemake.io import Wildcards

from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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
            target_coverage_report:
              path_target_interval_list_mapping:
              - pattern: "Agilent SureSelect Human All Exon V6.*"
                name: Agilent_SureSelect_Human_All_Exon_V6
                path: path/to/SureSelect_Human_All_Exon_V6_r2.bed
            bwa:
              path_index: /path/to/bwa/index.fasta.amb
            bwa_mem2:
              path_index: /path/to/bwa_mem2/index.fasta.amb
            minimap2:
              mapping_threads: 16
            star:
              path_index: /path/to/star/index
              transcriptome: false
              out_filter_intron_motifs: ""
              out_sam_strand_field: ""
            mbcs:
              mapping_tool: bwa
              use_barcodes: True
              recalibrate: True
            bqsr:
              common_variants: /path/to/common/variants
            agent:
              prepare:
                path: /path/to/trimmer
                lib_prep_type: v2
              mark_duplicates:
                path: /path/to/creak
                path_baits: /path/to/baits
                consensus_mode: HYBRID
            bam_collect_doc:
              enabled: true

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

    # Define expected
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
    config = ngs_mapping_workflow.config_model_class(**minimal_config_dict)

    # Create germline sample sheet
    germline_sheet_io = io.StringIO(germline_sheet_tsv)
    germline_sheet = read_germline_tsv_sheet(germline_sheet_io)

    # Create rna sample sheet
    rna_sheet_io = io.StringIO(generic_rna_sheet_tsv)
    rna_sheet = read_generic_tsv_sheet(rna_sheet_io)

    # Method returns None without exception, cause DNA sample sheet and DNA tool defined in config
    out = ngs_mapping_workflow.validate_project(config=config, sample_sheets_list=[germline_sheet])
    assert out is None, "No exception expected: DNA sample sheet and DNA tool defined in config."

    # Exception raised cause no RNA mapper defined in config
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(config=config, sample_sheets_list=[rna_sheet])
    error_msg = "RNA sample provided, but config only contains DNA mapper."
    assert exec_info.value.args[0] is not None, error_msg

    # Exception raised cause only DNA mapper defined in config
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(
            config=config, sample_sheets_list=[germline_sheet, rna_sheet]
        )
    error_msg = "DNA and RNA sample provided, but config only contains DNA mapper."
    assert exec_info.value.args[0] is not None, error_msg

    # Update config and remove RNA exception
    config.tools.rna = ["rna_mapper"]
    out = ngs_mapping_workflow.validate_project(
        config=config, sample_sheets_list=[germline_sheet, rna_sheet]
    )
    error_msg = (
        "No exception expected: DNA, RNA sample sheet and respective tools defined in config."
    )
    assert out is None, error_msg

    # Update config and introduce DNA exception
    config.tools.dna = []
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(
            config=config, sample_sheets_list=[germline_sheet, rna_sheet]
        )
    error_msg = "DNA and RNA sample provided, but config only contains RNA mapper."
    assert exec_info.value.args[0] is not None, error_msg

    # Exception raised cause no DNA mapper defined in config
    with pytest.raises(Exception) as exec_info:
        ngs_mapping_workflow.validate_project(config=config, sample_sheets_list=[germline_sheet])
    error_msg = "DNA and RNA sample provided, but config only contains RNA mapper."
    assert exec_info.value.args[0] is not None, error_msg


# Tests for BwaStepPart, BwaMem2StepPart & MBCsStepPart -----------------------


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
    expected_dict = {
        "bwa": {"threads": 16, "time": "3-00:00:00", "memory": "73728M", "partition": "medium"},
        "bwa_mem2": {
            "threads": 16,
            "time": "3-00:00:00",
            "memory": "73728M",
            "partition": "medium",
        },
        "mbcs": {"threads": 1, "time": "72:00:00", "memory": "4G", "partition": "medium"},
    }
    # Evaluate
    for tool, v in expected_dict.items():
        for resource, expected in v.items():
            msg_error = f"Assertion error for tool '{tool}' & resource '{resource}'."
            actual = ngs_mapping_workflow.get_resource(tool, "run", resource)()
            assert actual == expected, msg_error


def test_bwa_step_part_check_config(ngs_mapping_workflow):
    """Tests BaseStepPart.check_config()"""
    # Define expected
    for tool in ("bwa", "bwa_mem2", "mbcs"):
        ngs_mapping_workflow.sub_steps[tool].check_config()


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
    log_base_out = "work/star.{library_name}/log/star.{library_name}"
    expected = get_expected_output_files_dict(bam_base_out, report_base_out, log_base_out)
    expected["gene_counts"] = "work/star.{library_name}/out/star.{library_name}.GeneCounts.tab"
    expected["gene_counts_md5"] = expected["gene_counts"] + ".md5"
    expected["junctions"] = "work/star.{library_name}/out/star.{library_name}.Junctions.tab"
    expected["junctions_md5"] = expected["junctions"] + ".md5"
    expected["output_links"].extend(
        [
            "output/star.{library_name}/out/star.{library_name}.Junctions.tab",
            "output/star.{library_name}/out/star.{library_name}.Junctions.tab.md5",
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
        "json": "work/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz",
        "json_md5": "work/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz.md5",
    }
    assert ngs_mapping_workflow.get_output_files("target_coverage_report", "run") == expected


def test_target_coverage_report_step_part_run_get_log_file(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_log_file() - run"""
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


def test_target_coverage_report_step_part_run_get_params(ngs_mapping_workflow):
    """Tests TargetCoverageReportStepPart.get_params() - action 'run'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "path_targets_bed": "path/to/SureSelect_Human_All_Exon_V6_r2.bed",
    }
    actual = ngs_mapping_workflow.get_params("target_coverage_report", "run")(wildcards)
    assert actual == expected


# Tests for BamCollectDocStepPart -----------------------------------------------------------------


def test_generate_doc_files_step_part_run_get_input_files(ngs_mapping_workflow):
    """Tests BamCollectDocStepPart.get_input_files() - run"""
    # Define expected
    expected = {
        "bam": "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam",
        "bai": "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam.bai",
    }
    # Get actual
    actual = ngs_mapping_workflow.get_input_files("bam_collect_doc", "run")()
    assert actual == expected


def test_generate_doc_files_step_part_get_output_files(ngs_mapping_workflow):
    """Tests BamCollectDocStepPart.get_output_files() - run"""
    expected = {
        "cov_bw": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.bw",
        "cov_bw_md5": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.bw.md5",
        "mq_bw": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.mq.bw",
        "mq_bw_md5": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.mq.bw.md5",
        "vcf": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz",
        "vcf_md5": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.md5",
        "vcf_tbi": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.tbi",
        "vcf_tbi_md5": "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.tbi.md5",
        "output_links": [
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.md5",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.tbi",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.tbi.md5",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.bw",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.bw.md5",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.mq.bw",
            "output/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.mq.bw.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.log",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.log.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_info.txt",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_info.txt.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_list.txt",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_list.txt.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.wrapper.py",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.wrapper.py.md5",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.environment.yaml",
            "output/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.environment.yaml.md5",
        ],
    }
    assert ngs_mapping_workflow.get_output_files("bam_collect_doc", "run") == expected


def test_generate_doc_files_step_part_run_get_log_file(ngs_mapping_workflow):
    """Tests BamCollectDocStepPart.get_log_file() - run"""
    expected = {
        "log": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.log",
        "log_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.log.md5",
        "conda_info": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_info.txt",
        "conda_info_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_info.txt.md5",
        "conda_list": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_list.txt",
        "conda_list_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.conda_list.txt.md5",
        "wrapper": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.wrapper.py",
        "wrapper_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.wrapper.py.md5",
        "env_yaml": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.environment.yaml",
        "env_yaml_md5": "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc.environment.yaml.md5",
    }
    assert ngs_mapping_workflow.get_log_file("bam_collect_doc", "run") == expected


def test_generate_doc_files_step_part_get_resource(ngs_mapping_workflow):
    """Tests BamCollectDocStepPart.get_resource()"""
    expected_dict = {"threads": 1, "time": "24:00:00", "memory": "2G", "partition": "medium"}
    for resource, expected in expected_dict.items():
        actual = ngs_mapping_workflow.get_resource("bam_collect_doc", "run", resource)()
        assert actual == expected


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
    to the expected results from the tools, namely: bwa, external, link_in, link_out,
    link_out_bam, minimap2, star, target_coverage_report.
    """
    # Check result file construction
    expected = [
        "output/bwa.P00{i}-N1-DNA1-WGS1/out/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(i=i, ext=ext)
        for i in range(1, 7)
        for ext in ("bam", "bam.bai", "bam.md5", "bam.bai.md5")
    ]
    for infix in ("bam_collect_doc", "mapping", "target_cov_report", "ngs_chew_fingerprint"):
        expected += [
            "output/bwa.P00{i}-N1-DNA1-WGS1/log/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(i=i, ext=ext)
            for i in range(1, 7)
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
    expected += [
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/cov/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(
            i=i, ext=ext
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
        for i in range(1, 7)
    ]
    expected += [
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/fingerprint/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(
            i=i, ext=ext
        )
        for ext in ("npz", "npz.md5")
        for i in range(1, 7)
    ]
    expected += [
        "output/bwa.P00{i}-N1-DNA1-WGS1/report/alfred_qc/bwa.P00{i}-N1-DNA1-WGS1.{ext}".format(
            i=i, ext=ext
        )
        for ext in ("alfred.json.gz", "alfred.json.gz.md5")
        for i in range(1, 7)
    ]
    assert sorted(ngs_mapping_workflow.get_result_files()) == sorted(expected)
