# -*- coding: utf-8 -*-
"""Tests for the variant_export_external workflow module code"""
from copy import deepcopy
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.workflows.variant_export_external import VariantExportExternalWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
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
          variant_export_external:
            bam_available_flag: true
            merge_vcf_flag: true
            search_paths: [/search_path]
            search_patterns: [{"vcf": "*.vcf.gz", "bam": "*.bam", "bai": "*.bai"}]
            external_tool: dragen
            path_refseq_ser: /data/refseq_ser
            path_ensembl_ser: /data/ensembl_ser
            path_db: /data/db
            target_coverage_report:
              path_targets_bed: /path/to/targets.bed

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
def variant_export_external_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return VariantExportExternalWorkflow object pre-configured with germline sheet"""
    # Create annotation files
    for fk_file in ("refseq_ser", "ensembl_ser", "db"):
        germline_sheet_fake_fs.fs.create_file(
            "/data/" + fk_file,
            contents="",
            create_missing_dirs=True,
        )
    # Create search path
    germline_sheet_fake_fs.fs.makedirs("/search_path")
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_export_external", germline_sheet_fake_fs, mocker
    )
    # Construct the workflow object
    return VariantExportExternalWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_workflow_check_config_invalid_annotator_files(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Tests VariantExportExternalWorkflow.check_config() - invalid varfish-annotator files"""
    # Create search path
    germline_sheet_fake_fs.fs.makedirs("/search_path")
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_export_external", germline_sheet_fake_fs, mocker
    )
    # Construct the workflow object
    with pytest.raises(MissingConfiguration) as exec_info:
        VariantExportExternalWorkflow(
            dummy_workflow,
            minimal_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )
    assert "path_refseq_ser" in exec_info.value.args[0]
    assert "path_ensembl_ser" in exec_info.value.args[0]
    assert "path_db" in exec_info.value.args[0]


def test_workflow_check_config_invalid_search_directory(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Tests VariantExportExternalWorkflow.check_config() - no search directory"""
    # Create annotation files
    for fk_file in ("refseq_ser", "ensembl_ser", "db"):
        germline_sheet_fake_fs.fs.create_file(
            "/data/" + fk_file,
            contents="",
            create_missing_dirs=True,
        )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_export_external", germline_sheet_fake_fs, mocker
    )
    # Construct the workflow object
    with pytest.raises(MissingConfiguration) as exec_info:
        VariantExportExternalWorkflow(
            dummy_workflow,
            minimal_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )
    assert " is not a directory: /search_path" in exec_info.value.args[0]


def test_workflow_check_config_invalid_search_pattern(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return VariantExportExternalWorkflow object pre-configured with germline sheet"""
    # Create annotation files
    for fk_file in ("refseq_ser", "ensembl_ser", "db"):
        germline_sheet_fake_fs.fs.create_file(
            "/data/" + fk_file,
            contents="",
            create_missing_dirs=True,
        )
    # Create search path
    germline_sheet_fake_fs.fs.makedirs("/search_path")
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_export_external", germline_sheet_fake_fs, mocker
    )
    # Change search patterns to invalid
    modified_config = deepcopy(minimal_config)
    modified_config["step_config"]["variant_export_external"]["search_patterns"] = [
        "vcf",
        "*/*.vcf.gz",
    ]
    # Construct the workflow object
    with pytest.raises(MissingConfiguration) as exec_info:
        VariantExportExternalWorkflow(
            dummy_workflow,
            modified_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )
    assert "Value in 'search_patterns' is not a dictionary" in exec_info.value.args[0]


# Tests for BamReportsExternalStepPart (bam_qc)   --------------------------------------------------


def test_bam_reports_step_part_call_get_input_files_bam_qc(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_input_files_bam_qc()"""
    wildcards = Wildcards(fromdict={"mapper_lib": "dragen.P001-N1-DNA1-WGS1"})
    expected = [
        "work/input_links/P001-N1-DNA1-WGS1/.done_bai_external",
        "work/input_links/P001-N1-DNA1-WGS1/.done_bam_external",
    ]
    actual = sorted(
        list(variant_export_external_workflow.get_input_files("bam_reports", "bam_qc")(wildcards))
    )
    assert actual == expected


def test_bam_reports_step_part_call_get_output_files_bam_qc(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_output_files_bam_qc()"""
    expected = {
        "bamstats": "work/{mapper_lib}/report/bam_qc/{mapper_lib}.bam.bamstats.txt",
        "bamstats_md5": "work/{mapper_lib}/report/bam_qc/{mapper_lib}.bam.bamstats.txt.md5",
        "flagstats": "work/{mapper_lib}/report/bam_qc/{mapper_lib}.bam.flagstats.txt",
        "flagstats_md5": "work/{mapper_lib}/report/bam_qc/{mapper_lib}.bam.flagstats.txt.md5",
        "idxstats": "work/{mapper_lib}/report/bam_qc/{mapper_lib}.bam.idxstats.txt",
        "idxstats_md5": "work/{mapper_lib}/report/bam_qc/{mapper_lib}.bam.idxstats.txt.md5",
    }
    actual = variant_export_external_workflow.get_output_files("bam_reports", "bam_qc")
    assert actual == expected


def test_bam_reports_step_part_call_get_log_file_bam_qc(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_log_file_bam_qc()"""
    base_out = "work/{mapper_lib}/log/{mapper_lib}.bam_qc"
    expected = get_expected_log_files_dict(base_out)
    actual = variant_export_external_workflow.get_log_file("bam_reports", "bam_qc")
    assert actual == expected


def test_bam_reports_step_part_call_get_params_bam_qc(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_params_bam_qc()"""
    wildcards = Wildcards(fromdict={"mapper_lib": "dragen.P001-N1-DNA1-WGS1"})
    expected = {"bam": [], "bam_count": 0}
    actual = variant_export_external_workflow.get_params("bam_reports", "bam_qc")(wildcards)
    assert actual == expected


# Tests for BamReportsExternalStepPart (collect)   -------------------------------------------------


def test_bam_reports_step_part_call_get_input_files_collect(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_input_files_collect()"""
    expected = ["work/{mapper_lib}/report/cov_qc/{mapper_lib}.txt"]
    actual = variant_export_external_workflow.get_input_files("bam_reports", "collect")(None)
    assert actual == expected


# Tests for BamReportsExternalStepPart (run)    ----------------------------------------------------


def test_bam_reports_step_part_call_get_input_files_run(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_input_files_run()"""
    wildcards = Wildcards(fromdict={"mapper_lib": "dragen.P001-N1-DNA1-WGS1"})
    expected = [
        "work/input_links/P001-N1-DNA1-WGS1/.done_bai_external",
        "work/input_links/P001-N1-DNA1-WGS1/.done_bam_external",
    ]
    actual = sorted(
        list(variant_export_external_workflow.get_input_files("bam_reports", "run")(wildcards))
    )
    assert actual == expected


def test_bam_reports_step_part_call_get_params_run(
    variant_export_external_workflow,
):
    """Tests BamReportsExternalStepPart._get_params_run()"""
    wildcards = Wildcards(fromdict={"mapper_lib": "dragen.P001-N1-DNA1-WGS1"})
    expected = {
        "bam": [],
        "bam_count": 0,
        "path_targets_bed": "/path/to/targets.bed",
        "max_coverage": 200,
        "min_cov_warning": 20,
        "min_cov_ok": 50,
        "detailed_reporting": False,
    }
    actual = variant_export_external_workflow.get_params("bam_reports", "run")(wildcards)
    assert actual == expected


# Tests for VarfishAnnotatorExternalStepPart (gvcf_to_vcf)   -----------------------------------------


def test_varfish_annotator_step_part_call_get_input_files_gvcf_to_vcf(
    variant_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_input_files_gvcf_to_vcf()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = ["work/input_links/P001-N1-DNA1-WGS1/.done"]
    actual = variant_export_external_workflow.get_input_files(
        "varfish_annotator_external", "gvcf_to_vcf"
    )(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_call_get_output_files_gvcf_to_vcf(
    variant_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_output_files_gvcf_to_vcf()"""
    base_name = "work/dragen.{index_ngs_library}/out/dragen.{index_ngs_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = variant_export_external_workflow.get_output_files(
        "varfish_annotator_external", "gvcf_to_vcf"
    )
    assert actual == expected


def test_varfish_annotator_step_part_call_get_log_file_gvcf_to_vcf(
    variant_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_log_file_gvcf_to_vcf()"""
    base_name = "work/dragen.{index_ngs_library}/log/dragen.{index_ngs_library}.gvcf_to_vcf"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = variant_export_external_workflow.get_log_file(
        "varfish_annotator_external", "gvcf_to_vcf"
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage_gvcf_to_vcf(
    variant_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart.get_resource_usage() - action 'gvcf_to_vcf'"""
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "14336M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'gvcf_to_vcf'."
        actual = variant_export_external_workflow.get_resource(
            "varfish_annotator_external", "gvcf_to_vcf", resource
        )
        assert actual == expected, msg_error


def test_varfish_annotator_step_part_get_params_gvcf_to_vcf(variant_export_external_workflow):
    """Tests TargetCoverageReportStepPart._get_params_gvcf_to_vcf()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": [],
        "sample_names": ["P001", "P002", "P003"],
    }
    actual = variant_export_external_workflow.get_params(
        "varfish_annotator_external", "gvcf_to_vcf"
    )(wildcards)
    assert actual == expected


# Tests for VarfishAnnotatorExternalStepPart (merge_vcf)   -----------------------------------------


def test_varfish_annotator_step_part_call_get_input_files_merge_vcf(
    variant_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_input_files_merge_vcf()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = [
        "work/input_links/P001-N1-DNA1-WGS1/.done",
        "work/input_links/P002-N1-DNA1-WGS1/.done",
        "work/input_links/P003-N1-DNA1-WGS1/.done",
    ]
    actual = variant_export_external_workflow.get_input_files(
        "varfish_annotator_external", "merge_vcf"
    )(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_call_get_output_files_merge_vcf(
    variant_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_output_files_merge_vcf()"""
    base_name = "work/dragen.{index_ngs_library}/out/dragen.{index_ngs_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = variant_export_external_workflow.get_output_files(
        "varfish_annotator_external", "merge_vcf"
    )
    assert actual == expected


def test_varfish_annotator_step_part_call_get_log_file_merge_vcf(variant_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart._get_log_file_merge_vcf()"""
    base_name = "work/dragen.{index_ngs_library}/log/dragen.{index_ngs_library}.merge_vcf"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = variant_export_external_workflow.get_log_file(
        "varfish_annotator_external", "merge_vcf"
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage_merge_vcf(variant_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart.get_resource_usage() - action 'merge_vcf'"""
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "14336M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'merge_vcf'."
        actual = variant_export_external_workflow.get_resource(
            "varfish_annotator_external", "merge_vcf", resource
        )
        assert actual == expected, msg_error


def test_varfish_annotator_step_part_get_params_merge_vcf(variant_export_external_workflow):
    """Tests TargetCoverageReportStepPart.get_params() - action 'merge_vcf'"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": [],
        "sample_names": ["P001", "P002", "P003"],
        "merge_option": None,
    }
    actual = variant_export_external_workflow.get_params("varfish_annotator_external", "merge_vcf")(
        wildcards
    )
    assert actual == expected


# Tests for VarfishAnnotatorExternalStepPart (annotate)   ------------------------------------------


def test_varfish_annotator_step_part_get_input_files_annotate(variant_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart._get_input_files_annotate()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    base_name = "work/dragen.P001-N1-DNA1-WGS1/out/dragen.P001-N1-DNA1-WGS1"
    expected = {
        "ped": "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped",
        "vcf": base_name + ".vcf.gz",
        "vcf_md5": base_name + ".vcf.gz.md5",
        "tbi": base_name + ".vcf.gz.tbi",
        "tbi_md5": base_name + ".vcf.gz.tbi.md5",
    }
    actual = variant_export_external_workflow.get_input_files(
        "varfish_annotator_external", "annotate"
    )(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_get_output_files_annotate(variant_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_output_files_annotate()"""
    # Define expected
    base_name_out = (
        "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
    )
    expected = {
        "gts": base_name_out + ".gts.tsv.gz",
        "gts_md5": base_name_out + ".gts.tsv.gz.md5",
        "db_infos": base_name_out + ".db-infos.tsv.gz",
        "db_infos_md5": base_name_out + ".db-infos.tsv.gz.md5",
    }
    # Get actual
    actual = variant_export_external_workflow.get_output_files(
        "varfish_annotator_external", "annotate"
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_log_file_annotate(variant_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_log_file()_annotate"""
    # Define expected
    base_name_log = (
        "work/varfish_annotated.{index_ngs_library}/log/"
        "varfish_annotated.annotate.{index_ngs_library}"
    )
    base_name_wrapper = (
        "work/varfish_annotated.{index_ngs_library}/log/"
        "varfish_annotated.annotate.{index_ngs_library}"
    )
    wrapper_dict = {
        "wrapper": base_name_wrapper + ".wrapper.py",
        "wrapper_md5": base_name_wrapper + ".wrapper.py.md5",
    }
    log_dict = get_expected_log_files_dict(base_out=base_name_log)
    expected = {**wrapper_dict, **log_dict}
    # Get actual
    actual = variant_export_external_workflow.get_log_file("varfish_annotator_external", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_params_annotate(variant_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_params_annotate()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {"is_wgs": True, "step_name": "variant_export_external"}
    actual = variant_export_external_workflow.get_params("varfish_annotator_external", "annotate")(
        wildcards
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage_annotate(variant_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart.get_resource_usage() - action 'annotate'"""
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "14336M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'annotate'."
        actual = variant_export_external_workflow.get_resource(
            "varfish_annotator_external", "annotate", resource
        )
        assert actual == expected, msg_error


# Tests for VarfishAnnotatorExternalStepPart (bam_qc)   --------------------------------------------


def test_varfish_annotator_step_part_get_input_files_bam_qc(variant_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart._get_input_files_bam_qc()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    # Define expected
    donor_indices = (1, 2, 3)
    base_name_bam = "work/dragen.P00{i}-N1-DNA1-WGS1/report/bam_qc/dragen.P00{i}-N1-DNA1-WGS1.{ext}"
    base_name_cov = "work/dragen.P00{i}-N1-DNA1-WGS1/report/cov_qc/dragen.P00{i}-N1-DNA1-WGS1.txt"
    expected = {
        "bamstats": [base_name_bam.format(i=i, ext="bam.bamstats.txt") for i in donor_indices],
        "flagstats": [base_name_bam.format(i=i, ext="bam.flagstats.txt") for i in donor_indices],
        "idxstats": [base_name_bam.format(i=i, ext="bam.idxstats.txt") for i in donor_indices],
        "cov_qc": [base_name_cov.format(i=i) for i in donor_indices],
    }
    # Get actual
    actual = variant_export_external_workflow.get_input_files(
        "varfish_annotator_external", "bam_qc"
    )(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_get_output_files_bam_qc(variant_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_output_files_bam_qc()"""
    # Define expected
    base_name_out = (
        "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
    )
    expected = {
        "bam_qc": base_name_out + ".bam-qc.tsv.gz",
        "bam_qc_md5": base_name_out + ".bam-qc.tsv.gz.md5",
    }
    # Get actual
    actual = variant_export_external_workflow.get_output_files(
        "varfish_annotator_external", "bam_qc"
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_log_file_bam_qc(variant_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_log_file()_bam_qc"""
    # Define expected
    base_name_log = (
        "work/varfish_annotated.{index_ngs_library}/log/"
        "varfish_annotated.bam_qc.{index_ngs_library}"
    )
    base_name_wrapper = (
        "work/varfish_annotated.{index_ngs_library}/log/"
        "varfish_annotated.bam_qc.{index_ngs_library}"
    )
    wrapper_dict = {
        "wrapper": base_name_wrapper + ".wrapper.py",
        "wrapper_md5": base_name_wrapper + ".wrapper.py.md5",
    }
    log_dict = get_expected_log_files_dict(base_out=base_name_log)
    expected = {**wrapper_dict, **log_dict}
    # Get actual
    actual = variant_export_external_workflow.get_log_file("varfish_annotator_external", "bam_qc")
    assert actual == expected


def test_varfish_annotator_step_part_get_params_bam_qc(variant_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_params_bam_qc()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "P001-N1-DNA1-WGS1": "P001",
        "P002-N1-DNA1-WGS1": "P002",
        "P003-N1-DNA1-WGS1": "P003",
    }
    actual = variant_export_external_workflow.get_params("varfish_annotator_external", "bam_qc")(
        wildcards
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage_bam_qc(variant_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart.get_resource_usage() - action 'bam_qc'"""
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "14336M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'bam_qc'."
        actual = variant_export_external_workflow.get_resource(
            "varfish_annotator_external", "bam_qc", resource
        )
        assert actual == expected, msg_error


# Tests for VariantExportExternalWorkflow  ---------------------------------------------------------


def test_variant_export_external_workflow(variant_export_external_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "bam_reports",
        "link_in_bai_external",
        "link_in_bam_external",
        "link_in_vcf_external",
        "link_out",
        "varfish_annotator_external",
        "write_pedigree_with_sample_name",
    ]
    actual = list(sorted(variant_export_external_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/varfish_annotated.P00{i}-N1-DNA1-WGS1/{dir_}/"
        "varfish_annotated{type_}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    # Files in `out` directory
    expected = [
        tpl.format(dir_="out", i=i, ext=ext, type_="")
        for i in (1, 4)  # only for indices
        for ext in (
            "gts.tsv.gz",
            "gts.tsv.gz.md5",
            "db-infos.tsv.gz",
            "db-infos.tsv.gz.md5",
            "bam-qc.tsv.gz",
            "bam-qc.tsv.gz.md5",
        )
    ]
    # Files in `log` directory
    expected += [
        tpl.format(dir_="log", i=i, ext=ext, type_=type_)
        for i in (1, 4)  # only for indices
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
        )
        for type_ in (".annotate", ".bam_qc")
    ]
    expected = sorted(expected)
    actual = sorted(variant_export_external_workflow.get_result_files())
    assert actual == expected
