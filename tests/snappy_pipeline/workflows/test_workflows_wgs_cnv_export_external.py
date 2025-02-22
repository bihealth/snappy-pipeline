# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_export_external workflow module code"""

import textwrap
from copy import deepcopy

import pytest
import ruamel.yaml as ruamel_yaml
from pydantic import ValidationError
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_cnv_export_external import WgsCnvExportExternalWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for germline configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          wgs_cnv_export_external:
            merge_vcf_flag: true
            search_paths: [/search_path]
            search_patterns: [{"vcf": "*/*.vcf.gz"}]
            tool_ngs_mapping: null
            tool_wgs_cnv_calling: dragen
            path_refseq_ser: /data/refseq_ser
            path_ensembl_ser: /data/ensembl_ser
            path_db: /data/db

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
def wgs_cnv_export_external_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsCnvExportExternalWorkflow object pre-configured with germline sheet"""
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
    patch_module_fs("pathlib", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.wgs_cnv_export_external", germline_sheet_fake_fs, mocker
    )
    # Construct the workflow object
    return WgsCnvExportExternalWorkflow(
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
    """Tests WgsCnvExportExternalWorkflow.check_config() - invalid varfish-annotator files"""
    # Create search path
    germline_sheet_fake_fs.fs.makedirs("/search_path")
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("pathlib", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.wgs_cnv_export_external", germline_sheet_fake_fs, mocker
    )
    # Construct the workflow object
    with pytest.raises(ValidationError) as exec_info:
        WgsCnvExportExternalWorkflow(
            dummy_workflow,
            minimal_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )
    errors = exec_info.value.errors()
    assert len(errors) == 3
    assert (
        len(
            {"path_refseq_ser", "path_ensembl_ser", "path_db"}
            & set(s for e in errors for s in e["loc"])
        )
        == 3
    )


def test_workflow_check_config_invalid_search_directory(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Tests WgsCnvExportExternalWorkflow.check_config() - no search directory"""
    # Create annotation files
    for fk_file in ("refseq_ser", "ensembl_ser", "db"):
        germline_sheet_fake_fs.fs.create_file(
            "/data/" + fk_file,
            contents="",
            create_missing_dirs=True,
        )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("pathlib", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.wgs_cnv_export_external", germline_sheet_fake_fs, mocker
    )
    # Construct the workflow object
    with pytest.raises(ValidationError) as exec_info:
        WgsCnvExportExternalWorkflow(
            dummy_workflow,
            minimal_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )

    errors = exec_info.value.errors()
    assert len(errors) == 1

    assert "path_not_directory" in {e["type"] for e in errors}


def test_workflow_check_config_invalid_search_pattern(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsCnvExportExternalWorkflow object pre-configured with germline sheet"""
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
    patch_module_fs("pathlib", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.wgs_cnv_export_external", germline_sheet_fake_fs, mocker
    )
    # Change search patterns to invalid
    modified_config = deepcopy(minimal_config)
    modified_config["step_config"]["wgs_cnv_export_external"]["search_patterns"] = [
        "vcf",
        "*/*.vcf.gz",
    ]
    # Construct the workflow object
    with pytest.raises(ValidationError) as exec_info:
        WgsCnvExportExternalWorkflow(
            dummy_workflow,
            modified_config,
            config_lookup_paths,
            config_paths,
            work_dir,
        )

    errors = exec_info.value.errors()

    # there is 1 incorrectly defined search_patterns entry
    # which is incorrectly defined as a list instead of a dict/key-value pairs,
    # so pydantic tries to parse *2* dicts from the list and fails
    assert len(errors) == 2

    import pydantic
    pydantic_version = pydantic.version.version_short()
    expected_errors = [
        {
            "input": input_str,
            "loc": ("step_config", "wgs_cnv_export_external", "search_patterns", i),
            "msg": "Input should be a valid dictionary",
            "type": "dict_type",
            "url": f"https://errors.pydantic.dev/{pydantic_version}/v/dict_type",
        }
        for i, input_str in enumerate(["vcf", "*/*.vcf.gz"])
    ]
    assert expected_errors == errors


# Tests for VarfishAnnotatorExternalStepPart (merge_vcf)   -----------------------------------------


def test_varfish_annotator_step_part_call_get_input_files_merge_vcf(
    wgs_cnv_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_input_files_merge_vcf()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = [
        "work/input_links/P001-N1-DNA1-WGS1/.done",
        "work/input_links/P002-N1-DNA1-WGS1/.done",
        "work/input_links/P003-N1-DNA1-WGS1/.done",
    ]
    actual = wgs_cnv_export_external_workflow.get_input_files(
        "varfish_annotator_external", "merge_vcf"
    )(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_call_get_output_files_merge_vcf(
    wgs_cnv_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_output_files_merge_vcf()"""
    base_name = "work/dragen.{index_ngs_library}/out/dragen.{index_ngs_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = wgs_cnv_export_external_workflow.get_output_files(
        "varfish_annotator_external", "merge_vcf"
    )
    assert actual == expected


def test_varfish_annotator_step_part_call_get_log_file_merge_vcf(wgs_cnv_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart._get_log_file_merge_vcf()"""
    base_name = "work/dragen.{index_ngs_library}/log/dragen.{index_ngs_library}.merge_vcf"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = wgs_cnv_export_external_workflow.get_log_file(
        "varfish_annotator_external", "merge_vcf"
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_params_merge_vcf(wgs_cnv_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart._get_params_merge_vcf()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "input": [],
        "sample_names": ["P001", "P002", "P003"],
        "merge_option": "id",
        "gvcf_option": False,
    }
    actual = wgs_cnv_export_external_workflow.get_params("varfish_annotator_external", "merge_vcf")(
        wildcards
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage_merge_vcf(wgs_cnv_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart.get_resource_usage() - action 'annotate'"""
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "14336M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'merge_vcf'."
        actual = wgs_cnv_export_external_workflow.get_resource(
            "varfish_annotator_external", "merge_vcf", resource
        )()
        assert actual == expected, msg_error


# Tests for VarfishAnnotatorExternalStepPart (annotate)   ------------------------------------------


def test_varfish_annotator_step_part_call_get_input_files_annotate(
    wgs_cnv_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_input_files_annotate()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    base_name = "work/dragen.P001-N1-DNA1-WGS1/out/dragen.P001-N1-DNA1-WGS1"
    vcf_dict = get_expected_output_vcf_files_dict(base_out=base_name)
    ped_dict = {"ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped"}
    expected = {**ped_dict, **vcf_dict}
    # Get actual
    actual = wgs_cnv_export_external_workflow.get_input_files(
        "varfish_annotator_external", "annotate"
    )(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_call_get_output_files_annotate(
    wgs_cnv_export_external_workflow,
):
    """Tests VarfishAnnotatorExternalStepPart._get_output_files_annotate()"""
    base_name_out = (
        "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
    )
    expected = {
        "gts": base_name_out + ".gts.tsv.gz",
        "gts_md5": base_name_out + ".gts.tsv.gz.md5",
        "feature_effects": base_name_out + ".feature-effects.tsv.gz",
        "feature_effects_md5": base_name_out + ".feature-effects.tsv.gz.md5",
        "db_infos": base_name_out + ".db-infos.tsv.gz",
        "db_infos_md5": base_name_out + ".db-infos.tsv.gz.md5",
        "output_links": [],
    }
    actual = wgs_cnv_export_external_workflow.get_output_files(
        "varfish_annotator_external", "annotate"
    )
    assert actual == expected


def test_varfish_annotator_step_part_call_get_log_file_annotate(wgs_cnv_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart._get_log_file_annotate()"""
    base_name = (
        "work/varfish_annotated.{index_ngs_library}/log/varfish_annotated.{index_ngs_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name, extended=True)
    actual = wgs_cnv_export_external_workflow.get_log_file("varfish_annotator_external", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_params_annotate(wgs_cnv_export_external_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart._get_params_annotate()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "step_name": "wgs_cnv_export_external",
        "varfish_server_compatibility": False,
    }
    actual = wgs_cnv_export_external_workflow.get_params("varfish_annotator_external", "annotate")(
        wildcards
    )
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage_annotate(wgs_cnv_export_external_workflow):
    """Tests VarfishAnnotatorExternalStepPart.get_resource_usage() - action 'annotate'"""
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "14336M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'annotate'."
        actual = wgs_cnv_export_external_workflow.get_resource(
            "varfish_annotator_external", "annotate", resource
        )()
        assert actual == expected, msg_error


# Tests for WgsCnvExportExternalWorkflow   ----------------------------------------------------------


def test_wgs_sv_annotation_workflow(wgs_cnv_export_external_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = [
        "link_in_vcf_external",
        "link_out",
        "varfish_annotator_external",
        "write_pedigree_with_sample_name",
    ]
    actual = list(sorted(wgs_cnv_export_external_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/varfish_annotated.P00{i}-N1-DNA1-WGS1/{dir_}/"
        "varfish_annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    # Expected files in `out`
    expected = [
        tpl.format(i=i, dir_="out", ext=ext)
        for i in (1, 4)
        for ext in (
            "gts.tsv.gz",
            "gts.tsv.gz.md5",
            "db-infos.tsv.gz",
            "db-infos.tsv.gz.md5",
            "feature-effects.tsv.gz",
            "feature-effects.tsv.gz.md5",
        )
    ]
    # Expected files in `log`
    expected += [
        tpl.format(i=i, dir_="log", ext=ext)
        for i in (1, 4)
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
            "wrapper.py",
            "wrapper.py.md5",
            "environment.yaml",
            "environment.yaml.md5",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_cnv_export_external_workflow.get_result_files()))
    assert actual == expected
