# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_export workflow module code"""


import textwrap

import pytest
from ruamel import yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_cnv_export import WgsCnvExportWorkflow

from .common import get_expected_log_files_dict
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
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          wgs_cnv_calling:
            variant_calling_tool: gatk_ug
            tools:
            - erds_sv2

          wgs_cnv_annotation:
            path_ngs_mapping: ../ngs_mapping
            path_wgs_cnv_calling: ../wgs_cnv_calling
            tools_ngs_mapping: [bwa]
            tools_wgs_cnv_calling: [erds_sv2]

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
def wgs_cnv_export_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsCnvAnnotationWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "wgs_cnv_annotation": lambda x: "WGS_CNV_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return WgsCnvExportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VarfishAnnotatorAnnotateStepPart  ------------------------------------------------------


def test_varfish_annotator_step_part_get_input_files(wgs_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_input_files()"""
    wgs_cnv_annotation_path = (
        "WGS_CNV_ANNOTATION/output/{mapper}.{var_caller}.annotated.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.annotated.{index_ngs_library}"
    )
    expected = {
        "ped": "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped",
        "vcf": wgs_cnv_annotation_path + ".vcf.gz",
        "tbi": wgs_cnv_annotation_path + ".vcf.gz.tbi",
    }
    actual = wgs_cnv_export_workflow.get_input_files("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_output_files_call(wgs_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_output_files()"""
    base_file_name = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    expected = {
        "gts": base_file_name + ".gts.tsv.gz",
        "gts_md5": base_file_name + ".gts.tsv.gz.md5",
        "feature_effects": base_file_name + ".feature-effects.tsv.gz",
        "feature_effects_md5": base_file_name + ".feature-effects.tsv.gz.md5",
        "db_infos": base_file_name + ".db-infos.tsv.gz",
        "db_infos_md5": base_file_name + ".db-infos.tsv.gz.md5",
    }
    actual = wgs_cnv_export_workflow.get_output_files("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_log_file(wgs_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_log_file()"""
    # Define expected
    base_name_wrapper = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    base_name_log = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    expected_wrapper_dict = {
        "wrapper": base_name_wrapper + ".wrapper.py",
        "wrapper_md5": base_name_wrapper + ".wrapper.py.md5",
    }
    expected_log_dict = get_expected_log_files_dict(base_out=base_name_log)
    expected = {**expected_wrapper_dict, **expected_log_dict}
    # Get actual
    actual = wgs_cnv_export_workflow.get_log_file("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_params(wgs_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_params()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {"is_wgs": True, "step_name": "wgs_cnv_export"}
    actual = wgs_cnv_export_workflow.get_params("varfish_annotator", "annotate")(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage(wgs_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "14336M", "partition": None}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_cnv_export_workflow.get_resource("varfish_annotator", "annotate", resource)
        assert actual == expected, msg_error


# Tests for WgsCnvExportWorkflow  ------------------------------------------------------------------


def test_wgs_cnv_annotation_workflow(wgs_cnv_export_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "varfish_annotator", "write_pedigree"]
    actual = list(sorted(wgs_cnv_export_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{cnv_caller}.varfish_annotated.P00{i}-N1-DNA1-WGS1/{dir_}/"
        "{mapper}.{cnv_caller}.varfish_annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    # Expected files in `out`
    expected = [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, dir_="out", ext=ext)
        for i in (1, 4)
        for ext in (
            "gts.tsv.gz",
            "gts.tsv.gz.md5",
            "db-infos.tsv.gz",
            "db-infos.tsv.gz.md5",
            "feature-effects.tsv.gz",
            "feature-effects.tsv.gz.md5",
        )
        for mapper in ("bwa",)
        for cnv_caller in ("erds_sv2",)
    ]
    # Expected files in `log`
    expected += [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, dir_="log", ext=ext)
        for i in (1, 4)
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
        )
        for mapper in ("bwa",)
        for cnv_caller in ("erds_sv2",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_cnv_export_workflow.get_result_files()))
    assert actual == expected
