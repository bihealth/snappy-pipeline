# -*- coding: utf-8 -*-
"""Tests for the targeted_seq_cnv_export workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.targeted_seq_cnv_export import TargetedSeqCnvExportWorkflow

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

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          targeted_seq_cnv_calling:
            tools:
              - gcnv
            gcnv:
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed

          targeted_seq_cnv_export:
            path_ngs_mapping: ../ngs_mapping
            path_targeted_seq_cnv_annotation: ../CNV_ANNOTATION

        data_sets:
          first_batch:
            file: sheet_large_cohort_trio.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: germline_variants
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def targeted_seq_cnv_export_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2,
    mocker,
):
    """Return TargetedSeqCnvExportWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "targeted_seq_cnv_annotation": lambda x: "CNV_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return TargetedSeqCnvExportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VarfishAnnotatorAnnotateStepPart -------------------------------------------------------


def test_varfish_annotator_step_part_get_input_files(targeted_seq_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_input_files()"""
    cnv_base_name = (
        "CNV_ANNOTATION/output/{mapper}.{var_caller}.annotated.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.annotated.{index_ngs_library}"
    )
    expected = {
        "ped": "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped",
        "vcf": cnv_base_name + ".vcf.gz",
        "tbi": cnv_base_name + ".vcf.gz.tbi",
    }
    actual = targeted_seq_cnv_export_workflow.get_input_files("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_output_files(targeted_seq_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    expected = {
        "gts": base_name_out + ".gts.tsv.gz",
        "gts_md5": base_name_out + ".gts.tsv.gz.md5",
        "feature_effects": base_name_out + ".feature-effects.tsv.gz",
        "feature_effects_md5": base_name_out + ".feature-effects.tsv.gz.md5",
        "db_infos": base_name_out + ".db-infos.tsv.gz",
        "db_infos_md5": base_name_out + ".db-infos.tsv.gz.md5",
    }
    actual = targeted_seq_cnv_export_workflow.get_output_files("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_log_file(targeted_seq_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_log_file()"""
    base_name_wrapper = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    base_name_log = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    wrapper_dict = {
        "wrapper": base_name_wrapper + ".wrapper.py",
        "wrapper_md5": base_name_wrapper + ".wrapper.py.md5",
    }
    log_dict = get_expected_log_files_dict(base_out=base_name_log)
    expected = {**wrapper_dict, **log_dict}
    actual = targeted_seq_cnv_export_workflow.get_log_file("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage(targeted_seq_cnv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "14336M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = targeted_seq_cnv_export_workflow.get_resource(
            "varfish_annotator", "annotate", resource
        )
        assert actual == expected, msg_error


# Tests for TargetedSeqCnvExportWorkflow -----------------------------------------------------------


def test_targeted_seq_cnv_export_workflow(targeted_seq_cnv_export_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "varfish_annotator", "write_pedigree"]
    actual = list(sorted(targeted_seq_cnv_export_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.varfish_annotated.P{i}-N1-DNA1-WGS1/log/"
        "{mapper}.{var_caller}.varfish_annotated.P{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, var_caller=var_caller, i=i, ext=ext)
        for i in [str(index_i).zfill(3) for index_i in range(1, 501, 3)]  # only for indices
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
        )
        for mapper in ("bwa",)
        for var_caller in ("gcnv",)
    ]
    actual = targeted_seq_cnv_export_workflow.get_result_files()
    assert actual == expected
