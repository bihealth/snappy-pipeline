# -*- coding: utf-8 -*-
"""Tests for the wgs_mei_annotation workflow module code"""

from pathlib import Path
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.wgs_mei_annotation import WgsMeiAnnotationWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


def get_project_root():
    """Get project root

    :return: Return path to project root as a string.
    """
    return str(Path(__file__).parent.parent.parent.parent)


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (somatic) configuration"""
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
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          wgs_mei_calling:
            tools:
            - melt
            melt:
              path_genes_bed: /path/to/genes.bed
              path_me_refs: /path/to/me/refs

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
def wgs_mei_annotation_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsMeiAnnotationWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really a NGSMappingPipelineStep here
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "wgs_mei_calling": lambda x: "WGS_MEI_CALLING/" + x,
    }
    # Construct the workflow object
    return WgsMeiAnnotationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VcfMeiFilterStepPart   -----------------------------------------------------------------


def test_vcf_mei_filter_step_part_get_input_files(wgs_mei_annotation_workflow):
    """Tests VcfMeiFilterStepPart.get_input_files()"""
    # Define expected
    root = get_project_root()
    wgs_mei_base_name = (
        "WGS_MEI_CALLING/work/{mapper}.{caller}.merge_vcf/out/{mapper}.{caller}.merge_vcf"
    )
    expected = {
        "ped": root + "/work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped",
        "vcf": wgs_mei_base_name + ".vcf.gz",
        "tbi": wgs_mei_base_name + ".vcf.gz.tbi",
    }
    # Get actual
    actual = wgs_mei_annotation_workflow.get_input_files("vcf_mei_filter", "run")
    assert actual == expected


def test_vcf_mei_filter_step_part_get_output_files(wgs_mei_annotation_workflow):
    """Tests VcfMeiFilterStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/"
        "{mapper}.{caller}.annotated.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_mei_annotation_workflow.get_output_files("vcf_mei_filter", "run")
    assert actual == expected


def test_vcf_mei_filter_step_part_get_log_file(wgs_mei_annotation_workflow):
    """Tests VcfMeiFilterStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{caller}.annotated.{index_ngs_library}/log/snakemake.wgs_mei_filter.log"
    )
    actual = wgs_mei_annotation_workflow.get_log_file("vcf_mei_filter", "run")
    assert actual == expected


def test_vcf_mei_filter_step_part_get_resource_usage(wgs_mei_annotation_workflow):
    """Tests VcfMeiFilterStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "10240M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_mei_annotation_workflow.get_resource("vcf_mei_filter", "run", resource)
        assert actual == expected, msg_error


# Tests for WgsMeiAnnotationWorkflow   -------------------------------------------------------------


def test_sv_calling_workflow(wgs_mei_annotation_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "vcf_mei_filter", "write_pedigree"]
    actual = list(sorted(wgs_mei_annotation_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{mei_caller}.annotated.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{mei_caller}.annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, mei_caller=mei_caller, i=i, ext=ext)
        for i in [1, 4]  # only indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for mei_caller in ("melt",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_mei_annotation_workflow.get_result_files()))
    assert actual == expected
