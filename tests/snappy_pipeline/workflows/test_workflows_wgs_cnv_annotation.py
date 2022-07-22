# -*- coding: utf-8 -*-
"""Tests for the wgs_cnv_annotation workflow module code"""


import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_cnv_annotation import WgsCnvAnnotationWorkflow

from .common import get_expected_output_vcf_files_dict
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
def wgs_cnv_annotation_workflow(
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
        "wgs_cnv_calling": lambda x: "WGS_CNV_CALLING/" + x,
    }
    # Construct the workflow object
    return WgsCnvAnnotationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VcfCnvFilterStepPart  ------------------------------------------------------------------


def test_vcf_cnv_filter_step_part_get_input_files(wgs_cnv_annotation_workflow):
    """Tests VcfCnvFilterStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1", "caller": "erds_sv2"}
    )
    # Define expected
    wgs_cnv_calling_path = (
        "WGS_CNV_CALLING/work/bwa.erds_sv2.merge_genotypes/out/bwa.erds_sv2.merge_genotypes"
    )
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "vcf": wgs_cnv_calling_path + ".vcf.gz",
        "tbi": wgs_cnv_calling_path + ".vcf.gz.tbi",
    }
    # Get actual
    actual = wgs_cnv_annotation_workflow.get_input_files("vcf_cnv_filter", "run")(wildcards)
    assert actual == expected


def test_vcf_cnv_filter_step_part_get_output_files_call(wgs_cnv_annotation_workflow):
    """Tests VcfCnvFilterStepPart.get_output_files()"""
    # Define expected
    base_file_name = (
        "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/"
        "{mapper}.{caller}.annotated.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_file_name)
    # Get actual
    actual = wgs_cnv_annotation_workflow.get_output_files("vcf_cnv_filter", "run")
    assert actual == expected


def test_vcf_cnv_filter_step_part_get_log_file(wgs_cnv_annotation_workflow):
    """Tests VcfCnvFilterStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{caller}.annotated.{index_ngs_library}/log/snakemake.wgs_cnv_filter.log"
    )
    actual = wgs_cnv_annotation_workflow.get_log_file("vcf_cnv_filter", "run")
    assert actual == expected


def test_vcf_cnv_filter_step_part_get_resource_usage(wgs_cnv_annotation_workflow):
    """Tests VcfCnvFilterStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "10240M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_cnv_annotation_workflow.get_resource("vcf_cnv_filter", "run", resource)
        assert actual == expected, msg_error


# Tests for WgsCnvAnnotationWorkflow  --------------------------------------------------------------


def test_wgs_cnv_annotation_workflow(wgs_cnv_annotation_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "vcf_cnv_filter", "write_pedigree"]
    actual = list(sorted(wgs_cnv_annotation_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{cnv_caller}.annotated.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{cnv_caller}.annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, ext=ext)
        for i in (1, 4)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for cnv_caller in ("erds_sv2",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_cnv_annotation_workflow.get_result_files()))
    assert actual == expected
