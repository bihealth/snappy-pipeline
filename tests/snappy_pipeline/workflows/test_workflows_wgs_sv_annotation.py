# -*- coding: utf-8 -*-
"""Tests for the wgs_sv_annotation workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_sv_annotation import WgsSvAnnotationWorkflow

from .common import get_expected_output_vcf_files_dict
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
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          wgs_sv_annotation:
            path_ngs_mapping: ../NGS_MAPPING
            path_variant_calling: ../VAR_CALLING
            path_wgs_sv_calling: ../WGS_SV_CALLING
            tool_ngs_mapping_variant_calling: bwa
            tool_variant_calling: gatk_hc
            tools_ngs_mapping: [bwa]
            tools_wgs_sv_calling: [delly2]

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
def wgs_sv_annotation_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
    fai_file_content,
):
    """Return WgsSvAnnotationWorkflow object pre-configured with germline sheet"""
    # Add Fasta file
    # Create FASTA files
    germline_sheet_fake_fs.fs.create_file(
        "/path/to/ref.fa.fai", contents=fai_file_content, create_missing_dirs=True
    )
    germline_sheet_fake_fs.fs.create_file("/path/to/ref.fa", create_missing_dirs=True)
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.wgs_sv_calling", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really a NGSMappingPipelineStep here
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VAR_CALLING/" + x,
        "wgs_sv_calling": lambda x: "WGS_SV_CALLING/" + x,
    }
    # Construct the workflow object
    return WgsSvAnnotationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VcfSvFilterStepPart  -------------------------------------------------------------------


def test_vcf_sv_filter_step_part_call_get_input_files_delly2(wgs_sv_annotation_workflow):
    """Tests VcfSvFilterStepPart.get_input_files() - caller ``delly2``"""
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1", "caller": "delly2"}
    )
    # Define expected
    sv_base_name = (
        "WGS_SV_CALLING/output/bwa.delly2.P001-N1-DNA1-WGS1/out/bwa.delly2.P001-N1-DNA1-WGS1"
    )
    v_base_name = "VAR_CALLING/work/bwa.gatk_hc.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.P001-N1-DNA1-WGS1"
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "sv_bcf": sv_base_name + ".vcf.gz",
        "sv_csi": sv_base_name + ".vcf.gz.tbi",
        "var_vcf": v_base_name + ".vcf.gz",
        "var_tbi": v_base_name + ".vcf.gz.tbi",
    }
    # Get actual
    actual = wgs_sv_annotation_workflow.get_input_files("vcf_sv_filter", "run")(wildcards)
    assert actual == expected


def test_vcf_sv_filter_step_part_call_get_input_files_popdel(wgs_sv_annotation_workflow):
    """Tests VcfSvFilterStepPart.get_input_files() - caller ``popdel``"""
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1", "caller": "popdel"}
    )
    # Define expected
    sv_base_name = (
        "WGS_SV_CALLING/work/bwa.popdel.internal.concat_calls/out/bwa.popdel.internal.concat_calls"
    )
    v_base_name = "VAR_CALLING/work/bwa.gatk_hc.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.P001-N1-DNA1-WGS1"
    expected = {
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "sv_bcf": sv_base_name + ".vcf.gz",
        "sv_csi": sv_base_name + ".vcf.gz.tbi",
        "var_vcf": v_base_name + ".vcf.gz",
        "var_tbi": v_base_name + ".vcf.gz.tbi",
    }
    # Get actual
    actual = wgs_sv_annotation_workflow.get_input_files("vcf_sv_filter", "run")(wildcards)
    assert actual == expected


def test_vcf_sv_filter_step_part_call_get_output_files(wgs_sv_annotation_workflow):
    """Tests VcfSvFilterStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/"
        "{mapper}.{caller}.annotated.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_sv_annotation_workflow.get_output_files("vcf_sv_filter", "run")
    assert actual == expected


def test_vcf_sv_filter_step_part_call_get_log_file(wgs_sv_annotation_workflow):
    """Tests VcfSvFilterStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{caller}.annotated.{index_ngs_library}/log/snakemake.wgs_sv_filter.log"
    )
    actual = wgs_sv_annotation_workflow.get_log_file("vcf_sv_filter", "run")
    assert actual == expected


def test_vcf_sv_filter_step_part_get_resource_usage(wgs_sv_annotation_workflow):
    """Tests VcfSvFilterStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "10240M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_sv_annotation_workflow.get_resource("vcf_sv_filter", "run", resource)
        assert actual == expected, msg_error


# Tests for WgsSvAnnotationWorkflow   --------------------------------------------------------------


def test_wgs_sv_annotation_workflow(wgs_sv_annotation_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = ["link_out", "vcf_sv_filter", "write_pedigree"]
    actual = list(sorted(wgs_sv_annotation_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{sv_caller}.annotated.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{sv_caller}.annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, sv_caller=sv_caller, i=i, ext=ext)
        for mapper in ("bwa",)
        for sv_caller in ("delly2",)
        for i in [1, 4]  # only indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_sv_annotation_workflow.get_result_files()))
    assert actual == expected
