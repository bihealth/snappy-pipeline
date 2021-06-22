# -*- coding: utf-8 -*-
"""Tests for the wgs_sv_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_sv_calling import WgsSvCallingWorkflow

from .common import get_expected_output_bcf_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (somatic) configuration"""
    config_str = textwrap.dedent(
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

          wgs_sv_calling:
            tools:
              dna:
              - manta
              - delly2
              long_dna:
              - pb_honey_spots

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
    yaml_ = yaml.YAML()
    return yaml_.load(config_str)


@pytest.fixture
def wgs_sv_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsSvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return WgsSvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for Delly2StepPart (call) -----------------------------------------------------------------


def test_delly2_step_part_call_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "call")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_delly2_step_part_call_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.call.{library_name,[^\.]+}/out/"
        r"{mapper}.delly2.call.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "call")
    assert actual == expected


def test_delly_step_part_call_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.delly2.call.{library_name}/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("delly2", "call") == expected


def test_delly_step_part_call_update_cluster_config(wgs_sv_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["wgs_sv_calling_delly2_call"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for Delly2StepPart (merge_calls) ----------------------------------------------------------


def test_delly2_step_part_merge_calls_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "merge_calls")(wildcards)
    expected = [
        "work/bwa.delly2.call.P001-N1-DNA1-WGS1/out/bwa.delly2.call.P001-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P002-N1-DNA1-WGS1/out/bwa.delly2.call.P002-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P003-N1-DNA1-WGS1/out/bwa.delly2.call.P003-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P004-N1-DNA1-WGS1/out/bwa.delly2.call.P004-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P005-N1-DNA1-WGS1/out/bwa.delly2.call.P005-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P006-N1-DNA1-WGS1/out/bwa.delly2.call.P006-N1-DNA1-WGS1.bcf",
    ]
    assert actual == expected


def test_delly2_step_part_merge_calls_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = r"work/{mapper,[^\.]+}.delly2.merge_calls/out/{mapper}.delly2.merge_calls"
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "merge_calls")
    assert actual == expected


def test_delly_step_part_merge_calls_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.delly2.merge_calls/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("delly2", "merge_calls") == expected


def test_delly_step_part_merge_calls_update_cluster_config(
    wgs_sv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_sv_calling_delly2_merge_calls"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for Delly2StepPart (genotype) -------------------------------------------------------------


def test_delly2_step_part_genotype_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "genotype")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bcf": "work/bwa.delly2.merge_calls/out/bwa.delly2.merge_calls.bcf",
    }
    assert actual == expected


def test_delly2_step_part_genotype_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.genotype.{library_name,[^\.]+}/out/"
        r"{mapper}.delly2.genotype.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "genotype")
    assert actual == expected


def test_delly_step_part_genotype_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.delly2.genotype.{library_name}/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("delly2", "genotype") == expected


def test_delly_step_part_genotype_update_cluster_config(
    wgs_sv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_sv_calling_delly2_genotype"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for Delly2StepPart (merge_genotypes) ------------------------------------------------------


def test_delly2_step_part_merge_genotypes_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "merge_genotypes")(wildcards)
    expected = [
        "work/bwa.delly2.genotype.P001-N1-DNA1-WGS1/out/bwa.delly2.genotype.P001-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P002-N1-DNA1-WGS1/out/bwa.delly2.genotype.P002-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P003-N1-DNA1-WGS1/out/bwa.delly2.genotype.P003-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P004-N1-DNA1-WGS1/out/bwa.delly2.genotype.P004-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P005-N1-DNA1-WGS1/out/bwa.delly2.genotype.P005-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P006-N1-DNA1-WGS1/out/bwa.delly2.genotype.P006-N1-DNA1-WGS1.bcf",
    ]
    assert actual == expected


def test_delly2_step_part_merge_genotypes_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.merge_genotypes/out/{mapper}.delly2.merge_genotypes"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "merge_genotypes")
    assert actual == expected


def test_delly_step_part_merge_genotypes_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.delly2.merge_genotypes/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("delly2", "merge_genotypes") == expected


def test_delly_step_part_merge_genotypes_update_cluster_config(
    wgs_sv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_sv_calling_delly2_merge_genotypes"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for Delly2StepPart (reorder_vcf) ----------------------------------------------------------


def test_delly2_step_part_reorder_vcf_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "reorder_vcf")(wildcards)
    expected = {"bcf": "work/bwa.delly2.merge_genotypes/out/bwa.delly2.merge_genotypes.bcf"}
    assert actual == expected


def test_delly2_step_part_reorder_vcf_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.{index_ngs_library,[^\.]+}/out/"
        r"{mapper}.delly2.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "reorder_vcf")
    assert actual == expected


def test_delly_step_part_reorder_vcf_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.delly2.{index_ngs_library}/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("delly2", "reorder_vcf") == expected


def test_delly_step_part_reorder_vcf_update_cluster_config(
    wgs_sv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_sv_calling_delly2_reorder_vcf"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MantaStepPart -------------------------------------------------------------------------


def test_manta_step_part_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("manta", "run")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_manta_step_part_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.manta.{index_ngs_library}/out/{mapper}.manta.{index_ngs_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("manta", "run")
    assert actual == expected


def test_manta_step_part_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.manta.{index_ngs_library}/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("manta", "run") == expected


def test_manta_step_part_update_cluster_config(wgs_sv_calling_workflow, dummy_cluster_config):
    actual = set(dummy_cluster_config["wgs_sv_calling_manta_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for WgsSvCallingWorkflow -----------------------------------------------------------


def test_wgs_sv_calling_workflow(wgs_sv_calling_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = ["delly2", "link_out", "manta", "pb_honey_spots", "popdel", "sniffles", "svtk"]
    actual = list(sorted(wgs_sv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{sv_caller}.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{sv_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, sv_caller=sv_caller, i=i, ext=ext)
        for mapper in ("bwa",)
        for sv_caller in ("delly2", "manta")
        for i in [1, 4]  # only indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_sv_calling_workflow.get_result_files()))
    assert actual == expected
