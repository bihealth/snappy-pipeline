# -*- coding: utf-8 -*-
"""Tests for the wgs_mei_calling workflow module code"""


import pytest
import ruamel.yaml as yaml
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_mei_calling import WgsMeiCallingWorkflow

from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (somatic) configuration"""
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
def wgs_mei_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsMeiCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return WgsMeiCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for MeltStepPart (preprocess) -------------------------------------------------------------


def test_melt_step_part_get_input_files_preprocess(wgs_mei_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_mei_calling_workflow.get_input_files("melt", "preprocess")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_melt_step_part_get_output_files_preprocess(wgs_mei_calling_workflow):
    expected = {
        "disc_bai": "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam.disc.bai",
        "disc_bam": "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam.disc",
        "disc_fq": "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam.fq",
        "orig_bai": "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam.bai",
        "orig_bam": "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam",
    }
    assert wgs_mei_calling_workflow.get_output_files("melt", "preprocess") == expected


def test_melt_step_part_get_log_file_preprocess(wgs_mei_calling_workflow):
    expected = "work/{mapper}.melt.preprocess.{library_name}/log/snakemake.wgs_mei_calling.log"
    assert wgs_mei_calling_workflow.get_log_file("melt", "preprocess") == expected


def test_melt_step_part_update_cluster_config_preprocess(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_preprocess"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MeltStepPart (indiv_analysis) ---------------------------------------------------------


def test_melt_step_part_get_input_files_indiv_analysis(wgs_mei_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_mei_calling_workflow.get_input_files("melt", "indiv_analysis")(wildcards)
    expected = {
        "disc_bai": "work/bwa.melt.preprocess.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.bam.disc.bai",
        "disc_bam": "work/bwa.melt.preprocess.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.bam.disc",
        "orig_bai": "work/bwa.melt.preprocess.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.bam.bai",
        "orig_bam": "work/bwa.melt.preprocess.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_melt_step_part_get_output_files_indiv_analysis(wgs_mei_calling_workflow):
    expected = {"done": "work/{mapper}.melt.indiv_analysis.{me_type}/out/.done.{library_name}"}
    assert wgs_mei_calling_workflow.get_output_files("melt", "indiv_analysis") == expected


def test_melt_step_part_get_log_file_indiv_analysis(wgs_mei_calling_workflow):
    expected = "work/{mapper}.melt.indiv_analysis.{me_type}/log/snakemake.wgs_mei_calling.{library_name}.log"
    actual = wgs_mei_calling_workflow.get_log_file("melt", "indiv_analysis")
    assert actual == expected


def test_melt_step_part_update_cluster_config_indiv_analysis(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_indiv_analysis"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MeltStepPart (group_analysis) ---------------------------------------------------------


def test_melt_step_part_get_input_files_group_analysis(wgs_mei_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "me_type": "ALU"})
    actual = wgs_mei_calling_workflow.get_input_files("melt", "group_analysis")(wildcards)
    expected = [
        "work/bwa.melt.indiv_analysis.ALU/out/.done.P001-N1-DNA1-WGS1",
        "work/bwa.melt.indiv_analysis.ALU/out/.done.P002-N1-DNA1-WGS1",
        "work/bwa.melt.indiv_analysis.ALU/out/.done.P003-N1-DNA1-WGS1",
        "work/bwa.melt.indiv_analysis.ALU/out/.done.P004-N1-DNA1-WGS1",
        "work/bwa.melt.indiv_analysis.ALU/out/.done.P005-N1-DNA1-WGS1",
        "work/bwa.melt.indiv_analysis.ALU/out/.done.P006-N1-DNA1-WGS1",
    ]
    assert actual == expected


def test_melt_step_part_get_output_files_group_analysis(wgs_mei_calling_workflow):
    expected = {"done": "work/{mapper}.melt.group_analysis.{me_type}/out/.done"}
    assert wgs_mei_calling_workflow.get_output_files("melt", "group_analysis") == expected


def test_melt_step_part_get_log_file_group_analysis(wgs_mei_calling_workflow):
    expected = "work/{mapper}.melt.group_analysis.{me_type}/log/snakemake.wgs_mei_calling.log"
    assert wgs_mei_calling_workflow.get_log_file("melt", "group_analysis") == expected


def test_melt_step_part_update_cluster_config_group_analysis(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_group_analysis"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MeltStepPart (genotype) ---------------------------------------------------------------


def test_melt_step_part_get_input_files_genotype(wgs_mei_calling_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "me_type": "ALU", "library_name": "P001-N1-DNA1-WGS1"}
    )
    actual = wgs_mei_calling_workflow.get_input_files("melt", "genotype")(wildcards)
    expected = {
        "bam": "work/bwa.melt.preprocess.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.bam",
        "done": "work/bwa.melt.group_analysis.ALU/out/.done",
    }
    assert actual == expected


def test_melt_step_part_get_output_files_genotype(wgs_mei_calling_workflow):
    expected = {"done": "work/{mapper}.melt.genotype.{me_type}/out/.done.{library_name}"}
    assert wgs_mei_calling_workflow.get_output_files("melt", "genotype") == expected


def test_melt_step_part_get_log_file_genotype(wgs_mei_calling_workflow):
    expected = (
        "work/{mapper}.melt.genotype.{me_type}/log/snakemake.wgs_mei_calling.{library_name}.log"
    )
    actual = wgs_mei_calling_workflow.get_log_file("melt", "genotype")
    assert actual == expected


def test_melt_step_part_update_cluster_config_genotype(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_genotype"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MeltStepPart (make_vcf) ---------------------------------------------------------------


def test_melt_step_part_get_input_files_make_vcf(wgs_mei_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "me_type": "ALU"})
    actual = wgs_mei_calling_workflow.get_input_files("melt", "make_vcf")(wildcards)
    expected = [
        "work/bwa.melt.group_analysis.ALU/out/.done",
        "work/bwa.melt.genotype.ALU/out/.done.P001-N1-DNA1-WGS1",
        "work/bwa.melt.genotype.ALU/out/.done.P002-N1-DNA1-WGS1",
        "work/bwa.melt.genotype.ALU/out/.done.P003-N1-DNA1-WGS1",
        "work/bwa.melt.genotype.ALU/out/.done.P004-N1-DNA1-WGS1",
        "work/bwa.melt.genotype.ALU/out/.done.P005-N1-DNA1-WGS1",
        "work/bwa.melt.genotype.ALU/out/.done.P006-N1-DNA1-WGS1",
    ]
    assert actual == expected


def test_melt_step_part_get_output_files_make_vcf(wgs_mei_calling_workflow):
    expected = {
        "done": "work/{mapper}.melt.make_vcf.{me_type}/out/.done",
        "list_txt": "work/{mapper}.melt.genotype.{me_type}/out/list.txt",
        "tbi": "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz.tbi",
        "vcf": "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz",
    }
    assert wgs_mei_calling_workflow.get_output_files("melt", "make_vcf") == expected


def test_melt_step_part_get_log_file_make_vcf(wgs_mei_calling_workflow):
    expected = "work/{mapper}.melt.make_vcf.{me_type}/log/snakemake.wgs_mei_calling.log"
    assert wgs_mei_calling_workflow.get_log_file("melt", "make_vcf") == expected


def test_melt_step_part_update_cluster_config_make_vcf(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_make_vcf"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MeltStepPart (merge_vcf) --------------------------------------------------------------


def test_melt_step_part_get_input_files_merge_vcf(wgs_mei_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    actual = wgs_mei_calling_workflow.get_input_files("melt", "merge_vcf")(wildcards)
    expected = [
        "work/bwa.melt.merge_vcf.ALU/out/ALU.final_comp.vcf.gz",
        "work/bwa.melt.merge_vcf.LINE1/out/LINE1.final_comp.vcf.gz",
        "work/bwa.melt.merge_vcf.SVA/out/SVA.final_comp.vcf.gz",
    ]
    assert actual == expected


def test_melt_step_part_get_output_files_merge_vcf(wgs_mei_calling_workflow):
    expected = {
        "tbi": "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz",
        "vcf_md5": "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.md5",
    }
    assert wgs_mei_calling_workflow.get_output_files("melt", "merge_vcf") == expected


def test_melt_step_part_get_log_file_merge_vcf(wgs_mei_calling_workflow):
    expected = "work/{mapper}.melt.merge_vcf/log/snakemake.wgs_mei_calling.log"
    assert wgs_mei_calling_workflow.get_log_file("melt", "merge_vcf") == expected


def test_melt_step_part_update_cluster_config_merge_vcf(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_merge_vcf"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for MeltStepPart (reorder_vcf) ------------------------------------------------------------


def test_melt_step_part_get_input_files_reorder_vcf(wgs_mei_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_mei_calling_workflow.get_input_files("melt", "reorder_vcf")(wildcards)
    expected = {
        "tbi": "work/bwa.melt.merge_vcf/out/bwa.melt.merge_vcf.vcf.gz.tbi",
        "vcf": "work/bwa.melt.merge_vcf/out/bwa.melt.merge_vcf.vcf.gz",
    }
    assert actual == expected


def test_melt_step_part_get_output_files_reorder_vcf(wgs_mei_calling_workflow):
    expected = {
        "tbi": "work/{mapper}.melt.{index_library_name}/out/{mapper}.melt.{index_library_name}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.melt.{index_library_name}/out/{mapper}.melt.{index_library_name}.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.melt.{index_library_name}/out/{mapper}.melt.{index_library_name}.vcf.gz",
        "vcf_md5": "work/{mapper}.melt.{index_library_name}/out/{mapper}.melt.{index_library_name}.vcf.gz.md5",
    }
    assert wgs_mei_calling_workflow.get_output_files("melt", "reorder_vcf") == expected


def test_melt_step_part_get_log_file_reorder_vcf(wgs_mei_calling_workflow):
    expected = (
        "work/{mapper}.melt.reorder_vcf.{index_library_name}/log/snakemake.wgs_mei_calling.log"
    )
    assert wgs_mei_calling_workflow.get_log_file("melt", "reorder_vcf") == expected


def test_melt_step_part_update_cluster_config_reorder_vcf(
    wgs_mei_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["wgs_mei_calling_melt_reorder_vcf"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for SomaticWgsSvCallingWorkflow -----------------------------------------------------------


def test_sv_calling_workflow(wgs_mei_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["link_out", "melt"]
    assert list(sorted(wgs_mei_calling_workflow.sub_steps.keys())) == expected
    # Check result file construction
    tpl = (
        "output/{mapper}.{mei_caller}.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{mei_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, mei_caller=mei_caller, i=i, ext=ext)
        for i in [1, 4]  # only indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for mei_caller in ("melt",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_mei_calling_workflow.get_result_files()))
    assert expected == actual
