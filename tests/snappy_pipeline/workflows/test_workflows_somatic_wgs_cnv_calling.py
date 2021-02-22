# -*- coding: utf-8 -*-
"""Tests for the somatic_wgs_cnv_calling workflow module code"""


import pytest
import ruamel.yaml as yaml
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_wgs_cnv_calling import SomaticWgsCnvCallingWorkflow

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
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa
            star:
              path_index: /path/to/star/index
          somatic_wgs_cnv_calling:
            somatic_variant_calling_tool: mutect
            tools:
            - canvas
            - cnvetti
            - control_freec
            - cnvkit
            tools_ngs_mapping:
                - bwa
            canvas:
              reference: /path/to/reference.fasta
              filter_bed: /path/to/filter.bed
              genome_folder: /path/to/genome/folder

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def somatic_wgs_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticWgsCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had NGSMappingPipelineStep etc. here
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_calling": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return SomaticWgsCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CanvasSomaticWgsStepPart --------------------------------------------------------------


def test_canvas_somatic_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "cancer_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("canvas", "run")(wildcards)
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        #         "somatic_tbi": "SOMATIC_VARIANT_CALLING/output/bwa.mutect.P001-T1-DNA1-WGS1/out/bwa.mutect.P001-T1-DNA1-WGS1.vcf.gz.tbi",
        #         "somatic_vcf": "SOMATIC_VARIANT_CALLING/output/bwa.mutect.P001-T1-DNA1-WGS1/out/bwa.mutect.P001-T1-DNA1-WGS1.vcf.gz",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_canvas_somatic_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    expected = {
        "tbi": "work/{mapper}.canvas.{cancer_library}/out/{mapper}.canvas.{cancer_library}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.canvas.{cancer_library}/out/{mapper}.canvas.{cancer_library}.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.canvas.{cancer_library}/out/{mapper}.canvas.{cancer_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.canvas.{cancer_library}/out/{mapper}.canvas.{cancer_library}.vcf.gz.md5",
    }
    assert somatic_wgs_cnv_calling_workflow.get_output_files("canvas", "run") == expected


def test_canvas_somatic_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    expected = {
        "conda_info": "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}.conda_info.txt",
        "conda_info_md5": "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}.conda_info.txt.md5",
        "conda_list": "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}.conda_list.txt",
        "conda_list_md5": "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}.conda_list.txt.md5",
        "log": "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}.log",
        "log_md5": "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}.log.md5",
    }
    assert somatic_wgs_cnv_calling_workflow.get_log_file("canvas", "run") == expected


def test_canvas_somatic_step_part_update_cluster_config(
    somatic_wgs_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_wgs_cnv_calling_canvas_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for ControlFreecStepPart --------------------------------------------------------------


def test_control_freec_somatic_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "cancer_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("control_freec", "run")(wildcards)
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_control_freec_somatic_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    expected = {
        "ratio": "work/{mapper}.control_freec.{cancer_library}/out/{mapper}.control_freec.{cancer_library}.ratio.txt",
        "ratio_md5": "work/{mapper}.control_freec.{cancer_library}/out/{mapper}.control_freec.{cancer_library}.ratio.txt.md5",
    }
    assert somatic_wgs_cnv_calling_workflow.get_output_files("control_freec", "run") == expected


def test_control_freec_somatic_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    expected = {
        "conda_info": "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}.conda_info.txt",
        "conda_info_md5": "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}.conda_info.txt.md5",
        "conda_list": "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}.conda_list.txt",
        "conda_list_md5": "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}.conda_list.txt.md5",
        "log": "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}.log",
        "log_md5": "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}.log.md5",
    }
    assert somatic_wgs_cnv_calling_workflow.get_log_file("control_freec", "run") == expected


def test_control_freec_somatic_step_part_update_cluster_config(
    somatic_wgs_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_wgs_cnv_calling_control_freec_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvkitSomaticWgsStepPart ----------------------------------------------------------


def test_cnvkit_somatic_wgs_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "cancer_library": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "run")(wildcards)
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_cnvkit_somatic_wgs_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    expected = {
        "segment": "work/{mapper}.cnvkit.{cancer_library}/out/{mapper}.cnvkit.{cancer_library}.cns",
        "segment_md5": "work/{mapper}.cnvkit.{cancer_library}/out/{mapper}.cnvkit.{cancer_library}.cns.md5",
        "bins": "work/{mapper}.cnvkit.{cancer_library}/out/{mapper}.cnvkit.{cancer_library}.cnr",
        "bins_md5": "work/{mapper}.cnvkit.{cancer_library}/out/{mapper}.cnvkit.{cancer_library}.cnr.md5",
        "scatter": "work/{mapper}.cnvkit.{cancer_library}/out/{mapper}.cnvkit.{cancer_library}.scatter.png",
        "scatter_md5": "work/{mapper}.cnvkit.{cancer_library}/out/{mapper}.cnvkit.{cancer_library}.scatter.png.md5",
    }
    assert somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "run") == expected


def test_cnvkit_somatic_wgs_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    expected = {
        "conda_info": "work/{mapper}.cnvkit.{cancer_library}/log/{mapper}.cnvkit.{cancer_library}.conda_info.txt",
        "conda_info_md5": "work/{mapper}.cnvkit.{cancer_library}/log/{mapper}.cnvkit.{cancer_library}.conda_info.txt.md5",
        "conda_list": "work/{mapper}.cnvkit.{cancer_library}/log/{mapper}.cnvkit.{cancer_library}.conda_list.txt",
        "conda_list_md5": "work/{mapper}.cnvkit.{cancer_library}/log/{mapper}.cnvkit.{cancer_library}.conda_list.txt.md5",
        "log": "work/{mapper}.cnvkit.{cancer_library}/log/{mapper}.cnvkit.{cancer_library}.log",
        "log_md5": "work/{mapper}.cnvkit.{cancer_library}/log/{mapper}.cnvkit.{cancer_library}.log.md5",
    }
    assert somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "run") == expected


def test_cnvkit_somatic_wgs_step_part_update_cluster_config(
    somatic_wgs_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_wgs_cnv_calling_cnvkit_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for SomaticWgsCnvCallingWorkflow ----------------------------------------------------------


def test_somatic_cnv_calling_workflow(somatic_wgs_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["canvas", "cnvetti", "cnvkit", "control_freec", "link_out"]
    assert list(sorted(somatic_wgs_cnv_calling_workflow.sub_steps.keys())) == expected
    # Check result file construction
    tpl = (
        "output/{mapper}.{cnv_caller}.P00{i}-T{t}-DNA1-WGS1/out/"
        "{mapper}.{cnv_caller}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = []
    # -- add files from canvas
    expected += [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for cnv_caller in ("canvas",)
    ]
    # -- add files from cnvetti
    expected += [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("bcf", "bcf.md5", "bcf.csi", "bcf.csi.md5")
        for mapper in ("bwa",)
        for cnv_caller in ("cnvetti",)
    ]
    tpl2 = "output/{mapper}.cnvetti_plot.P00{i}/out/" "{mapper}.cnvetti_plot.P00{i}_{chrom}.{ext}"
    expected += [
        tpl2.format(mapper=mapper, i=i, ext=ext, chrom=chrom)
        for i in (1, 2)
        for ext in ("png", "png.md5")
        for mapper in ("bwa",)
        for chrom in (["chr{}".format(x) for x in (list(range(1, 23)) + ["X", "Y"])] + ["genome"])
    ]
    # -- add files from control_freec
    expected += [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "ratio.txt",
            "ratio.txt.md5",
            "gene_log2.txt",
            "gene_call.txt",
            "segments.txt",
            "heatmap.png",
            "scatter.png",
            "diagram.pdf",
        )
        for mapper in ("bwa",)
        for cnv_caller in ("control_freec",)
    ]
    # -- add files from cnvkit
    expected += [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("cnr", "cnr.md5", "cns", "cns.md5", "scatter.png", "scatter.png.md5")
        for mapper in ("bwa",)
        for cnv_caller in ("cnvkit",)
    ]
    # Perform the comparison
    expected = list(sorted(expected))
    actual = list(sorted(somatic_wgs_cnv_calling_workflow.get_result_files()))
    assert expected == actual
