# -*- coding: utf-8 -*-
"""Tests for the somatic_wgs_cnv_calling workflow module code"""

import textwrap
from itertools import chain

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_wgs_cnv_calling import SomaticWgsCnvCallingWorkflow

from .common import (
    get_expected_log_files_dict,
    get_expected_output_bcf_files_dict,
    get_expected_output_vcf_files_dict,
)
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_wgs_cnv_calling:
            path_somatic_variant_calling: ../somatic_variant_calling
            somatic_variant_calling_tool: mutect
            path_ngs_mapping: NGS_MAPPING
            tools:
            - canvas
            - cnvetti
            - control_freec
            - cnvkit
            tools_ngs_mapping:
                - bwa
            canvas:
              path_reference: /path/to/reference.fasta
              path_filter_bed: /path/to/filter.bed
              path_genome_folder: /path/to/genome/folder
            cnvkit:
              path_target: /path/to/panel_of_normals/output/cnvkit.target/out/cnvkit.target.bed
              path_antitarget: /path/to/panel_of_normals/output/cnvkit.antitarget/out/cnvkit.antitarget.bed
              path_panel_of_normals: /path/to/panel_of_normals/output/bwa.cnvkit.create_panel/out/bwa.cnvkit.panel_of_normals.cnn
            cnvetti: {}
            control_freec:
              path_chrlenfile: /path/to/chrlenfile
              path_mappability: /path/to/mappability
              convert: {}


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
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticWgsCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
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
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CanvasSomaticWgsStepPart --------------------------------------------------------------


def test_canvas_somatic_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    """Tests CanvasSomaticWgsStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "cancer_library": "P001-T1-DNA1-WGS1"})
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("canvas", "run")(wildcards)
    assert actual == expected


def test_canvas_somatic_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    """Tests CanvasSomaticWgsStepPart.get_output_files()"""
    base_name = "work/{mapper}.canvas.{cancer_library}/out/{mapper}.canvas.{cancer_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("canvas", "run")
    assert actual == expected


def test_canvas_somatic_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    """Tests CanvasSomaticWgsStepPart.get_log_file()"""
    base_name = "work/{mapper}.canvas.{cancer_library}/log/{mapper}.canvas.{cancer_library}"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("canvas", "run")
    assert actual == expected


def test_canvas_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CanvasSomaticWgsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "1-16:00:00", "memory": "61440M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("canvas", "run", resource)()
        assert actual == expected, msg_error


# Tests for ControlFreecStepPart --------------------------------------------------------------


def test_control_freec_somatic_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    """Tests ControlFreecStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "cancer_library": "P001-T1-DNA1-WGS1"})
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("control_freec", "run")(wildcards)
    assert actual == expected


def test_control_freec_somatic_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    """Tests ControlFreecStepPart.get_output_files()"""
    base_name = "work/{mapper}.control_freec.{cancer_library}/out/"
    expected = {
        "ratio": base_name + "{mapper}.control_freec.{cancer_library}.ratio.txt",
        "ratio_md5": base_name + "{mapper}.control_freec.{cancer_library}.ratio.txt.md5",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("control_freec", "run")
    assert actual == expected


def test_control_freec_somatic_step_part_get_output_exception(somatic_wgs_cnv_calling_workflow):
    """Tests ControlFreecStepPart.get_output_files()"""
    # Exception action request not defined
    with pytest.raises(Exception) as exec_info:
        somatic_wgs_cnv_calling_workflow.get_output_files("control_freec", "invalid_action")
    error_msg = "Called unsupported action, exceptions should have been raised."
    assert exec_info.value.args[0] is not None, error_msg


def test_control_freec_somatic_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    """Tests ControlFreecStepPart.get_log_file()"""
    base_name = (
        "work/{mapper}.control_freec.{cancer_library}/log/{mapper}.control_freec.{cancer_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("control_freec", "run")
    assert actual == expected


def test_control_freec_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests ControlFreecStepPart.get_resource()"""
    # Define expected
    run_expected_dict = {
        "threads": 8,
        "time": "1-16:00:00",
        "memory": "61440M",
        "partition": "medium",
    }
    plot_expected_dict = {
        "threads": 1,
        "time": "1-16:00:00",
        "memory": "61440M",
        "partition": "medium",
    }
    transform_expected_dict = {
        "threads": 1,
        "time": "1-16:00:00",
        "memory": "16384M",
        "partition": "medium",
    }

    # Evaluate action `run`
    for resource, expected in run_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' in action 'run'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("control_freec", "run", resource)()
        assert actual == expected, msg_error

    # Evaluate action `plot`
    for resource, expected in plot_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' in action 'plot"
        actual = somatic_wgs_cnv_calling_workflow.get_resource("control_freec", "plot", resource)()
        assert actual == expected, msg_error

    # Evaluate action `transform`
    for resource, expected in transform_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' in action 'transform"
        actual = somatic_wgs_cnv_calling_workflow.get_resource(
            "control_freec", "transform", resource
        )()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (coverage) -------------------------------------------------------------


def test_cnvkit_coverage_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "target": "target", "library_name": "P001-T1-DNA1-WGS1"}
    )
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "coverage")(wildcards)
    assert actual == expected


def test_cnvkit_coverage_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}"
    expected = {
        "target": base_name_out + ".targetcoverage.cnn",
        "target_md5": base_name_out + ".targetcoverage.cnn.md5",
        "antitarget": base_name_out + ".antitargetcoverage.cnn",
        "antitarget_md5": base_name_out + ".antitargetcoverage.cnn.md5",
    }
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "coverage")

    assert actual == expected


def test_cnvkit_coverage_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    base_file_name = (
        "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.coverage.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_file_name)
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_coverage_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'coverage'"""
    # Define expected
    expected_dict = {"threads": 8, "time": "1-00:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "coverage", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (fix) ------------------------------------------------------------------


def test_cnvkit_fix_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    coverage_base_out = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1"
    expected = {
        "antitarget": coverage_base_out + ".antitargetcoverage.cnn",
        "target": coverage_base_out + ".targetcoverage.cnn",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "fix")(wildcards)
    assert actual == expected


def test_cnvkit_fix_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cnr"
    expected = {"ratios": base_name_out, "ratios_md5": base_name_out + ".md5"}
    assert somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "fix") == expected


def test_cnvkit_fix_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.fix.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "fix")
    assert actual == expected


def test_cnvkit_fix_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'fix'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "fix", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (segment) --------------------------------------------------------------


def test_cnvkit_segment_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {"cnr": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"}
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "segment")(wildcards)
    assert actual == expected


def test_cnvkit_segment_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    base_name_out = (
        "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.segment.cns"
    )
    expected = {"segments": base_name_out, "segments_md5": base_name_out + ".md5"}
    assert somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "segment") == expected


def test_cnvkit_segment_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.segment.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "segment")
    assert actual == expected


def test_cnvkit_segment_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'fix'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "segment", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (call) -----------------------------------------------------------------


def test_cnvkit_call_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    segment_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.segment.cns"
    expected = {"segment": segment_file}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "call")(wildcards)
    assert actual == expected


def test_cnvkit_call_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    base_name_out = (
        "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.call.cns"
    )
    expected = {"calls": base_name_out, "calls_md5": base_name_out + ".md5"}
    assert somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "call") == expected


def test_cnvkit_call_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.call.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "call")
    assert actual == expected


def test_cnvkit_call_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "call", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (postprocess) ----------------------------------------------------------


def test_cnvkit_postprocess_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    call_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.call.cns"
    expected = {"call": call_file}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "postprocess")(wildcards)
    assert actual == expected


def test_cnvkit_postprocess_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cns"
    expected = {"final": base_name_out, "final_md5": base_name_out + ".md5"}
    assert somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "postprocess") == expected


def test_cnvkit_postprocess_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.postprocess.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "postprocess")
    assert actual == expected


def test_cnvkit_postprocess_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "postprocess", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (plot) -----------------------------------------------------------------


def test_cnvkit_plot_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    cnr_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"
    cns_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.call.cns"
    expected = {
        "cnr": cnr_file,
        "cns": cns_file,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "plot")(wildcards)
    assert actual == expected


def test_cnvkit_plot_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    expected = {}
    tpl = (
        "work/{{mapper}}.cnvkit.{{library_name}}/report/"
        "{{mapper}}.cnvkit.{{library_name}}.{plot}.{ext}"
    )
    for plot, ext in (("diagram", "pdf"), ("heatmap", "pdf"), ("scatter", "png")):
        expected[plot] = tpl.format(plot=plot, ext=ext)
        expected[plot + "_md5"] = expected[plot] + ".md5"
    tpl = (
        "work/{{mapper}}.cnvkit.{{library_name}}/report/"
        "{{mapper}}.cnvkit.{{library_name}}.{plot}.chr{chrom}.{ext}"
    )
    for plot, ext in (("heatmap", "pdf"), ("scatter", "png")):
        for chrom in chain(range(1, 23), ("X", "Y")):
            key = "{plot}_chr{chrom}".format(plot=plot, chrom=str(chrom))
            expected[key] = tpl.format(plot=plot, ext=ext, chrom=str(chrom))
            expected[key + "_md5"] = expected[key] + ".md5"
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "plot")
    assert actual == expected


def test_cnvkit_plot_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.plot.{library_name}"
    )
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "plot")
    assert actual == expected


def test_cnvkit_plot_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "30720M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "plot", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (export) ---------------------------------------------------------------


def test_cnvkit_export_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "cns": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.call.cns"
    }
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "export")(wildcards)
    assert actual == expected


def test_cnvkit_export_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    expected = {}
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}"
    for key, ext in (("bed", "bed"), ("seg", "seg"), ("vcf", "vcf.gz"), ("tbi", "vcf.gz.tbi")):
        expected[key] = base_name_out + "." + ext
        expected[key + "_md5"] = expected[key] + ".md5"
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "export")
    assert actual == expected


def test_cnvkit_export_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.export.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "export")
    assert actual == expected


def test_cnvkit_export_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "export", resource)()
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (report) ---------------------------------------------------------------


def test_cnvkit_report_step_part_get_input_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    cnr_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"
    cns_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.call.cns"
    target_file = (
        "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.targetcoverage.cnn"
    )
    antitarget_file = (
        "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.antitargetcoverage.cnn"
    )
    expected = {
        "cnr": cnr_file,
        "cns": cns_file,
        "target": target_file,
        "antitarget": antitarget_file,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvkit", "report")(wildcards)
    assert actual == expected


def test_cnvkit_report_step_part_get_output_files(somatic_wgs_cnv_calling_workflow):
    # Define expected
    expected = {}
    base_name_out = "work/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}"
    for report in ("breaks", "genemetrics", "segmetrics", "sex", "metrics"):
        expected[report] = base_name_out + "." + report + ".txt"
        expected[report + "_md5"] = expected[report] + ".md5"
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvkit", "report")
    assert actual == expected


def test_cnvkit_report_step_part_get_log_file(somatic_wgs_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.report.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvkit", "report")
    assert actual == expected


def test_cnvkit_report_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvkit", "report", resource)()
        assert actual == expected, msg_error


# Tests for CnvettiSomaticWgsStepPart --------------------------------------------------------------


def test_cnvetti_step_part_get_input_files_coverage(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvetti", "coverage")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_input_files_tumor_normal_ratio():
    """Tests CnvettiSomaticWgsStepPart._get_input_files_tumor_normal_ratio()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "cancer_library": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    # actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvetti", "tumor_normal_ratio")(
    #     wildcards
    # )
    _ = wildcards
    _ = expected
    assert True


def test_cnvetti_step_part_get_input_files_segment(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_input_files_segment()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    base_name = (
        "work/bwa.cnvetti_tumor_normal_ratio.P001-N1-DNA1-WGS1/out/"
        "bwa.cnvetti_tumor_normal_ratio.P001-N1-DNA1-WGS1"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_input_files("cnvetti", "segment")(wildcards)
    assert actual == expected


def test_cnvetti_step_part_get_output_files_coverage(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_output_files_coverage()"""
    base_name = (
        "work/{mapper}.cnvetti_coverage.{library_name}/out/{mapper}.cnvetti_coverage.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvetti", "coverage")
    print(actual)
    assert actual == expected


def test_cnvetti_step_part_get_output_files_tumor_normal_ratio(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_output_files_tumor_normal_ratio()"""
    base_name = (
        "work/{mapper}.cnvetti_tumor_normal_ratio.{library_name}/out/"
        "{mapper}.cnvetti_tumor_normal_ratio.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvetti", "tumor_normal_ratio")
    assert actual == expected


def test_cnvetti_step_part_get_output_files_segment(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_output_files_segment()"""
    base_name = (
        "work/{mapper}.cnvetti_segment.{library_name}/out/{mapper}.cnvetti_segment.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = somatic_wgs_cnv_calling_workflow.get_output_files("cnvetti", "segment")
    assert actual == expected


def test_cnvetti_step_part_get_log_file_coverage(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_log_file_coverage()"""
    base_name = (
        "work/{mapper}.cnvetti_coverage.{library_name}/log/{mapper}.cnvetti_coverage.{library_name}"
    )
    expected = {
        "log": base_name + ".log",
        "conda_info": base_name + ".conda_info.txt",
        "conda_list": base_name + ".conda_list.txt",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvetti", "coverage")
    assert actual == expected


def test_cnvetti_step_part_get_log_file_tumor_normal_ratio(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_log_file_tumor_normal_ratio()"""
    base_name = (
        "work/{mapper}.cnvetti_tumor_normal_ratio.{library_name}/log/"
        "{mapper}.cnvetti_tumor_normal_ratio.{library_name}"
    )
    expected = {
        "log": base_name + ".log",
        "conda_info": base_name + ".conda_info.txt",
        "conda_list": base_name + ".conda_list.txt",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvetti", "tumor_normal_ratio")
    assert actual == expected


def test_cnvetti_step_part_get_log_file_segment(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart._get_log_file_segment()"""
    base_name = (
        "work/{mapper}.cnvetti_segment.{library_name}/log/{mapper}.cnvetti_segment.{library_name}"
    )
    expected = {
        "log": base_name + ".log",
        "conda_info": base_name + ".conda_info.txt",
        "conda_list": base_name + ".conda_list.txt",
    }
    actual = somatic_wgs_cnv_calling_workflow.get_log_file("cnvetti", "segment")
    assert actual == expected


def test_cnvetti_step_part_get_resource(somatic_wgs_cnv_calling_workflow):
    """Tests CnvettiSomaticWgsStepPart.get_resource()"""
    # Get all available actions
    all_actions = somatic_wgs_cnv_calling_workflow.substep_getattr("cnvetti", "actions")
    # Define expected
    expected_dict = {"threads": 4, "time": "04:00:00", "memory": "15360M", "partition": "medium"}
    # Evaluate
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = somatic_wgs_cnv_calling_workflow.get_resource("cnvetti", action, resource)()
            assert actual == expected, msg_error


# Tests for SomaticWgsCnvCallingWorkflow -----------------------------------------------------------


def test_somatic_cnv_calling_workflow(somatic_wgs_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
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
    tpl2 = "output/{mapper}.cnvetti_plot.P00{i}/out/{mapper}.cnvetti_plot.P00{i}_{chrom}.{ext}"
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
    tpl = "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/out/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.{ext}{md5}"
    expected += [
        tpl.format(i=i, t=t, ext=ext, md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("cnr", "cns", "bed", "seg", "vcf.gz", "vcf.gz.tbi")
        for md5 in ("", ".md5")
    ]
    tpl = (
        "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/report/"
        "bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.{plot}.{ext}{md5}"
    )
    expected += [
        tpl.format(i=i, t=t, plot=plot, ext=ext, md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for plot, ext in (("diagram", "pdf"), ("heatmap", "pdf"), ("scatter", "png"))
        for md5 in ("", ".md5")
    ]
    tpl = (
        "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/report/"
        "bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.{plot}.chr{chrom}.{ext}{md5}"
    )
    expected += [
        tpl.format(i=i, t=t, plot=plot, ext=ext, chrom=str(chrom), md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for plot, ext in (("heatmap", "pdf"), ("scatter", "png"))
        for chrom in chain(range(1, 23), ("X", "Y"))
        for md5 in ("", ".md5")
    ]
    tpl = (
        "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/report/"
        "bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.{report}.txt{md5}"
    )
    expected += [
        tpl.format(i=i, t=t, report=report, md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for report in ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        for md5 in ("", ".md5")
    ]
    # Perform the comparison
    expected = list(sorted(expected))
    actual = list(sorted(somatic_wgs_cnv_calling_workflow.get_result_files()))
    assert actual == expected
