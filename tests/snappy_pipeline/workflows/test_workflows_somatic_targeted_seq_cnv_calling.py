# -*- coding: utf-8 -*-
"""Tests for the somatic_targeted_seq_cnv_calling workflow module code"""


from itertools import chain
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    SomaticTargetedSeqCnvCallingWorkflow,
)

from .common import get_expected_log_files_dict, get_expected_output_bcf_files_dict
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
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_targeted_seq_cnv_calling:
            tools:
            - cnvetti_on_target
            - cnvkit
            - copywriter
            cnvkit:
              path_target: /path/to/panel_of_normals/output/cnvkit.target/out/cnvkit.target.bed
              path_antitarget: /path/to/panel_of_normals/output/cnvkit.antitarget/out/cnvkit.antitarget.bed
              path_panel_of_normals: /path/to/panel_of_normals/output/bwa.cnvkit.create_panel/out/bwa.cnvkit.panel_of_normals.cnn

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
def somatic_targeted_seq_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticTargetedSeqCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SomaticTargetedSeqCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CnvettiOnTargetStepPart  ---------------------------------------------------------------


def test_cnvetti_on_target_step_part_get_input_files_coverage(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files(
        "cnvetti_on_target", "coverage"
    )(wildcards)
    assert actual == expected


def test_cnvetti_on_target_step_part_get_input_files_segment(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart._get_input_files_segment()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    base_name = (
        "work/bwa.cnvetti_on_target_coverage.P001-T1-DNA1-WGS1/out/"
        "bwa.cnvetti_on_target_coverage.P001-T1-DNA1-WGS1"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files(
        "cnvetti_on_target", "segment"
    )(wildcards)
    assert actual == expected


def test_cnvetti_on_target_step_part_get_input_files_postprocess(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart._get_input_files_postprocess()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    base_name = (
        "work/bwa.cnvetti_on_target_segment.P001-T1-DNA1-WGS1/out/"
        "bwa.cnvetti_on_target_segment.P001-T1-DNA1-WGS1"
    )
    expected = {
        "targets_bcf": base_name + ".targets.bcf",
        "targets_bcf_md5": base_name + ".targets.bcf.md5",
        "targets_bcf_csi": base_name + ".targets.bcf.csi",
        "targets_bcf_csi_md5": base_name + ".targets.bcf.csi.md5",
        "segments_bcf": base_name + ".segments.bcf",
        "segments_bcf_md5": base_name + ".segments.bcf.md5",
        "segments_bcf_csi": base_name + ".segments.bcf.csi",
        "segments_bcf_csi_md5": base_name + ".segments.bcf.csi.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files(
        "cnvetti_on_target", "postprocess"
    )(wildcards)
    assert actual == expected


def test_cnvetti_on_target_step_part_get_output_files_coverage(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart._get_output_files_coverage()"""
    base_name = (
        "work/{mapper}.cnvetti_on_target_coverage.{library_name}/out/"
        "{mapper}.cnvetti_on_target_coverage.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files(
        "cnvetti_on_target", "coverage"
    )
    assert actual == expected


def test_cnvetti_on_target_step_part_get_output_files_segment(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart._get_output_files_segment()"""
    base_name = (
        "work/{mapper}.cnvetti_on_target_segment.{library_name}/out/"
        "{mapper}.cnvetti_on_target_segment.{library_name}"
    )
    expected = {
        "targets_bcf": base_name + ".targets.bcf",
        "targets_bcf_md5": base_name + ".targets.bcf.md5",
        "targets_bcf_csi": base_name + ".targets.bcf.csi",
        "targets_bcf_csi_md5": base_name + ".targets.bcf.csi.md5",
        "segments_bcf": base_name + ".segments.bcf",
        "segments_bcf_md5": base_name + ".segments.bcf.md5",
        "segments_bcf_csi": base_name + ".segments.bcf.csi",
        "segments_bcf_csi_md5": base_name + ".segments.bcf.csi.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files(
        "cnvetti_on_target", "segment"
    )
    assert actual == expected


def test_cnvetti_on_target_step_part_get_output_files_postprocess(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart._get_output_files_postprocess()"""
    base_name = (
        "work/{mapper}.cnvetti_on_target_postprocess.{library_name}/out/"
        "{mapper}.cnvetti_on_target_postprocess.{library_name}"
    )
    expected = {
        "targets_txt": base_name + "_targets.txt",
        "targets_md5": base_name + "_targets.txt.md5",
        "targets_segmented_txt": base_name + "_targets_segmented.txt",
        "targets_segmented_md5": base_name + "_targets_segmented.txt.md5",
        "segments_txt": base_name + "_segments.txt",
        "segments_md5": base_name + "_segments.txt.md5",
        "gene_call_txt": base_name + "_gene_call.txt",
        "gene_call_md5": base_name + "_gene_call.txt.md5",
        "gene_log2_txt": base_name + "_gene_log2.txt",
        "gene_log2_md5": base_name + "_gene_log2.txt.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files(
        "cnvetti_on_target", "postprocess"
    )
    assert actual == expected


def test_cnvetti_on_target_step_part_get_log_file_coverage(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart.get_log_file() - action 'coverage'"""
    base_name = (
        "work/{mapper}.cnvetti_on_target_coverage.{library_name}/log/"
        "{mapper}.cnvetti_on_target_coverage.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvetti_on_target", "coverage")
    assert actual == expected


def test_cnvetti_on_target_step_part_get_log_file_segment(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart.get_log_file() - action 'segment'"""
    base_name = (
        "work/{mapper}.cnvetti_on_target_segment.{library_name}/log/"
        "{mapper}.cnvetti_on_target_segment.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvetti_on_target", "segment")
    assert actual == expected


def test_cnvetti_on_target_step_part_get_log_file_postprocess(
    somatic_targeted_seq_cnv_calling_workflow,
):
    """Tests CnvettiOnTargetStepPart.get_log_file() - action 'segment'"""
    base_name = (
        "work/{mapper}.cnvetti_on_target_postprocess.{library_name}/log/"
        "{mapper}.cnvetti_on_target_postprocess.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file(
        "cnvetti_on_target", "postprocess"
    )
    assert actual == expected


def test_cnvetti_on_target_step_part_get_resource_usage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvettiOnTargetStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "7500M", "partition": "medium"}
    # Evaluate
    all_actions = somatic_targeted_seq_cnv_calling_workflow.substep_getattr(
        "cnvetti_on_target", "actions"
    )
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action {action}."
            actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
                "cnvetti_on_target", action, resource
            )
            assert actual == expected, msg_error


# Tests for CnvKitStepPart (coverage) -------------------------------------------------------------


def test_cnvkit_coverage_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "target": "target", "library_name": "P001-T1-DNA1-WGS1"}
    )
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "coverage")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_coverage_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}"
    expected = {
        "target": base_name_out + ".targetcoverage.cnn",
        "target_md5": base_name_out + ".targetcoverage.cnn.md5",
        "antitarget": base_name_out + ".antitargetcoverage.cnn",
        "antitarget_md5": base_name_out + ".antitargetcoverage.cnn.md5",
    }
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "coverage")

    assert actual == expected


def test_cnvkit_coverage_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    base_file_name = (
        "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.coverage.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_file_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_coverage_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'coverage'"""
    # Define expected
    expected_dict = {"threads": 8, "time": "08:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "cnvkit", "coverage", resource
        )
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (fix) ------------------------------------------------------------------


def test_cnvkit_fix_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    coverage_base_out = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1"
    expected = {
        "antitarget": coverage_base_out + ".antitargetcoverage.cnn",
        "target": coverage_base_out + ".targetcoverage.cnn",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "fix")(wildcards)
    assert actual == expected


def test_cnvkit_fix_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cnr"
    expected = {"ratios": base_name_out, "ratios_md5": base_name_out + ".md5"}
    assert somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "fix") == expected


def test_cnvkit_fix_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.fix.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "fix")
    assert actual == expected


def test_cnvkit_fix_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'fix'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "03:59:59", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource("cnvkit", "fix", resource)
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (segment) --------------------------------------------------------------


def test_cnvkit_segment_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {"cnr": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "segment")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_segment_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    base_name_out = (
        "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.segment.cns"
    )
    expected = {"segments": base_name_out, "segments_md5": base_name_out + ".md5"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "segment")
    assert actual == expected


def test_cnvkit_segment_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.segment.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "segment")
    assert actual == expected


def test_cnvkit_segment_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'fix'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "03:59:59", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "cnvkit", "segment", resource
        )
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (call) -----------------------------------------------------------------


def test_cnvkit_call_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    segment_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.segment.cns"
    expected = {"segment": segment_file}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "call")(wildcards)
    assert actual == expected


def test_cnvkit_call_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cns"
    expected = {"calls": base_name_out, "calls_md5": base_name_out + ".md5"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "call")
    assert actual == expected


def test_cnvkit_call_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.call.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "call")
    assert actual == expected


def test_cnvkit_call_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "03:59:59", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource("cnvkit", "call", resource)
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (plot) -----------------------------------------------------------------


def test_cnvkit_plot_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    cnr_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"
    cns_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns"
    expected = {
        "cnr": cnr_file,
        "cns": cns_file,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "plot")(wildcards)
    assert actual == expected


def test_cnvkit_plot_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
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
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "plot")
    assert actual == expected


def test_cnvkit_plot_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.plot.{library_name}"
    )
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "plot")
    assert actual == expected


def test_cnvkit_plot_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "08:00:00", "memory": "30720M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource("cnvkit", "plot", resource)
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (export) ---------------------------------------------------------------


def test_cnvkit_export_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {"cns": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "export")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_export_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    expected = {}
    base_name_out = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}"
    for key, ext in (("bed", "bed"), ("seg", "seg"), ("vcf", "vcf.gz"), ("tbi", "vcf.gz.tbi")):
        expected[key] = base_name_out + "." + ext
        expected[key + "_md5"] = expected[key] + ".md5"
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "export")
    assert actual == expected


def test_cnvkit_export_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.export.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "export")
    assert actual == expected


def test_cnvkit_export_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "03:59:59", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "cnvkit", "export", resource
        )
        assert actual == expected, msg_error


# Tests for CnvKitStepPart (report) ---------------------------------------------------------------


def test_cnvkit_report_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    cnr_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"
    cns_file = "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns"
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
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "report")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_report_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    expected = {}
    base_name_out = "work/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}"
    for report in ("breaks", "genemetrics", "segmetrics", "sex", "metrics"):
        expected[report] = base_name_out + "." + report + ".txt"
        expected[report + "_md5"] = expected[report] + ".md5"
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "report")
    assert actual == expected


def test_cnvkit_report_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.report.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "report")
    assert actual == expected


def test_cnvkit_report_step_part_get_resource(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CnvKitStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "03:59:59", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "cnvkit", "report", resource
        )
        assert actual == expected, msg_error


# Tests for CopywriterStepPart   -------------------------------------------------------------------


def test_copywriter_step_part_get_input_files_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_input_files() - action 'run'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("copywriter", "run")(
        wildcards
    )
    assert actual == expected


def test_copywriter_step_part_get_input_files_call(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_input_files() - action 'call'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "input": "work/bwa.copywriter.P001-T1-DNA1-WGS1/CNAprofiles/input.Rdata",
        "segment": "work/bwa.copywriter.P001-T1-DNA1-WGS1/CNAprofiles/segment.Rdata",
        "counts": "work/bwa.copywriter.P001-T1-DNA1-WGS1/CNAprofiles/read_counts.txt",
        "log2": "work/bwa.copywriter.P001-T1-DNA1-WGS1/CNAprofiles/log2_read_counts.igv",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("copywriter", "call")(
        wildcards
    )
    assert actual == expected


def test_copywriter_step_part_get_output_files_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_output_files() - action 'run'"""
    expected = {
        "input": "work/{mapper}.copywriter.{library_name}/CNAprofiles/input.Rdata",
        "input_md5": "work/{mapper}.copywriter.{library_name}/CNAprofiles/input.Rdata.md5",
        "segment": "work/{mapper}.copywriter.{library_name}/CNAprofiles/segment.Rdata",
        "segment_md5": "work/{mapper}.copywriter.{library_name}/CNAprofiles/segment.Rdata.md5",
        "counts": "work/{mapper}.copywriter.{library_name}/CNAprofiles/read_counts.txt",
        "counts_md5": "work/{mapper}.copywriter.{library_name}/CNAprofiles/read_counts.txt.md5",
        "log2": "work/{mapper}.copywriter.{library_name}/CNAprofiles/log2_read_counts.igv",
        "log2_md5": "work/{mapper}.copywriter.{library_name}/CNAprofiles/log2_read_counts.igv.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("copywriter", "run")
    assert actual == expected


def test_copywriter_step_part_get_output_files_call(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_output_files() - action 'call'"""
    base_name = "work/{mapper}.copywriter.{library_name}/out/{mapper}.copywriter.{library_name}"
    expected = {
        "bins_txt": base_name + "_bins.txt",
        "bins_txt_md5": base_name + "_bins.txt.md5",
        "gene_call_txt": base_name + "_gene_call.txt",
        "gene_call_txt_md5": base_name + "_gene_call.txt.md5",
        "gene_log2_txt": base_name + "_gene_log2.txt",
        "gene_log2_txt_md5": base_name + "_gene_log2.txt.md5",
        "segments_txt": base_name + "_segments.txt",
        "segments_txt_md5": base_name + "_segments.txt.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("copywriter", "call")
    assert actual == expected


def test_copywriter_step_part_get_log_file_prepare(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_log_file() - action 'prepare'"""
    expected = {
        "log": "work/copywriter.prepare/log/snakemake.log",
        "log_md5": "work/copywriter.prepare/log/snakemake.log.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("copywriter", "prepare")
    assert actual == expected


def test_copywriter_step_part_get_log_file_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_log_file() - action 'run'"""
    base_name = "work/{mapper}.copywriter.{library_name}/log/{mapper}.copywriter.{library_name}.run"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("copywriter", "run")
    assert actual == expected


def test_copywriter_step_part_get_log_file_call(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_log_file() - action 'run'"""
    base_name = (
        "work/{mapper}.copywriter.{library_name}/log/{mapper}.copywriter.{library_name}.call"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("copywriter", "call")
    assert actual == expected


def test_copywriter_step_part_get_resource_usage_prepare(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_resource_usage() - action 'prepare'"""
    # Define expected
    expected_dict = {"threads": 1, "time": "02:00:00", "memory": "4000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "copywriter", "prepare", resource
        )
        assert actual == expected, msg_error


def test_copywriter_step_part_get_resource_usage_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_resource_usage() - action 'run'"""
    # Define expected
    expected_dict = {"threads": 2, "time": "16:00:00", "memory": "80000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "copywriter", "run", resource
        )
        assert actual == expected, msg_error


def test_copywriter_step_part_get_resource_usage_call(somatic_targeted_seq_cnv_calling_workflow):
    """Tests CopywriterStepPart.get_resource_usage() - action 'call'"""
    # Define expected
    expected_dict = {"threads": 8, "time": "03:59:00", "memory": "8000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
            "copywriter", "call", resource
        )
        assert actual == expected, msg_error


# Tests for SomaticTargetedSeqCnvCallingWorkflow --------------------------------------------------


def test_somatic_targeted_seq_cnv_calling_workflow(somatic_targeted_seq_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnvetti_off_target", "cnvetti_on_target", "cnvkit", "copywriter", "link_out"]
    actual = list(sorted(somatic_targeted_seq_cnv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    # cnvetti
    tpl = (
        "output/bwa.cnvetti_on_target_coverage.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvetti_on_target_coverage.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("bcf", "bcf.md5", "bcf.csi", "bcf.csi.md5")
    ]
    tpl = (
        "output/bwa.cnvetti_on_target_segment.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvetti_on_target_segment.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "segments.bcf",
            "segments.bcf.csi",
            "segments.bcf.csi.md5",
            "segments.bcf.md5",
            "targets.bcf",
            "targets.bcf.csi",
            "targets.bcf.csi.md5",
            "targets.bcf.md5",
        )
    ]
    tpl = (
        "output/bwa.cnvetti_on_target_postprocess.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvetti_on_target_postprocess.P00{i}-T{t}-DNA1-WGS1_{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "gene_call.txt",
            "gene_call.txt.md5",
            "gene_log2.txt",
            "gene_log2.txt.md5",
            "segments.txt",
            "segments.txt.md5",
            "targets.txt",
            "targets.txt.md5",
            "targets_segmented.txt",
            "targets_segmented.txt.md5",
        )
    ]
    # cnvkit
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
    # copywriter
    tpl = (
        "output/bwa.copywriter.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.copywriter.P00{i}-T{t}-DNA1-WGS1_{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "bins.txt",
            "gene_call.txt",
            "gene_log2.txt",
            "segments.txt",
            "bins.txt.md5",
            "gene_call.txt.md5",
            "gene_log2.txt.md5",
            "segments.txt.md5",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_targeted_seq_cnv_calling_workflow.get_result_files()))
    # HACK TODO
    actual = [f for f in actual if "/log/" not in f]
    assert expected == actual
