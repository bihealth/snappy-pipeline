# -*- coding: utf-8 -*-
"""Tests for the somatic_targeted_seq_cnv_calling workflow module code"""


from itertools import chain
import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    SomaticTargetedSeqCnvCallingWorkflow,
)

from .common import get_expected_log_files_dict
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

          somatic_targeted_seq_cnv_calling:
            tools:
            - cnvetti_on_target
            - cnvkit
            - copywriter

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
    yaml_ = yaml.YAML()
    return yaml_.load(config_str)


@pytest.fixture
def somatic_targeted_seq_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
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
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CnvKitStepPart (access) ---------------------------------------------------------------


def test_cnvkit_access_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "access")
    assert actual is None


def test_cnvkit_access_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = "work/cnvkit.access/out/access.bed"
    assert (
        somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "access") == expected
    )


def test_cnvkit_access_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    expected = get_expected_log_files_dict(base_out="work/cnvkit.access/log/cnvkit.access")
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "access")
    assert actual == expected


def test_cnvkit_access_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_access"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (target) ---------------------------------------------------------------


def test_cnvkit_target_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "target")(
        wildcards
    )
    assert actual == {"access": "work/cnvkit.access/out/access.bed"}


def test_cnvkit_target_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = "work/cnvkit.target/out/target.bed"
    assert (
        somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "target") == expected
    )


def test_cnvkit_target_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    expected = get_expected_log_files_dict(base_out="work/cnvkit.target/log/cnvkit.target")
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "target")
    assert actual == expected


def test_cnvkit_target_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_target"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (antitarget) -----------------------------------------------------------


def test_cnvkit_antitarget_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "access": "work/cnvkit.access/out/access.bed",
        "target": "work/cnvkit.target/out/target.bed",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "antitarget")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_antitarget_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = "work/cnvkit.antitarget/out/antitarget.bed"
    assert (
        somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "antitarget")
        == expected
    )


def test_cnvkit_antitarget_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    expected = get_expected_log_files_dict(base_out="work/cnvkit.antitarget/log/cnvkit.antitarget")
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "antitarget")
    assert actual == expected


def test_cnvkit_antitarget_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_antitarget"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (coverage) -------------------------------------------------------------


def test_cnvkit_coverage_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "target": "target", "library_name": "P001-T1-DNA1-WGS1"}
    )
    expected = {
        "antitarget": "work/cnvkit.antitarget/out/antitarget.bed",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "target": "work/cnvkit.target/out/target.bed",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "coverage")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_coverage_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.coverage.{library_name}/out/{mapper}.cnvkit.coverage.{library_name}"
    )
    expected = {
        "target": base_name_out + ".targetcoverage.cnn",
        "antitarget": base_name_out + ".antitargetcoverage.cnn",
    }
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "coverage")

    assert actual == expected


def test_cnvkit_coverage_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    base_file_name = (
        "work/{mapper}.cnvkit.coverage.{library_name}/log/{mapper}.cnvkit.coverage.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_file_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_coverage_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_coverage"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (reference) ------------------------------------------------------------


def test_cnvkit_reference_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "reference")(
        wildcards
    )
    expected = {
        "antitarget": "work/cnvkit.antitarget/out/antitarget.bed",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "target": "work/cnvkit.target/out/target.bed",
    }
    assert actual == expected


def test_cnvkit_reference_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = (
        "work/{mapper}.cnvkit.reference.{library_name}/out/"
        "{mapper}.cnvkit.reference.{library_name}.cnn"
    )
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "reference")
    assert actual == expected


def test_cnvkit_reference_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.reference.{library_name}/log/{mapper}.cnvkit.reference.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "reference")

    assert actual == expected


def test_cnvkit_reference_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_reference"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (fix) ------------------------------------------------------------------


def test_cnvkit_fix_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    coverage_base_out = (
        "work/bwa.cnvkit.coverage.P001-T1-DNA1-WGS1/out/bwa.cnvkit.coverage.P001-T1-DNA1-WGS1"
    )
    reference_base_out = (
        "work/bwa.cnvkit.reference.P001-T1-DNA1-WGS1/out/bwa.cnvkit.reference.P001-T1-DNA1-WGS1"
    )
    expected = {
        "antitarget": coverage_base_out + ".antitargetcoverage.cnn",
        "ref": reference_base_out + ".cnn",
        "target": coverage_base_out + ".targetcoverage.cnn",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "fix")(wildcards)
    assert actual == expected


def test_cnvkit_fix_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = "work/{mapper}.cnvkit.fix.{library_name}/out/{mapper}.cnvkit.fix.{library_name}.cnr"
    assert somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "fix") == expected


def test_cnvkit_fix_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.cnvkit.fix.{library_name}/log/{mapper}.cnvkit.fix.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "fix")

    assert actual == expected


def test_cnvkit_fix_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_fix"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (segment) --------------------------------------------------------------


def test_cnvkit_segment_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "cnr": "work/bwa.cnvkit.fix.P001-T1-DNA1-WGS1/out/bwa.cnvkit.fix.P001-T1-DNA1-WGS1.cnr"
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "segment")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_segment_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = (
        "work/{mapper}.cnvkit.segment.{library_name}/out/{mapper}.cnvkit.segment.{library_name}.cns"
    )
    assert (
        somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "segment") == expected
    )


def test_cnvkit_segment_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.segment.{library_name}/log/{mapper}.cnvkit.segment.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "segment")

    assert actual == expected


def test_cnvkit_segment_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_segment"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (call) -----------------------------------------------------------------


def test_cnvkit_call_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    segment_file = (
        "work/bwa.cnvkit.segment.P001-T1-DNA1-WGS1/out/bwa.cnvkit.segment.P001-T1-DNA1-WGS1.cns"
    )
    expected = {"segment": segment_file}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "call")(wildcards)
    assert actual == expected


def test_cnvkit_call_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    expected = (
        "work/{mapper}.cnvkit.call.{library_name}/out/{mapper}.cnvkit.call.{library_name}.cns"
    )
    assert somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "call") == expected


def test_cnvkit_call_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.call.{library_name}/log/{mapper}.cnvkit.call.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "call")
    assert actual == expected


def test_cnvkit_call_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_call"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (plot) -----------------------------------------------------------------


def test_cnvkit_plot_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    cnr_file = "work/bwa.cnvkit.fix.P001-T1-DNA1-WGS1/out/bwa.cnvkit.fix.P001-T1-DNA1-WGS1.cnr"
    cns_file = (
        "work/bwa.cnvkit.segment.P001-T1-DNA1-WGS1/out/bwa.cnvkit.segment.P001-T1-DNA1-WGS1.cns"
    )
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
    base_out = "work/{mapper}.cnvkit.plot.{library_name}/out/{mapper}.cnvkit.plot.{library_name}"
    expected = {
        "diagram": base_out + ".diagram.pdf",
        "heatmap": base_out + ".heatmap.pdf",
        "scatter": base_out + ".scatter.pdf",
    }
    for chrom in chain(range(1, 23), ("X", "Y")):
        for diagram in ("heatmap", "scatter"):
            tpl = base_out + ".%s.chr%s.pdf"
            expected["{}_chr{}".format(diagram, chrom)] = tpl % (diagram, chrom)

    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "plot")
    assert actual == expected


def test_cnvkit_plot_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.cnvkit.plot.{library_name}/log/{mapper}.cnvkit.plot.{library_name}"
    )
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "plot")
    assert actual == expected


def test_cnvkit_plot_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_plot"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (export) ---------------------------------------------------------------


def test_cnvkit_export_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "cns": "work/bwa.cnvkit.call.P001-T1-DNA1-WGS1/out/bwa.cnvkit.call.P001-T1-DNA1-WGS1.cns"
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "export")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_export_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.export.{library_name}/out/{mapper}.cnvkit.export.{library_name}"
    )
    expected = {
        "bed": base_name_out + ".bed",
        "seg": base_name_out + ".seg",
        "tbi": base_name_out + ".vcf.gz.tbi",
        "vcf": base_name_out + ".vcf.gz",
    }
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "export")
    assert actual == expected


def test_cnvkit_export_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.export.{library_name}/log/{mapper}.cnvkit.export.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "export")
    assert actual == expected


def test_cnvkit_export_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_export"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for CnvKitStepPart (report) ---------------------------------------------------------------


def test_cnvkit_report_step_part_get_input_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    cnr_file = "work/bwa.cnvkit.fix.P001-T1-DNA1-WGS1/out/bwa.cnvkit.fix.P001-T1-DNA1-WGS1.cnr"
    cns_file = (
        "work/bwa.cnvkit.segment.P001-T1-DNA1-WGS1/out/bwa.cnvkit.segment.P001-T1-DNA1-WGS1.cns"
    )
    expected = {
        "cnr": cnr_file,
        "cns": cns_file,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("cnvkit", "report")(
        wildcards
    )
    assert actual == expected


def test_cnvkit_report_step_part_get_output_files(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.report.{library_name}/out/{mapper}.cnvkit.report.{library_name}"
    )
    expected = {
        "breaks": base_name_out + ".breaks.txt",
        "gainloss": base_name_out + ".gainloss.txt",
        "gender": base_name_out + ".gender.txt",
        "metrics": base_name_out + ".metrics.txt",
        "segmetrics": base_name_out + ".segmetrics.txt",
    }
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("cnvkit", "report")
    assert actual == expected


def test_cnvkit_report_step_part_get_log_file(somatic_targeted_seq_cnv_calling_workflow):
    # Define expected
    base_name_out = (
        "work/{mapper}.cnvkit.report.{library_name}/log/{mapper}.cnvkit.report.{library_name}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out)
    # Get actual
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("cnvkit", "report")
    assert actual == expected


def test_cnvkit_report_step_part_update_cluster_config(
    somatic_targeted_seq_cnv_calling_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["somatic_targeted_seq_cnv_calling_cnvkit_report"].keys())
    expected = {"mem", "time", "ntasks"}
    assert actual == expected


# Tests for SomaticTargetedSeqCnvCallingWorkflow --------------------------------------------------


def test_somatic_targeted_seq_cnv_calling_workflow(somatic_targeted_seq_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = ["cnvetti_off_target", "cnvetti_on_target", "cnvkit", "copywriter", "link_out"]
    assert list(sorted(somatic_targeted_seq_cnv_calling_workflow.sub_steps.keys())) == expected
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
    tpl = (
        "output/bwa.cnvkit.call.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvkit.call.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext) for i, t in ((1, 1), (2, 1), (2, 2)) for ext in ("cns",)
    ]
    tpl = (
        "output/bwa.cnvkit.export.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvkit.export.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("bed", "seg", "vcf.gz", "vcf.gz.tbi")
    ]
    tpl = (
        "output/bwa.cnvkit.plot.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvkit.plot.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("diagram.pdf", "heatmap.pdf", "scatter.pdf")
    ]
    tpl = (
        "output/bwa.cnvkit.plot.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvkit.plot.P00{i}-T{t}-DNA1-WGS1.{diagram}.chr{chrom}.pdf"
    )
    expected += [
        tpl.format(i=i, t=t, diagram=diagram, chrom=chrom)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for diagram in ("heatmap", "scatter")
        for chrom in chain(range(1, 23), ("X", "Y"))
    ]
    tpl = (
        "output/bwa.cnvkit.report.P00{i}-T{t}-DNA1-WGS1/out/"
        "bwa.cnvkit.report.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("breaks.txt", "gainloss.txt", "gender.txt", "metrics.txt", "segmetrics.txt")
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
