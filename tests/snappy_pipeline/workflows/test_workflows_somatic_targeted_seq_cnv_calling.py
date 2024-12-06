# -*- coding: utf-8 -*-
"""Tests for the somatic_targeted_seq_cnv_calling workflow module code"""

import textwrap
from itertools import chain

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
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_targeted_seq_cnv_calling:
            tools:
            - cnvetti_on_target
            - copywriter
            - sequenza
            - purecn
            purecn:
              path_container: /path/to/purecn/container
              path_intervals: /path/to/interval/list
              path_panel_of_normals: /path/to/purecn/pon
              path_mapping_bias: /path/to/mapping/bias
              path_somatic_variants: /path/to/somatic/variants
            cnvetti_on_target:
              path_target_regions: /path/to/target/regions
            copywriter:
              path_target_regions: /path/to/target/regions
              plot_genes: "/path/to/civic/annotation??"
            sequenza: {}  # use defaults, no required fields.

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
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticTargetedSeqCnvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variants": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
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
            )()
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
        )()
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
        )()
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
        )()
        assert actual == expected, msg_error


# Tests for SequenzaStepPart ----------------------------------------------------------------------


def test_sequenza_step_part_get_output_files_install(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_output_files() - action 'install'"""
    expected = {"done": "work/R_packages/out/sequenza.done"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("sequenza", "install")
    assert actual == expected


def test_sequenza_step_part_get_log_file_install(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_log_file() - action 'install'"""
    base_name = "work/R_packages/log/sequenza"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("sequenza", "install")
    assert actual == expected


def test_sequenza_step_part_get_output_files_gcreference(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_output_files() - action 'gcreference'"""
    expected = {"gc": "work/static_data/out/sequenza.50.wig.gz"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("sequenza", "gcreference")
    assert actual == expected


def test_sequenza_step_part_get_log_file_gcreference(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_log_file() - action 'gcreference'"""
    base_name = "work/static_data/log/sequenza.50"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("sequenza", "gcreference")
    assert actual == expected


def test_sequenza_step_part_get_input_files_coverage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_input_files() - action 'coverage'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "gc": "work/static_data/out/sequenza.50.wig.gz",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("sequenza", "coverage")(
        wildcards
    )
    assert actual == expected


def test_sequenza_step_part_get_output_files_coverage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_output_files() - action 'coverage'"""
    base_name = "{mapper}.sequenza.{library_name}"
    expected = {
        "seqz": f"work/{base_name}/out/{base_name}.seqz.gz",
        "seqz_md5": f"work/{base_name}/out/{base_name}.seqz.gz.md5",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("sequenza", "coverage")
    assert actual == expected


def test_sequenza_step_part_get_log_file_coverage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_log_file() - action 'coverage'"""
    base_name = (
        "work/{mapper}.sequenza.{library_name}/log/{mapper}.sequenza.{library_name}.coverage"
    )
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("sequenza", "coverage")
    assert actual == expected


def test_sequenza_step_part_get_input_files_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_input_files() - action 'run'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "packages": "work/R_packages/out/sequenza.done",
        "seqz": "work/{mapper}.sequenza.{library_name}/out/{mapper}.sequenza.{library_name}.seqz.gz",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("sequenza", "run")(wildcards)
    assert actual == expected


def test_sequenza_step_part_get_output_files_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_output_files() - action 'run'"""
    base_name = "{mapper}.sequenza.{library_name}"
    expected = {
        "seg": f"work/{base_name}/out/{base_name}_dnacopy.seg",
        "seg_md5": f"work/{base_name}/out/{base_name}_dnacopy.seg.md5",
        "done": f"work/{base_name}/report/.done",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("sequenza", "run")
    assert actual == expected


def test_sequenza_step_part_get_log_file_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_log_file() - action 'run'"""
    base_name = "work/{mapper}.sequenza.{library_name}/log/{mapper}.sequenza.{library_name}.run"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("sequenza", "run")
    assert actual == expected


def test_sequenza_step_part_get_resource_usage_call(somatic_targeted_seq_cnv_calling_workflow):
    """Tests SequenzaStepPart.get_resource_usage()"""
    # Define expected
    expected_dicts = {
        "coverage": {"threads": 1, "time": "24:00:00", "memory": "24G", "partition": "medium"},
        "run": {"threads": 4, "time": "24:00:00", "memory": "64G", "partition": "medium"},
    }
    # Evaluate
    for action, resources in expected_dicts.items():
        for resource, expected in resources.items():
            msg_error = f"Assertion error for resource '{resource}' in '{action}'."
            actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
                "sequenza", action, resource
            )()
            assert actual == expected, msg_error


# Tests for PureCNStepPart ----------------------------------------------------------------------


def test_purecn_step_part_get_input_files_coverage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_input_files() - action 'coverage'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("purecn", "coverage")(
        wildcards
    )
    assert actual == expected


def test_purecn_step_part_get_output_files_coverage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_output_files() - action 'coverage'"""
    name_pattern = "{mapper}.purecn.{library_name}"
    expected = {"coverage": f"work/{name_pattern}/out/{name_pattern}_coverage_loess.txt.gz"}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("purecn", "coverage")
    assert actual == expected


def test_purecn_step_part_get_log_file_coverage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_log_file() - action 'coverage'"""
    base_name = "work/{mapper}.purecn.{library_name}/log/{mapper}.purecn.{library_name}.coverage"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("purecn", "coverage")
    assert actual == expected


def test_purecn_step_part_get_input_files_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_input_files() - action 'run'"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "tumor": "work/bwa.purecn.P001-T1-DNA1-WGS1/out/bwa.purecn.P001-T1-DNA1-WGS1_coverage_loess.txt.gz",
        "vcf": "SOMATIC_VARIANT_CALLING/output/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1.full.vcf.gz",
    }
    actual = somatic_targeted_seq_cnv_calling_workflow.get_input_files("purecn", "run")(wildcards)
    assert actual == expected


def test_purecn_step_part_get_output_files_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_output_files() - action 'run'"""
    base_name = "work/{mapper}.purecn.{library_name}/out/{mapper}.purecn.{library_name}"
    expected = {
        "segments": base_name + "_dnacopy.seg",
        "ploidy": base_name + ".csv",
        "pvalues": base_name + "_amplification_pvalues.csv",
        "vcf": base_name + ".vcf.gz",
        "vcf_tbi": base_name + ".vcf.gz.tbi",
        "loh": base_name + "_loh.csv",
    }
    expected = {**expected, **{k + "_md5": v + ".md5" for k, v in expected.items()}}
    actual = somatic_targeted_seq_cnv_calling_workflow.get_output_files("purecn", "run")
    assert actual == expected


def test_purecn_step_part_get_log_file_run(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_log_file() - action 'run'"""
    base_name = "work/{mapper}.purecn.{library_name}/log/{mapper}.purecn.{library_name}.run"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_targeted_seq_cnv_calling_workflow.get_log_file("purecn", "run")
    assert actual == expected


def test_purecn_step_part_get_resource_usage(somatic_targeted_seq_cnv_calling_workflow):
    """Tests PureCNStepPart.get_resource_usage()"""
    # Define expected
    expected_dicts = {
        "coverage": {"threads": 1, "time": "04:00:00", "memory": "24G", "partition": "medium"},
        "run": {"threads": 4, "time": "24:00:00", "memory": "96G", "partition": "medium"},
    }
    # Evaluate
    for action, resources in expected_dicts.items():
        for resource, expected in resources.items():
            msg_error = f"Assertion error for resource '{resource}' in '{action}'."
            actual = somatic_targeted_seq_cnv_calling_workflow.get_resource(
                "purecn", action, resource
            )()
            assert actual == expected, msg_error


# Tests for SomaticTargetedSeqCnvCallingWorkflow --------------------------------------------------


def test_somatic_targeted_seq_cnv_calling_workflow(somatic_targeted_seq_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "cnvetti_off_target",
        "cnvetti_on_target",
        "copywriter",
        "link_out",
        "purecn",
        "sequenza",
    ]
    actual = list(sorted(somatic_targeted_seq_cnv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    name_pattern = "bwa.{method}.P00{{i}}-T{{t}}-DNA1-WGS1"

    # cnvetti
    tpl = (
        f"output/{name_pattern}/out/{name_pattern}".format(method="cnvetti_on_target_coverage")
        + ".{ext}"
    )
    expected = [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("bcf", "bcf.md5", "bcf.csi", "bcf.csi.md5")
    ]
    tpl = (
        f"output/{name_pattern}/out/{name_pattern}".format(method="cnvetti_on_target_segment")
        + ".{ext}"
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
        f"output/{name_pattern}/out/{name_pattern}".format(method="cnvetti_on_target_postprocess")
        + "_{ext}"
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
    # copywriter
    tpl = f"output/{name_pattern}/out/{name_pattern}".format(method="copywriter") + "_{ext}{md5}"
    expected += [
        tpl.format(i=i, t=t, ext=ext, md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "bins.txt",
            "gene_call.txt",
            "gene_log2.txt",
            "segments.txt",
        )
        for md5 in ("", ".md5")
    ]
    # purecn
    tpl = f"output/{name_pattern}/out/{name_pattern}".format(method="purecn") + "{ext}{md5}"
    expected += [
        tpl.format(i=i, t=t, ext=ext, md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            ".csv",
            ".vcf.gz",
            ".vcf.gz.tbi",
            "_dnacopy.seg",
            "_amplification_pvalues.csv",
            "_loh.csv",
        )
        for md5 in ("", ".md5")
    ]
    # sequenza
    tpl = f"output/{name_pattern}/out/{name_pattern}".format(method="sequenza") + "{ext}{md5}"
    expected += [
        tpl.format(i=i, t=t, ext=ext, md5=md5)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("_dnacopy.seg", ".seqz.gz")
        for md5 in ("", ".md5")
    ]
    tpl = f"output/{name_pattern}/report/.done".format(method="sequenza")
    expected += [tpl.format(i=i, t=t) for i, t in ((1, 1), (2, 1), (2, 2))]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_targeted_seq_cnv_calling_workflow.get_result_files()))
    # HACK TODO
    actual = [f for f in actual if "/log/" not in f]
    assert expected == actual
