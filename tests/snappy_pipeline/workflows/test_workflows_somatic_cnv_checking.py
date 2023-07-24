# -*- coding: utf-8 -*-
"""Tests for the somatic_cnv_checking workflow module code"""


import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_cnv_checking import (
    SomaticCnvCheckingWorkflow,
)

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
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

          somatic_targeted_seq_cnv_calling:
            tools: ["cnvkit"]

          somatic_cnv_checking:
            path_ngs_mapping: ../ngs_mapping
            path_cnv_calling: ../somatic_targeted_seq_cnv_calling
            cnv_calling_tool: cnvkit

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
def somatic_cnv_checking_workflow(
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
    return SomaticCnvCheckingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CnvCheckingPileupStepPart --------------------------------------------------------------


def test_pileup_normal_step_part_get_input_files(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_input_files() - action normal"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library": "P001-T1-DNA1-WGS1"})
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_cnv_checking_workflow.get_input_files("pileup", "normal")(wildcards)
    assert actual == expected


def test_pileup_normal_step_part_get_output_files(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_output_files() - action normal"""
    base_name = "work/{mapper}.{library}/out/{mapper}.{library}.normal"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_output_files("pileup", "normal")
    assert actual == expected


def test_pileup_normal_step_part_get_log_file(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_log_file() - action normal"""
    base_name = "work/{mapper}.{library}/log/{mapper}.{library}.normal"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_log_file("pileup", "normal")
    assert actual == expected


def test_pileup_normal_step_part_get_resource(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_resource() - action normal"""
    expected_dict = {"threads": 2, "time": "12:00:00", "memory": "7577M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_cnv_checking_workflow.get_resource("pileup", "normal", resource)
        assert actual == expected, msg_error


def test_pileup_tumor_step_part_get_input_files(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_input_files() - action tumor"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library": "P001-T1-DNA1-WGS1"})
    expected = {
        "vcf": "work/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.normal.vcf.gz",
        "tbi": "work/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.normal.vcf.gz.tbi",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_cnv_checking_workflow.get_input_files("pileup", "tumor")(wildcards)
    assert actual == expected


def test_pileup_tumor_step_part_get_output_files(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_output_files() - action tumor"""
    base_name = "work/{mapper}.{library}/out/{mapper}.{library}.tumor"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_output_files("pileup", "tumor")
    assert actual == expected


def test_pileup_tumor_step_part_get_log_file(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_log_file() - action tumor"""
    base_name = "work/{mapper}.{library}/log/{mapper}.{library}.tumor"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_log_file("pileup", "tumor")
    assert actual == expected


def test_pileup_tumor_step_part_get_resource(somatic_cnv_checking_workflow):
    """Tests CnvCheckingPileupStepPart.get_resource() - action tumor"""
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7577M", "partition": "medium"}
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_cnv_checking_workflow.get_resource("pileup", "tumor", resource)
        assert actual == expected, msg_error


def test_cnv_run_step_part_get_input_files(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingCnvStepPart.get_input_files()"""
    expected = {
        "cnv": "../somatic_targeted_seq_cnv_calling/output/{mapper}.cnvkit.{library}/out/{mapper}.cnvkit.{library}.bed.gz",
        "cnv_tbi": "../somatic_targeted_seq_cnv_calling/output/{mapper}.cnvkit.{library}/out/{mapper}.cnvkit.{library}.bed.gz.tbi",
        "normal": "work/{mapper}.{library}/out/{mapper}.{library}.normal.vcf.gz",
        "normal_tbi": "work/{mapper}.{library}/out/{mapper}.{library}.normal.vcf.gz.tbi",
        "tumor": "work/{mapper}.{library}/out/{mapper}.{library}.tumor.vcf.gz",
        "tumor_tbi": "work/{mapper}.{library}/out/{mapper}.{library}.tumor.vcf.gz.tbi",
    }
    actual = somatic_cnv_checking_workflow.get_input_files("cnv", "run")
    assert actual == expected


def test_cnv_run_step_part_get_output_files(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingCnvStepPart.get_output_files()"""
    base_name = "work/{mapper}.cnvkit.{library}/out/{mapper}.cnvkit.{library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_output_files("cnv", "run")
    assert actual == expected


def test_cnv_run_step_part_get_log_file(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingCnvStepPart.get_log_file()"""
    base_name = "work/{mapper}.cnvkit.{library}/log/{mapper}.cnvkit.{library}"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_log_file("cnv", "run")
    assert actual == expected


def test_report_run_step_part_get_input_files(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingReportStepPart.get_input_files()"""
    expected = {"vcf": "work/{mapper}.cnvkit.{library}/out/{mapper}.cnvkit.{library}.vcf.gz"}
    actual = somatic_cnv_checking_workflow.get_input_files("report", "run")
    assert actual == expected


def test_report_run_step_part_get_output_files(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingReportStepPart.get_output_files()"""
    base_name = "work/{mapper}.cnvkit.{library}/report/{mapper}.cnvkit.{library}"
    expected = {
        "cnv": base_name + ".cnv.pdf",
        "cnv_md5": base_name + ".cnv.pdf.md5",
        "locus": base_name + ".locus.pdf",
        "locus_md5": base_name + ".locus.pdf.md5",
    }
    actual = somatic_cnv_checking_workflow.get_output_files("report", "run")
    assert actual == expected


def test_report_run_step_part_get_log_file(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingReportStepPart.get_log_file()"""
    base_name = "work/{mapper}.cnvkit.{library}/log/{mapper}.cnvkit.{library}.report"
    expected = get_expected_log_files_dict(base_out=base_name)
    actual = somatic_cnv_checking_workflow.get_log_file("report", "run")
    assert actual == expected


def test_report_run_step_part_get_params(somatic_cnv_checking_workflow):
    """Tests SomaticCnvCheckingReportStepPart.get_params()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library": "P001-T1-DNA1-WGS1"})
    expected = {"library_name": "P001-T1-DNA1-WGS1"}
    actual = somatic_cnv_checking_workflow.get_params("report", "run")(wildcards)
    assert actual == expected


# Tests for SomaticCnvCheckingWorkflow --------------------------------------------------


def test_somatic_cnv_checking_workflow(somatic_cnv_checking_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnv", "link_out", "pileup", "report"]
    actual = list(sorted(somatic_cnv_checking_workflow.sub_steps.keys()))
    assert actual == expected

    # main output
    tpl = "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/out/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.{ext}"
    expected = [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "vcf.gz",
            "vcf.gz.tbi",
            "vcf.gz.md5",
            "vcf.gz.tbi.md5",
        )
    ]

    # report
    tpl = "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/report/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.{ext}"
    expected += [
        tpl.format(i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "cnv.pdf",
            "cnv.pdf.md5",
            "locus.pdf",
            "locus.pdf.md5",
        )
    ]

    # logs
    tpl = "output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/log/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1{part}.{ext}{chksum}"
    expected += [
        tpl.format(i=i, t=t, part=part, ext=ext, chksum=chksum)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for part in ("", ".normal", ".tumor", ".report")
        for ext in ("log", "conda_info.txt", "conda_list.txt")
        for chksum in ("", ".md5")
    ]

    expected = list(sorted(expected))
    actual = list(sorted(somatic_cnv_checking_workflow.get_result_files()))
    assert expected == actual
