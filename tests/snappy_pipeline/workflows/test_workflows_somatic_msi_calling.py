# -*- coding: utf-8 -*-
"""Tests for the somatic_msi_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_msi_calling import SomaticMsiCallingWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs

__author__ = "Clemens Messerschmidt"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fasta
          somatic_msi_calling:
            tools: ["mantis_msi2"]
            path_ngs_mapping: ../ngs_mapping  # REQUIRED
            loci_bed: /path/to/hg19/loci.bed  # REQUIRED

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
def somatic_msi_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticMsiCallingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SomaticMsiCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for Mantis_msi2 ------------------------------------------------------------------


def test_mantis_msi2_step_part_get_input_files(somatic_msi_calling_workflow):
    """Tests mantis_msi2.get_input_files()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_msi_calling_workflow.get_input_files("mantis_msi2", "run")(wildcards)
    assert actual == expected


def test_mantis_msi2_step_part_get_output_files(somatic_msi_calling_workflow):
    """Test mantis_msi2.get_output_files"""
    base_out = "work/{mapper}.mantis_msi2.{tumor_library}/out/{mapper}.mantis_msi2.{tumor_library}"
    expected = {
        "result": base_out + ".results.txt",
        "status": base_out + ".results.txt.status",
        "result_md5": base_out + ".results.txt.md5",
        "status_md5": base_out + ".results.txt.status.md5",
    }
    actual = somatic_msi_calling_workflow.get_output_files("mantis_msi2", "run")
    assert actual == expected


def test_mantis_msi2_step_part_get_log_files(somatic_msi_calling_workflow):
    """Tests mantis_msi2.get_log_files()"""
    base_out = (
        "work/{mapper}.mantis_msi2.{tumor_library}/log/" "{mapper}.mantis_msi2.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_msi_calling_workflow.get_log_file("mantis_msi2", "run")
    assert actual == expected


def test_mantis_msi2_step_part_get_resource_usage(
    somatic_msi_calling_workflow,
):
    """Tests mantis_msi2.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 3, "time": "02:00:00", "memory": f"{30 * 1024 * 3}M"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_msi_calling_workflow.get_resource("mantis_msi2", "run", resource)
        assert actual == expected, msg_error


# Tests for somatic msi calling workflow------------------------------------------------------------------


def test_somatic_msi_calling_workflow(somatic_msi_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "mantis_msi2"]
    assert set(expected).issubset(list(sorted(somatic_msi_calling_workflow.sub_steps.keys())))

    # Check result file construction
    tpl = (
        "output/{mapper}.{msi_caller}.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "{mapper}.{msi_caller}.P00{i}-T{t}-DNA1-WGS1{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, msi_caller=msi_caller, i=i, t=t, ext=ext, dir_="out")
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            ".results.txt",
            ".results.txt.status",
            ".results.txt.md5",
            ".results.txt.status.md5",
        )
        for mapper in ("bwa",)
        for msi_caller in ("mantis_msi2",)
    ]
    expected += [
        tpl.format(mapper=mapper, msi_caller=msi_caller, i=i, t=t, ext=ext, dir_="log")
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            ".conda_info.txt",
            ".conda_list.txt",
            ".log",
            ".conda_info.txt.md5",
            ".conda_list.txt.md5",
            ".log.md5",
        )
        for mapper in ("bwa",)
        for msi_caller in ("mantis_msi2",)
    ]
    expected = list(sorted(expected))

    actual = list(sorted(somatic_msi_calling_workflow.get_result_files()))
    assert expected == actual
