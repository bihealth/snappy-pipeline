# -*- coding: utf-8 -*-
"""Tests for the somatic_wgs_cnv_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_wgs_sv_calling import SomaticWgsSvCallingWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


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

          somatic_wgs_sv_calling:
              path_ngs_mapping: ../ngs_mapping
              tools: ['manta']

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
def somatic_wgs_sv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticWgsSvCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)

    # Construct the workflow object
    return SomaticWgsSvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for MantaStepPart --------------------------------------------------------------------------


def test_manta_somatic_step_part_get_input_files(somatic_wgs_sv_calling_workflow):
    """Tests MantaStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "cancer_library": "P001-T1-DNA1-WGS1"})
    expected = {
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = somatic_wgs_sv_calling_workflow.get_input_files("manta", "run")(wildcards)
    assert actual == expected


def test_manta_somatic_step_part_get_output_files(somatic_wgs_sv_calling_workflow):
    """Tests MantaStepPart.get_output_files()"""
    base_name = "work/{mapper}.manta.{cancer_library}/out/{mapper}.manta.{cancer_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name)
    actual = somatic_wgs_sv_calling_workflow.get_output_files("manta", "run")
    assert actual == expected


def test_manta_somatic_step_part_get_log_file(somatic_wgs_sv_calling_workflow):
    """Tests MantaStepPart.get_log_file()"""
    expected = "work/{mapper}.manta.{cancer_library}/log/snakemake.somatic_wgs_sv_calling.log"
    actual = somatic_wgs_sv_calling_workflow.get_log_file("manta", "run")
    assert actual == expected


def test_manta_somatic_step_part_get_resource(somatic_wgs_sv_calling_workflow):
    """Tests MantaStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "1-16:00:00", "memory": "61440M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_wgs_sv_calling_workflow.get_resource("manta", "run", resource)
        assert actual == expected, msg_error


# Tests for Delly2StepPart -------------------------------------------------------------------------


def test_delly2_somatic_step_part_get_resource(somatic_wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_resource()"""
    # Get all available actions
    all_actions = somatic_wgs_sv_calling_workflow.substep_getattr("delly2", "actions")
    # Define expected
    expected_dict = {"threads": 4, "time": "6-06:00:00", "memory": "15360M", "partition": "medium"}
    # Evaluate
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = somatic_wgs_sv_calling_workflow.get_resource("delly2", action, resource)
            assert actual == expected, msg_error


# Tests for SomaticWgsSvCallingWorkflow ------------------------------------------------------------


def test_somatic_wgs_sv_calling_workflow(somatic_wgs_sv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["delly2", "link_out", "manta"]
    actual = list(sorted(somatic_wgs_sv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{cnv_caller}.P00{i}-T{t}-DNA1-WGS1/out/"
        "{mapper}.{cnv_caller}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for cnv_caller in ("manta",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_wgs_sv_calling_workflow.get_result_files()))
    assert actual == expected
