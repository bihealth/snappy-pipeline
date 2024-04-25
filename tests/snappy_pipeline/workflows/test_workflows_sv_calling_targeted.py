# -*- coding: utf-8 -*-
"""Tests for the sv_calling_targeted workflow module code"""

import copy
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.base import InvalidConfiguration, UnsupportedActionException
from snappy_pipeline.workflows.sv_calling_targeted import SvCallingTargetedWorkflow

from .common import get_expected_gcnv_log_file, get_expected_output_vcf_files_dict
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
            bwa:
              path_index: /path/to/bwa/index.fa

          sv_calling_targeted:
            path_ngs_mapping: NGS_MAPPING
            tools:
              - delly2
              - manta
              - gcnv
            delly2: {}  # use defaults
            manta: {}   # use defaults
            gcnv:
              # path_uniquely_mapable_bed: /path/to/uniquely/mappable/variable/GRCh37/file.bed.gz
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed
              precomputed_model_paths:
                - library: "Agilent SureSelect Human All Exon V6"
                  contig_ploidy: /path/to/ploidy-model
                  model_pattern: "/data/model_*"

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
def minimal_config_large_cohort(minimal_config):
    """Returns minimum configuration file for large trio cohort."""
    minimal_config_adjusted = copy.deepcopy(minimal_config)
    minimal_config_adjusted["data_sets"]["first_batch"]["file"] = "sheet_large_cohort_trio.tsv"
    return minimal_config_adjusted


@pytest.fixture
def minimal_config_large_cohort_background(minimal_config_large_cohort):
    """Returns minimum configuration file for large trio cohort."""
    minimal_config_adjusted = copy.deepcopy(minimal_config_large_cohort)
    # Set sample sheet as background
    minimal_config_adjusted["data_sets"]["first_batch"]["is_background"] = True
    return minimal_config_adjusted


@pytest.fixture
def minimal_config_large_cohort_without_gcnv_model(minimal_config_large_cohort):
    """Returns minimum configuration file for large trio cohort, no precomputed models defined."""
    minimal_config_adjusted = copy.deepcopy(minimal_config_large_cohort)
    minimal_config_adjusted["step_config"]["sv_calling_targeted"]["gcnv"][
        "precomputed_model_paths"
    ] = []
    return minimal_config_adjusted


@pytest.fixture
def sv_calling_targeted_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2_gcnv_model,
    aligner_indices_fake_fs,
    mocker,
):
    """
    Return SvCallingTargetedWorkflow object pre-configured with germline sheet - small cohort
    """
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2_gcnv_model, mocker
    )
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.common.gcnv.gcnv_run",
        germline_sheet_fake_fs2_gcnv_model,
        mocker,
    )
    # Patch glob with expected model directories
    mocker.patch(
        "snappy_pipeline.workflows.common.gcnv.gcnv_run.glob",
        return_value=["/data/model_01", "/data/model_02", "/data/model_03"],
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SvCallingTargetedWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


@pytest.fixture
def sv_calling_targeted_workflow_large_cohort(
    dummy_workflow,
    minimal_config_large_cohort,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2,
    aligner_indices_fake_fs,
    mocker,
):
    """
    Return SvCallingTargetedWorkflow object pre-configured with germline sheet -
    large trio cohort.
    """
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SvCallingTargetedWorkflow(
        dummy_workflow,
        minimal_config_large_cohort,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


@pytest.fixture
def sv_calling_targeted_workflow_large_cohort_background(
    dummy_workflow,
    minimal_config_large_cohort_background,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SvCallingTargetedWorkflow object pre-configured with germline sheet -
    large trio cohort as background."""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SvCallingTargetedWorkflow(
        dummy_workflow,
        minimal_config_large_cohort_background,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Global tests -------------------------------------------------------------------------------------


def test_validate_request(
    dummy_workflow,
    minimal_config_large_cohort_without_gcnv_model,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2_gcnv_model,
    aligner_indices_fake_fs,
    mocker,
):
    """Tests SvCallingTargetedWorkflow.validate_request()"""
    # Patch out file-system
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2_gcnv_model, mocker
    )
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.common.gcnv.gcnv_run",
        germline_sheet_fake_fs2_gcnv_model,
        mocker,
    )
    # Patch glob with expected model directories
    mocker.patch(
        "snappy_pipeline.workflows.common.gcnv.gcnv_run.glob",
        return_value=["/data/model_01", "/data/model_02", "/data/model_03"],
    )

    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Empty precomputed models, should fail
    with pytest.raises(InvalidConfiguration):
        SvCallingTargetedWorkflow(
            dummy_workflow,
            minimal_config_large_cohort_without_gcnv_model,
            config_lookup_paths,
            config_paths,
            work_dir,
        )


def test_target_seq_cnv_calling_workflow_get_result_files(sv_calling_targeted_workflow):
    """Tests SvCallingTargetedWorkflow.get_result_files()

    Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools, namely: gCNV.
    """
    # Define expected
    pattern_out = "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/out/bwa.{tool}.P00{i}-N1-DNA1-WGS1{ext}"
    expected = [
        pattern_out.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv", "manta", "delly2")
        for ext in (
            ".vcf.gz",
            ".vcf.gz.md5",
            ".vcf.gz.tbi",
            ".vcf.gz.tbi.md5",
        )
    ]
    pattern_log = (
        "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/log/"
        "bwa.{tool}.P00{i}-N1-DNA1-WGS1.{step_name}{ext}"
    )
    expected += [
        pattern_log.format(i=i, tool=tool, step_name=step_name, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool, step_name in (
            ("gcnv", "merge_multikit_families"),
            ("manta", "sv_calling"),
            ("delly2", "sv_calling"),
        )
        for ext in (
            ".log",
            ".log.md5",
            ".conda_info.txt",
            ".conda_info.txt.md5",
            ".conda_list.txt",
            ".conda_list.txt.md5",
            ".wrapper.py",
            ".wrapper.py.md5",
            ".environment.yaml",
            ".environment.yaml.md5",
        )
    ]
    expected = sorted(expected)
    # Get actual
    actual = sorted(sv_calling_targeted_workflow.get_result_files())
    assert actual == expected


def test_target_seq_cnv_calling_workflow_all_donors(
    sv_calling_targeted_workflow, sv_calling_targeted_workflow_large_cohort
):
    """Tests SvCallingTargetedWorkflow.all_donors()"""
    # ----------------------- #
    # Test small sample sheet #
    # ----------------------- #
    # Define expected
    expected = ["P00{i}-N1-DNA1-WGS1".format(i=i) for i in range(1, 7)]
    # Get actual
    actual = sv_calling_targeted_workflow.all_donors()
    assert len(actual) == 6, "Small sample sheet should contain only 6 donors."
    assert_message_tpl = "Value {{donor_name}} not in expected list: {expected}".format(
        expected=", ".join(expected)
    )
    for donor in actual:
        msg = assert_message_tpl.format(donor_name=donor.dna_ngs_library.name)
        assert donor.dna_ngs_library.name in expected, msg
    # ----------------------- #
    # Test large sample sheet #
    # ----------------------- #
    # Define expected
    expected = ["P{i}-N1-DNA1-WGS1".format(i=str(i).zfill(3)) for i in range(1, 502)]
    # Get actual
    actual = sv_calling_targeted_workflow_large_cohort.all_donors()
    assert len(actual) == 501, "Large sample sheet should contain 501 donors."
    assert_message_tpl = (
        "Value {donor_name} not in expected list: P001-N1-DNA1-WGS1 up to P501-N1-DNA1-WGS1."
    )
    for donor in actual:
        msg = assert_message_tpl.format(donor_name=donor.dna_ngs_library.name)
        assert donor.dna_ngs_library.name in expected, msg


def test_target_seq_cnv_calling_workflow_get_library_count(sv_calling_targeted_workflow):
    """Tests SvCallingTargetedWorkflow.get_library_count()"""
    # Test undefined library kit
    expected = 0
    actual = sv_calling_targeted_workflow.get_library_count("_not_a_library_kit_name_")
    assert actual == expected, "It should return zero as the library kit name is not defined."
    # Test defined library kit - foreground
    expected = 6
    actual = sv_calling_targeted_workflow.get_library_count("Agilent_SureSelect_Human_All_Exon_V6")
    assert actual == expected


def test_pick_kits_and_donors(
    sv_calling_targeted_workflow, sv_calling_targeted_workflow_large_cohort
):
    """Tests SvCallingTargetedWorkflow.pick_kits_and_donors()"""

    # Test small cohort - 6 individuals
    expected_library_kits = ["Agilent_SureSelect_Human_All_Exon_V6"]
    expected_kit_counts = {"Agilent_SureSelect_Human_All_Exon_V6": 6}
    library_kits, _, kit_counts = sv_calling_targeted_workflow.pick_kits_and_donors()
    assert library_kits == expected_library_kits
    assert expected_kit_counts == kit_counts

    # Test large trio cohort - 501 individuals
    expected_library_kits = ["Agilent_SureSelect_Human_All_Exon_V6"]
    expected_kit_counts = {"Agilent_SureSelect_Human_All_Exon_V6": 501}
    (
        library_kits,
        _,
        kit_counts,
    ) = sv_calling_targeted_workflow_large_cohort.pick_kits_and_donors()
    assert library_kits == expected_library_kits
    assert expected_kit_counts == kit_counts


# Global RunGcnvTargetSeqStepPart Tests ------------------------------------------------------------


def test_gcnv_call_assertion(sv_calling_targeted_workflow):
    """Tests raise UnsupportedActionException"""
    with pytest.raises(UnsupportedActionException):
        sv_calling_targeted_workflow.get_input_files("gcnv", "_undefined_action_")


def test_gcnv_step_part_get_resource_usage(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_resource()"""
    # Define tested actions
    high_resource_actions = (
        "call_cnvs",
        "post_germline_calls",
        "joint_germline_cnv_segmentation",
    )
    all_actions = sv_calling_targeted_workflow.substep_getattr("gcnv", "actions")
    default_actions = [action for action in all_actions if action not in high_resource_actions]
    # Define expected
    high_res_expected_dict = {
        "threads": 16,
        "time": "4-00:00:00",
        "memory": "46080M",
        "partition": "medium",
    }
    default_expected_dict = {
        "threads": 1,
        "time": "1-00:00:00",
        "memory": "7680M",
        "partition": "medium",
    }
    # Evaluate - high resource actions
    for action in high_resource_actions:
        for resource, expected in high_res_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = sv_calling_targeted_workflow.get_resource("gcnv", action, resource)()
            assert actual == expected, msg_error

    # Evaluate - all other actions
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = sv_calling_targeted_workflow.get_resource("gcnv", action, resource)()
            assert actual == expected, msg_error


def test_gcnv_get_params(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_params for all actions"""
    all_actions = (
        "preprocess_intervals",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "merge_cohort_vcfs",
        "joint_germline_cnv_segmentation",
    )
    actions_w_params = ("call_cnvs", "contig_ploidy", "post_germline_calls")
    for action in all_actions:
        if action in actions_w_params:
            sv_calling_targeted_workflow.get_params("gcnv", action)
        else:
            with pytest.raises(UnsupportedActionException):
                sv_calling_targeted_workflow.get_params("gcnv", action)


def test_gcnv_validate_ploidy_model_directory(
    fake_fs, mocker, sv_calling_targeted_workflow, ploidy_model_files
):
    """Tests RunGcnvTargetSeqStepPart.validate_ploidy_model_directory()"""
    # Create data directories
    fake_fs.fs.makedirs("/data", exist_ok=True)
    fake_fs.fs.makedirs("/empty", exist_ok=True)
    # Create required files
    tpl = "/data/{file_}"
    for file_ in ploidy_model_files:
        fake_fs.fs.create_file(tpl.format(file_=file_))
    # Patch out file-system
    patch_module_fs("snappy_pipeline.workflows.common.gcnv.gcnv_run", fake_fs, mocker)

    # Should return True as it is a directory and it contains the expected files
    assert sv_calling_targeted_workflow.substep_getattr("gcnv", "validate_ploidy_model_directory")(
        "/data"
    )
    # Should return False as empty directory
    assert not sv_calling_targeted_workflow.substep_getattr(
        "gcnv", "validate_ploidy_model_directory"
    )("/empty")
    # Should return False not a directory
    assert not sv_calling_targeted_workflow.substep_getattr(
        "gcnv", "validate_ploidy_model_directory"
    )("__not_a_directory__")


def test_gcnv_validate_call_model_directory(
    fake_fs, mocker, sv_calling_targeted_workflow, call_model_files
):
    """Tests RunGcnvTargetSeqStepPart.validate_call_model_directory()"""
    # Create data directories
    fake_fs.fs.makedirs("/call_data", exist_ok=True)
    fake_fs.fs.makedirs("/empty", exist_ok=True)
    # Create required files
    tpl = "/call_data/{file_}"
    for file_ in call_model_files:
        fake_fs.fs.create_file(tpl.format(file_=file_))
    # Patch out file-system
    patch_module_fs("snappy_pipeline.workflows.common.gcnv.gcnv_run", fake_fs, mocker)

    # Should return True as it is a directory and it contains the expected files
    assert sv_calling_targeted_workflow.substep_getattr("gcnv", "validate_call_model_directory")(
        "/call_data"
    )
    # Should return False as empty directory
    assert not sv_calling_targeted_workflow.substep_getattr(
        "gcnv", "validate_call_model_directory"
    )("/empty")
    # Should return False not a directory
    assert not sv_calling_targeted_workflow.substep_getattr(
        "gcnv", "validate_call_model_directory"
    )("__not_a_directory__")


def test_gcnv_get_result_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_result_files()"""
    expected = [
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/out/bwa.gcnv.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.conda_info.txt",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.conda_info.txt.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.conda_list.txt",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.conda_list.txt.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.wrapper.py",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.wrapper.py.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.environment.yaml",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.environment.yaml.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.log",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.merge_multikit_families.log.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.conda_info.txt",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.conda_info.txt.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.conda_list.txt",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.conda_list.txt.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.wrapper.py",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.wrapper.py.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.environment.yaml",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.environment.yaml.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.log",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.merge_multikit_families.log.md5",
    ]
    actual = sv_calling_targeted_workflow.substep_getattr("gcnv", "get_result_files")()
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (preprocess_intervals) ----------------------------------------


def test_gcnv_preprocess_intervals_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_preprocess_intervals()"""
    # Define expected - empty dictionary for all
    expected = {}
    # Get actual
    actual = sv_calling_targeted_workflow.get_input_files("gcnv", "preprocess_intervals")(None)
    assert actual == expected


def test_gcnv_preprocess_intervals_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_preprocess_intervals()"""
    # Define expected
    output_path = (
        "work/gcnv_preprocess_intervals.{library_kit}/out/"
        "gcnv_preprocess_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("gcnv", "preprocess_intervals")
    assert actual == expected


def test_gcnv_target_step_part_get_log_file(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_log_file for 'preprocess_intervals' step"""
    # Define expected
    expected = (
        "work/gcnv_preprocess_intervals.{library_kit}/log/"
        "gcnv_preprocess_intervals.{library_kit}.log"
    )
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "preprocess_intervals")
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (coverage) ----------------------------------------------------


def test_gcnv_coverage_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_coverage()"""
    # Define expected
    interval_list_out = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    bam_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1"
    expected = {
        "interval_list": interval_list_out,
        "bam": bam_out + ".bam",
        "bai": bam_out + ".bam.bai",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = sv_calling_targeted_workflow.get_input_files("gcnv", "coverage")(wildcards)
    assert actual == expected


def test_gcnv_coverage_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_coverage()"""
    # Define expected
    tsv_out = (
        "work/{mapper}.gcnv_coverage.{library_name}/out/{mapper}.gcnv_coverage.{library_name}.tsv"
    )
    expected = {"tsv": tsv_out}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("gcnv", "coverage")
    assert actual == expected


def test_gcnv_coverage_step_part_get_log_file(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_log_file for 'coverage' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_coverage.{library_name}/log/{mapper}.gcnv_coverage.{library_name}.log"
    )
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "coverage")
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (contig_ploidy) -----------------------------------------------


def test_gcnv_contig_ploidy_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_contig_ploidy()"""
    # Define expected
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    expected = {
        "tsv": tsv_list_out,
        "ped": [
            "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
            "work/write_pedigree.P004-N1-DNA1-WGS1/out/P004-N1-DNA1-WGS1.ped",
        ],
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = sv_calling_targeted_workflow.get_input_files("gcnv", "contig_ploidy")(wildcards)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_contig_ploidy()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_contig_ploidy.{library_kit}/out/"
        "{mapper}.gcnv_contig_ploidy.{library_kit}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("gcnv", "contig_ploidy")
    assert actual == expected


def test_gcnv_get_params_ploidy_model(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_params_ploidy_model()"""
    # Initialise wildcard
    wildcards = Wildcards(fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6"})
    wildcards_fake = Wildcards(fromdict={"library_kit": "__not_a_library_kit__"})
    # Test large cohort - model defined in config
    expected = {"model": "/path/to/ploidy-model"}
    actual = sv_calling_targeted_workflow.get_params("gcnv", "contig_ploidy")(wildcards)
    assert actual == expected
    # Test large cohort - model not defined in config
    expected = {"model": "__no_ploidy_model_for_library_in_config__"}
    actual = sv_calling_targeted_workflow.get_params("gcnv", "contig_ploidy")(wildcards_fake)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_log_file(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_log_file for 'contig_ploidy' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="contig_ploidy")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "contig_ploidy")
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (call_cnvs) ---------------------------------------------------


def test_gcnv_call_cnvs_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_call_cnvs()"""
    # Define expected
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    ploidy_out = (
        "work/bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    expected = {
        "tsv": tsv_list_out,
        "ploidy": ploidy_out,
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = sv_calling_targeted_workflow.get_input_files("gcnv", "call_cnvs")(wildcards)
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_call_cnvs()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/out/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("gcnv", "call_cnvs")
    assert actual == expected


def test_gcnv_get_params_model(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_params_model()"""
    # Initialise wildcard
    wildcards_01 = Wildcards(
        fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6", "shard": "01"}
    )
    wildcards_02 = Wildcards(
        fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6", "shard": "02"}
    )
    wildcards_fake = Wildcards(fromdict={"library_kit": "__not_a_library_kit__"})
    # Test large cohort - model defined in config - shard 01
    expected = {"model": "/data/model_01"}
    actual = sv_calling_targeted_workflow.get_params("gcnv", "call_cnvs")(wildcards_01)
    assert actual == expected
    # Test large cohort - model defined in config - shard 02
    expected = {"model": "/data/model_02"}
    actual = sv_calling_targeted_workflow.get_params("gcnv", "call_cnvs")(wildcards_02)
    assert actual == expected
    # Test large cohort - model not defined in config
    expected = {"model": "__no_model_for_library_in_config__"}
    actual = sv_calling_targeted_workflow.get_params("gcnv", "call_cnvs")(wildcards_fake)
    assert actual == expected


def test_gcnv_call_cnvs_step_part_get_log_file(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.get_log_file for 'call_cnvs' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/log/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}.log"
    )
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "call_cnvs")
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (post_germline_calls) -----------------------------------------


def test_gcnv_post_germline_calls_step_part_get_input_files(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_post_germline_calls()"""
    # Define expected
    call_pattern = (
        "work/bwa.gcnv_call_cnvs.Agilent_SureSelect_Human_All_Exon_V6.0{i}/out/"
        "bwa.gcnv_call_cnvs.Agilent_SureSelect_Human_All_Exon_V6.0{i}/.done"
    )
    call_list_out = [call_pattern.format(i=i) for i in range(1, 4)]  # model 01 - 03
    ploidy_out = (
        "work/bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    expected = {
        "calls": call_list_out,
        "ploidy": ploidy_out,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = sv_calling_targeted_workflow.get_input_files("gcnv", "post_germline_calls")(wildcards)
    assert actual == expected


def test_gcnv_get_params_post_germline_calls(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart._get_params_post_germline_calls()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = {"model": ["/data/model_01", "/data/model_02", "/data/model_03"]}
    actual = sv_calling_targeted_workflow.get_params("gcnv", "post_germline_calls")(wildcards)
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (joint_germline_cnv_segmentation) -----------------------------


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_input_files(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_joint_germline_cnv_segmentation()"""
    # Define expected
    pattern_out = (
        "work/bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1.vcf.gz"
    )
    expected = {
        "interval_list": (
            "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
            "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
        ),
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "vcf": [pattern_out.format(i=i) for i in range(1, 4)],  # P001 - P003
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files(
        "gcnv",
        "joint_germline_cnv_segmentation",
    )(wildcards)
    assert actual == expected


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_output_files(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_joint_germline_cnv_segmentation()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.gcnv_joint_segmentation.{kit}.{library_name}/out/"
        "{mapper}.gcnv_joint_segmentation.{kit}.{library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files(
        "gcnv", "joint_germline_cnv_segmentation"
    )
    assert actual == expected


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_log_file(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart.get_log_file for 'joint_germline_cnv_segmentation' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="gcnv_joint_segmentation", extended=True)
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "joint_germline_cnv_segmentation")
    assert actual == expected


# Tests for RunGcnvTargetSeqStepPart (merge_multikit_families) -------------------------------------


def test_gcnv_merge_multikit_families_step_part_get_input_files(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart._get_input_files_merge_multikit_families()"""
    # Define expected
    expected = {
        "vcf": [
            (
                "work/bwa.gcnv_joint_segmentation.Agilent_SureSelect_Human_All_Exon_V6.P001-N1-DNA1-WGS1/out/"
                "bwa.gcnv_joint_segmentation.Agilent_SureSelect_Human_All_Exon_V6.P001-N1-DNA1-WGS1.vcf.gz"
            )
        ],
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files(
        "gcnv",
        "merge_multikit_families",
    )(wildcards)
    assert actual == expected


def test_gcnv_merge_multikit_families_step_part_get_output_files(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_merge_multikit_families()"""
    # Define expected
    pattern_out = "work/{mapper}.gcnv.{library_name}/out/" "{mapper}.gcnv.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    expected["output_links"] = [
        "output/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}.vcf.gz",
        "output/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}.vcf.gz.md5",
        "output/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}.vcf.gz.tbi",
        "output/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}.vcf.gz.tbi.md5",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.conda_info.txt",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.conda_info.txt.md5",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.conda_list.txt",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.conda_list.txt.md5",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.wrapper.py",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.wrapper.py.md5",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.environment.yaml",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.environment.yaml.md5",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.log",
        "output/{mapper}.gcnv.{library_name}/log/{mapper}.gcnv.{library_name}.merge_multikit_families.log.md5",
    ]
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("gcnv", "merge_multikit_families")
    assert actual == expected


def test_gcnv_merge_multikit_families_step_part_get_log_file(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart.get_log_file for 'merge_multikit_families' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="merge_multikit_families", extended=True)
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "merge_multikit_families")
    assert actual == expected
