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

# List of valid actions - XHMM
XHMM_ACTIONS = [
    "coverage",
    "merge_cov",
    "ref_stats",
    "filter_center",
    "pca",
    "normalize",
    "zscore_center",
    "refilter",
    "discover",
    "genotype",
    "extract_ped",
]


def get_expected_xhmm_log_file(step_name):
    """
    :param step_name: Step name.
    :type step_name: str

    :return: Returns expected log file path for basic steps in XHMM.
    """
    expected_log = (
        "work/{mapper}.xhmm_" + step_name + ".{library_kit}/log/snakemake.sv_calling_targeted.log"
    )
    return expected_log


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

          sv_calling_targeted:
            tools:
              - xhmm
              - gcnv
            xhmm:
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed
            gcnv:
              path_target_interval_list_mapping:
                - pattern: "Agilent SureSelect Human All Exon V6.*"
                  name: "Agilent_SureSelect_Human_All_Exon_V6"
                  path: /path/to/Agilent/SureSelect_Human_All_Exon_V6_r2/GRCh37/Exons.bed
                  path_uniquely_mapable_bed: /path/to/uniquely/mappable/variable/GRCh37/file.bed.gz
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
    mocker,
):
    """
    Return SvCallingTargetedWorkflow object pre-configured with germline sheet - small cohort
    """
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2_gcnv_model, mocker
    )
    patch_module_fs(
        "snappy_pipeline.workflows.gcnv.gcnv_run",
        germline_sheet_fake_fs2_gcnv_model,
        mocker,
    )
    # Patch glob with expected model directories
    mocker.patch(
        "snappy_pipeline.workflows.gcnv.gcnv_run.glob",
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
    mocker,
):
    """
    Return SvCallingTargetedWorkflow object pre-configured with germline sheet -
    large trio cohort.
    """
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
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
    mocker,
):
    """Return SvCallingTargetedWorkflow object pre-configured with germline sheet -
    large trio cohort as background."""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
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
    mocker,
):
    """Tests SvCallingTargetedWorkflow.validate_request()"""
    # Patch out file-system
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2_gcnv_model, mocker
    )
    patch_module_fs(
        "snappy_pipeline.workflows.gcnv.gcnv_run",
        germline_sheet_fake_fs2_gcnv_model,
        mocker,
    )
    # Patch glob with expected model directories
    mocker.patch(
        "snappy_pipeline.workflows.gcnv.gcnv_run.glob",
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
    to the expected results from the tools, namely: gCNV and XHMM.
    """
    # Define expected
    pattern_out = "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/out/bwa.{tool}.P00{i}-N1-DNA1-WGS1{ext}"
    expected = [
        pattern_out.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv", "xhmm")
        for ext in (
            ".vcf.gz",
            ".vcf.gz.md5",
            ".vcf.gz.tbi",
            ".vcf.gz.tbi.md5",
        )
    ]
    pattern_log = (
        "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/log/"
        "bwa.{tool}.P00{i}-N1-DNA1-WGS1.joint_germline_segmentation{ext}"
    )
    expected += [
        pattern_log.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv",)
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


def test_target_seq_cnv_calling_workflow_all_background_donors(
    sv_calling_targeted_workflow, sv_calling_targeted_workflow_large_cohort_background
):
    """Tests SvCallingTargetedWorkflow.all_background_donors()"""
    # Test small foreground sample sheet
    actual = sv_calling_targeted_workflow.all_background_donors()
    assert len(actual) == 0, "Small sample sheet should contain zero background donors."
    # Test large background sample sheet
    actual = sv_calling_targeted_workflow_large_cohort_background.all_background_donors()
    assert len(actual) == 501, "Large sample sheet should contain 501 background donors."


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
    )
    all_actions = sv_calling_targeted_workflow.substep_getattr("gcnv", "actions")
    default_actions = [action for action in all_actions if action not in high_resource_actions]
    # Define expected
    high_res_expected_dict = {
        "threads": 16,
        "time": "2-00:00:00",
        "memory": "46080M",
        "partition": "medium",
    }
    default_expected_dict = {
        "threads": 1,
        "time": "04:00:00",
        "memory": "7680M",
        "partition": "medium",
    }
    # Evaluate - high resource actions
    for action in high_resource_actions:
        for resource, expected in high_res_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = sv_calling_targeted_workflow.get_resource("gcnv", action, resource)
            assert actual == expected, msg_error

    # Evaluate - all other actions
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = sv_calling_targeted_workflow.get_resource("gcnv", action, resource)
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


def test_gcnv_validate_precomputed_model_paths_config(sv_calling_targeted_workflow):
    """Tests RunGcnvTargetSeqStepPart.validate_model_requirements()"""
    # Initialise input
    valid_dict = {
        "library": "library",
        "contig_ploidy": "/path/to/ploidy-model",
        "model_pattern": "/path/to/model_*",
    }
    typo_dict = {
        "library_n": "library",
        "contig_ploidy": "/path/to/ploidy-model",
        "model_pattern": "/path/to/model_*",
    }
    missing_key_dict = {"model_pattern": "/path/to/model_*"}

    # Sanity check
    sv_calling_targeted_workflow.substep_getattr("gcnv", "validate_precomputed_model_paths_config")(
        config=[valid_dict]
    )
    # Test key typo
    with pytest.raises(InvalidConfiguration):
        sv_calling_targeted_workflow.substep_getattr(
            "gcnv", "validate_precomputed_model_paths_config"
        )(config=[valid_dict, typo_dict])
    # Test key missing
    with pytest.raises(InvalidConfiguration):
        sv_calling_targeted_workflow.substep_getattr(
            "gcnv", "validate_precomputed_model_paths_config"
        )(config=[valid_dict, missing_key_dict])


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
    patch_module_fs("snappy_pipeline.workflows.gcnv.gcnv_run", fake_fs, mocker)

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
    patch_module_fs("snappy_pipeline.workflows.gcnv.gcnv_run", fake_fs, mocker)

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
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml.md5",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.log",
        "output/bwa.gcnv.P001-N1-DNA1-WGS1/log/bwa.gcnv.P001-N1-DNA1-WGS1.joint_germline_segmentation.log.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/out/bwa.gcnv.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_info.txt.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.conda_list.txt.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.wrapper.py.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.environment.yaml.md5",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.log",
        "output/bwa.gcnv.P004-N1-DNA1-WGS1/log/bwa.gcnv.P004-N1-DNA1-WGS1.joint_germline_segmentation.log.md5",
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
    expected = {"tsv": tsv_list_out}
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
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = sv_calling_targeted_workflow.get_input_files(
        "gcnv", "joint_germline_cnv_segmentation"
    )(wildcards)
    assert actual == expected


def test_gcnv_joint_germline_cnv_segmentation_step_part_get_output_files(
    sv_calling_targeted_workflow,
):
    """Tests RunGcnvTargetSeqStepPart._get_output_files_joint_germline_cnv_segmentation()"""
    # Define expected
    pattern_out = "work/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}"
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
    expected = get_expected_gcnv_log_file(step_name="gcnv", extended=True)
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("gcnv", "joint_germline_cnv_segmentation")
    assert actual == expected


# Global XhmmStepPart Tests ------------------------------------------------------------------------


def test_xhmm_call_assertion(sv_calling_targeted_workflow):
    """Tests raise UnsupportedActionException"""
    with pytest.raises(UnsupportedActionException):
        sv_calling_targeted_workflow.get_input_files("xhmm", "_undefined_action_")


def test_xhmm_get_params(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_params for all actions"""
    for action in XHMM_ACTIONS:
        if action == "coverage":
            sv_calling_targeted_workflow.get_params("xhmm", action)
        else:
            with pytest.raises(UnsupportedActionException):
                sv_calling_targeted_workflow.get_params("xhmm", action)


# Tests for XhmmStepPart (coverage) ----------------------------------------------------------------


def test_xhmm_coverage_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_coverage()"""
    # Define expected
    bam_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1"
    expected = {"bam": bam_out + ".bam", "bai": bam_out + ".bam.bai"}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "coverage")(wildcards)
    assert actual == expected


def test_xhmm_coverage_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_coverage()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.xhmm_coverage.{library_name}/out/{mapper}.xhmm_coverage.{library_name}.DATA"
    )
    expected = {
        "sample_interval_statistics": pattern_out + ".sample_interval_statistics",
        "sample_interval_summary": pattern_out + ".sample_interval_summary",
        "sample_statistics": pattern_out + ".sample_statistics",
        "sample_summary": pattern_out + ".sample_summary",
    }
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "coverage")
    assert actual == expected


def test_xhmm_coverage_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'coverage' step"""
    # Define expected
    expected = "work/{mapper}.xhmm_coverage.{library_name}/log/snakemake.sv_calling_targeted.log"
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "coverage")
    assert actual == expected


# Tests for XhmmStepPart (merge_cov) ---------------------------------------------------------------


def test_xhmm_merge_cov_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_merge_cov()"""
    # Define expected
    base_out = (
        "work/bwa.xhmm_coverage.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.xhmm_coverage.P00{i}-N1-DNA1-WGS1.DATA.sample_interval_summary"
    )
    expected = [base_out.format(i=i) for i in range(1, 7)]  # P001 - P006
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "merge_cov")(wildcards)
    assert actual == expected


def test_xhmm_merge_cov_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_merge_cov()"""
    # Define expected
    name_out = (
        "work/{mapper}.xhmm_merge_cov.{library_kit}/out/"
        "{mapper}.xhmm_merge_cov.{library_kit}.RD.txt"
    )
    expected = [name_out]
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "merge_cov")
    assert actual == expected


def test_xhmm_merge_cov_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'merge_cov' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="merge_cov")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "merge_cov")
    assert actual == expected


# Tests for XhmmStepPart (ref_stats) ---------------------------------------------------------------


def test_xhmm_ref_stats_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_ref_stats()"""
    # Define expected
    name_out = (
        "work/{mapper}.xhmm_ref_stats.{library_kit}/out/"
        "{mapper}.xhmm_ref_stats.{library_kit}.extreme_gc_targets.txt"
    )
    expected = {"extreme_gc_targets": name_out}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "ref_stats")
    assert actual == expected


def test_xhmm_ref_stats_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'ref_stats' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="ref_stats")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "ref_stats")
    assert actual == expected


# Tests for XhmmStepPart (filter_center) -----------------------------------------------------------


def test_xhmm_filter_center_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_filter_center()"""
    # Define expected
    base_out = (
        "work/bwa.{step}.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.{step}.Agilent_SureSelect_Human_All_Exon_V6.{type}.txt"
    )
    expected = {
        "merge_cov": base_out.format(step="xhmm_merge_cov", type="RD"),
        "extreme_gc": base_out.format(step="xhmm_ref_stats", type="extreme_gc_targets"),
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "filter_center")(wildcards)
    assert actual == expected


def test_xhmm_filter_center_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_filter_center()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.xhmm_filter_center.{library_kit}/out/"
        "{mapper}.xhmm_filter_center.{library_kit}"
    )
    expected = {
        "centered": pattern_out + ".centered.txt",
        "filtered_targets": pattern_out + ".filtered_targets.txt",
        "filtered_samples": pattern_out + ".filtered_samples.txt",
    }
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "filter_center")
    assert actual == expected


def test_xhmm_filter_center_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'filter_center' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="filter_center")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "filter_center")
    assert actual == expected


# Tests for XhmmStepPart (pca) ---------------------------------------------------------------------


def test_xhmm_pca_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_pca()"""
    # Define expected
    base_out = (
        "work/bwa.xhmm_filter_center.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_filter_center.Agilent_SureSelect_Human_All_Exon_V6.centered.txt"
    )
    expected = [base_out]
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "pca")(wildcards)
    assert actual == expected


def test_xhmm_pca_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_pca()"""
    # Define expected
    pattern_out = "work/{mapper}.xhmm_pca.{library_kit}/out/{mapper}.xhmm_pca.{library_kit}"
    expected = {
        "pc_loading": pattern_out + ".PC_LOADINGS.txt",
        "pc_sd": pattern_out + ".PC_SD.txt",
        "pc": pattern_out + ".PC.txt",
    }
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "pca")
    assert actual == expected


def test_xhmm_pca_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'pca' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="pca")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "pca")
    assert actual == expected


# Tests for XhmmStepPart (normalize) ---------------------------------------------------------------


def test_xhmm_normalize_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_normalize()"""
    # Define expected
    centered_out = (
        "work/bwa.xhmm_filter_center.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_filter_center.Agilent_SureSelect_Human_All_Exon_V6.centered.txt"
    )
    pca_out = (
        "work/bwa.xhmm_pca.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_pca.Agilent_SureSelect_Human_All_Exon_V6.PC.txt"
    )
    expected = {"centered": centered_out, "pca": pca_out}
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "normalize")(wildcards)
    assert actual == expected


def test_xhmm_normalize_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_normalize()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.xhmm_normalize.{library_kit}/out/{mapper}.xhmm_normalize.{library_kit}"
    )
    expected = {"normalized": pattern_out, "num_removed": pattern_out + ".num_removed_PC.txt"}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "normalize")
    assert actual == expected


def test_xhmm_normalize_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'normalize' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="normalize")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "normalize")
    assert actual == expected


# Tests for XhmmStepPart (zscore_center) -----------------------------------------------------------


def test_xhmm_zscore_center_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_zscore_center()"""
    # Define expected
    base_out = (
        "work/bwa.xhmm_normalize.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_normalize.Agilent_SureSelect_Human_All_Exon_V6"
    )
    expected = [base_out]
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "zscore_center")(wildcards)
    assert actual == expected


def test_xhmm_zscore_center_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_zscore_center()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.xhmm_zscore_center.{library_kit}/out/"
        "{mapper}.xhmm_zscore_center.{library_kit}"
    )
    expected = {
        "zscore_center": pattern_out + "",
        "filtered_samples": pattern_out + ".filtered_samples.txt",
        "filtered_targets": pattern_out + ".filtered_targets.txt",
    }
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "zscore_center")
    assert actual == expected


def test_xhmm_zscore_center_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'zscore_center' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="zscore_center")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "zscore_center")
    assert actual == expected


# Tests for XhmmStepPart (zscore_center) -----------------------------------------------------------


def test_xhmm_refilter_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_refilter()"""
    # Define expected
    original_out = (
        "work/bwa.xhmm_merge_cov.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_merge_cov.Agilent_SureSelect_Human_All_Exon_V6.RD.txt"
    )
    pattern_out = (
        "work/bwa.xhmm_{step}_center.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_{step}_center.Agilent_SureSelect_Human_All_Exon_V6.filtered_{used}.txt"
    )
    expected = {
        "original": original_out,
        "filtered_samples_filter_center": pattern_out.format(step="filter", used="samples"),
        "filtered_targets_filter_center": pattern_out.format(step="filter", used="targets"),
        "filtered_samples_zscore_center": pattern_out.format(step="zscore", used="samples"),
        "filtered_targets_zscore_center": pattern_out.format(step="zscore", used="targets"),
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "refilter")(wildcards)
    assert actual == expected


def test_xhmm_refilter_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_refilter()"""
    # Define expected
    base_out = (
        "work/{mapper}.xhmm_refilter.{library_kit}/out/{mapper}.xhmm_refilter.{library_kit}.RD.txt"
    )
    expected = [base_out]
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "refilter")
    assert actual == expected


def test_xhmm_refilter_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'refilter' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="refilter")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "refilter")
    assert actual == expected


# Tests for XhmmStepPart (discover) ----------------------------------------------------------------


def test_xhmm_discover_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_discover()"""
    # Define expected
    center_zscore_out = (
        "work/bwa.xhmm_zscore_center.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_zscore_center.Agilent_SureSelect_Human_All_Exon_V6"
    )
    refilter_original_out = (
        "work/bwa.xhmm_refilter.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_refilter.Agilent_SureSelect_Human_All_Exon_V6.RD.txt"
    )
    expected = {"center_zscore": center_zscore_out, "refilter_original": refilter_original_out}
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "discover")(wildcards)
    assert actual == expected


def test_xhmm_discover_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_discover()"""
    # Define expected
    base_out = "work/{mapper}.xhmm_discover.{library_kit}/out/{mapper}.xhmm_discover.{library_kit}"
    expected = {"xcnv": base_out + ".xcnv", "aux_xcnv": base_out + ".aux_xcnv"}
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "discover")
    assert actual == expected


def test_xhmm_discover_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'discover' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="discover")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "discover")
    assert actual == expected


# Tests for XhmmStepPart (genotype) ----------------------------------------------------------------


def test_xhmm_genotype_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_genotype()"""
    # Define expected
    pattern_out = (
        "work/bwa.xhmm_{action}.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_{action}.Agilent_SureSelect_Human_All_Exon_V6{ext}"
    )
    expected = {
        "center_zscore": pattern_out.format(action="zscore_center", ext=""),
        "refilter_original": pattern_out.format(action="refilter", ext=".RD.txt"),
        "discover_xcnv": pattern_out.format(action="discover", ext=".xcnv"),
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_kit": "Agilent_SureSelect_Human_All_Exon_V6",
        }
    )
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "genotype")(wildcards)
    assert actual == expected


def test_xhmm_genotype_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_genotype()"""
    # Define expected
    base_out = "work/{mapper}.xhmm_genotype.{library_kit}/out/{mapper}.xhmm_genotype.{library_kit}"
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "genotype")
    assert actual == expected


def test_xhmm_genotype_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'genotype' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="genotype")
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "genotype")
    assert actual == expected


# Tests for XhmmStepPart (extract_ped) -------------------------------------------------------------


def test_xhmm_extract_ped_step_part_get_input_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_input_files_extract_ped()"""
    # Define expected
    filtered_samples_out = (
        "work/bwa.xhmm_filter_center.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_filter_center.Agilent_SureSelect_Human_All_Exon_V6.filtered_samples.txt"
    )
    vcf_pattern_out = (
        "work/bwa.xhmm_genotype.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.xhmm_genotype.Agilent_SureSelect_Human_All_Exon_V6.vcf.gz"
    )
    expected = {
        "filtered_samples": filtered_samples_out,
        "vcf": vcf_pattern_out,
        "vcf_tbi": vcf_pattern_out + ".tbi",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = sv_calling_targeted_workflow.get_input_files("xhmm", "extract_ped")(wildcards)
    assert actual == expected


def test_xhmm_extract_ped_step_part_get_output_files(sv_calling_targeted_workflow):
    """Tests XhmmStepPart._get_output_files_genotype()"""
    # Define expected
    base_out = "work/{mapper}.xhmm.{library_name}/out/{mapper}.xhmm.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    # Get actual
    actual = sv_calling_targeted_workflow.get_output_files("xhmm", "extract_ped")
    assert actual == expected


def test_xhmm_extract_ped_part_get_log_file(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_log_file for 'genotype' step"""
    # Define expected
    expected = "work/{mapper}.xhmm.{library_name}/log/snakemake.sv_calling_targeted.log"
    # Get actual
    actual = sv_calling_targeted_workflow.get_log_file("xhmm", "extract_ped")
    assert actual == expected


def test_xhmm_step_part_get_resource_usage(sv_calling_targeted_workflow):
    """Tests XhmmStepPart.get_resource()"""
    # Define expected
    merge_cov_expected_dict = {
        "threads": 1,
        "time": "1-00:00:00",
        "memory": "12G",
        "partition": "medium",
    }
    default_expected_dict = {
        "threads": 1,
        "time": "08:00:00",
        "memory": "12G",
        "partition": "medium",
    }

    # Evaluate - merge_cov
    for resource, expected in merge_cov_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = sv_calling_targeted_workflow.get_resource("xhmm", "merge_cov", resource)
        assert actual == expected, msg_error

    # Evaluate - all other actions
    all_actions = sv_calling_targeted_workflow.substep_getattr("xhmm", "actions")
    default_actions = [action for action in all_actions if action != "merge_cov"]
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = sv_calling_targeted_workflow.get_resource("xhmm", action, resource)
            assert actual == expected, msg_error
