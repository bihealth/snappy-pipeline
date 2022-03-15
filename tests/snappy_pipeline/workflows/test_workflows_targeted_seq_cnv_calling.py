# -*- coding: utf-8 -*-
"""Tests for the targeted_seq_cnv_calling workflow module code"""

import copy
import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.base import InvalidConfiguration, UnsupportedActionException
from snappy_pipeline.workflows.targeted_seq_cnv_calling import TargetedSeqCnvCallingWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

# List of valid actions - gCNV
GCNV_ACTIONS = [
    "preprocess_intervals",
    "annotate_gc",
    "filter_intervals",
    "scatter_intervals",
    "coverage",
    "contig_ploidy",
    "call_cnvs_cohort_mode",
    "post_germline_calls_cohort_mode",
    "merge_cohort_vcfs",
    "extract_ped",
]

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

# Ploidy model required files
PLOIDY_MODEL_FILES = [
    "contig_ploidy_prior.tsv",
    "gcnvkernel_version.json",
    "interval_list.tsv",
    "mu_mean_bias_j_lowerbound__.tsv",
    "mu_psi_j_log__.tsv",
    "ploidy_config.json",
    "std_mean_bias_j_lowerbound__.tsv",
    "std_psi_j_log__.tsv",
]

# Call model required files
CALL_MODEL_FILES = [
    "calling_config.json",
    "gcnvkernel_version.json",
    "log_q_tau_tk.tsv",
    "mu_ard_u_log__.tsv",
    "mu_psi_t_log__.tsv",
    "std_ard_u_log__.tsv",
    "std_psi_t_log__.tsv",
    "denoising_config.json",
    "interval_list.tsv",
    "mu_W_tu.tsv",
    "mu_log_mean_bias_t.tsv",
    "std_W_tu.tsv",
    "std_log_mean_bias_t.tsv",
]


def get_expected_xhmm_log_file(step_name):
    """
    :param step_name: Step name.
    :type step_name: str

    :return: Returns expected log file path for basic steps in XHMM.
    """
    expected_log = (
        "work/{mapper}.xhmm_"
        + step_name
        + ".{library_kit}/log/snakemake.targeted_seq_cnv_calling.log"
    )
    return expected_log


def get_expected_gcnv_log_file(step_name):
    """
    :param step_name: Step name.
    :type step_name: str

    :return: Returns expected log file path for basic steps in gCNV.
    """
    expected_log = (
        "work/{mapper}.gcnv_"
        + step_name
        + ".{library_kit}/log/{mapper}.gcnv_"
        + step_name
        + ".{library_kit}.log"
    )
    return expected_log


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
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          targeted_seq_cnv_calling:
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
def minimal_config_large_cohort_with_gcnv_model(minimal_config_large_cohort):
    """Returns minimum configuration file for large trio cohort. Path defined in fixture
    ``germline_sheet_fake_fs2_gcnv_model``.
    """
    minimal_config_adjusted = copy.deepcopy(minimal_config_large_cohort)
    minimal_config_adjusted["step_config"]["targeted_seq_cnv_calling"]["gcnv"][
        "precomputed_model_paths"
    ] = [
        {
            "library": "Agilent SureSelect Human All Exon V6",
            "contig_ploidy": "/path/to/ploidy-model",
            "model_pattern": "/data/model_*",
        }
    ]
    return minimal_config_adjusted


@pytest.fixture
def germline_sheet_fake_fs2_gcnv_model(germline_sheet_fake_fs2):
    """Return fake file system setup with files for the germline_sheet_tsv and gCNV files."""
    # Create contig-ploidy model
    ploidy_dir = "/path/to/ploidy-model"
    germline_sheet_fake_fs2.fs.makedirs(ploidy_dir, exist_ok=True)
    # Create required files
    tpl = ploidy_dir + "/{file_}"
    for file_ in PLOIDY_MODEL_FILES:
        germline_sheet_fake_fs2.fs.create_file(tpl.format(file_=file_))

    # Create model directories
    for model_n in ("01", "02", "03"):
        model_path = "/data/model_{0}".format(model_n)
        germline_sheet_fake_fs2.fs.makedirs(model_path, exist_ok=True)
        # Create required files
        tpl = model_path + "/{file_}"
        for file_ in CALL_MODEL_FILES:
            germline_sheet_fake_fs2.fs.create_file(tpl.format(file_=file_))
    return germline_sheet_fake_fs2


@pytest.fixture
def targeted_seq_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """
    Return TargetedSeqCnvCallingWorkflow object pre-configured with germline sheet - small cohort
    """
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return TargetedSeqCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


@pytest.fixture
def targeted_seq_cnv_calling_workflow_large_cohort(
    dummy_workflow,
    minimal_config_large_cohort,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2,
    mocker,
):
    """
    Return TargetedSeqCnvCallingWorkflow object pre-configured with germline sheet -
    large trio cohort.
    """
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return TargetedSeqCnvCallingWorkflow(
        dummy_workflow,
        minimal_config_large_cohort,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


@pytest.fixture
def targeted_seq_cnv_calling_workflow_large_cohort_background(
    dummy_workflow,
    minimal_config_large_cohort_background,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2,
    mocker,
):
    """Return TargetedSeqCnvCallingWorkflow object pre-configured with germline sheet -
    large trio cohort as background."""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return TargetedSeqCnvCallingWorkflow(
        dummy_workflow,
        minimal_config_large_cohort_background,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


@pytest.fixture
def targeted_seq_cnv_calling_workflow_with_gcnv_model(
    dummy_workflow,
    minimal_config_large_cohort_with_gcnv_model,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs2_gcnv_model,
    mocker,
):
    """Return TargetedSeqCnvCallingWorkflow object pre-configured with germline sheet -
    large trio cohort as background."""
    # Patch out file-system
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2_gcnv_model, mocker
    )
    patch_module_fs(
        "snappy_pipeline.workflows.targeted_seq_cnv_calling",
        germline_sheet_fake_fs2_gcnv_model,
        mocker,
    )
    # Patch glob.glob with expected model directories
    mocker.patch(
        "snappy_pipeline.workflows.targeted_seq_cnv_calling.glob.glob",
        return_value=["/data/model_01", "/data/model_02", "/data/model_03"],
    )

    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return TargetedSeqCnvCallingWorkflow(
        dummy_workflow,
        minimal_config_large_cohort_with_gcnv_model,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Global tests -------------------------------------------------------------------------------------


def test_target_seq_cnv_calling_workflow_files(targeted_seq_cnv_calling_workflow):
    """Tests TargetedSeqCnvCallingWorkflow.get_result_files()

    Tests simple functionality of the workflow: checks if file structure is created according
    to the expected results from the tools, namely: gCNV and XHMM.
    """
    # Define expected
    pattern_out = "output/bwa.{tool}.P00{i}-N1-DNA1-WGS1/out/bwa.{tool}.P00{i}-N1-DNA1-WGS1.{ext}"
    expected = [
        pattern_out.format(i=i, tool=tool, ext=ext)
        for i in (1, 4)  # only index: P001, P004
        for tool in ("gcnv", "xhmm")
        for ext in (
            "vcf.gz",
            "vcf.gz.md5",
            "vcf.gz.tbi",
            "vcf.gz.tbi.md5",
        )
    ]
    expected = sorted(expected)
    # Get actual
    actual = sorted(targeted_seq_cnv_calling_workflow.get_result_files())
    assert actual == expected


def test_target_seq_cnv_calling_workflow_all_donors(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_large_cohort
):
    """Tests TargetedSeqCnvCallingWorkflow.all_donors()"""
    # ----------------------- #
    # Test small sample sheet #
    # ----------------------- #
    # Define expected
    expected = ["P00{i}-N1-DNA1-WGS1".format(i=i) for i in range(1, 7)]
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.all_donors()
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
    actual = targeted_seq_cnv_calling_workflow_large_cohort.all_donors()
    assert len(actual) == 501, "Large sample sheet should contain 501 donors."
    assert_message_tpl = (
        "Value {donor_name} not in expected list: P001-N1-DNA1-WGS1 up to P501-N1-DNA1-WGS1."
    )
    for donor in actual:
        msg = assert_message_tpl.format(donor_name=donor.dna_ngs_library.name)
        assert donor.dna_ngs_library.name in expected, msg


def test_target_seq_cnv_calling_workflow_all_background_donors(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_large_cohort_background
):
    """Tests TargetedSeqCnvCallingWorkflow.all_background_donors()"""
    # Test small foreground sample sheet
    actual = targeted_seq_cnv_calling_workflow.all_background_donors()
    assert len(actual) == 0, "Small sample sheet should contain zero background donors."
    # Test large background sample sheet
    actual = targeted_seq_cnv_calling_workflow_large_cohort_background.all_background_donors()
    assert len(actual) == 501, "Large sample sheet should contain 501 background donors."


def test_target_seq_cnv_calling_workflow_get_library_count(targeted_seq_cnv_calling_workflow):
    """Tests TargetedSeqCnvCallingWorkflow.get_library_count()"""
    # Test undefined library kit
    expected = 0
    actual = targeted_seq_cnv_calling_workflow.get_library_count("_not_a_library_kit_name_")
    assert actual == expected, "It should return zero as the library kit name is not defined."
    # Test defined library kit - foreground
    expected = 6
    actual = targeted_seq_cnv_calling_workflow.get_library_count(
        "Agilent_SureSelect_Human_All_Exon_V6"
    )
    assert actual == expected


def test_pick_kits_and_donors(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_large_cohort
):
    """Tests TargetedSeqCnvCallingWorkflow.pick_kits_and_donors()"""

    # Test small cohort - 6 individuals
    expected_library_kits = ["Agilent_SureSelect_Human_All_Exon_V6"]
    expected_kit_counts = {"Agilent_SureSelect_Human_All_Exon_V6": 6}
    library_kits, _, kit_counts = targeted_seq_cnv_calling_workflow.pick_kits_and_donors()
    assert library_kits == expected_library_kits
    assert expected_kit_counts == kit_counts

    # Test large trio cohort - 501 individuals
    expected_library_kits = ["Agilent_SureSelect_Human_All_Exon_V6"]
    expected_kit_counts = {"Agilent_SureSelect_Human_All_Exon_V6": 501}
    (
        library_kits,
        _,
        kit_counts,
    ) = targeted_seq_cnv_calling_workflow_large_cohort.pick_kits_and_donors()
    assert library_kits == expected_library_kits
    assert expected_kit_counts == kit_counts


# Global GcnvStepPart Tests ------------------------------------------------------------------------


def test_gcnv_call_assertion(targeted_seq_cnv_calling_workflow):
    """Tests raise UnsupportedActionException"""
    with pytest.raises(UnsupportedActionException):
        targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "_undefined_action_")


def test_gcnv_step_part_get_resource_usage(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_resource()"""
    # Define tested actions
    high_resource_actions = (
        "call_cnvs_cohort_mode",
        "call_cnvs_case_mode",
        "post_germline_calls",
        "post_germline_calls_cohort_mode",
        "post_germline_calls_case_mode",
    )
    all_actions = targeted_seq_cnv_calling_workflow.substep_getattr("gcnv", "actions")
    default_actions = [action for action in all_actions if action not in high_resource_actions]
    # Define expected
    high_res_expected_dict = {
        "threads": 16,
        "time": "2-00:00:00",
        "memory": "46080M",
        "partition": None,
    }
    default_expected_dict = {"threads": 1, "time": "04:00:00", "memory": "7680M", "partition": None}
    # Evaluate - high resource actions
    for action in high_resource_actions:
        for resource, expected in high_res_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = targeted_seq_cnv_calling_workflow.get_resource("gcnv", action, resource)
            assert actual == expected, msg_error

    # Evaluate - all other actions
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = targeted_seq_cnv_calling_workflow.get_resource("gcnv", action, resource)
            assert actual == expected, msg_error


def test_gcnv_get_params(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_params for all actions"""
    actions_w_params = ("model",)
    for action in GCNV_ACTIONS:
        if action in actions_w_params:
            targeted_seq_cnv_calling_workflow.get_params("gcnv", action)
        else:
            with pytest.raises(UnsupportedActionException):
                targeted_seq_cnv_calling_workflow.get_params("gcnv", action)


def test_gcnv_validate_request(
    targeted_seq_cnv_calling_workflow,
    targeted_seq_cnv_calling_workflow_with_gcnv_model,
):
    """Tests GcnvStepPart.validate_request()"""
    # Test small cohort - COHORT MODE
    expected = "cohort_mode"
    actual = targeted_seq_cnv_calling_workflow.substep_getattr("gcnv", "validate_request")()
    assert actual == expected

    # Test large background cohort - CASE MODE, gCNV model provided
    expected = "case_mode"
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.substep_getattr(
        "gcnv", "validate_request"
    )()
    assert actual == expected


def test_gcnv_get_analysis_type(
    targeted_seq_cnv_calling_workflow,
    targeted_seq_cnv_calling_workflow_with_gcnv_model,
):
    """Tests GcnvStepPart.get_analysis_type()"""
    # Test small cohort - COHORT MODE
    expected = "cohort_mode"
    actual = targeted_seq_cnv_calling_workflow.substep_getattr("gcnv", "get_analysis_type")()
    assert actual == expected

    # Test large background cohort - CASE MODE, gCNV model provided
    expected = "case_mode"
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.substep_getattr(
        "gcnv", "get_analysis_type"
    )()
    assert actual == expected


# validate_precomputed_model_paths_config
def test_gcnv_validate_precomputed_model_paths_config(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.validate_model_requirements()"""
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
    targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_precomputed_model_paths_config"
    )(config=[valid_dict])
    # Test key typo
    with pytest.raises(InvalidConfiguration):
        targeted_seq_cnv_calling_workflow.substep_getattr(
            "gcnv", "validate_precomputed_model_paths_config"
        )(config=[valid_dict, typo_dict])
    # Test key missing
    with pytest.raises(InvalidConfiguration):
        targeted_seq_cnv_calling_workflow.substep_getattr(
            "gcnv", "validate_precomputed_model_paths_config"
        )(config=[valid_dict, missing_key_dict])


def test_gcnv_validate_ploidy_model_directory(fake_fs, mocker, targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.validate_ploidy_model_directory()"""
    # Create data directories
    fake_fs.fs.makedirs("/data", exist_ok=True)
    fake_fs.fs.makedirs("/empty", exist_ok=True)
    # Create required files
    tpl = "/data/{file_}"
    for file_ in PLOIDY_MODEL_FILES:
        fake_fs.fs.create_file(tpl.format(file_=file_))
    # Patch out file-system
    patch_module_fs("snappy_pipeline.workflows.targeted_seq_cnv_calling", fake_fs, mocker)

    # Should return True as it is a directory and it contains the expected files
    assert targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_ploidy_model_directory"
    )("/data")
    # Should return False as empty directory
    assert not targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_ploidy_model_directory"
    )("/empty")
    # Should return False not a directory
    assert not targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_ploidy_model_directory"
    )("__not_a_directory__")


def test_gcnv_validate_call_model_directory(fake_fs, mocker, targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.validate_call_model_directory()"""
    # Create data directories
    fake_fs.fs.makedirs("/call_data", exist_ok=True)
    fake_fs.fs.makedirs("/empty", exist_ok=True)
    # Create required files
    tpl = "/call_data/{file_}"
    for file_ in CALL_MODEL_FILES:
        fake_fs.fs.create_file(tpl.format(file_=file_))
    # Patch out file-system
    patch_module_fs("snappy_pipeline.workflows.targeted_seq_cnv_calling", fake_fs, mocker)

    # Should return True as it is a directory and it contains the expected files
    assert targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_call_model_directory"
    )("/call_data")
    # Should return False as empty directory
    assert not targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_call_model_directory"
    )("/empty")
    # Should return False not a directory
    assert not targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "validate_call_model_directory"
    )("__not_a_directory__")


def test_gcnv_get_cnv_model_result_files(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_large_cohort
):
    """Tests GcnvStepPart.get_cnv_model_result_files()"""

    # Test small cohort - 6 individuals, not enough to build a model (<10)
    expected = []
    actual = targeted_seq_cnv_calling_workflow.substep_getattr(
        "gcnv", "get_cnv_model_result_files"
    )(None)
    actual = sorted(actual)
    assert actual == expected

    # Test large trio cohort - 501 individuals, all Agilent v6, enough for a model (>10)
    interval_file = (
        "work/bwa.gcnv_filter_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_filter_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    ploidy_file = (
        "work/bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    expected = sorted([interval_file, ploidy_file])
    actual = targeted_seq_cnv_calling_workflow_large_cohort.substep_getattr(
        "gcnv", "get_cnv_model_result_files"
    )(None)
    actual = sorted(actual)
    assert actual == expected


# Tests for GcnvStepPart (preprocess_intervals) ----------------------------------------------------


def test_gcnv_preprocess_intervals_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_preprocess_intervals()"""
    # Define expected - empty dictionary for all
    expected = {}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "preprocess_intervals")(None)
    assert actual == expected


def test_gcnv_preprocess_intervals_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_preprocess_intervals()"""
    # Define expected
    output_path = (
        "work/gcnv_preprocess_intervals.{library_kit}/out/"
        "gcnv_preprocess_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "preprocess_intervals")
    assert actual == expected


def test_gcnv_target_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'preprocess_intervals' step"""
    # Define expected
    expected = (
        "work/gcnv_preprocess_intervals.{library_kit}/log/"
        "gcnv_preprocess_intervals.{library_kit}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "preprocess_intervals")
    assert actual == expected


# Tests for GcnvStepPart (annotate_gc) -------------------------------------------------------------


def test_gcnv_annotate_gc_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_annotate_gc()"""
    # Define expected
    output_path = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual - Note: library kit defined in conftest: germline_sheet_tsv
    wildcards = Wildcards(fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6"})
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "annotate_gc")(wildcards)
    assert actual == expected


def test_gcnv_annotate_gc_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_annotate_gc()"""
    # Define expected
    expected = {"tsv": "work/gcnv_annotate_gc.{library_kit}/out/gcnv_annotate_gc.{library_kit}.tsv"}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "annotate_gc")
    assert actual == expected


def test_gcnv_annotate_gc_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'annotate_gc' step"""
    # Define expected
    expected = "work/gcnv_annotate_gc.{library_kit}/log/gcnv_annotate_gc.{library_kit}.log"
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "annotate_gc")
    assert actual == expected


# Tests for GcnvStepPart (filter_intervals) --------------------------------------------------------


def test_gcnv_filter_intervals_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_filter_intervals()"""
    # Define expected
    interval_list_out = (
        "work/gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_preprocess_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    tsv_out = (
        "work/gcnv_annotate_gc.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "gcnv_annotate_gc.Agilent_SureSelect_Human_All_Exon_V6.tsv"
    )
    csv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    csv_list_out = [csv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    expected = {"interval_list": interval_list_out, "tsv": tsv_out, "covs": csv_list_out}
    # Get actual. Notes:
    # - library kit defined in conftest: `germline_sheet_tsv`
    # - mapper defined in `minimal_config`
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "filter_intervals")(
        wildcards
    )
    assert actual == expected


def test_gcnv_filter_intervals_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_filter_intervals()"""
    # Define expected
    interval_list_out = (
        "work/{mapper}.gcnv_filter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_filter_intervals.{library_kit}.interval_list"
    )
    expected = {"interval_list": interval_list_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "filter_intervals")
    assert actual == expected


def test_gcnv_filter_intervals_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'filter_intervals' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="filter_intervals")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "filter_intervals")
    assert actual == expected


# Tests for GcnvStepPart (scatter_intervals) -------------------------------------------------------


def test_gcnv_scatter_intervals_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_scatter_intervals()"""
    # Define expected
    output_path = (
        "work/bwa.gcnv_filter_intervals.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_filter_intervals.Agilent_SureSelect_Human_All_Exon_V6.interval_list"
    )
    expected = {"interval_list": output_path}
    # Get actual - Note: library kit defined in conftest: germline_sheet_tsv
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "scatter_intervals")(
        wildcards
    )
    assert actual == expected


def test_gcnv_scatter_intervals_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_scatter_intervals()"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_scatter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_scatter_intervals.{library_kit}"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "scatter_intervals")
    assert actual == expected


def test_gcnv_scatter_intervals_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'scatter_intervals' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="scatter_intervals")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "scatter_intervals")
    assert actual == expected


# Tests for GcnvStepPart (coverage) ----------------------------------------------------------------


def test_gcnv_coverage_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_coverage()"""
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "coverage")(wildcards)
    assert actual == expected


def test_gcnv_coverage_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_coverage()"""
    # Define expected
    tsv_out = (
        "work/{mapper}.gcnv_coverage.{library_name}/out/{mapper}.gcnv_coverage.{library_name}.tsv"
    )
    expected = {"tsv": tsv_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "coverage")
    assert actual == expected


def test_gcnv_coverage_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'coverage' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_coverage.{library_name}/log/{mapper}.gcnv_coverage.{library_name}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "coverage")
    assert actual == expected


# Tests for GcnvStepPart (contig_ploidy) -----------------------------------------------------------


def test_gcnv_contig_ploidy_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_contig_ploidy()"""
    # Define expected
    interval_list_out = (
        "work/{mapper}.gcnv_filter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_filter_intervals.{library_kit}.interval_list"
    )
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    expected = {"interval_list": interval_list_out, "tsv": tsv_list_out}
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "contig_ploidy")(wildcards)
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_contig_ploidy()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_contig_ploidy.{library_kit}/out/"
        "{mapper}.gcnv_contig_ploidy.{library_kit}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "contig_ploidy")
    assert actual == expected


def test_gcnv_contig_ploidy_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'contig_ploidy' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="contig_ploidy")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "contig_ploidy")
    assert actual == expected


# Tests for GcnvStepPart (contig_ploidy_case_mode) -------------------------------------------------


def test_gcnv_contig_ploidy_case_mode_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_contig_ploidy_case_mode()"""
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "contig_ploidy_case_mode")(
        wildcards
    )
    assert actual == expected


def test_gcnv_contig_ploidy_case_mode_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_contig_ploidy_case_mode()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_contig_ploidy_case_mode.{library_kit}/out/"
        "{mapper}.gcnv_contig_ploidy_case_mode.{library_kit}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "contig_ploidy_case_mode")
    assert actual == expected


def test_gcnv_get_params_ploidy_model(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_with_gcnv_model
):
    """Tests GcnvStepPart._get_params_ploidy_model()"""
    # Initialise wildcard
    wildcards = Wildcards(fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6"})
    wildcards_fake = Wildcards(fromdict={"library_kit": "__not_a_library_kit__"})
    # Test small cohort - undefined model and wrong mode
    expected = {}
    actual = targeted_seq_cnv_calling_workflow.get_params("gcnv", "ploidy_model")(wildcards)
    assert actual == expected
    # Test large cohort - model defined in config
    expected = {"model": "/path/to/ploidy-model"}
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_params("gcnv", "ploidy_model")(
        wildcards
    )
    assert actual == expected
    # Test large cohort - model not defined in config
    expected = {"model": "__no_ploidy_model_for_library_in_config__"}
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_params("gcnv", "ploidy_model")(
        wildcards_fake
    )
    assert actual == expected


def test_gcnv_contig_ploidy_case_mode_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'contig_ploidy' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="contig_ploidy_case_mode")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "contig_ploidy_case_mode")
    assert actual == expected


# Tests for GcnvStepPart (call_cnvs_cohort_mode) ---------------------------------------------------


def test_gcnv_call_cnvs_cohort_mode_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_call_cnvs_cohort_mode()"""
    # Define expected
    interval_list_shard_out = (
        "work/{mapper}.gcnv_scatter_intervals.{library_kit}/out/"
        "{mapper}.gcnv_scatter_intervals.{library_kit}/temp_{shard}/scattered.interval_list"
    )
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    ploidy_out = (
        "work/bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    intervals_out = "work/gcnv_annotate_gc.{library_kit}/out/gcnv_annotate_gc.{library_kit}.tsv"
    expected = {
        "interval_list_shard": interval_list_shard_out,
        "tsv": tsv_list_out,
        "ploidy": ploidy_out,
        "intervals": intervals_out,
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "call_cnvs_cohort_mode")(
        wildcards
    )
    assert actual == expected


def test_gcnv_call_cnvs_cohort_mode_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_call_cnvs_cohort_mode()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/out/"
        "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "call_cnvs_cohort_mode")
    assert actual == expected


def test_gcnv_call_cnvs_cohort_mode_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'call_cnvs_cohort_mode' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_call_cnvs_cohort_mode.{library_kit}.{shard}/log/"
        "{mapper}.gcnv_call_cnvs_cohort_mode.{library_kit}.{shard}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "call_cnvs_cohort_mode")
    assert actual == expected


# Tests for GcnvStepPart (call_cnvs_case_mode) -----------------------------------------------------


def test_gcnv_call_cnvs_case_mode_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_call_cnvs_case_mode()"""
    # Define expected
    tsv_pattern = (
        "work/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1/out/bwa.gcnv_coverage.P00{i}-N1-DNA1-WGS1.tsv"
    )
    tsv_list_out = [tsv_pattern.format(i=i) for i in range(1, 7)]  # P001 - P006
    ploidy_out = (
        "work/bwa.gcnv_contig_ploidy_case_mode.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy_case_mode.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    expected = {
        "tsv": tsv_list_out,
        "ploidy": ploidy_out,
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "call_cnvs_case_mode")(
        wildcards
    )
    assert actual == expected


def test_gcnv_call_cnvs_case_mode_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_call_cnvs_case_mode()"""
    # Define expected
    done_out = (
        "work/{mapper}.gcnv_call_cnvs_case_mode.{library_kit}.{shard}/out/"
        "{mapper}.gcnv_call_cnvs_case_mode.{library_kit}.{shard}/.done"
    )
    expected = {"done": done_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "call_cnvs_case_mode")
    assert actual == expected


def test_gcnv_get_params_model(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_with_gcnv_model
):
    """Tests GcnvStepPart._get_params_model()"""
    # Initialise wildcard
    wildcards_01 = Wildcards(
        fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6", "shard": "01"}
    )
    wildcards_02 = Wildcards(
        fromdict={"library_kit": "Agilent_SureSelect_Human_All_Exon_V6", "shard": "02"}
    )
    wildcards_fake = Wildcards(fromdict={"library_kit": "__not_a_library_kit__"})
    # Test small cohort - undefined model
    expected = {}
    actual = targeted_seq_cnv_calling_workflow.get_params("gcnv", "model")(wildcards_01)
    assert actual == expected
    # Test large cohort - model defined in config - shard 01
    expected = {"model": "/data/model_01"}
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_params("gcnv", "model")(
        wildcards_01
    )
    assert actual == expected
    # Test large cohort - model defined in config - shard 02
    expected = {"model": "/data/model_02"}
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_params("gcnv", "model")(
        wildcards_02
    )
    assert actual == expected
    # Test large cohort - model not defined in config
    expected = {"model": "__no_model_for_library_in_config__"}
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_params("gcnv", "model")(
        wildcards_fake
    )
    assert actual == expected


def test_gcnv_call_cnvs_case_mode_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'call_cnvs_case_mode' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_call_cnvs_case_mode.{library_kit}.{shard}/log/"
        "{mapper}.gcnv_call_cnvs_case_mode.{library_kit}.{shard}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "call_cnvs_case_mode")
    assert actual == expected


# Tests for GcnvStepPart (post_germline_calls_case_mode) -------------------------------------------


def test_gcnv_post_germline_calls_case_mode_step_part_get_input_files(
    targeted_seq_cnv_calling_workflow_with_gcnv_model,
):
    """Tests GcnvStepPart._get_input_files_post_germline_calls_case_mode()"""
    # Define expected
    call_pattern = (
        "work/bwa.gcnv_call_cnvs_case_mode.Agilent_SureSelect_Human_All_Exon_V6.0{i}/out/"
        "bwa.gcnv_call_cnvs_case_mode.Agilent_SureSelect_Human_All_Exon_V6.0{i}/.done"
    )
    call_list_out = [call_pattern.format(i=i) for i in range(1, 4)]  # model 01 - 03
    ploidy_out = (
        "work/bwa.gcnv_contig_ploidy_case_mode.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_contig_ploidy_case_mode.Agilent_SureSelect_Human_All_Exon_V6/.done"
    )
    expected = {
        "calls": call_list_out,
        "ploidy": ploidy_out,
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_input_files(
        "gcnv", "post_germline_calls_case_mode"
    )(wildcards)
    assert actual == expected


def test_gcnv_get_params_postgermline_models(
    targeted_seq_cnv_calling_workflow, targeted_seq_cnv_calling_workflow_with_gcnv_model
):
    """Tests GcnvStepPart._get_params_postgermline_models()"""
    # Initialise wildcard
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    # Test small cohort - undefined model and wrong mode
    expected = {}
    actual = targeted_seq_cnv_calling_workflow.get_params("gcnv", "postgermline_models")(wildcards)
    assert actual == expected
    # Test large cohort - model defined in config
    expected = {"model": ["/data/model_01", "/data/model_02", "/data/model_03"]}
    actual = targeted_seq_cnv_calling_workflow_with_gcnv_model.get_params(
        "gcnv", "postgermline_models"
    )(wildcards)
    assert actual == expected


# Tests for GcnvStepPart (merge_cohort_vcfs) -------------------------------------------------------


def test_gcnv_merge_cohort_vcfs_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_merge_cohort_vcfs()"""
    # Define expected
    pattern_out = (
        "work/bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gcnv_post_germline_calls.P00{i}-N1-DNA1-WGS1.vcf.gz"
    )
    expected = [pattern_out.format(i=i) for i in range(1, 7)]  # P001 - P006
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "library_kit": "Agilent_SureSelect_Human_All_Exon_V6"}
    )
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "merge_cohort_vcfs")(
        wildcards
    )
    assert actual == expected


def test_gcnv_merge_cohort_vcfs_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_merge_cohort_vcfs()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.gcnv_merge_cohort_vcfs.{library_kit}/out/"
        "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "merge_cohort_vcfs")
    assert actual == expected


def test_gcnv_merge_cohort_vcfs_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'merge_cohort_vcfs' step"""
    # Define expected
    expected = get_expected_gcnv_log_file(step_name="merge_cohort_vcfs")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "merge_cohort_vcfs")
    assert actual == expected


# Tests for GcnvStepPart (extract_ped) -------------------------------------------------------------


def test_gcnv_extract_ped_vcfs_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_input_files_extract_ped()"""
    # Define expected
    pattern_out = (
        "work/bwa.gcnv_merge_cohort_vcfs.Agilent_SureSelect_Human_All_Exon_V6/out/"
        "bwa.gcnv_merge_cohort_vcfs.Agilent_SureSelect_Human_All_Exon_V6"
    )
    expected = {"vcf": pattern_out + ".vcf.gz", "tbi": pattern_out + ".vcf.gz.tbi"}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = targeted_seq_cnv_calling_workflow.get_input_files("gcnv", "extract_ped")(wildcards)
    assert actual == expected


def test_gcnv_extract_ped_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart._get_output_files_extract_ped()"""
    # Define expected
    pattern_out = "work/{mapper}.gcnv.{library_name}/out/{mapper}.gcnv.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=pattern_out)
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("gcnv", "extract_ped")
    assert actual == expected


def test_gcnv_extract_ped_step_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests GcnvStepPart.get_log_file for 'extract_ped' step"""
    # Define expected
    expected = (
        "work/{mapper}.gcnv_extract_ped.{library_name}/log/"
        "{mapper}.gcnv_extract_ped.{library_name}.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("gcnv", "extract_ped")
    assert actual == expected


# Global XhmmStepPart Tests ------------------------------------------------------------------------


def test_xhmm_call_assertion(targeted_seq_cnv_calling_workflow):
    """Tests raise UnsupportedActionException"""
    with pytest.raises(UnsupportedActionException):
        targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "_undefined_action_")


def test_xhmm_get_params(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_params for all actions"""
    for action in XHMM_ACTIONS:
        if action == "coverage":
            targeted_seq_cnv_calling_workflow.get_params("xhmm", action)
        else:
            with pytest.raises(UnsupportedActionException):
                targeted_seq_cnv_calling_workflow.get_params("xhmm", action)


# Tests for XhmmStepPart (coverage) ----------------------------------------------------------------


def test_xhmm_coverage_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_input_files_coverage()"""
    # Define expected
    bam_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1"
    expected = {"bam": bam_out + ".bam", "bai": bam_out + ".bam.bai"}
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "coverage")(wildcards)
    assert actual == expected


def test_xhmm_coverage_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "coverage")
    assert actual == expected


def test_xhmm_coverage_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'coverage' step"""
    # Define expected
    expected = (
        "work/{mapper}.xhmm_coverage.{library_name}/log/snakemake.targeted_seq_cnv_calling.log"
    )
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "coverage")
    assert actual == expected


# Tests for XhmmStepPart (merge_cov) ---------------------------------------------------------------


def test_xhmm_merge_cov_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "merge_cov")(wildcards)
    assert actual == expected


def test_xhmm_merge_cov_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_merge_cov()"""
    # Define expected
    name_out = (
        "work/{mapper}.xhmm_merge_cov.{library_kit}/out/"
        "{mapper}.xhmm_merge_cov.{library_kit}.RD.txt"
    )
    expected = [name_out]
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "merge_cov")
    assert actual == expected


def test_xhmm_merge_cov_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'merge_cov' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="merge_cov")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "merge_cov")
    assert actual == expected


# Tests for XhmmStepPart (ref_stats) ---------------------------------------------------------------


def test_xhmm_ref_stats_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_ref_stats()"""
    # Define expected
    name_out = (
        "work/{mapper}.xhmm_ref_stats.{library_kit}/out/"
        "{mapper}.xhmm_ref_stats.{library_kit}.extreme_gc_targets.txt"
    )
    expected = {"extreme_gc_targets": name_out}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "ref_stats")
    assert actual == expected


def test_xhmm_ref_stats_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'ref_stats' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="ref_stats")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "ref_stats")
    assert actual == expected


# Tests for XhmmStepPart (filter_center) -----------------------------------------------------------


def test_xhmm_filter_center_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "filter_center")(wildcards)
    assert actual == expected


def test_xhmm_filter_center_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "filter_center")
    assert actual == expected


def test_xhmm_filter_center_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'filter_center' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="filter_center")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "filter_center")
    assert actual == expected


# Tests for XhmmStepPart (pca) ---------------------------------------------------------------------


def test_xhmm_pca_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "pca")(wildcards)
    assert actual == expected


def test_xhmm_pca_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_pca()"""
    # Define expected
    pattern_out = "work/{mapper}.xhmm_pca.{library_kit}/out/{mapper}.xhmm_pca.{library_kit}"
    expected = {
        "pc_loading": pattern_out + ".PC_LOADINGS.txt",
        "pc_sd": pattern_out + ".PC_SD.txt",
        "pc": pattern_out + ".PC.txt",
    }
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "pca")
    assert actual == expected


def test_xhmm_pca_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'pca' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="pca")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "pca")
    assert actual == expected


# Tests for XhmmStepPart (normalize) ---------------------------------------------------------------


def test_xhmm_normalize_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "normalize")(wildcards)
    assert actual == expected


def test_xhmm_normalize_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_normalize()"""
    # Define expected
    pattern_out = (
        "work/{mapper}.xhmm_normalize.{library_kit}/out/{mapper}.xhmm_normalize.{library_kit}"
    )
    expected = {"normalized": pattern_out, "num_removed": pattern_out + ".num_removed_PC.txt"}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "normalize")
    assert actual == expected


def test_xhmm_normalize_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'normalize' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="normalize")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "normalize")
    assert actual == expected


# Tests for XhmmStepPart (zscore_center) -----------------------------------------------------------


def test_xhmm_zscore_center_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "zscore_center")(wildcards)
    assert actual == expected


def test_xhmm_zscore_center_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "zscore_center")
    assert actual == expected


def test_xhmm_zscore_center_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'zscore_center' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="zscore_center")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "zscore_center")
    assert actual == expected


# Tests for XhmmStepPart (zscore_center) -----------------------------------------------------------


def test_xhmm_refilter_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "refilter")(wildcards)
    assert actual == expected


def test_xhmm_refilter_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_refilter()"""
    # Define expected
    base_out = (
        "work/{mapper}.xhmm_refilter.{library_kit}/out/{mapper}.xhmm_refilter.{library_kit}.RD.txt"
    )
    expected = [base_out]
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "refilter")
    assert actual == expected


def test_xhmm_refilter_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'refilter' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="refilter")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "refilter")
    assert actual == expected


# Tests for XhmmStepPart (discover) ----------------------------------------------------------------


def test_xhmm_discover_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "discover")(wildcards)
    assert actual == expected


def test_xhmm_discover_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_discover()"""
    # Define expected
    base_out = "work/{mapper}.xhmm_discover.{library_kit}/out/{mapper}.xhmm_discover.{library_kit}"
    expected = {"xcnv": base_out + ".xcnv", "aux_xcnv": base_out + ".aux_xcnv"}
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "discover")
    assert actual == expected


def test_xhmm_discover_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'discover' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="discover")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "discover")
    assert actual == expected


# Tests for XhmmStepPart (genotype) ----------------------------------------------------------------


def test_xhmm_genotype_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "genotype")(wildcards)
    assert actual == expected


def test_xhmm_genotype_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_genotype()"""
    # Define expected
    base_out = "work/{mapper}.xhmm_genotype.{library_kit}/out/{mapper}.xhmm_genotype.{library_kit}"
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "genotype")
    assert actual == expected


def test_xhmm_genotype_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'genotype' step"""
    # Define expected
    expected = get_expected_xhmm_log_file(step_name="genotype")
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "genotype")
    assert actual == expected


# Tests for XhmmStepPart (extract_ped) -------------------------------------------------------------


def test_xhmm_extract_ped_step_part_get_input_files(targeted_seq_cnv_calling_workflow):
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
        "tbi": vcf_pattern_out + ".tbi",
    }
    # Get actual
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = targeted_seq_cnv_calling_workflow.get_input_files("xhmm", "extract_ped")(wildcards)
    assert actual == expected


def test_xhmm_extract_ped_step_part_get_output_files(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart._get_output_files_genotype()"""
    # Define expected
    base_out = "work/{mapper}.xhmm.{library_name}/out/{mapper}.xhmm.{library_name}"
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_output_files("xhmm", "extract_ped")
    assert actual == expected


def test_xhmm_extract_ped_part_get_log_file(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_log_file for 'genotype' step"""
    # Define expected
    expected = "work/{mapper}.xhmm.{library_name}/log/snakemake.targeted_seq_cnv_calling.log"
    # Get actual
    actual = targeted_seq_cnv_calling_workflow.get_log_file("xhmm", "extract_ped")
    assert actual == expected


def test_xhmm_step_part_get_resource_usage(targeted_seq_cnv_calling_workflow):
    """Tests XhmmStepPart.get_resource()"""
    # Define expected
    merge_cov_expected_dict = {
        "threads": 1,
        "time": "1-00:00:00",
        "memory": "12G",
        "partition": None,
    }
    default_expected_dict = {"threads": 1, "time": "08:00:00", "memory": "12G", "partition": None}

    # Evaluate - merge_cov
    for resource, expected in merge_cov_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = targeted_seq_cnv_calling_workflow.get_resource("xhmm", "merge_cov", resource)
        assert actual == expected, msg_error

    # Evaluate - all other actions
    all_actions = targeted_seq_cnv_calling_workflow.substep_getattr("xhmm", "actions")
    default_actions = [action for action in all_actions if action != "merge_cov"]
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = targeted_seq_cnv_calling_workflow.get_resource("xhmm", action, resource)
            assert actual == expected, msg_error
