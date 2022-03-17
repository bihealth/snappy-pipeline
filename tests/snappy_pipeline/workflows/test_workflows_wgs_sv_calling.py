# -*- coding: utf-8 -*-
"""Tests for the wgs_sv_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_sv_calling import WgsSvCallingWorkflow

from .common import get_expected_output_bcf_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


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

          wgs_sv_calling:
            tools:
              dna:
              - manta
              - delly2
              long_dna:
              - pb_honey_spots

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
def wgs_sv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
    fai_file_content,
):
    """Return WgsSvCallingWorkflow object pre-configured with germline sheet"""
    # Add Fasta file
    # Create FASTA files
    germline_sheet_fake_fs.fs.create_file(
        "/path/to/ref.fa.fai", contents=fai_file_content, create_missing_dirs=True
    )
    germline_sheet_fake_fs.fs.create_file("/path/to/ref.fa", create_missing_dirs=True)
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.wgs_sv_calling", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really a NGSMappingPipelineStep here
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return WgsSvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for Delly2StepPart (all) -------------------------------------------------------------------


def test_delly2_step_part_get_ped_members(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_ped_members()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected = "P001-N1-DNA1-WGS1 P002-N1-DNA1-WGS1 P003-N1-DNA1-WGS1"
    # Get actual
    actual = wgs_sv_calling_workflow.substep_getattr("delly2", "get_ped_members")(wildcards)
    assert actual == expected


def test_delly2_step_part_get_resource_usage(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_resource_usage()"""
    # Set actions
    all_actions = wgs_sv_calling_workflow.substep_getattr("delly2", "actions")
    cheap_actions = ("merge_genotypes", "merge_calls", "reorder_vcf")
    default_actions = [action for action in all_actions if action not in cheap_actions]

    # Define expected
    cheap_expected_dict = {
        "threads": 2,
        "time": "4-00:00:00",
        "memory": "14336M",
        "partition": None,
    }
    default_expected_dict = {
        "threads": 2,
        "time": "7-00:00:00",
        "memory": "40960M",
        "partition": None,
    }

    # Evaluate - cheap actions
    for action in cheap_actions:
        for resource, expected in cheap_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' for action '{action}'."
            actual = wgs_sv_calling_workflow.get_resource("delly2", action, resource)
            assert actual == expected, msg_error

    # Evaluate - default actions
    for action in default_actions:
        for resource, expected in default_expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' for action '{action}'."
            actual = wgs_sv_calling_workflow.get_resource("delly2", action, resource)
            assert actual == expected, msg_error


# Tests for Delly2StepPart (call) ------------------------------------------------------------------


def test_delly2_step_part_call_get_input_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_call()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "call")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    assert actual == expected


def test_delly2_step_part_call_get_output_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - call"""
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.call.{library_name,[^\.]+}/out/"
        r"{mapper}.delly2.call.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "call")
    assert actual == expected


def test_delly_step_part_call_get_log_file(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - call"""
    expected = "work/{mapper}.delly2.call.{library_name}/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("delly2", "call")
    assert actual == expected


# Tests for Delly2StepPart (merge_calls) ----------------------------------------------------------


def test_delly2_step_part_merge_calls_get_input_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_merge_calls()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "merge_calls")(wildcards)
    expected = [
        "work/bwa.delly2.call.P001-N1-DNA1-WGS1/out/bwa.delly2.call.P001-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P002-N1-DNA1-WGS1/out/bwa.delly2.call.P002-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P003-N1-DNA1-WGS1/out/bwa.delly2.call.P003-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P004-N1-DNA1-WGS1/out/bwa.delly2.call.P004-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P005-N1-DNA1-WGS1/out/bwa.delly2.call.P005-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.call.P006-N1-DNA1-WGS1/out/bwa.delly2.call.P006-N1-DNA1-WGS1.bcf",
    ]
    assert actual == expected


def test_delly2_step_part_merge_calls_get_output_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - merge_calls"""
    # Define expected
    base_name_out = r"work/{mapper,[^\.]+}.delly2.merge_calls/out/{mapper}.delly2.merge_calls"
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "merge_calls")
    assert actual == expected


def test_delly_step_part_merge_calls_get_log_file(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - merge_calls"""
    expected = "work/{mapper}.delly2.merge_calls/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("delly2", "merge_calls")
    assert actual == expected


# Tests for Delly2StepPart (genotype) -------------------------------------------------------------


def test_delly2_step_part_genotype_get_input_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_genotype()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "genotype")(wildcards)
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bcf": "work/bwa.delly2.merge_calls/out/bwa.delly2.merge_calls.bcf",
    }
    assert actual == expected


def test_delly2_step_part_genotype_get_output_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - genotype"""
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.genotype.{library_name,[^\.]+}/out/"
        r"{mapper}.delly2.genotype.{library_name}"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "genotype")
    assert actual == expected


def test_delly_step_part_genotype_get_log_file(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - genotype"""
    expected = "work/{mapper}.delly2.genotype.{library_name}/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("delly2", "genotype")
    assert actual == expected


# Tests for Delly2StepPart (merge_genotypes) ------------------------------------------------------


def test_delly2_step_part_merge_genotypes_get_input_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_merge_genotypes()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "merge_genotypes")(wildcards)
    expected = [
        "work/bwa.delly2.genotype.P001-N1-DNA1-WGS1/out/bwa.delly2.genotype.P001-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P002-N1-DNA1-WGS1/out/bwa.delly2.genotype.P002-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P003-N1-DNA1-WGS1/out/bwa.delly2.genotype.P003-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P004-N1-DNA1-WGS1/out/bwa.delly2.genotype.P004-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P005-N1-DNA1-WGS1/out/bwa.delly2.genotype.P005-N1-DNA1-WGS1.bcf",
        "work/bwa.delly2.genotype.P006-N1-DNA1-WGS1/out/bwa.delly2.genotype.P006-N1-DNA1-WGS1.bcf",
    ]
    assert actual == expected


def test_delly2_step_part_merge_genotypes_get_output_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - merge_genotypes"""
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.merge_genotypes/out/{mapper}.delly2.merge_genotypes"
    )
    expected = get_expected_output_bcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "merge_genotypes")
    assert actual == expected


def test_delly_step_part_merge_genotypes_get_log_file(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - merge_genotypes"""
    expected = "work/{mapper}.delly2.merge_genotypes/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("delly2", "merge_genotypes")
    assert actual == expected


# Tests for Delly2StepPart (reorder_vcf) ----------------------------------------------------------


def test_delly2_step_part_reorder_vcf_get_input_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart._get_input_files_reorder_vcf()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("delly2", "reorder_vcf")(wildcards)
    expected = {"bcf": "work/bwa.delly2.merge_genotypes/out/bwa.delly2.merge_genotypes.bcf"}
    assert actual == expected


def test_delly2_step_part_reorder_vcf_get_output_files(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_output_files() - reorder_vcf"""
    # Define expected
    base_name_out = (
        r"work/{mapper,[^\.]+}.delly2.{index_ngs_library,[^\.]+}/out/"
        r"{mapper}.delly2.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("delly2", "reorder_vcf")
    assert actual == expected


def test_delly2_step_part_reorder_vcf_get_log_file(wgs_sv_calling_workflow):
    """Tests Delly2StepPart.get_log_file() - reorder_vcf"""
    expected = "work/{mapper}.delly2.{index_ngs_library}/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("delly2", "reorder_vcf")
    assert actual == expected


# Tests for MantaStepPart -------------------------------------------------------------------------


def test_manta_step_part_get_input_files(wgs_sv_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("manta", "run")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_manta_step_part_get_output_files(wgs_sv_calling_workflow):
    # Define expected
    base_name_out = "work/{mapper}.manta.{index_ngs_library}/out/{mapper}.manta.{index_ngs_library}"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("manta", "run")
    assert actual == expected


def test_manta_step_part_get_log_file(wgs_sv_calling_workflow):
    expected = "work/{mapper}.manta.{index_ngs_library}/log/snakemake.log"
    assert wgs_sv_calling_workflow.get_log_file("manta", "run") == expected


def test_manta_step_part_get_resource_usage(wgs_sv_calling_workflow):
    """Tests MantaStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "1-16:00:00", "memory": "61440M", "partition": None}
    # Evaluate
    # Note: only action available is 'run'
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'run'."
        actual = wgs_sv_calling_workflow.get_resource("manta", "run", resource)
        assert actual == expected, msg_error


# Tests for SvTkStepPart ---------------------------------------------------------------------------


def test_svtk_step_part_get_input_files(wgs_sv_calling_workflow):
    """Tests SvTkStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "delly2", "library_name": "P001-N1-DNA1-WGS1"}
    )
    expected = {
        "calls": "output/bwa.delly2.call.P001-N1-DNA1-WGS1/out/bwa.delly2.call.P001-N1-DNA1-WGS1.bcf"
    }
    actual = wgs_sv_calling_workflow.get_input_files("svtk", "standardize")(wildcards)
    assert actual == expected


def test_svtk_step_part_get_output_files(wgs_sv_calling_workflow):
    """Tests SvTkStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{caller}.svtk_standardize.{library_name}/out/"
        "{mapper}.{caller}.svtk_standardize.{library_name}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_sv_calling_workflow.get_output_files("svtk", "standardize")
    assert actual == expected


def test_svtk_step_part_get_log_file(wgs_sv_calling_workflow):
    """Tests SvTkStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{caller}.svtk_standardize.{library_name}/log/"
        "{mapper}.{caller}.svtk_standardize.{library_name}.log"
    )
    actual = wgs_sv_calling_workflow.get_log_file("svtk", "standardize")
    assert actual == expected


def test_svtk_step_part_get_resource_usage(wgs_sv_calling_workflow):
    """Tests SvTkStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "04:00:00", "memory": "4096M", "partition": None}
    # Evaluate
    # Note: only action available is 'standardize'
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'standardize'."
        actual = wgs_sv_calling_workflow.get_resource("svtk", "standardize", resource)
        assert actual == expected, msg_error


def test_svtk_step_part_get_ped_members(wgs_sv_calling_workflow):
    """Tests SvTkStepPart.get_ped_members()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected = "P001-N1-DNA1-WGS1 P002-N1-DNA1-WGS1 P003-N1-DNA1-WGS1"
    # Get actual
    actual = wgs_sv_calling_workflow.substep_getattr("svtk", "get_ped_members")(wildcards)
    assert actual == expected


# Tests for PopDelStepPart ---------------------------------------------------------------------------


def test_popdel_step_part_get_input_files_profile(wgs_sv_calling_workflow):
    """Tests PopDelStepPart._get_input_files_profile()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = wgs_sv_calling_workflow.get_input_files("popdel", "profile")(wildcards)
    assert actual == expected


def test_popdel_step_part_get_input_files_call(wgs_sv_calling_workflow):
    """Tests PopDelStepPart._get_input_files_call()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    patient_id = (1, 2, 3, 4, 5, 6)
    base_out = (
        "work/bwa.popdel.internal.profile.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.popdel.internal.profile.P00{i}-N1-DNA1-WGS1.profile"
    )
    expected = {
        "profile": [base_out.format(i=i) for i in patient_id],
    }
    actual = wgs_sv_calling_workflow.get_input_files("popdel", "call")(wildcards)
    assert actual == expected


def test_popdel_step_part_get_input_files_concat_calls(wgs_sv_calling_workflow):
    """Tests PopDelStepPart._get_input_files_concat_calls()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    # Based on chr1 from GRCh37
    regions = (
        "1-10020000",
        "9980000-20020000",
        "19980000-30020000",
        "29980000-40020000",
        "39980000-50020000",
        "49980000-60020000",
        "59980000-70020000",
        "69980000-80020000",
        "79980000-90020000",
        "89980000-100020000",
        "99980000-110020000",
        "109980000-120020000",
        "119980000-130020000",
        "129980000-140020000",
        "139980000-150020000",
        "149980000-160020000",
        "159980000-170020000",
        "169980000-180020000",
        "179980000-190020000",
        "189980000-200020000",
        "199980000-210020000",
        "209980000-220020000",
        "219980000-230020000",
        "229980000-240020000",
        "239980000-249250621",
    )
    base_out = (
        "work/bwa.popdel.internal.call.1-{region}/out/" "bwa.popdel.internal.call.1-{region}.vcf.gz"
    )
    expected = {"vcf": [base_out.format(region=region) for region in regions]}
    actual = wgs_sv_calling_workflow.get_input_files("popdel", "concat_calls")(wildcards)
    assert actual == expected


def test_popdel_step_part_get_input_files_reorder_vcf(wgs_sv_calling_workflow):
    """Tests PopDelStepPart._get_input_files_reorder_vcf()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "vcf": "work/bwa.popdel.internal.concat_calls/out/bwa.popdel.internal.concat_calls.vcf.gz"
    }

    actual = wgs_sv_calling_workflow.get_input_files("popdel", "reorder_vcf")(wildcards)
    assert actual == expected


def test_popdel_step_part_get_output_files_profile(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_output_files() - profile"""
    base_name_out = (
        "work/{mapper}.popdel.internal.profile.{index_ngs_library}/out/"
        "{mapper}.popdel.internal.profile.{index_ngs_library}"
    )
    expected = {
        "profile": base_name_out + ".profile",
        "profile_md5": base_name_out + ".profile.md5",
    }
    actual = wgs_sv_calling_workflow.get_output_files("popdel", "profile")
    assert actual == expected


def test_popdel_step_part_get_output_files_call(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_output_files() - call"""
    base_name_out = (
        "work/{mapper}.popdel.internal.call.{chrom}-{begin}-{end}/out/"
        "{mapper}.popdel.internal.call.{chrom}-{begin}-{end}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_sv_calling_workflow.get_output_files("popdel", "call")
    assert actual == expected


def test_popdel_step_part_get_output_files_concat_calls(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_output_files() - concat_calls"""
    base_name_out = (
        "work/{mapper}.popdel.internal.concat_calls/out/{mapper}.popdel.internal.concat_calls"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_sv_calling_workflow.get_output_files("popdel", "concat_calls")
    assert actual == expected


def test_popdel_step_part_get_output_files_reorder_vcf(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_output_files() - reorder_vcf"""
    base_name_out = (
        "work/{mapper}.popdel.{index_ngs_library}/out/{mapper}.popdel.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = wgs_sv_calling_workflow.get_output_files("popdel", "reorder_vcf")
    assert actual == expected


def test_popdel_step_part_get_log_file_profile(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_log_file() - profile"""
    expected = (
        "work/{mapper}.popdel.internal.profile.{index_ngs_library}/log/"
        "{mapper}.popdel.internal.profile.{index_ngs_library}.log"
    )
    actual = wgs_sv_calling_workflow.get_log_file("popdel", "profile")
    assert actual == expected


def test_popdel_step_part_get_log_file_call(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_log_file() - call"""
    expected = (
        "work/{mapper}.popdel.internal.call.{chrom}-{begin}-{end}/log/"
        "{mapper}.popdel.internal.call.{chrom}-{begin}-{end}.log"
    )
    actual = wgs_sv_calling_workflow.get_log_file("popdel", "call")
    assert actual == expected


def test_popdel_step_part_get_log_file_concat_calls(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_log_file() - concat_calls"""
    expected = (
        "work/{mapper}.popdel.internal.concat_calls/log/{mapper}.popdel.internal.concat_calls.log"
    )
    actual = wgs_sv_calling_workflow.get_log_file("popdel", "concat_calls")
    assert actual == expected


def test_popdel_step_part_get_log_file_reorder_vcf(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_log_file() - reorder_vcf"""
    expected = (
        "work/{mapper}.popdel.{index_ngs_library}/log/{mapper}.popdel.{index_ngs_library}.log"
    )
    actual = wgs_sv_calling_workflow.get_log_file("popdel", "reorder_vcf")
    assert actual == expected


def test_popdel_step_part_get_resource_usage(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_resource_usage()"""
    all_actions = wgs_sv_calling_workflow.substep_getattr("popdel", "actions")
    # Define expected
    expected_dict = {"threads": 2, "time": "4-00:00:00", "memory": "24576M", "partition": None}
    # Evaluate
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' for action '{action}'."
            actual = wgs_sv_calling_workflow.get_resource("popdel", action, resource)
            assert actual == expected, msg_error


def test_popdel_step_part_get_ped_members(wgs_sv_calling_workflow):
    """Tests PopDelStepPart.get_ped_members()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected = "P001-N1-DNA1-WGS1 P002-N1-DNA1-WGS1 P003-N1-DNA1-WGS1"
    # Get actual
    actual = wgs_sv_calling_workflow.substep_getattr("popdel", "get_ped_members")(wildcards)
    assert actual == expected


# Tests for PbHoneySpotsStepPart -------------------------------------------------------------------


def test_pb_honey_spots_step_part_get_input_files(wgs_sv_calling_workflow):
    """Tests PbHoneySpotsStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("pb_honey_spots", "run")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_pb_honey_spots_step_part_get_output_files(wgs_sv_calling_workflow):
    """Tests PbHoneySpotsStepPart.get_output_files()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.pb_honey_spots.{index_ngs_library}/out/"
        "{mapper}.pb_honey_spots.{index_ngs_library}"
    )
    expected = {
        "bed": base_name_out + ".bed",
        "bed_md5": base_name_out + ".bed.md5",
    }
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("pb_honey_spots", "run")
    assert actual == expected


def test_pb_honey_spots_step_part_get_log_file(wgs_sv_calling_workflow):
    """Tests PbHoneySpotsStepPart.get_log_file()"""
    expected = "work/{mapper}.pb_honey_spots.{index_ngs_library}/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("pb_honey_spots", "run")
    assert actual == expected


def test_pb_honey_spots_step_part_get_resource_usage(wgs_sv_calling_workflow):
    """Tests PbHoneySpotsStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "1-16:00:00", "memory": "61440M", "partition": None}
    # Evaluate
    # Note: only action available is 'run'
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'run'."
        actual = wgs_sv_calling_workflow.get_resource("pb_honey_spots", "run", resource)
        assert actual == expected, msg_error


# Tests for SnifflesStepPart -----------------------------------------------------------------------


def test_sniffles_spots_step_part_get_input_files(wgs_sv_calling_workflow):
    """Tests SnifflesStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "index_ngs_library": "P001-N1-DNA1-WGS1"})
    actual = wgs_sv_calling_workflow.get_input_files("sniffles", "run")(wildcards)
    expected = [
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam",
        "NGS_MAPPING/output/bwa.P003-N1-DNA1-WGS1/out/bwa.P003-N1-DNA1-WGS1.bam.bai",
    ]
    assert actual == expected


def test_sniffles_step_part_get_output_files(wgs_sv_calling_workflow):
    """Tests SnifflesStepPart.get_output_files()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.sniffles.{index_ngs_library}/out/{mapper}.sniffles.{index_ngs_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    # Get actual
    actual = wgs_sv_calling_workflow.get_output_files("sniffles", "run")
    assert actual == expected


def test_sniffles_step_part_get_log_file(wgs_sv_calling_workflow):
    """Tests SnifflesStepPart.get_log_file()"""
    expected = "work/{mapper}.sniffles.{index_ngs_library}/log/snakemake.log"
    actual = wgs_sv_calling_workflow.get_log_file("sniffles", "run")
    assert actual == expected


def test_sniffles_step_part_get_resource_usage(wgs_sv_calling_workflow):
    """Tests SnifflesStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 16, "time": "1-16:00:00", "memory": "61440M", "partition": None}
    # Evaluate
    # Note: only action available is 'run'
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'run'."
        actual = wgs_sv_calling_workflow.get_resource("sniffles", "run", resource)
        assert actual == expected, msg_error


# Tests for WgsSvCallingWorkflow -----------------------------------------------------------


def test_wgs_sv_calling_workflow(wgs_sv_calling_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = ["delly2", "link_out", "manta", "pb_honey_spots", "popdel", "sniffles", "svtk"]
    actual = list(sorted(wgs_sv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{sv_caller}.P00{i}-N1-DNA1-WGS1/out/"
        "{mapper}.{sv_caller}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, sv_caller=sv_caller, i=i, ext=ext)
        for mapper in ("bwa",)
        for sv_caller in ("delly2", "manta")
        for i in [1, 4]  # only indices
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_sv_calling_workflow.get_result_files()))
    assert actual == expected
