# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

from snappy_pipeline.workflows.tumor_mutational_burden import (
    TumorMutationalBurdenCalculationWorkflow,
)

from .common import get_expected_log_files_dict, get_expected_output_json_files_dict
from .conftest import patch_module_fs


# Test tumor mutational burden calculation with vcf file from somatic variant calling step
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
          cosmic:
            path: /path/to/cosmic.vcf.gz
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_variant_calling:
            tools:
            - mutect2
            - scalpel
            mutect2: {}
            scalpel:
              path_target_regions: /path/to/target/regions.bed

          somatic_variant_annotation:
            path_somatic_variant_calling: ../somatic_variant_calling
            tools: ["vep", "jannovar"]
            jannovar:
              path_jannovar_ser: /path/to/jannover.ser
              flag_off_target: False
              dbnsfp: {}
            vep:
              cache_dir: /path/to/dir/cache

          somatic_variant_filtration:
            filtration_schema: sets
            filter_sets:
              dkfz_only: ~
              dkfz_and_ebfilter: {}
              dkfz_and_oxog: {}
              dkfz_and_ebfilter_and_oxog: {}
            exon_lists: {}

          tumor_mutational_burden:
            path_somatic_variant: ../somatic_variant_filtration
            tools_ngs_mapping: []
            has_annotation: True # REQUIRED
            is_filtered: True
            tools_somatic_variant_calling: []
            target_regions: /path/to/regions.bed

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
def tumor_mutational_burden_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return TumorMutationalBurdenCalculationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)

    # Construct the workflow object
    return TumorMutationalBurdenCalculationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for TumorMutationalBurdenCalculationStepPart -----------------------------------------------------


def test_tumor_mutational_step_part_get_input_files(tumor_mutational_burden_workflow):
    """Test TumorMutationalBurdenCalculationStepPart.get_input_files()"""
    base_out = (
        "../somatic_variant_filtration/output/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter}.{region}/out/"
        "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter}.{region}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = tumor_mutational_burden_workflow.get_input_files("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_output_files(tumor_mutational_burden_workflow):
    """Tests TumorMutationalBurdenCalculationStepPart.get_output_files()"""
    base_out = (
        "output/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tmb.{tumor_library}.{filter}.{region}/out/"
        "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tmb.{tumor_library}.{filter}.{region}"
    )
    expected = get_expected_output_json_files_dict(base_out=base_out)
    actual = tumor_mutational_burden_workflow.get_output_files("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_log_files(tumor_mutational_burden_workflow):
    """Tests TumorMutationalBurdenCalculationStepPart.get_log_files()"""
    base_out = (
        "output/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tmb.{tumor_library}.{filter}.{region}/log/"
        "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tmb.{tumor_library}.{filter}.{region}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = tumor_mutational_burden_workflow.get_log_file("tmb_gathering", "run")
    assert actual == expected


def test_tumor_mutational_step_part_get_resource_usage(tumor_mutational_burden_workflow):
    """Tests TumorMutationalBurdenCalculationStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "1:00:00", "memory": "4096M"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = tumor_mutational_burden_workflow.get_resource("tmb_gathering", "run", resource)()
        assert actual == expected, msg_error


# Tests for TumorMutationalBurdenCalculationWorkflow -------------------------------------------------------


def test_tumor_mutational_burden_workflow(tumor_mutational_burden_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "tmb_gathering"]
    actual = list(sorted(tumor_mutational_burden_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tmb.P00{i}-T{t}-DNA1-WGS1.{filt}.{region}/{dir_}/"
        "{mapper}.{var_caller}.{anno_caller}.dkfz_bias_filter.eb_filter.tmb.P00{i}-T{t}-DNA1-WGS1.{filt}.{region}.{ext}"
    )
    expected = [
        tpl.format(
            mapper=mapper,
            var_caller=var_caller,
            anno_caller=anno_caller,
            filt=filt,
            region=region,
            i=i,
            t=t,
            ext=ext,
            dir_="out",
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("json", "json.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect2", "scalpel")
        for anno_caller in ("vep", "jannovar")
        for filt in (
            "no_filter",
            "dkfz_only",
            "dkfz_and_ebfilter",
            "dkfz_and_ebfilter_and_oxog",
            "dkfz_and_oxog",
        )
        for region in ("genome_wide",)
    ]
    expected += [
        tpl.format(
            mapper=mapper,
            var_caller=var_caller,
            anno_caller=anno_caller,
            filt=filt,
            region=region,
            i=i,
            t=t,
            ext=ext,
            dir_="log",
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "conda_info.txt",
            "conda_list.txt",
            "log",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
            "log.md5",
        )
        for mapper in ("bwa",)
        for var_caller in ("mutect2", "scalpel")
        for anno_caller in ("vep", "jannovar")
        for filt in (
            "no_filter",
            "dkfz_only",
            "dkfz_and_ebfilter",
            "dkfz_and_ebfilter_and_oxog",
            "dkfz_and_oxog",
        )
        for region in ("genome_wide",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(tumor_mutational_burden_workflow.get_result_files()))
    assert expected == actual
