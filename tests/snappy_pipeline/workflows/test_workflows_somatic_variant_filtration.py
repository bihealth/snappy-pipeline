# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_filtration workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_filtration import SomaticVariantFiltrationWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


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

          somatic_variant_filtration:
            tools_ngs_mapping: ['bwa']
            tools_somatic_variant_calling: ['mutect2']
            tools_somatic_variant_annotation: ['jannovar']

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
def somatic_variant_filtration_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticVariantFiltrationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return SomaticVariantFiltrationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for DkfzBiasFilterStepPart -----------------------------------------------------------------


def test_dkfz_bias_filter_step_part_get_input_files(somatic_variant_filtration_workflow):
    """Tests DkfzBiasFilterStepPart.get_input_files()"""
    somatic_base_out = (
        "SOMATIC_VARIANT_ANNOTATION/output/"
        "{mapper}.{var_caller}.{annotator}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{annotator}.{tumor_library}"
    )
    expected = {
        "vcf": somatic_base_out + ".vcf.gz",
        "vcf_tbi": somatic_base_out + ".vcf.gz.tbi",
        "bam": "NGS_MAPPING/output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}.bam",
        "bai": "NGS_MAPPING/output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}.bam.bai",
    }
    actual = somatic_variant_filtration_workflow.get_input_files("dkfz_bias_filter", "run")
    assert actual == expected


def test_dkfz_bias_filter_step_part_get_output_files(somatic_variant_filtration_workflow):
    """Tests DkfzBiasFilterStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter."
        "{tumor_library,[^\\.]+}/out/{mapper}.{var_caller}.{annotator}."
        "dkfz_bias_filter.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow.get_output_files("dkfz_bias_filter", "run")
    assert actual == expected


def test_dkfz_bias_filter_step_part_get_log_file(somatic_variant_filtration_workflow):
    """Tests DkfzBiasFilterStepPart.get_log_file()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.{tumor_library}/"
        "log/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow.get_log_file("dkfz_bias_filter", "run")
    assert actual == expected


def test_dkfz_bias_filter_step_part_get_resource_usage(somatic_variant_filtration_workflow):
    """Tests DkfzBiasFilterStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "3-00:00:00", "memory": "3072M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_filtration_workflow.get_resource(
            "dkfz_bias_filter", "run", resource
        )
        assert actual == expected, msg_error


# Tests for EbFilterStepPart -----------------------------------------------------------------


def test_eb_filter_step_part_get_input_files_run(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart._get_input_files_run()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "annotator": "jannovar",
            "tumor_library": "P001-T1-DNA1-WGS1",
        }
    )
    base_out = (
        "work/bwa.mutect2.jannovar.dkfz_bias_filter.P001-T1-DNA1-WGS1/out/"
        "bwa.mutect2.jannovar.dkfz_bias_filter.P001-T1-DNA1-WGS1"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "txt": "work/bwa.eb_filter.panel_of_normals/out/bwa.eb_filter.panel_of_normals.txt",
    }

    actual = somatic_variant_filtration_workflow.get_input_files("eb_filter", "run")(wildcards)
    assert actual == expected


def test_eb_filter_step_part_get_input_files_write_panel(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart._get_input_files_write_panel()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "bam": [
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        ],
        "bai": [
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        ],
    }
    actual = somatic_variant_filtration_workflow.get_input_files("eb_filter", "write_panel")(
        wildcards
    )
    assert actual == expected


def test_eb_filter_step_part_get_output_files_run(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart._get_output_files_run()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter."
        "eb_filter.{tumor_library,[^\\.]+}/out/{mapper}.{var_caller}.{annotator}."
        "dkfz_bias_filter.eb_filter.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow.get_output_files("eb_filter", "run")
    assert actual == expected


def test_eb_filter_step_part_get_output_files_write_panel(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart._get_output_files_run()"""
    expected = {
        "txt": (
            "work/{mapper}.eb_filter.panel_of_normals/out/{mapper}.eb_filter.panel_of_normals.txt"
        )
    }
    actual = somatic_variant_filtration_workflow.get_output_files("eb_filter", "write_panel")
    assert actual == expected


def test_eb_eb_filter_step_part_get_log_file_run(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart._get_log_file_run()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.eb_filter."
        "{tumor_library}/log/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter."
        "eb_filter.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow.get_log_file("eb_filter", "run")
    assert actual == expected


def test_eb_eb_filter_step_part_get_log_file_write_panel(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart._get_log_file_write_panel()"""
    expected = {}
    actual = somatic_variant_filtration_workflow.get_log_file("eb_filter", "write_panel")
    assert actual == expected


def test_eb_filter_step_part_get_resource_usage(somatic_variant_filtration_workflow):
    """Tests EbFilterStepPart.get_resource()"""
    # All actions
    actions = ("run", "write_panel")
    # Define expected
    expected_dict = {"threads": 1, "time": "04:00:00", "memory": "8192M", "partition": "medium"}
    # Evaluate
    for action in actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = somatic_variant_filtration_workflow.get_resource("eb_filter", action, resource)
            assert actual == expected, msg_error


# Tests for ApplyFiltersStepPart -------------------------------------------------------------------


def test_apply_filters_step_part_get_input_files(somatic_variant_filtration_workflow):
    """Tests ApplyFiltersStepPart.get_input_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.eb_filter."
        "{tumor_library}/out/{mapper}.{var_caller}.{annotator}."
        "dkfz_bias_filter.eb_filter.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = somatic_variant_filtration_workflow.get_input_files("apply_filters", "run")
    assert actual == expected


def test_apply_filters_step_part_get_output_files(somatic_variant_filtration_workflow):
    """Tests ApplyFiltersStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.eb_filter."
        "{tumor_library}.{filter_set}.genome_wide/out/{mapper}.{var_caller}."
        "{annotator}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
        "genome_wide"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow.get_output_files("apply_filters", "run")
    assert actual == expected


def test_apply_filters_step_part_get_log_file(somatic_variant_filtration_workflow):
    """Tests ApplyFiltersStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.eb_filter."
        "{tumor_library}.{filter_set}.genome_wide/log/{mapper}.{var_caller}."
        "{annotator}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
        "genome_wide.log"
    )
    actual = somatic_variant_filtration_workflow.get_log_file("apply_filters", "run")
    assert actual == expected


def test_apply_filters_step_part_get_resource_usage(somatic_variant_filtration_workflow):
    """Tests ApplyFiltersStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_filtration_workflow.get_resource("apply_filters", "run", resource)
        assert actual == expected, msg_error


# Tests for FilterToExonsStepPart ------------------------------------------------------------------


# def test_filter_to_exons_step_part_get_input_files(somatic_variant_filtration_workflow):
#     """Tests FilterToExonsStepPart.get_input_files()"""
#     # TODO: fix use of `tumor_library` in get_input_files()
#     wildcards = Wildcards(
#         fromdict={
#             "mapper": "bwa",
#             "var_caller": "mutect2",
#             "filter_set": "custom_exon_list",  # totally made up, not sure it should look like
#             "tumor_library": "P001-T1-DNA1-WGS1",
#         }
#     )
#     base_out = (
#         "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.eb_filter."
#         "{tumor_library}/out/{mapper}.{var_caller}.{annotator}."
#         "dkfz_bias_filter.eb_filter.{tumor_library}"
#     )
#     expected = {
#         "vcf": base_out + ".vcf.gz",
#         "vcf_tbi": base_out + ".vcf.gz.tbi",
#     }
#     # _ = somatic_variant_filtration_workflow.get_input_files("filter_to_exons", "run")(wildcards)
#     _ = expected
#     assert True


def test_filter_to_exons_step_part_get_output_files(somatic_variant_filtration_workflow):
    """Tests FilterToExonsStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter.eb_filter."
        "{tumor_library}.{filter_set}.{exon_list}/out/{mapper}.{var_caller}."
        "{annotator}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
        "{exon_list}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow.get_output_files("filter_to_exons", "run")
    assert actual == expected


def test_filter_to_exons_step_part_get_log_file(somatic_variant_filtration_workflow):
    """Tests FilterToExonsStepPart.get_log_file()"""
    expected = (
        "work/{mapper}.{var_caller}.{annotator}.dkfz_bias_filter."
        "eb_filter.{tumor_library}.{filter_set}.{exon_list}/log/{mapper}.{var_caller}."
        "{annotator}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
        "{exon_list}.log"
    )
    actual = somatic_variant_filtration_workflow.get_log_file("filter_to_exons", "run")
    assert actual == expected


def test_filter_to_exons_step_part_get_resource_usage(somatic_variant_filtration_workflow):
    """Tests FilterToExonsStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "7680M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_filtration_workflow.get_resource(
            "filter_to_exons", "run", resource
        )
        assert actual == expected, msg_error


# Tests for SomaticVariantFiltrationWorkflow -------------------------------------------------------


def test_somatic_variant_filtration_workflow(somatic_variant_filtration_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = [
        "apply_filters",
        "dkfz_bias_filter",
        "eb_filter",
        "filter_to_exons",
        "last_filter",
        "link_out",
        "one_bcftools",
        "one_dkfz",
        "one_ebfilter",
        "one_protected",
        "one_regions",
    ]
    actual = list(sorted(somatic_variant_filtration_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/bwa.mutect2.jannovar.dkfz_bias_filter."
        "eb_filter.P00{i}-T{t}-DNA1-WGS1.{filter}/out/"
        "bwa.mutect2.jannovar.dkfz_bias_filter.eb_filter."
        "P00{i}-T{t}-DNA1-WGS1.{filter}.{ext}"
    )
    expected = [
        tpl.format(mapper=mapper, filter=filter_, i=i, t=t, ext=ext)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for filter_ in (
            "dkfz_and_ebfilter.genome_wide",
            "dkfz_and_ebfilter_and_oxog.genome_wide",
            "dkfz_and_oxog.genome_wide",
            "dkfz_only.genome_wide",
            "no_filter.genome_wide",
        )
    ]
    tpl = (
        "output/bwa.mutect2.jannovar.dkfz_bias_filter."
        "eb_filter.P00{i}-T{t}-DNA1-WGS1.{filter}/log/"
        "bwa.mutect2.jannovar.dkfz_bias_filter.eb_filter."
        "P00{i}-T{t}-DNA1-WGS1.{filter}.{ext}{chksum}"
    )
    expected += [
        tpl.format(mapper=mapper, filter=filter_, i=i, t=t, ext=ext, chksum=chksum)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("log", "conda_list.txt", "conda_info.txt")
        for chksum in ("", ".md5")
        for mapper in ("bwa",)
        for filter_ in (
            "dkfz_and_ebfilter.genome_wide",
            "dkfz_and_ebfilter_and_oxog.genome_wide",
            "dkfz_and_oxog.genome_wide",
            "dkfz_only.genome_wide",
            "no_filter.genome_wide",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_filtration_workflow.get_result_files()))
    assert expected == actual


# Tests for filtration using filter_list -----------------------------------------------------------


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config_list():
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

          somatic_variant_filtration:
            tools_ngs_mapping: ['bwa']
            tools_somatic_variant_calling: ['mutect2']
            tools_somatic_variant_annotation: ['jannovar']
            filter_list:
            - bcftools:
                include: "include_statment"
            - dkfz:
            - ebfilter:
                threshold: 2.3
            - regions:
                path_bed: /path/to/regions.bed
            - protected:
                path_bed: /path/to/protected.bed

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
def somatic_variant_filtration_workflow_list(
    dummy_workflow,
    minimal_config_list,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return SomaticVariantFiltrationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return SomaticVariantFiltrationWorkflow(
        dummy_workflow,
        minimal_config_list,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for one_filter --------------------------------------------------------


def test_one_filter_step_part_get_input_files(somatic_variant_filtration_workflow_list):
    """Tests ApplyFiltersStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "annotator": "jannovar",
            "tumor_library": "P001-T1-DNA1-WGS1",
            "filter_nb": 1,
        }
    )
    expected = {
        "vcf": "../somatic_variant_annotation/output/{mapper}.{var_caller}.{annotator}.{tumor_library}/out/{mapper}.{var_caller}.{annotator}.{tumor_library}.vcf.gz",
        "bam": "../ngs_mapping/output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}.bam",
        "normal": "../ngs_mapping/output/{mapper}.P001-N1-DNA1-WGS1/out/{mapper}.P001-N1-DNA1-WGS1.bam",
    }
    actual = somatic_variant_filtration_workflow_list.get_input_files("one_bcftools", "run")(
        wildcards
    )
    assert actual == expected


def test_one_filter_step_part_get_output_files(somatic_variant_filtration_workflow_list):
    """Tests ApplyFiltersStepPart.get_output_files()"""
    base_out = "work/{mapper}.{var_caller}.{annotator}.{tumor_library}/out/{mapper}.{var_caller}.{annotator}.{tumor_library}.dkfz_{filter_nb}"
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow_list.get_output_files("one_dkfz", "run")
    assert actual == expected


def test_one_filter_step_part_get_log_file(somatic_variant_filtration_workflow_list):
    """Tests ApplyFiltersStepPart.get_log_file()"""
    base_out = "work/{mapper}.{var_caller}.{annotator}.{tumor_library}/log/{mapper}.{var_caller}.{annotator}.{tumor_library}.ebfilter_{filter_nb}"
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow_list.get_log_file("one_ebfilter", "run")
    assert actual == expected


def test_one_filter_step_part_get_resource_usage(somatic_variant_filtration_workflow_list):
    """Tests ApplyFiltersStepPart.get_resource()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "annotator": "jannovar",
            "tumor_library": "P001-T1-DNA1-WGS1",
            "filter_nb": 3,
        }
    )
    # Define expected
    expected_dict = {"threads": 1, "time": "24:00:00", "memory": "2048M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        if resource == "threads" or resource == "partition":
            actual = somatic_variant_filtration_workflow_list.get_resource(
                "one_regions", "run", resource
            )
        else:
            actual = somatic_variant_filtration_workflow_list.get_resource(
                "one_regions", "run", resource
            )(wildcards)
        assert actual == expected, msg_error


# Tests for last_filter -------------------------------------------------------


def test_last_filter_step_part_get_input_files(somatic_variant_filtration_workflow_list):
    """Tests ApplyFiltersStepPart.get_input_files()"""
    base_tpl = "{{mapper}}.{{var_caller}}.{{annotator}}.{{tumor_library}}"
    log_tpl = "work/" + base_tpl + "/log/" + base_tpl + ".{filt}.{ext}{chksum}"
    expected = {
        "logs": [
            log_tpl.format(filt=filt, ext=ext, chksum=chksum)
            for filt in ("bcftools_1", "dkfz_2", "ebfilter_3", "regions_4", "protected_5")
            for ext in ("log", "conda_list.txt", "conda_info.txt")
            for chksum in ("", ".md5")
        ],
        "vcf": "work/{mapper}.{var_caller}.{annotator}.{tumor_library}/out/{mapper}.{var_caller}.{annotator}.{tumor_library}.protected_5.vcf.gz",
    }
    actual = somatic_variant_filtration_workflow_list.get_input_files("last_filter", "run")
    assert actual == expected


def test_last_filter_step_part_get_output_files(somatic_variant_filtration_workflow_list):
    """Tests ApplyFiltersStepPart.get_output_files()"""
    base_name = "{mapper}.{var_caller}.{annotator}.filtered.{tumor_library}"
    base_out = "work/" + base_name + "/out/" + base_name
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    expected["full"] = base_out + ".full.vcf.gz"
    expected["full_tbi"] = expected["full"] + ".tbi"
    expected["full_md5"] = expected["full"] + ".md5"
    expected["full_tbi_md5"] = expected["full_tbi"] + ".md5"
    base_out = "work/" + base_name + "/log/" + base_name
    expected["log"] = base_out + ".merged.tar.gz"
    expected["log_md5"] = expected["log"] + ".md5"
    actual = somatic_variant_filtration_workflow_list.get_output_files("last_filter", "run")
    assert actual == expected


# Tests for SomaticVariantFiltrationWorkflow (filter_list) ---------------------


def test_somatic_variant_filtration_workflow_list(somatic_variant_filtration_workflow_list):
    """Test simple functionality of the workflow"""
    # Check result file construction
    tpl = "output/bwa.mutect2.jannovar.filtered.P00{i}-T{t}-DNA1-WGS1/out/bwa.mutect2.jannovar.filtered.P00{i}-T{t}-DNA1-WGS1.{ext}{chksum}"
    expected = [
        tpl.format(i=i, t=t, ext=ext, chksum=chksum)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.tbi", "full.vcf.gz", "full.vcf.gz.tbi")
        for chksum in ("", ".md5")
    ]
    tpl = "output/bwa.mutect2.jannovar.filtered.P00{i}-T{t}-DNA1-WGS1/log/bwa.mutect2.jannovar.filtered.P00{i}-T{t}-DNA1-WGS1.{ext}{chksum}"
    expected += [
        tpl.format(i=i, t=t, ext=ext, chksum=chksum)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("log", "merged.tar.gz", "conda_list.txt", "conda_info.txt")
        for chksum in ("", ".md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_filtration_workflow_list.get_result_files()))
    assert expected == actual
