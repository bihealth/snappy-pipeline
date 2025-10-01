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
            tools_somatic_variant_annotation: ['vep']
            path_somatic_variant: "/SOMATIC_VARIANT_ANNOTATION"
            filter_list:
            - dkfz: {}
            - ebfilter:
                ebfilter_threshold: 2.3
            - bcftools:
                include: "include_statment"
            - regions:
                exclude: /path/to/regions.bed
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
# @pytest.mark.filterwarnings("ignore:.*Could not find .*")
def somatic_variant_filtration_workflow_list(
    dummy_workflow,
    minimal_config_list,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticVariantFiltrationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "../ngs_mapping/" + x,
        "somatic_variant": lambda x: "/SOMATIC_VARIANT_ANNOTATION/" + x,
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
    """Tests OneFilterStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "annotator": "vep",
            "tumor_library": "P001-T1-DNA1-WGS1",
            "filter_nb": 1,
        }
    )
    expected = {
        "vcf": "/SOMATIC_VARIANT_ANNOTATION/output/bwa.mutect2.vep.P001-T1-DNA1-WGS1/out/bwa.mutect2.vep.P001-T1-DNA1-WGS1.vcf.gz",
        "bam": "../ngs_mapping/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "reference": "/path/to/ref.fa",
    }
    actual = somatic_variant_filtration_workflow_list.get_input_files("one_dkfz", "run")(wildcards)
    assert actual == expected

    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "annotator": "vep",
            "tumor_library": "P001-T1-DNA1-WGS1",
            "filter_nb": 2,
        }
    )
    expected = {
        "vcf": "work/bwa.mutect2.vep.P001-T1-DNA1-WGS1/out/bwa.mutect2.vep.P001-T1-DNA1-WGS1.dkfz_1.vcf.gz",
        "bam": "../ngs_mapping/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal": "../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "reference": "/path/to/ref.fa",
        "txt": "work/bwa.eb_filter.panel_of_normals/out/bwa.eb_filter.panel_of_normals.txt",
    }
    actual = somatic_variant_filtration_workflow_list.get_input_files("one_ebfilter", "run")(
        wildcards
    )
    assert actual == expected


def test_one_filter_step_part_get_output_files(somatic_variant_filtration_workflow_list):
    """Tests OneFilterStepPart.get_output_files()"""
    base_out = "work/{mapper}.{var_caller}.{annotator}.{tumor_library}/out/{mapper}.{var_caller}.{annotator}.{tumor_library}.dkfz_{filter_nb}"
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow_list.get_output_files("one_dkfz", "run")
    assert actual == expected

    expected = {
        "txt": "work/{mapper}.eb_filter.panel_of_normals/out/{mapper}.eb_filter.panel_of_normals.txt"
    }
    actual = somatic_variant_filtration_workflow_list.get_output_files(
        "one_ebfilter", "write_panel"
    )
    assert actual == expected


def test_one_filter_step_part_get_log_file(somatic_variant_filtration_workflow_list):
    """Tests OneFilterStepPart.get_log_file()"""
    base_out = "work/{mapper}.{var_caller}.{annotator}.{tumor_library}/log/{mapper}.{var_caller}.{annotator}.{tumor_library}.ebfilter_{filter_nb}"
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_filtration_workflow_list.get_log_file("one_ebfilter", "run")
    assert actual == expected


def test_one_filter_step_part_get_args(somatic_variant_filtration_workflow_list):
    """Tests OneFilterStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"filter_nb": 1})
    expected = {"filter_name": "dkfz_1"}
    actual = somatic_variant_filtration_workflow_list.get_args("one_dkfz", "run")(wildcards)
    assert actual == expected

    wildcards = Wildcards(fromdict={"filter_nb": 2})
    expected = {
        "filter_name": "ebfilter_2",
        "ebfilter_threshold": 2.3,
        "has_annotation": True,
        "shuffle_seed": 1,
        "panel_of_normals_size": 25,
        "min_mapq": 20,
        "min_baseq": 15,
        "path_panel_of_normals_sample_list": "",
    }
    actual = somatic_variant_filtration_workflow_list.get_args("one_ebfilter", "run")(wildcards)
    assert actual == expected

    wildcards = Wildcards(fromdict={"filter_nb": 3})
    expected = {"filter_name": "bcftools_3", "include": "include_statment", "exclude": ""}
    actual = somatic_variant_filtration_workflow_list.get_args("one_bcftools", "run")(wildcards)
    assert actual == expected

    wildcards = Wildcards(fromdict={"filter_nb": 4})
    expected = {
        "filter_name": "regions_4",
        "exclude": "/path/to/regions.bed",
        "include": "",
        "path_bed": "",
    }
    actual = somatic_variant_filtration_workflow_list.get_args("one_regions", "run")(wildcards)
    assert actual == expected

    wildcards = Wildcards(fromdict={"filter_nb": 5})
    expected = {"filter_name": "protected_5", "path_bed": "/path/to/protected.bed"}
    actual = somatic_variant_filtration_workflow_list.get_args("one_protected", "run")(wildcards)
    assert actual == expected


def test_one_filter_step_part_get_resource_usage(somatic_variant_filtration_workflow_list):
    """Tests OneFilterStepPart.get_resource()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "24:00:00", "memory": "2048M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_filtration_workflow_list.get_resource(
            "one_ebfilter", "run", resource
        )()
        assert actual == expected, msg_error


# Tests for last_filter -------------------------------------------------------


def test_last_filter_step_part_get_input_files(somatic_variant_filtration_workflow_list):
    """Tests LastFilterStepPart.get_input_files()"""
    base_tpl = "{{mapper}}.{{var_caller}}.{{annotator}}.{{tumor_library}}"
    log_tpl = "work/" + base_tpl + "/log/" + base_tpl + ".{filt}.{ext}{chksum}"
    expected = {
        "logs": [
            log_tpl.format(filt=filt, ext=ext, chksum=chksum)
            for filt in ("dkfz_1", "ebfilter_2", "bcftools_3", "regions_4", "protected_5")
            for ext in ("log", "conda_list.txt", "conda_info.txt")
            for chksum in ("", ".md5")
        ],
        "vcf": "work/{mapper}.{var_caller}.{annotator}.{tumor_library}/out/{mapper}.{var_caller}.{annotator}.{tumor_library}.protected_5.vcf.gz",
    }
    actual = somatic_variant_filtration_workflow_list.get_input_files("last_filter", "run")
    assert actual == expected


def test_last_filter_step_part_get_output_files(somatic_variant_filtration_workflow_list):
    """Tests LastFilterStepPart.get_output_files()"""
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
    tpl = "output/bwa.mutect2.vep.filtered.P00{i}-T{t}-DNA1-WGS1/out/bwa.mutect2.vep.filtered.P00{i}-T{t}-DNA1-WGS1.{ext}{chksum}"
    expected = [
        tpl.format(i=i, t=t, ext=ext, chksum=chksum)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.tbi", "full.vcf.gz", "full.vcf.gz.tbi")
        for chksum in ("", ".md5")
    ]
    tpl = "output/bwa.mutect2.vep.filtered.P00{i}-T{t}-DNA1-WGS1/log/bwa.mutect2.vep.filtered.P00{i}-T{t}-DNA1-WGS1.{ext}{chksum}"
    expected += [
        tpl.format(i=i, t=t, ext=ext, chksum=chksum)
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("log", "merged.tar.gz", "conda_list.txt", "conda_info.txt")
        for chksum in ("", ".md5")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_filtration_workflow_list.get_result_files()))
    assert expected == actual
