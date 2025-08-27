# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
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
            - mutect
            - scalpel
            scalpel:
              path_target_regions: /path/to/target/regions.bed
            mutect: {}

          somatic_variant_annotation:
            path_somatic_variant: /path/to/somatic_variant_calling
            tools: ["jannovar", "vep"]
            jannovar:
              dbnsfp: {}
              flag_off_target: true
              path_jannovar_ser: /path/to/jannover.ser
            vep:
              cache_dir: /path/to/dir/cache

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
def somatic_variant_annotation_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticVariantAnnotationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return SomaticVariantAnnotationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for JannovarAnnotateSomaticVcfStepPart -----------------------------------------------------


def test_jannovar_step_part_get_input_files(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_input_files()"""
    base_out = (
        "SOMATIC_VARIANT_CALLING/output/{mapper}.{var_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
    }
    actual = somatic_variant_annotation_workflow.get_input_files("jannovar", "annotate_somatic_vcf")
    assert actual == expected


def test_jannovar_step_part_get_output_files(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.jannovar.{tumor_library}/out/"
        "{mapper}.{var_caller}.jannovar.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_output_files(
        "jannovar", "annotate_somatic_vcf"
    )
    assert actual == expected


def test_jannovar_step_part_get_log_file(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.jannovar.{tumor_library}/log/"
        "{mapper}.{var_caller}.jannovar.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_log_file("jannovar", "annotate_somatic_vcf")
    assert actual == expected


def test_jannovar_step_part_get_args(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {
        "static_data_config": {
            "reference": {"path": "/path/to/ref.fa"},
            "cosmic": {"path": "/path/to/cosmic.vcf.gz"},
            "dbsnp": {"path": "/path/to/dbsnp.vcf.gz"},
            "dbnsfp": None,
            "features": None,
        },
        "jannovar": {
            "dbnsfp": {
                "col_contig": 1,
                "col_pos": 2,
                "columns": [],
            },
            "flag_off_target": True,
            "path_jannovar_ser": "/path/to/jannover.ser",
            "annotation_tracks_bed": [],
            "annotation_tracks_tsv": [],
            "annotation_tracks_vcf": [],
            "window_length": 50000000,
            "num_jobs": 100,
            "use_profile": True,
            "restart_times": 5,
            "max_jobs_per_second": 10,
            "max_status_checks_per_second": 10,
            "ignore_chroms": ["NC_007605", "hs37d5", "chrEBV", "GL*", "*_decoy", "HLA-*"],
        },
    }
    actual = somatic_variant_annotation_workflow.get_args("jannovar", "annotate_somatic_vcf")
    assert actual == expected


def test_jannovar_step_part_get_resource_usage(somatic_variant_annotation_workflow):
    """Tests JannovarAnnotateSomaticVcfStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_annotation_workflow.get_resource(
            "jannovar", "annotate_somatic_vcf", resource
        )()
        assert actual == expected, msg_error


# Tests for VepAnnotateSomaticVcfStepPart ----------------------------------------------------------


def test_vep_step_part_get_input_files(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_input_files()"""
    base_out = (
        "SOMATIC_VARIANT_CALLING/output/{mapper}.{var_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{tumor_library}"
    )
    expected = {
        "vcf": base_out + ".vcf.gz",
        "vcf_tbi": base_out + ".vcf.gz.tbi",
        "reference": "/path/to/ref.fa",
    }
    actual = somatic_variant_annotation_workflow.get_input_files("vep", "run")
    assert actual == expected


def test_vep_step_part_get_output_files(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.vep.{tumor_library}/out/"
        "{mapper}.{var_caller}.vep.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_out)
    full = {k.replace("vcf", "full"): v.replace(".vcf", ".full.vcf") for k, v in expected.items()}
    expected = {**expected, **full}
    actual = somatic_variant_annotation_workflow.get_output_files("vep", "run")
    assert actual == expected


def test_vep_step_part_get_log_file(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_output_files()"""
    base_out = (
        "work/{mapper}.{var_caller}.vep.{tumor_library}/log/"
        "{mapper}.{var_caller}.vep.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_variant_annotation_workflow.get_log_file("vep", "run")
    assert actual == expected


def test_vep_step_part_get_args(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {
        "config": {
            "cache_dir": "/path/to/dir/cache",
            "species": "homo_sapiens",
            "assembly": "GRCh38",
            "cache_version": "102",
            "tx_flag": "gencode_basic",
            "output_options": ["everything"],
            "buffer_size": 1000,
            "num_threads": 8,
            "pick_order": [
                "biotype",
                "mane_select",
                "mane_plus_clinical",
                "appris",
                "tsl",
                "ccds",
                "canonical",
                "rank",
                "length",
            ]
        },
    }
    actual = somatic_variant_annotation_workflow.get_args("vep", "run")
    assert actual == expected


def test_vep_step_part_get_resource_usage(somatic_variant_annotation_workflow):
    """Tests VepAnnotateSomaticVcfStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 8, "time": "24:00:00", "memory": "16384M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_variant_annotation_workflow.get_resource("vep", "run", resource)()
        assert actual == expected, msg_error


# Tests for SomaticVariantAnnotationWorkflow -------------------------------------------------------


def test_somatic_variant_annotation_workflow(somatic_variant_annotation_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["jannovar", "link_out", "vep"]
    actual = list(sorted(somatic_variant_annotation_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.{annotator}.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "{mapper}.{var_caller}.{annotator}.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(
            mapper=mapper, var_caller=var_caller, annotator=annotator, i=i, t=t, ext=ext, dir_="out"
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect", "scalpel")
        for annotator in ("jannovar", "vep")
    ]
    expected += [
        tpl.format(
            mapper=mapper, var_caller=var_caller, annotator=annotator, i=i, t=t, ext=ext, dir_="out"
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("full.vcf.gz", "full.vcf.gz.md5", "full.vcf.gz.tbi", "full.vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect", "scalpel")
        for annotator in ("vep",)
    ]
    expected += [
        tpl.format(
            mapper=mapper, var_caller=var_caller, annotator=annotator, i=i, t=t, ext=ext, dir_="log"
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
        for var_caller in ("mutect", "scalpel")
        for annotator in ("jannovar", "vep")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_variant_annotation_workflow.get_result_files()))
    assert expected == actual
