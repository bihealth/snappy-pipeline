# -*- coding: utf-8 -*-
"""Tests for the somatic_purity_ploidy_estimate workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_purity_ploidy_estimate import (
    SomaticPurityPloidyEstimateWorkflow,
)

from .conftest import patch_module_fs
from .common import get_expected_log_files_dict


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

          somatic_purity_ploidy_estimate:
            tools: ['ascat']
            ascat:
              genomeVersion: hg38
              sex:
                source: auto
              allele_counter:
                loci_prefix: /path/to/locii_chr
                allele_prefix: /path/to/alleles_chr
              path_gc_content: /path/to/gc_content.txt
              path_reptiming: /path/to/reptiming.txt

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
def somatic_purity_ploidy_estimate_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    guess_sex_result_fake_fs,
    mocker,
):
    """Return SomaticPurityPloidyEstimateWorkflow object pre-configured with cancer sheet"""
    guess_sex_result_fake_fs.fs.create_file(
        file_path="work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.sex.txt",
        contents="female\n",
        create_missing_dirs=True,
    )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.somatic_purity_ploidy_estimate", guess_sex_result_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return SomaticPurityPloidyEstimateWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for AscatStepPart --------------------------------------------------------------------------


def test_ascat_step_part_get_input_files_build_baf(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_build_baf()"""
    wildcards = Wildcards(
        fromdict={
            "library_name": "P001-T1-DNA1-WGS1",
            "mapper": "bwa",
            "chrom_name": "1"
        }
    )
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "reference": "/path/to/ref.fa",
        "locii": "/path/to/locii_chr1.txt"
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "build_baf")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_prepare_hts(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_prepare_hts()"""
    chromosome_names = list(map(str, range(1, 23))) + ["X"]
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "tumorAlleleCounts": [
            f"work/bwa.ascat.P001-T1-DNA1-WGS1/out/AlleleCounts/P001-T1-DNA1-WGS1_alleleFrequencies_chr{chr}.txt"
            for chr in chromosome_names
        ],
        "normalAlleleCounts": [
            f"work/bwa.ascat.P001-N1-DNA1-WGS1/out/AlleleCounts/P001-N1-DNA1-WGS1_alleleFrequencies_chr{chr}.txt"
            for chr in chromosome_names
        ],
        "alleles": [f"/path/to/alleles_chr{chr}.txt" for chr in chromosome_names],
        "sex": "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.sex.txt",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "prepare_hts")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_run(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_run()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "tumor_logr": "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.Tumor_LogR.txt",
        "tumor_baf": "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.Tumor_BAF.txt",
        "normal_logr": "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.Germline_LogR.txt",
        "normal_baf": "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.Germline_BAF.txt",
        "sex": "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.sex.txt",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "run")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_input_files_guess_sex(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_input_files_guess_sex()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_input_files("ascat", "guess_sex")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_output_files_build_baf(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_build_baf()"""
    expected = {
        "vcf": "work/{mapper}.ascat.{library_name}/out/AlleleCounts/{library_name}_alleleFrequencies_chr{chrom_name}.vcf.gz",
        "txt": "work/{mapper}.ascat.{library_name}/out/AlleleCounts/{library_name}_alleleFrequencies_chr{chrom_name}.txt",
        "vcf_md5": "work/{mapper}.ascat.{library_name}/out/AlleleCounts/{library_name}_alleleFrequencies_chr{chrom_name}.vcf.gz.md5",
        "txt_md5": "work/{mapper}.ascat.{library_name}/out/AlleleCounts/{library_name}_alleleFrequencies_chr{chrom_name}.txt.md5",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "build_baf")
    assert actual == expected


def test_ascat_step_part_get_output_files_prepare_hts(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_prepare_hts()"""
    expected = {
        "tumor_logr": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Tumor_LogR.txt",
        "tumor_baf": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Tumor_BAF.txt",
        "normal_logr": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Germline_LogR.txt",
        "normal_baf": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Germline_BAF.txt",
        "tumor_logr_md5": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Tumor_LogR.txt.md5",
        "tumor_baf_md5": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Tumor_BAF.txt.md5",
        "normal_logr_md5": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Germline_LogR.txt.md5",
        "normal_baf_md5": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.Germline_BAF.txt.md5",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "prepare_hts")
    assert actual == expected


def test_ascat_step_part_get_output_files_run(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_run()"""
    tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}."
    expected = {
        "RData": tpl + "RData",
        "circos": tpl + "circos.txt",
        "goodness_of_fit": tpl + "goodness_of_fit.txt",
        "na": tpl + "na.txt",
        "nb": tpl + "nb.txt",
        "ploidy": tpl + "ploidy.txt",
        "purity": tpl + "purity.txt",
        "segments": tpl + "segments.txt",
        "segments_raw": tpl + "segments_raw.txt",
        "RData_md5": tpl + "RData.md5",
        "circos_md5": tpl + "circos.txt.md5",
        "goodness_of_fit_md5": tpl + "goodness_of_fit.txt.md5",
        "na_md5": tpl + "na.txt.md5",
        "nb_md5": tpl + "nb.txt.md5",
        "ploidy_md5": tpl + "ploidy.txt.md5",
        "purity_md5": tpl + "purity.txt.md5",
        "segments_md5": tpl + "segments.txt.md5",
        "segments_raw_md5": tpl + "segments_raw.txt.md5",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "run")
    assert actual == expected


def test_ascat_step_part_get_output_files_guess_sex(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_output_files_guess_sex()"""
    expected = {
        "table": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.sex.tsv",
        "decision": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.sex.txt",
        "table_md5": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.sex.tsv.md5",
        "decision_md5": "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.sex.txt.md5",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_output_files("ascat", "guess_sex")
    assert actual == expected


def test_ascat_step_part_get_log_file(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart.get_log_file()"""
    # Set test cases
    base_out = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}."
    actions = ("prepare_hts", "run", "guess_sex")
    # Evaluate all actions
    for action in actions:
        err_msg = f"Assert error for action '{action}'."
        expected = get_expected_log_files_dict(base_out=base_out + action)
        if action == "guess_sex":
            expected["script"] = base_out + action + ".sh"
        else:
            expected["script"] = base_out + action + ".R"
        expected["script_md5"] = expected["script"] + ".md5"
        actual = somatic_purity_ploidy_estimate_workflow.get_log_file("ascat", action)
        assert actual == expected, err_msg
    expected = get_expected_log_files_dict(base_out=base_out + "chr{chrom_name}.build_baf")
    expected["script"] = base_out + "chr{chrom_name}.build_baf.sh"
    expected["script_md5"] = expected["script"] + ".md5"
    actual = somatic_purity_ploidy_estimate_workflow.get_log_file("ascat", "build_baf")
    assert actual == expected


def test_ascat_step_part_get_args_build_baf(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_args_build_baf()"""
    expected = {
        "minCounts": 10,
        "min_base_qual": 20,
        "min_map_qual": 35,
        "max_coverage": 8000,
        "exclude_flags": 3852,
        "include_flags": 3,
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_args("ascat", "build_baf")({}, [])
    assert actual == expected


def test_ascat_step_part_get_args_prepare_hts(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_args_prepare_hts()"""
    chromosome_names = list(map(str, range(1, 23))) + ["X"]
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "tumorname": "P001-T1-DNA1-WGS1",
        "normalname": "P001-N1-DNA1-WGS1",
        "genomeVersion": "hg38",
        "minCounts": 10,
        "min_base_qual": 20,
        "min_map_qual": 35,
        "seed": 1234567,
        "chrom_names": chromosome_names,
        "gender": "XX",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_args("ascat", "prepare_hts")(wildcards, [])
    assert actual == expected


def test_ascat_step_part_get_args_run(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_args_run()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "gamma": 1.0,
        "penalty": 70.0,
        "genomeVersion": "hg38",
        "min_ploidy": 1.5,
        "max_ploidy": 5.5,
        "min_purity": 0.1,
        "max_purity": 1.05,
        "seed": 1234567,
        "rho_manual": "NA",
        "psi_manual": "NA",
        "gender": "XX",
    }
    actual = somatic_purity_ploidy_estimate_workflow.get_args("ascat", "run")(wildcards, [])
    assert actual == expected


def test_ascat_step_part_get_args_guess_sex(somatic_purity_ploidy_estimate_workflow):
    """Tests AscatStepPart._get_args_guess_sex()"""
    expected = {}
    actual = somatic_purity_ploidy_estimate_workflow.get_args("ascat", "guess_sex")({}, [])
    assert actual == expected


# Tests for SomaticPurityPloidyEstimateWorkflow ----------------------------------------------------


def test_somatic_purity_ploidy_estimate_workflow(somatic_purity_ploidy_estimate_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["ascat", "link_out"]
    actual = list(sorted(somatic_purity_ploidy_estimate_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction for ascat
    out_exts = ("goodness_of_fit", "na", "nb", "purity", "ploidy", "segments")
    actions = ("prepare_hts", "run")
    log_exts = ("conda_info.txt", "conda_list.txt", "log", "R")
    mapper = "bwa"
    expected = []
    for lib in ((1, 1), (2, 1), (2,2)):
        tpl = f"{mapper}.ascat.P00{lib[0]}-T{lib[1]}-DNA1-WGS1"
        for out_ext in out_exts:
            expected.append(f"output/{tpl}/out/{tpl}.{out_ext}.txt")
            expected.append(f"output/{tpl}/out/{tpl}.{out_ext}.txt.md5")
        for action in actions:
            for log_ext in log_exts:
                expected.append(f"output/{tpl}/log/{tpl}.{action}.{log_ext}")
                expected.append(f"output/{tpl}/log/{tpl}.{action}.{log_ext}.md5")

    expected = set(expected)
    actual = set(somatic_purity_ploidy_estimate_workflow.get_result_files())
    assert actual == expected
