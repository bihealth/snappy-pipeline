# -*- coding: utf-8 -*-
"""Tests for the somatic_purity_ploidy_estimate workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_cnv_calling import (
    SomaticCnvCallingWorkflow,
)
from snappy_pipeline.models.bcftools import BcftoolsBamFlag

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

          guess_sex:
            tools: [samtools]
              
          somatic_cnv_calling:
            tools: ['ascat']
            path_target_interval_list_mapping:
              - name: "Agilent SureSelect Human All Exon V6"
                pattern: pattern
                path: /path/to/baits.bed
            ignore_chroms: ["hs37d5", "MT"]
            sex:
              source: auto
            ascat:
              ignore_chroms: ["X"]
              genomeVersion: hg19
              allele_counter:
                loci_prefix: /path/to/loci_chr
                allele_prefix: /path/to/alleles_chr
                skip_any_set: [UNMAP, SECONDARY, DUP, QCFAIL, 2048]
                skip_any_unset: []
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
def somatic_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticCnvCallingWorkflow object pre-configured with cancer sheet"""
    cancer_sheet_fake_fs.fs.create_file(
        "/GUESS_SEX/output/bwa.samtools.P001-T1-DNA1-WGS1/out/bwa.samtools.P001-T1-DNA1-WGS1.txt",
        contents="female",
        create_missing_dirs=True,
    )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.somatic_cnv_calling", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "guess_sex": lambda x: "/GUESS_SEX/" + x,
    }
    # Construct the workflow object
    return SomaticCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for AscatStepPart --------------------------------------------------------------------------


def test_ascat_step_part_get_input_files_alleleCounter(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_input_files_alleleCounter()"""
    wildcards = Wildcards(
        fromdict={
            "library_name": "P001-T1-DNA1-WGS1",
            "mapper": "bwa",
            "chrom_name": "21",
        }
    )
    expected = {
        "reference": "/path/to/ref.fa",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "loci": "/path/to/loci_chr21.txt",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("ascat", "alleleCounter")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_args_alleleCounter(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_args_alleleCounter()"""
    wildcards = Wildcards(
        fromdict={
            "library_name": "P001-T1-DNA1-WGS1",
            "mapper": "bwa",
            "chrom_name": "21",
        }
    )
    skip_any_set = (
        BcftoolsBamFlag.UNMAP
        | BcftoolsBamFlag.SECONDARY
        | BcftoolsBamFlag.DUP
        | BcftoolsBamFlag.QCFAIL
        | BcftoolsBamFlag.SUPPLEMENTARY
    )
    skip_any_unset = BcftoolsBamFlag.NONE
    expected = {
            "seed": 1234567,
            "max-depth": 8000,
            "min-MQ": 35,
            "min-BQ": 20,
            "skip-any-set": skip_any_set.value,
            "skip-any-unset": skip_any_unset.value,
    }
    actual = somatic_cnv_calling_workflow.get_args("ascat", "alleleCounter")(
        wildcards, []
    )
    assert actual == expected


def test_ascat_step_part_get_output_files_alleleCounter(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_output_files_alleleCounter()"""
    tpl = "work/{mapper}.ascat.{tumor_name}/tmp/{library_name}_alleleFrequencies_chr{chrom_name}"
    expected = {"freq": tpl + ".txt", "vcf": tpl + ".vcf.gz"}
    actual = somatic_cnv_calling_workflow.get_output_files("ascat", "alleleCounter")
    assert actual == expected


def test_ascat_step_part_get_input_files_prepareHTS(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_input_files_prepareHTS()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "reference": "/path/to/ref.fa",
        "tumor_alleleFrequencies": [
            f"work/bwa.ascat.P001-T1-DNA1-WGS1/tmp/P001-T1-DNA1-WGS1_alleleFrequencies_chr{chrom_name}.txt"
            for chrom_name in list(map(str, range(1, 23))) 
        ],
        "normal_alleleFrequencies": [
            f"work/bwa.ascat.P001-T1-DNA1-WGS1/tmp/P001-N1-DNA1-WGS1_alleleFrequencies_chr{chrom_name}.txt"
            for chrom_name in list(map(str, range(1, 23)))
        ],
        "alleles": [
            f"/path/to/alleles_chr{chrom_name}.txt"
            for chrom_name in list(map(str, range(1, 23)))
        ],
        "sex": "/GUESS_SEX/output/bwa.samtools.P001-T1-DNA1-WGS1/out/bwa.samtools.P001-T1-DNA1-WGS1.txt",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("ascat", "prepareHTS")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_args_prepareHTS(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_args_prepareHTS()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "genomeVersion": "hg19",
        "minCounts": 10,
        "seed": 1234567,
        "chrom_names": list(map(str, range(1, 23))),
        "gender": "XX",
        "tumorname": "P001-T1-DNA1-WGS1",
        "normalname": "P001-N1-DNA1-WGS1",
    }
    actual = somatic_cnv_calling_workflow.get_args("ascat", "prepareHTS")(wildcards, [])
    assert actual == expected


def test_ascat_step_part_get_output_files_prepareHTS(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_output_files_prepareHTS()"""
    tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}."
    expected = {
        "tumor_logr": tpl + "Tumor_LogR.txt",
        "tumor_baf": tpl + "Tumor_BAF.txt",
        "normal_logr": tpl + "Germline_LogR.txt",
        "normal_baf": tpl + "Germline_BAF.txt",
        "tumor_logr_md5": tpl + "Tumor_LogR.txt.md5",
        "tumor_baf_md5": tpl + "Tumor_BAF.txt.md5",
        "normal_logr_md5": tpl + "Germline_LogR.txt.md5",
        "normal_baf_md5": tpl + "Germline_BAF.txt.md5",
    }
    actual = somatic_cnv_calling_workflow.get_output_files("ascat", "prepareHTS")
    assert actual == expected


def test_ascat_step_part_get_input_files_run(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_input_files_run()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    tpl = "work/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1."
    expected = {
        "GCcontent": "/path/to/gc_content.txt",
        "reptiming": "/path/to/reptiming.txt",
        "tumor_logr": tpl + "Tumor_LogR.txt",
        "tumor_baf": tpl + "Tumor_BAF.txt",
        "normal_logr": tpl + "Germline_LogR.txt",
        "normal_baf": tpl + "Germline_BAF.txt",
        "sex": "/GUESS_SEX/output/bwa.samtools.P001-T1-DNA1-WGS1/out/bwa.samtools.P001-T1-DNA1-WGS1.txt",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("ascat", "run")(
        wildcards
    )
    assert actual == expected


def test_ascat_step_part_get_args_run(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_args_run()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "genomeVersion": "hg19",
        "seed": 1234567,
        "gender": "XX",
        "y_limit": 5.0,
        "penalty": 70.0,
        "gamma": 1.0,
        "min_purity": 0.1,
        "max_purity": 1.05,
        "min_ploidy": 1.5,
        "max_ploidy": 5.5,
        "rho_manual": "NA",
        "psi_manual": "NA",
        "chrom_names": list(map(str, range(1, 23))),
    }
    actual = somatic_cnv_calling_workflow.get_args("ascat", "run")(wildcards, [])
    assert actual == expected


def test_ascat_step_part_get_output_files_run(somatic_cnv_calling_workflow):
    """Tests AscatStepPart._get_output_files_run()"""
    tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}."
    expected = {
        "RData": tpl + "RData",
        "goodness_of_fit": tpl + "goodness_of_fit.txt",
        "na": tpl + "na.txt",
        "nb": tpl + "nb.txt",
        "ploidy": tpl + "ploidy.txt",
        "purity": tpl + "purity.txt",
        "segments": tpl + "segments.txt",
        "segments_raw": tpl + "segments_raw.txt",
        "RData_md5": tpl + "RData.md5",
        "goodness_of_fit_md5": tpl + "goodness_of_fit.txt.md5",
        "na_md5": tpl + "na.txt.md5",
        "nb_md5": tpl + "nb.txt.md5",
        "ploidy_md5": tpl + "ploidy.txt.md5",
        "purity_md5": tpl + "purity.txt.md5",
        "segments_md5": tpl + "segments.txt.md5",
        "segments_raw_md5": tpl + "segments_raw.txt.md5",
    }
    actual = somatic_cnv_calling_workflow.get_output_files("ascat", "run")
    assert actual == expected


def test_ascat_step_part_get_log_file(somatic_cnv_calling_workflow):
    """Tests AscatStepPart.get_log_file()"""
    tpl = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}."
    for action in ("prepareHTS", "run"):
        err_msg = f"Assert error for action '{action}'."
        expected = {}
        for k, v in (
            ("conda_info", "conda_info.txt"),
            ("conda_list", "conda_list.txt"),
            ("log", "log"),
            ("script", "R"),
        ):
            expected[k] = tpl + action + "." + v
            expected[k + "_md5"] = tpl + action + "." + v + ".md5"
        actual = somatic_cnv_calling_workflow.get_log_file("ascat", action)
        assert actual == expected, err_msg
    
    tpl = "work/{mapper}.ascat.{tumor_name}/log/alleleCounter.{library_name}.chr_{chrom_name}."
    for k, v in (
        ("conda_info", "conda_info.txt"),
        ("conda_list", "conda_list.txt"),
        ("log", "conda_info.txt"),
        ("script", "sh"),
    ):
        expected[v] = tpl + v
        expected[k + "_md5"] = tpl + v + ".md5"
    actual = somatic_cnv_calling_workflow.get_log_file("ascat", "alleleCounter")


# Tests for SomaticCnvCallingWorkflow ----------------------------------------------------


def test_somatic_cnv_calling_workflow(somatic_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["ascat", "link_out"]
    actual = list(sorted(somatic_cnv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # ASCAT
    expected = []
    tpl = "output/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.{ext}.txt{md5}"
    for mapper in ("bwa",):
        for lib in ("P001-T1-DNA1-WGS1", "P002-T1-DNA1-WGS1", "P002-T2-DNA1-WGS1"):
            for ext in ("goodness_of_fit", "na", "nb", "ploidy", "purity", "segments"):
                for md5 in ("", ".md5"):
                    expected.append(tpl.format(mapper=mapper, library_name=lib, ext=ext, md5=md5))

    tpl = "output/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}.{action}.{ext}{md5}"
    for mapper in ("bwa",):
        for lib in ("P001-T1-DNA1-WGS1", "P002-T1-DNA1-WGS1", "P002-T2-DNA1-WGS1"):
            for action in ("prepareHTS", "run"):
                for ext in ("conda_info.txt", "conda_list.txt", "log", "R"):
                    for md5 in ("", ".md5"):
                        expected.append(
                            tpl.format(
                                mapper=mapper, library_name=lib, action=action, ext=ext, md5=md5
                            )
                        )

    actual = somatic_cnv_calling_workflow.get_result_files()
    assert set(actual) == set(expected)
