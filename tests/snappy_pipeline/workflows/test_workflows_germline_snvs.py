# -*- coding: utf-8 -*-
"""Tests for the somatic_purity_ploidy_estimate workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.germline_snvs import GermlineSnvsWorkflow

from .conftest import patch_module_fs
from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict


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

          germline_snvs:
            tools: ['bcftools']
            ignore_chroms: [ignored]
            bcftools:
              skip_indels: False
              pileup:
                min_BQ: 20
                min_MQ: 35
                skip_any_unset: [PROPER_PAIR, PAIRED]
                extra_args:
                  "illumina1.3+": True
                  "full_BAQ": False
                  "no_BAQ": True
                  "redo_BAQ": False
                  seed: 1234567
              call:
                caller: multiallelic
                extra_args:
                  "prior-freqs": ["AN", "AC"]
              annotate:
              - path_annotation: /path/to/gnomad.vcf
                columns:
                - INFO/POP_AF:=INFO/AF
                - =TAG
                extra_args:
                  merge_logic: TAG:append-missing
                mark_sites: +MARK
              - path_annotation: /path/to/dbsnp.vcf
                columns:
                - ^INFO/AD

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
def germline_snvs_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return GermlineSnvsWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return GermlineSnvsWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for SamtoolsStepPart --------------------------------------------------------------------------


def test_bcftools_step_part_get_input_files_pileup(germline_snvs_workflow):
    """Tests SamtoolsStepPart._get_input_files_pileup()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bams": ["NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam"],
        "bais": ["NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai"],
        "fasta_ref": "/path/to/ref.fa",
        "regions_file": "work/bcftools/out/regions.bed.gz",
    }
    actual = germline_snvs_workflow.get_input_files("bcftools", "pileup")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files_annotate(germline_snvs_workflow):
    """Tests SamtoolsStepPart._get_input_files_annotate()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa", "n": "1"})
    expected = {
        "vcf": "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.annotate_0.P001-T1-DNA1-WGS1.vcf.gz",
        "annotations": "/path/to/dbsnp.vcf",
    }
    actual = germline_snvs_workflow.get_input_files("bcftools", "annotate")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files_remove_unseen(germline_snvs_workflow):
    """Tests SamtoolsStepPart._get_input_files_remove_unseen()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "vcf": "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.annotate_1.P001-T1-DNA1-WGS1.vcf.gz",
    }
    actual = germline_snvs_workflow.get_input_files("bcftools", "remove_unseen")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files(germline_snvs_workflow):
    """Test BcfToolsStepPart._get_input_files (not pileup, annotate or remove_unseen)"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    actions = ("pileup", "call", "filter")
    tpl = "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.{action}.P001-T1-DNA1-WGS1.vcf"
    for i in range(1, len(actions)):
        expected = {"vcf": tpl.format(action=actions[i - 1])}
        actual = germline_snvs_workflow.get_input_files("bcftools", actions[i])(wildcards)
        assert actual == expected


def test_bcftools_step_part_get_output_files(germline_snvs_workflow):
    """Tests SamtoolsStepPart.get_output_files()"""
    for action in ("pileup", "call"):
        tpl = f"work/{{mapper}}.bcftools.{{library_name}}/out/{{mapper}}.{action}.{{library_name}}.vcf"
        expected = {"vcf": tpl}
        actual = germline_snvs_workflow.get_output_files("bcftools", action)
        assert actual == expected
    for action, label in (("filter", "filter"), ("annotate", "annotate_{n}"), ("remove_unseen", "bcftools")):
        base_out = "work/{mapper}.bcftools.{library_name}/out/{mapper}." + label + ".{library_name}"
        expected = get_expected_output_vcf_files_dict(base_out)
        actual = germline_snvs_workflow.get_output_files("bcftools", action)
        assert actual == expected


def test_bcftools_step_part_get_log_file(germline_snvs_workflow):
    """Tests SamtoolsStepPart.get_log_file()"""
    for action in ("pileup", "call", "filter", "annotate", "remove_unseen"):
        if action == "annotate":
            label = "annotate_{n}"
        elif action == "remove_unseen":
            label = "bcftools"
        else:
            label = action
        base_out = f"work/{{mapper}}.bcftools.{{library_name}}/log/{{mapper}}.{label}.{{library_name}}"
        expected = get_expected_log_files_dict(base_out=base_out)
        expected["script"] = base_out + ".sh"
        expected["script_md5"] = expected["script"] + ".md5"
        actual = germline_snvs_workflow.get_log_file("bcftools", action)
        assert actual == expected


def test_bcftools_step_part_get_args_pileup(germline_snvs_workflow):
    """Test BcfToolsStepPart._get_args_pileup()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "extra_args": [
            "--illumina1.3+",
            "--max-depth 5000",
            "--min-BQ 20",
            "--min-MQ 35",
            "--no-BAQ",
            "--seed 1234567",
            "--skip-any-set 3852",
            "--skip-any-unset 3",
        ]
    }
    actual = germline_snvs_workflow.get_args("bcftools", "pileup")(wildcards, [])
    assert actual == expected


def test_bcftools_step_part_get_args_call(germline_snvs_workflow):
    """Test BcfToolsStepPart._get_args_call()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "extra_args": [
            "--multiallelic-caller",
            "--ploidy '2'",
            "--prior-freqs 'AN,AC'",
            "--variants-only",
        ]
    }
    actual = germline_snvs_workflow.get_args("bcftools", "call")(wildcards, [])
    assert actual == expected


def test_bcftools_step_part_get_args_filter(germline_snvs_workflow):
    """Test BcfToolsStepPart._get_args_filter()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "index": True,
        "extra_args": [
            (
                "--include '"
                '(FORMAT/AD[0:0] + FORMAT/AD[0:1] >= 50) & '
                '((0.8181818181818181 * FORMAT/AD[0:0] <= FORMAT/AD[0:1]) & (FORMAT/AD[0:1] <= 1.2222222222222223 * FORMAT/AD[0:0])) & '
                '((FORMAT/GT[0] == "0/1") | (FORMAT/GT[0] == "0|1")) & '
                '(QUAL >= 20.0)'
                "'"
            ),
        ]
    }
    actual = germline_snvs_workflow.get_args("bcftools", "filter")(wildcards, [])
    assert actual == expected


def test_bcftools_step_part_get_args_annotate(germline_snvs_workflow):
    """Test BcfToolsStepPart._get_args_annotate()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa", "n": "0"})
    expected = {
        "index": True,
        "extra_args": [
            "--columns 'INFO/POP_AF:=INFO/AF,=TAG'",
            "--mark-sites '+MARK'",
            "--merge-logic 'TAG:append-missing'",
        ]
    }
    actual = germline_snvs_workflow.get_args("bcftools", "annotate")(wildcards, [])
    assert actual == expected


def test_bcftools_step_part_get_args_remove_unseen(germline_snvs_workflow):
    """Test BcfToolsStepPart._get_args_remove_unseen()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {"index": True,"extra_args": ["--trim-alt-alleles", "--trim-unseen-allele"]}
    actual = germline_snvs_workflow.get_args("bcftools", "remove_unseen")(wildcards, [])
    assert actual == expected


# Tests for GermlineSnvsWorkflow ----------------------------------------------------


def test_germline_snvs_workflow(germline_snvs_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["bcftools", "link_out"]
    actual = list(sorted(germline_snvs_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction for ascat
    log_exts = ("conda_info.txt", "conda_list.txt", "log", "sh")
    mapper = "bwa"
    expected = []
    for lib in ((1, "N1"), (1, "T1"), (2, "N1"), (2, "T1"), (2, "T2")):
        tpl = f"{mapper}.bcftools.P00{lib[0]}-{lib[1]}-DNA1-WGS1"
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz")
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz.md5")
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz.tbi")
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz.tbi.md5")

    expected = set(expected)
    actual = set(germline_snvs_workflow.get_result_files())
    assert actual == expected
