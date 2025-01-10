# -*- coding: utf-8 -*-
"""Tests for the somatic_purity_ploidy_estimate workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_variants_for_cnv import SomaticVariantsForCnvWorkflow

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

          somatic_variants_for_cnv:
            tools: ['bcftools']
            ignore_chroms: [ignored]
            bcftools:
              annotate:
              - path_annotation: /path/to/gnomad.vcf
                columns: 
                - INFO/POP_AF:=INFO/AF
                mark_sites: +DB

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
def somatic_variants_for_cnv_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticVariantsForCnvWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "germline_snvs": lambda x: "GERMLINE_SNVS/" + x,
        "somatic_variant_calling": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return SomaticVariantsForCnvWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for SamtoolsStepPart --------------------------------------------------------------------------


def test_bcftools_step_part_get_input_files_merge(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart._get_input_files_merge()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "germline": "GERMLINE_SNVS/output/bwa.bcftools.P001-N1-DNA1-WGS1/out/bwa.bcftools.P001-N1-DNA1-WGS1.vcf.gz",
        "somatic": "SOMATIC_VARIANT_CALLING/output/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1.vcf.gz",
        "regions": "work/bcftools/out/regions.bed.gz",
    }
    actual = somatic_variants_for_cnv_workflow.get_input_files("bcftools", "merge")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files_pileup(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart._get_input_files_pileup()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bams": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        ],
        "bais": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
            "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        ],
        "fasta_ref": "/path/to/ref.fa",
        "regions_file": "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.merge.P001-T1-DNA1-WGS1.bed.gz",
    }
    actual = somatic_variants_for_cnv_workflow.get_input_files("bcftools", "pileup")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files_annotate(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart._get_input_files_annotate()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa", "n": "0"})
    expected = {
        "vcf": "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.pileup.P001-T1-DNA1-WGS1.vcf.gz",
        "annotations": "/path/to/gnomad.vcf",
    }
    actual = somatic_variants_for_cnv_workflow.get_input_files("bcftools", "annotate")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files_flag(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart._get_input_files_flag()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "vcf": "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.annotate_0.P001-T1-DNA1-WGS1.vcf.gz",
        "annotations": "SOMATIC_VARIANT_CALLING/output/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1.vcf.gz",
    }
    actual = somatic_variants_for_cnv_workflow.get_input_files("bcftools", "flag")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_input_files_last(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart._get_input_files_last()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "vcf": "work/bwa.bcftools.P001-T1-DNA1-WGS1/out/bwa.flag.P001-T1-DNA1-WGS1.vcf",
    }
    actual = somatic_variants_for_cnv_workflow.get_input_files("bcftools", "last")(wildcards)
    assert actual == expected


def test_bcftools_step_part_get_output_files(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart.get_output_files()"""
    base_out = "work/{mapper}.bcftools.{library_name}/out/{mapper}."
    actions = {
        "pileup": "pileup.{library_name}",
        "annotate": "annotate_{n}.{library_name}",
        "flag": "flag.{library_name}.vcf",
        "last": "bcftools.{library_name}",
    }
    for action, ext in actions.items():
        if ext.endswith(".vcf"):
            expected = {"vcf": base_out + ext}
        else:
            expected = get_expected_output_vcf_files_dict(base_out + ext)
        actual = somatic_variants_for_cnv_workflow.get_output_files("bcftools", action)
        assert actual == expected


def test_bcftools_step_part_get_log_file(somatic_variants_for_cnv_workflow):
    """Tests SamtoolsStepPart.get_log_file()"""
    for action in ("pileup", "annotate", "flag", "last"):
        if action == "annotate":
            label = "annotate_{n}"
        elif action == "last":
            label = "bcftools"
        else:
            label = action
        base_out = f"work/{{mapper}}.bcftools.{{library_name}}/log/{{mapper}}.{label}.{{library_name}}"
        expected = get_expected_log_files_dict(base_out=base_out)
        expected["script"] = base_out + ".sh"
        expected["script_md5"] = expected["script"] + ".md5"
        actual = somatic_variants_for_cnv_workflow.get_log_file("bcftools", action)
        assert actual == expected


def test_bcftools_step_part_get_args_pileup(somatic_variants_for_cnv_workflow):
    """Test BcfToolsStepPart._get_args_pileup()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "extra_args": [
            "--full-BAQ",
            "--max-depth 5000",
            "--min-BQ 20",
            "--min-MQ 35",
            "--redo-BAQ",
            "--skip-any-set 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY'",
            "--skip-any-unset 'PAIRED,PROPER_PAIR'"
        ],
        "index": True,
    }
    actual = somatic_variants_for_cnv_workflow.get_args("bcftools", "pileup")(wildcards, [])
    assert actual == expected


def test_bcftools_step_part_get_args_annotate(somatic_variants_for_cnv_workflow):
    """Test BcfToolsStepPart._get_args_annotate()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa", "n": "0"})
    expected = {
        "extra_args": [
            "--columns 'INFO/POP_AF:=INFO/AF'",
            "--mark-sites '+DB'",
            "--regions-overlap 2",
        ],
        "index": True,
    }
    actual = somatic_variants_for_cnv_workflow.get_args("bcftools", "annotate")(wildcards, [])
    assert actual == expected


def test_bcftools_step_part_get_args_flag(somatic_variants_for_cnv_workflow):
    """Test BcfToolsStepPart._get_args_flag()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {"extra_args": ["--mark-sites '+SOMATIC'", "--regions-overlap 2"]}
    actual = somatic_variants_for_cnv_workflow.get_args("bcftools", "flag")(wildcards, [])
    assert actual == expected


def test_bcftools_steep_part_get_args_last(somatic_variants_for_cnv_workflow):
    """Test BcfToolsStepPart._get_args_last()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {"extra_args": ["--trim-alt-alleles", "--trim-unseen-allele"], "index": True}
    actual = somatic_variants_for_cnv_workflow.get_args("bcftools", "last")(wildcards, [])
    assert actual == expected


# Tests for GermlineSnvsWorkflow ----------------------------------------------------


def test_somatic_variants_for_cnv_workflow(somatic_variants_for_cnv_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["bcftools", "link_out"]
    actual = list(sorted(somatic_variants_for_cnv_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction for ascat
    mapper = "bwa"
    expected = []
    for lib in ((1, "T1"), (2, "T1"), (2, "T2")):
        tpl = f"{mapper}.bcftools.P00{lib[0]}-{lib[1]}-DNA1-WGS1"
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz")
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz.md5")
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz.tbi")
        expected.append(f"output/{tpl}/out/{tpl}.vcf.gz.tbi.md5")

    expected = set(expected)
    actual = set(somatic_variants_for_cnv_workflow.get_result_files())
    assert actual == expected
