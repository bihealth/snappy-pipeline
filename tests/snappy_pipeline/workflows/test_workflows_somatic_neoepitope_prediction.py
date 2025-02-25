# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_neoepitope_prediction import (
    SomaticNeoepitopePredictionWorkflow,
)

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
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
          features:
            path: /path/to/genecode.gtf

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
              rna: ['star']
            star:
              path_index: /path/to/star/index
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_variant_calling:
            tools:
            - mutect2
            - scalpel
            mutect2: {}
            scalpel:
              path_target_regions: /path/to/target/regions.bed

          hla_typing:
            path_link_in: ''  # OPTIONAL Override data set configuration search paths for FASTQ files
            tools: [optitype]   # REQUIRED - available: 'optitype' and 'arcashla'
            optitype:
                max_reads: 5000
                num_mapping_threads: 4

          somatic_variant_annotation:
            path_somatic_variant_calling: ../somatic_variant_calling
            tools: ["vep"]
            vep:
                cache_dir: /path/to/dir/cache
                species: homo_sapiens_merged
                assembly: GRCh38
                cache_version: '102'
                tx_flag: gencode_basic  
                num_threads: 8
                buffer_size: 500
                output_options:
                - everything
          somatic_neoepitope_prediction:
            path_somatic_variant_annotation: ../somatic_variant_annotation  # REQUIRED
            path_container: ../somatic_neoepitope_prediction/work/containers/out/pvactools.simg
            path_rna_ngs_mapping: ../ngs_mapping
            path_hla_typing: ../hla_typing
            tools_somatic_variant_annotation: [vep]
            preparation:
                format: snappy_custom    # REQUIRED - The file format of the expression file to process. (stringtie, kallisto, cufflinks, snappy_custom, custom)
                mode: gene   # REQUIRED - Determine whether the expression file contains gene or transcript TPM values.
            prediction:
                algorithms: ['MHCflurry','MHCnuggetsI']
                
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
def somatic_neoepitope_prediction_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping.model", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there

    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_annotation": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
        "hla_typing": lambda x: "HLA_TYPING/" + x,
    }
    # Construct the workflow object
    return SomaticNeoepitopePredictionWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_somatic_neoepitope_preparation_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "anno_caller": "vep",
            "tumor_library": "P001-T1-DNA1-WGS1",
        }
    )
    # Define expected
    vcf_base_out = (
        "SOMATIC_VARIANT_ANNOTATION/output/{mapper}.{var_caller}.{anno_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_caller}.{tumor_library}"
    )
    # vcf_base_out = (
    #     "SOMATIC_VARIANT_ANNOTATION/output/bwa.mutect2.vep.P001-T1-DNA1-WES1/out/"
    #     "{mapper}.{var_caller}.{anno_caller}.{tumor_library}"
    # )
    ngs_mapping_base_out = "NGS_MAPPING/output/star.P001-T1-RNA1-mRNA_seq1/out/"
    expected = {
        "vcf": vcf_base_out + ".full.vcf.gz",
        "vcf_tbi": vcf_base_out + ".full.vcf.gz.tbi",
        "expression": ngs_mapping_base_out + "star.P001-T1-RNA1-mRNA_seq1.GeneCounts.tab",
        "bam": ngs_mapping_base_out + "star.P001-T1-RNA1-mRNA_seq1.bam",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "prepare")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "mutect2",
            "anno_caller": "vep",
            "tumor_library": "P001-T1-DNA1-WGS1",
        }
    )
    vcf_base_out = (
        "work/prepare/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}"
    )
    hla_typing_base_out = "HLA_TYPING/output/"
    expected = {
        "container": "work/containers/out/pvactools.simg",
        "combine_vcf": vcf_base_out + ".vcf.gz",
        "hla_normal_dna": hla_typing_base_out
        + "optitype.P001-N1-DNA1-WGS1/out/optitype.P001-N1-DNA1-WGS1.txt",
        "hla_tumor_dna": hla_typing_base_out
        + "optitype.P001-T1-DNA1-WGS1/out/optitype.P001-T1-DNA1-WGS1.txt",
        "hla_tumor_rna": hla_typing_base_out
        + "optitype.P001-T1-RNA1-mRNA_seq1/out/optitype.P001-T1-RNA1-mRNA_seq1.txt",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "predict")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_preparation_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    base_out = (
        "work/prepare/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out)
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "prepare")
    assert actual == expected


def test_somatic_neoepitope_preparation_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    base_out = (
        "work/prepare/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/log/"
        "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "prepare")
    assert actual == expected


def test_somatic_neoepitope_preparation_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "6G"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_neoepitope_prediction_workflow.get_resource(
            "pvacseq", "prepare", resource
        )()
        assert actual == expected, msg_error


def test_somatic_neoepitope_prediction_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    base_out = "work/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}.epitopes/out/MHC_Class_I/{tumor_library}"
    expected = {
        "all_epitopes": base_out + ".all_epitopes.tsv",
        "all_epitopes_md5": base_out + ".all_epitopes.tsv.md5",
        "filtered_epitopes": base_out + ".filtered.tsv",
        "filtered_epitopes_md5": base_out + ".filtered.tsv.md5",
    }
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "predict")
    assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    base_out = (
        "work/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}.epitopes/"
        "log/{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "predict")
    assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {"threads": 4, "time": "4:00:00", "memory": "30G"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_neoepitope_prediction_workflow.get_resource(
            "pvacseq", "predict", resource
        )()
        assert actual == expected, msg_error


def test_somatic_neoepitope_prediction_workflow(somatic_neoepitope_prediction_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "pvacseq"]
    actual = list(sorted(somatic_neoepitope_prediction_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.P00{i}-T{t}-DNA1-WGS1.epitopes/"
        "out/MHC_Class_I/P00{i}-T{t}-DNA1-WGS1{ext}"
    )
    log_tpl = (
        "output/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.P00{i}-T{t}-DNA1-WGS1.epitopes/"
        "log/P00{i}-T{t}-DNA1-WGS1{ext}"
    )
    expected = [
        log_tpl.format(
            mapper="bwa", var_caller=var_caller, anno_caller="vep", i=i, t=t, ext=ext, mode="GX"
        )
        for (i, t) in ((1, 1), (2, 2))
        for var_caller in ("mutect2", "scalpel")
        for ext in (
            ".log",
            ".log.md5",
            ".conda_info.txt",
            ".conda_info.txt.md5",
            ".conda_list.txt",
            ".conda_list.txt.md5",
        )
    ]

    expected += [
        tpl.format(
            mapper="bwa", var_caller=var_caller, anno_caller="vep", i=i, t=t, ext=ext, mode="GX"
        )
        for (i, t) in ((1, 1), (2, 2))
        for var_caller in ("mutect2", "scalpel")
        for ext in (
            ".all_epitopes.tsv",
            ".filtered.tsv",
            ".all_epitopes.tsv.md5",
            ".filtered.tsv.md5",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_neoepitope_prediction_workflow.get_result_files()))
    assert expected == actual
