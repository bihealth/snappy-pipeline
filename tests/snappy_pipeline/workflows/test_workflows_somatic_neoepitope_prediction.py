# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml

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
              path_index: /path/to/index
              path_features:
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_variant_calling:
            tools:
            - mutect2
            - scalpel
            scalpel:
              path_target_regions: /path/to/target/regions.bed

          somatic_variant_annotation:
            tools: ["vep"]
            vep:
              path_dir_cache: /path/to/dir/cache

          somatic_neoepitope_prediction:
            path_somatic_variant_annotation: ../somatic_variant_annotation
            path_rna_ngs_mapping: ../ngs_mapping
            tools_somatic_variant_annotation: []
            tools_rna_mapping: []
            tools_ngs_mapping: []
            tools_somatic_variant_calling: []
            max_depth: "4000"
            preparation:
                format: 'snappy_custom'
                mode: 'gene'
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
    mocker,
):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_annotation": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
    }
    # Construct the workflow object
    return SomaticNeoepitopePredictionWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_neoepitope_preparation_step_part_get_input_files(somatic_neoepitope_prediction_workflow):
    # Define expected
    vcf_base_out = (
        "SOMATIC_VARIANT_ANNOTATION/output/{mapper}.{var_caller}.{anno_caller}.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_caller}.{tumor_library}"
    )
    ngs_mapping_base_out = "NGS_MAPPING/output/star.P002-T2-RNA1-mRNA_seq1/out/"
    expected = {
        "vcf": vcf_base_out + ".full.vcf.gz",
        "vcf_tbi": vcf_base_out + ".full.vcf.gz.tbi",
        "expression": ngs_mapping_base_out + "star.P002-T2-RNA1-mRNA_seq1.GeneCounts.tab",
        "bam": ngs_mapping_base_out + "star.P002-T2-RNA1-mRNA_seq1.bam",
        "bai": ngs_mapping_base_out + "star.P002-T2-RNA1-mRNA_seq1.bam.bai",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("neoepitope_preparation", "run")
    assert actual == expected


def test_neoepitope_preparation_step_part_get_output_files(somatic_neoepitope_prediction_workflow):
    base_out = (
        "work/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/out/"
        "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out)
    actual = somatic_neoepitope_prediction_workflow.get_output_files(
        "neoepitope_preparation", "run"
    )
    assert actual == expected


def test_neoepitope_preparation_step_part_get_log_files(somatic_neoepitope_prediction_workflow):
    base_out = (
        "work/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/log/"
        "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_neoepitope_prediction_workflow.get_log_file("neoepitope_preparation", "run")
    assert actual == expected


def test_neoepitope_preparation_step_part_get_resource_usage(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {"threads": 2, "time": "1:00:00", "memory": "6144M"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = somatic_neoepitope_prediction_workflow.get_resource(
            "neoepitope_preparation", "run", resource
        )
        assert actual == expected, msg_error


def test_somatic_neoepitope_prediction_workflow(somatic_neoepitope_prediction_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "neoepitope_preparation"]
    actual = list(sorted(somatic_neoepitope_prediction_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{var_caller}.{anno_caller}.GX.P00{i}-T{t}-DNA1-WGS1/{dir_}/"
        "{mapper}.{var_caller}.{anno_caller}.GX.P00{i}-T{t}-DNA1-WGS1.{ext}"
    )
    expected = [
        tpl.format(
            mapper=mapper,
            var_caller=var_caller,
            anno_caller=anno_caller,
            i=i,
            t=t,
            ext=ext,
            dir_="out",
        )
        for i, t in ((1, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.tbi", "vcf.gz.md5", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for var_caller in ("mutect2", "scalpel")
        for anno_caller in ("vep",)
    ]

    expected += [
        tpl.format(
            mapper=mapper,
            var_caller=var_caller,
            anno_caller=anno_caller,
            i=i,
            t=t,
            ext=ext,
            dir_="log",
        )
        for i, t in ((1, 1), (2, 2))
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
        for anno_caller in ("vep",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_neoepitope_prediction_workflow.get_result_files()))
    assert expected == actual
