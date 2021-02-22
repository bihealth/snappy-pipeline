# -*- coding: utf-8 -*-
"""Tests for the variant_denovo_filtration workflow starting off variant_phasing.
"""

import pytest
import ruamel.yaml as yaml
import textwrap

from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_denovo_filtration import VariantDeNovoFiltrationWorkflow

from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    return yaml.round_trip_load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa
            star:
              path_index: /path/to/star/index
          variant_calling:
            tools:
            - gatk_hc
          variant_denovo_filtration:
            path_variant_phasing: ../variant_phasing

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: germline_variants
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def variant_de_novo_filtration_workflow(
    dummy_workflow,
    minimal_config,
    dummy_cluster_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return VariantCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    germline_sheet_fake_fs.fs.create_file(
        file_path="/path/to/ref.fa.fai",
        contents="1\t249250621\t52\t60\t61\n2\t243199373\t253404903\t60\t61\n",
        create_missing_dirs=True,
    )
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_calling", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_annotation", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_phasing", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_denovo_filtration", germline_sheet_fake_fs, mocker
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VARIANT_CALLING/" + x,
        "variant_annotation": lambda x: "VARIANT_ANNOTATION/" + x,
        "variant_phasing": lambda x: "VARIANT_PHASING/" + x,
    }
    # Construct the workflow object
    return VariantDeNovoFiltrationWorkflow(
        dummy_workflow,
        minimal_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FilterDeNovosStepPart ------------------------------------------------------------------


def test_filter_de_novo_from_variant_phasing_step_part_get_input_files(
    variant_de_novo_filtration_workflow,
):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    expected = {
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
        "tbi": "VARIANT_PHASING/output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "tbi_md5": "VARIANT_PHASING/output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "vcf": "VARIANT_PHASING/output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz",
        "vcf_md5": "VARIANT_PHASING/output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.P001-N1-DNA1-WGS1.vcf.gz.md5",
    }
    actual = variant_de_novo_filtration_workflow.get_input_files("filter_denovo", "run")(wildcards)
    assert actual == expected


def test_filter_de_novo_from_variant_phasing_step_part_get_output_files(
    variant_de_novo_filtration_workflow,
):
    expected = {
        "tbi": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.vcf.gz.md5",
    }
    assert expected == variant_de_novo_filtration_workflow.get_output_files("filter_denovo", "run")


def test_filter_de_novo_from_variant_phasing_step_part_get_log_file(
    variant_de_novo_filtration_workflow,
):
    expected = "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\.]+}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.log"
    assert expected == variant_de_novo_filtration_workflow.get_log_file("filter_denovo", "run")


def test_filter_de_novo_from_variant_phasing_step_part_update_cluster_config(
    variant_de_novo_filtration_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_denovo_filtration_filter_denovo_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for FilterDeNovosHardStepPart --------------------------------------------------------------


def test_filter_de_novo_from_variant_annotationhard_step_part_get_input_files(
    variant_de_novo_filtration_workflow,
):
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = variant_de_novo_filtration_workflow.get_input_files("filter_denovo_hard", "run")
    expected = {
        "tbi": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.vcf.gz.tbi",
        "vcf": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos.{index_library}.vcf.gz",
    }
    assert expected == actual


def test_filter_de_novo_from_variant_annotationhard_step_part_get_output_files(
    variant_de_novo_filtration_workflow,
):
    expected = {
        "summary": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.summary.txt",
        "summary_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.summary.txt.md5",
        "tbi": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.vcf.gz.tbi",
        "tbi_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.vcf.gz.tbi.md5",
        "vcf": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.vcf.gz",
        "vcf_md5": "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\\.]+}/out/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.vcf.gz.md5",
    }
    assert expected == variant_de_novo_filtration_workflow.get_output_files(
        "filter_denovo_hard", "run"
    )


def test_filter_de_novo_from_variant_annotationhard_step_part_get_log_file(
    variant_de_novo_filtration_workflow,
):
    expected = "work/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library,[^\.]+}/log/{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.{index_library}.log"
    assert expected == variant_de_novo_filtration_workflow.get_log_file("filter_denovo_hard", "run")


def test_filter_de_novo_from_variant_annotationhard_step_part_update_cluster_config(
    variant_de_novo_filtration_workflow, dummy_cluster_config
):
    actual = set(dummy_cluster_config["variant_denovo_filtration_filter_denovo_hard_run"].keys())
    expected = {"mem", "time", "ntasks"}
    assert expected == actual


# Tests for VariantDeNovoFiltrationWorkflow ----------------------------------------------------------------


def test_de_novo_filtration_workflow(variant_de_novo_filtration_workflow):
    """Test simple functionality of the workflow"""
    # Perform the tests
    #
    # Check created sub steps
    expected = [
        "collect_msdn",
        "filter_denovo",
        "filter_denovo_hard",
        "link_out",
        "summarize_counts",
        "write_pedigree",
    ]
    assert expected == list(sorted(variant_de_novo_filtration_workflow.sub_steps.keys()))
    # Check result file construction
    expected = [
        "output/bwa.denovo_count_summary/out/bwa.denovo_count_summary.txt",
        "output/bwa.denovo_count_summary/out/bwa.denovo_count_summary.txt.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1.summary.txt",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1.summary.txt.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P001-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1.summary.txt",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1.summary.txt.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1.vcf.gz",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1.vcf.gz.md5",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1.vcf.gz.tbi",
        "output/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1/out/bwa.gatk_hc.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.de_novos_hard.P004-N1-DNA1-WGS1.vcf.gz.tbi.md5",
        "output/bwa.multisite_de_novo/out/bwa.multisite_de_novo.txt",
        "output/bwa.multisite_de_novo/out/bwa.multisite_de_novo.txt.md5",
    ]
    expected = list(sorted(expected))
    actual = list(sorted(variant_de_novo_filtration_workflow.get_result_files()))
    assert expected == actual
