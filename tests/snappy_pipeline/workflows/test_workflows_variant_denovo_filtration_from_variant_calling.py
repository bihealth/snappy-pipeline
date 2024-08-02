# -*- coding: utf-8 -*-
"""Tests for the variant_denovo_filtration workflow starting off variant_calling."""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.variant_denovo_filtration import VariantDeNovoFiltrationWorkflow

from .common import get_expected_output_vcf_files_dict
from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa

          variant_calling:
            tools:
            - gatk3_hc
            gatk3_hc: {}
          variant_denovo_filtration:
            path_variant_calling: ../variant_calling

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
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
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
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_calling", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.variant_denovo_filtration", germline_sheet_fake_fs, mocker
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return VariantDeNovoFiltrationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for FilterDeNovosStepPart ------------------------------------------------------------------


def test_filter_de_novo_from_variant_calling_step_part_get_input_files(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosStepPart.get_input_files()"""
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk3_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    # Define expected
    ngs_mapping_base_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/"
    bam_ped_dict = {
        "bai": ngs_mapping_base_out + "bwa.P001-N1-DNA1-WGS1.bam.bai",
        "bam": ngs_mapping_base_out + "bwa.P001-N1-DNA1-WGS1.bam",
        "ped": "work/write_pedigree.P001-N1-DNA1-WGS1/out/P001-N1-DNA1-WGS1.ped",
    }
    variant_calling_base_out = (
        "VARIANT_CALLING/output/bwa.gatk3_hc.P001-N1-DNA1-WGS1/out/bwa.gatk3_hc.P001-N1-DNA1-WGS1"
    )
    vcf_dict = get_expected_output_vcf_files_dict(base_out=variant_calling_base_out)
    expected = {**bam_ped_dict, **vcf_dict}
    # Get actual
    actual = variant_de_novo_filtration_workflow.get_input_files("filter_denovo", "run")(wildcards)
    assert actual == expected


def test_filter_de_novo_from_variant_calling_step_part_get_output_files(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosStepPart.get_output_files()"""
    base_name_out = (
        r"work/{mapper}.{caller}.de_novos.{index_library,[^\.]+}/out/"
        r"{mapper}.{caller}.de_novos.{index_library}"
    )
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    actual = variant_de_novo_filtration_workflow.get_output_files("filter_denovo", "run")
    assert actual == expected


def test_filter_de_novo_from_variant_calling_step_part_get_log_file(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.de_novos.{index_library,[^\.]+}/log/"
        r"{mapper}.{caller}.de_novos.{index_library}.log"
    )
    actual = variant_de_novo_filtration_workflow.get_log_file("filter_denovo", "run")
    assert actual == expected


def test_filter_de_novo_from_variant_calling_step_part_get_resource(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "14336M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_de_novo_filtration_workflow.get_resource(
            "filter_denovo", "run", resource
        )()
        assert actual == expected, msg_error


# Tests for FilterDeNovosHardStepPart --------------------------------------------------------------


def test_filter_de_novo_from_variant_annotationhard_step_part_get_input_files(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosHardStepPart.get_input_files()"""
    # Define expected
    base_name_out = (
        r"work/{mapper}.{caller}.de_novos.{index_library,[^\.]+}/out/"
        r"{mapper}.{caller}.de_novos.{index_library}"
    )
    expected = {
        "vcf_tbi": base_name_out + ".vcf.gz.tbi",
        "vcf": base_name_out + ".vcf.gz",
    }
    # Get actual
    actual = variant_de_novo_filtration_workflow.get_input_files("filter_denovo_hard", "run")
    assert actual == expected


def test_filter_de_novo_from_variant_annotationhard_step_part_get_output_files(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosHardStepPart.get_output_files()"""
    # Define expected
    base_out = (
        r"work/{mapper}.{caller}.de_novos_hard.{index_library,[^\.]+}/out/"
        r"{mapper}.{caller}.de_novos_hard.{index_library}"
    )
    summary_dict = {
        "summary": base_out + ".summary.txt",
        "summary_md5": base_out + ".summary.txt.md5",
    }
    vcf_dict = get_expected_output_vcf_files_dict(base_out=base_out)
    expected = {**summary_dict, **vcf_dict}
    # Get actual
    actual = variant_de_novo_filtration_workflow.get_output_files("filter_denovo_hard", "run")
    assert actual == expected


def test_filter_de_novo_from_variant_annotationhard_step_part_get_log_file(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosHardStepPart.get_log_file()"""
    expected = (
        r"work/{mapper}.{caller}.de_novos_hard.{index_library,[^\.]+}/log/"
        r"{mapper}.{caller}.de_novos_hard.{index_library}.log"
    )
    actual = variant_de_novo_filtration_workflow.get_log_file("filter_denovo_hard", "run")
    assert actual == expected


def test_filter_de_novo_from_variant_annotationhard_step_part_get_resource(
    variant_de_novo_filtration_workflow,
):
    """Tests FilterDeNovosHardStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 1, "time": "1-00:00:00", "memory": "14336M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = variant_de_novo_filtration_workflow.get_resource(
            "filter_denovo", "run", resource
        )()
        assert actual == expected, msg_error


# Tests for VariantDeNovoFiltrationWorkflow --------------------------------------------------------


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
    actual = list(sorted(variant_de_novo_filtration_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = [
        "output/bwa.denovo_count_summary/out/bwa.denovo_count_summary.txt",
        "output/bwa.denovo_count_summary/out/bwa.denovo_count_summary.txt.md5",
        "output/bwa.multisite_de_novo/out/bwa.multisite_de_novo.txt",
        "output/bwa.multisite_de_novo/out/bwa.multisite_de_novo.txt.md5",
    ]
    base_name_out = (
        "output/bwa.gatk3_hc.de_novos_hard.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gatk3_hc.de_novos_hard.P00{i}-N1-DNA1-WGS1{ext}"
    )
    expected += [
        base_name_out.format(i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in (
            ".summary.txt",
            ".summary.txt.md5",
            ".vcf.gz",
            ".vcf.gz.md5",
            ".vcf.gz.tbi",
            ".vcf.gz.tbi.md5",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(variant_de_novo_filtration_workflow.get_result_files()))
    assert expected == actual
