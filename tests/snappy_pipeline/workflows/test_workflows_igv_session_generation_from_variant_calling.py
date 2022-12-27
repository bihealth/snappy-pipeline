# -*- coding: utf-8 -*-
"""Tests for the igv_session_generation workflow starting off variant_calling."""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.igv_session_generation import IgvSessionGenerationWorkflow

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
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          variant_calling:
            tools:
            - gatk3_hc

          igv_session_generation:
            path_ngs_mapping: ../ngs_mapping
            path_variant_calling: ../variant_calling
            tools_variant_calling: ['gatk3_hc']

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
def igv_session_generation(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return VariantCallingWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.variant_calling", germline_sheet_fake_fs, mocker)
    patch_module_fs(
        "snappy_pipeline.workflows.igv_session_generation", germline_sheet_fake_fs, mocker
    )
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VARIANT_CALLING/" + x,
    }
    # Construct the workflow object
    return IgvSessionGenerationWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for WriteIgvSessionFileStepPart ------------------------------------------------------------


def test_igv_session_generation_from_variant_calling_step_part_get_input_files(
    igv_session_generation,
):
    """Tests WriteIgvSessionFileStepPart.get_input_files()"""
    # Define expected
    ngs_mapping_base_out = "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/"
    variant_calling_base_out = (
        "VARIANT_CALLING/output/bwa.gatk3_hc.P001-N1-DNA1-WGS1/out/bwa.gatk3_hc.P001-N1-DNA1-WGS1"
    )
    expected = {
        # TODO: underline method is probably wrong, hence the weirdness of this test (`[n] * 3`)
        "bam": [ngs_mapping_base_out + "bwa.P001-N1-DNA1-WGS1.bam"] * 3,
        "vcf": variant_calling_base_out + ".vcf.gz",
    }
    # Get actual
    wildcards = Wildcards(
        fromdict={"mapper": "bwa", "caller": "gatk3_hc", "index_library": "P001-N1-DNA1-WGS1"}
    )
    actual = igv_session_generation.get_input_files("write_igv_session_file", "run")(wildcards)
    assert actual == expected


def test_igv_session_generation_from_variant_calling_step_part_get_output_files(
    igv_session_generation,
):
    """Tests WriteIgvSessionFileStepPart.get_output_files()"""

    # Define expected
    base_name_out = "work/{mapper}.{caller}.{index_library}/out/{mapper}.{caller}.{index_library}"
    expected = {
        "xml": base_name_out + ".igv_session.xml",
        "xml_md5": base_name_out + ".igv_session.xml.md5",
    }

    # Get actual
    actual = igv_session_generation.get_output_files("write_igv_session_file", "run")
    assert actual == expected


def test_igv_session_generation_from_variant_calling_step_part_get_resource_usage(
    igv_session_generation,
):
    """Tests WriteIgvSessionFileStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = igv_session_generation.get_resource("write_igv_session_file", "run", resource)
        assert actual == expected, msg_error


# Tests for IgvSessionGenerationWorkflow   ---------------------------------------------------------


def test_igv_session_generation_workflow(igv_session_generation):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "write_igv_session_file"]
    actual = list(sorted(igv_session_generation.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    base_name_out = (
        "output/bwa.gatk3_hc.P00{i}-N1-DNA1-WGS1/out/"
        "bwa.gatk3_hc.P00{i}-N1-DNA1-WGS1.igv_session{ext}"
    )
    expected = [
        base_name_out.format(i=i, ext=ext)
        for i in (1, 4)  # only for indices
        for ext in (
            ".xml",
            ".xml.md5",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(igv_session_generation.get_result_files()))
    print(actual)
    assert expected == actual
