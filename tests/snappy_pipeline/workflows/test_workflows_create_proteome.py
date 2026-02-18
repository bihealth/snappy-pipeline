# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.create_proteome import CreateProteomeWorkflow

from .common import get_expected_log_files_dict
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
          features:
            path: /path/to/features.gtf

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa
          
          germline_variant_calling:
            tools: [gatk4_hc]
            gatk4_hc: {}

          germline_variant_filtration:
            path_variant: ../germline_variant_calling
            has_annotation: false
            filter_list:
            - bcftools:
                include: 'depth > min'
        
          create_proteome:
            variant_type: filtration
            path_variant: GERMLINE_VARIANT_FILTRATION
            is_filtered: true
            path_proteome: /path/to/aa.fa

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
def create_proteome_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return CreateProteomeWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "variant": lambda x: "GERMLINE_VARIANT_FILTRATION/" + x,
    }
    # Construct the workflow object
    return CreateProteomeWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CreateProteomeStepPart ----------------------------------------------------------


def test_combine_step_part_get_input_files_run(create_proteome_workflow):
    """Tests CreateProteomeStepPart.get_input_files()"""
    wildcards: Wildcards = Wildcards(fromdict={"library": "P001-N1-DNA1-WGS1"})
    base_out = "{dir}/output/bwa.{caller}.{library}/out/bwa.{caller}.{library}.vcf.gz"
    expected = {
        "reference": "/path/to/ref.fa",
        "features": "/path/to/features.gtf",
        "proteome": "/path/to/aa.fa",
        "vcf": base_out.format(dir="GERMLINE_VARIANT_FILTRATION", caller="gatk4_hc.filtered", library="P001-N1-DNA1-WGS1"),
    }
    actual = create_proteome_workflow.get_input_files("create_proteome", "run")(wildcards)
    assert actual == expected


def test_combine_step_part_get_output_files_run(create_proteome_workflow):
    """Tests CreateProteomeStepPart.get_output_files()"""
    base_out = "bwa.gatk4_hc.filtered.{library}"
    expected = {"vcf": f"work/{base_out}/out/{base_out}.vcf.gz"}
    actual = create_proteome_workflow.get_output_files("create_proteome", "run")
    assert actual == expected


def test_combine_step_part_get_log_file_run(create_proteome_workflow):
    """Tests CreateProteomeStepPart.get_output_files()"""
    base_out = "bwa.gatk4_hc.filtered.{library}"
    expected = get_expected_log_files_dict(base_out=f"work/{base_out}/log/{base_out}")
    actual = create_proteome_workflow.get_log_file("create_proteome", "run")
    assert actual == expected


def test_combine_step_part_get_args_run(create_proteome_workflow):
    """Tests CreateProteomeStepPart.get_args()"""
    wildcards: Wildcards = Wildcards(fromdict={"library": "P001-T1-DNA1-WGS1"})
    expected = {"add_reference": False}
    actual = create_proteome_workflow.get_args("create_proteome", "run")(wildcards)
    assert actual == expected


def test_combine_step_part_get_resource_usage(create_proteome_workflow):
    """Tests CreateProteomeStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {
        "run": {"threads": 1, "time": "03:59:59", "memory": "4G", "partition": "medium"}
    }
    # Evaluate
    for action in ("run",):
        for resource, expected in expected_dict[action].items():
            msg_error = f"Assertion error for resource '{resource}'."
            actual = create_proteome_workflow.get_resource("create_proteome", "run", resource)()
            assert actual == expected, msg_error


# Tests for CreateProteomeWorkflow -------------------------------------------------------


def test_create_proteome_workflow(create_proteome_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["create_proteome", "link_out"]
    actual = list(sorted(create_proteome_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = "output/bwa.gatk4_hc.filtered.{sample}-DNA1-WGS1/{dir_}/bwa.gatk4_hc.filtered.{sample}-DNA1-WGS1.{ext}"

    expected = [
        tpl.format(sample=sample, ext=ext, dir_="out")
        for sample in ("P001-N1", "P001-T1", "P002-N1", "P002-T1", "P002-T2")
        for ext in ("fa.gz", "fa.gz.md5")
    ]
    expected += [
        tpl.format(sample=sample, ext=ext, dir_="log")
        for sample in ("P001-N1", "P001-T1", "P002-N1", "P002-T1", "P002-T2")
        for ext in (
            "conda_info.txt",
            "conda_list.txt",
            "log",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
            "log.md5",
        )
    ]
    expected = list(sorted(expected))
    actual = list(sorted(create_proteome_workflow.get_result_files()))
    assert expected == actual
