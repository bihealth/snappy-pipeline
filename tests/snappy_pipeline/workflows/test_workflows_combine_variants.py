# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.combine_variants import CombineVariantsWorkflow

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
          cosmic:
            path: /path/to/cosmic.vcf.gz
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa

          somatic_variant_calling:
            tools:
            - mutect2
            mutect2:
              contamination: {}

          somatic_variant_annotation:
            path_somatic_variant: /path/to/somatic_variant_filtration
            tools: ["vep"]
            vep:
              cache_dir: /path/to/dir/cache
              plugins:
                - {'name': 'local_plugin', 'path': '/path/to/local_plugin.pm'}
                - {'name': 'remote_plugin', 'url': 'https://to.com/remote_plugin.pm'}
          
          germline_variant_calling:
            tools: [gatk4_hc]
            gatk4_hc: {}

          germline_variant_filtration:
            path_variant: ../germline_variant_calling
            has_annotation: false
            filter_list:
            - bcftools:
                include: 'depth > min'
        
          combine_variants:
            somatic_variant_type: annotation
            path_somatic_variant: SOMATIC_VARIANT_ANNOTATION
            tool_somatic_variant_annotation: vep
            germline_variant_type: filtration
            path_germline_variant: GERMLINE_VARIANT_FILTRATION
            is_germline_variant_filtered: true
            rename_combined: tumor

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
def combine_variants_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return CombineVariantsWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "somatic_variant": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
        "germline_variant": lambda x: "GERMLINE_VARIANT_FILTRATION/" + x,
    }
    # Construct the workflow object
    return CombineVariantsWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CombineVariantsStepPart ----------------------------------------------------------


def test_combine_step_part_get_input_files_run(combine_variants_workflow):
    """Tests CombineVariantsStepPart.get_input_files()"""
    wildcards: Wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    base_out = "{dir}/output/bwa.{caller}.{library}/out/bwa.{caller}.{library}.vcf.gz"
    expected = {
        "somatic_vcf": base_out.format(dir="SOMATIC_VARIANT_ANNOTATION", caller="mutect2.vep", library="P001-T1-DNA1-WGS1"),
        "germline_vcf": base_out.format(dir="GERMLINE_VARIANT_FILTRATION", caller="gatk4_hc.filtered", library="P001-N1-DNA1-WGS1"),
        "reference": "/path/to/ref.fa"
    }
    actual = combine_variants_workflow.get_input_files("combine", "run")(wildcards)
    assert actual == expected


def test_combine_step_part_get_output_files_run(combine_variants_workflow):
    """Tests CombineVariantsStepPart.get_output_files()"""
    base_out = "bwa.combined.{tumor_library}"
    expected = {"vcf": f"work/{base_out}/out/{base_out}.vcf.gz"}
    actual = combine_variants_workflow.get_output_files("combine", "run")
    assert actual == expected


def test_combine_step_part_get_log_file_run(combine_variants_workflow):
    """Tests CombineVariantsStepPart.get_output_files()"""
    base_out = "bwa.combined.{tumor_library}"
    expected = get_expected_log_files_dict(base_out=f"work/{base_out}/log/{base_out}")
    actual = combine_variants_workflow.get_log_file("combine", "run")
    assert actual == expected


def test_combine_step_part_get_args_run(combine_variants_workflow):
    """Tests CombineVariantsStepPart.get_args()"""
    wildcards: Wildcards = Wildcards(fromdict={"tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {"sample_name": "P001-T1-DNA1-WGS1"}
    actual = combine_variants_workflow.get_args("combine", "run")(wildcards)
    assert actual == expected


def test_combine_step_part_get_resource_usage(combine_variants_workflow):
    """Tests CombineVariantsStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {
        "run": {"threads": 1, "time": "03:59:59", "memory": "4G", "partition": "medium"}
    }
    # Evaluate
    for action in ("run",):
        for resource, expected in expected_dict[action].items():
            msg_error = f"Assertion error for resource '{resource}'."
            actual = combine_variants_workflow.get_resource("combine", "run", resource)()
            assert actual == expected, msg_error


# Tests for CombineVariantsWorkflow -------------------------------------------------------


def test_combine_variants_workflow(combine_variants_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["combine", "link_out"]
    actual = list(sorted(combine_variants_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = "output/bwa.combined.P00{i}-T{t}-DNA1-WGS1/{dir_}/bwa.combined.P00{i}-T{t}-DNA1-WGS1.{ext}"

    expected = [
        tpl.format(
            mapper=mapper, i=i, t=t, ext=ext, dir_="out"
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
    ]
    expected += [
        tpl.format(
            mapper=mapper, i=i, t=t, ext=ext, dir_="log"
        )
        for i, t in ((1, 1), (2, 1), (2, 2))
        for ext in (
            "conda_info.txt",
            "conda_list.txt",
            "log",
            "conda_info.txt.md5",
            "conda_list.txt.md5",
            "log.md5",
        )
        for mapper in ("bwa",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(combine_variants_workflow.get_result_files()))
    assert expected == actual
