# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.germline_variant_calling import GermlineVariantCallingWorkflow

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

          germline_variant_calling:
            tools: ["gatk4_hc"]
            gatk4_hc: {}

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: generic
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def germline_variant_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    generic_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return SomaticVariantAnnotationWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", generic_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return GermlineVariantCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for GATK4_hc_step_part ----------------------------------------------------------

def test_gatk4_hc_discover_step_part_get_input_files(germline_variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    expected = {
        "reference": "/path/to/ref.fa",
        "dbsnp": "/path/to/dbsnp.vcf.gz",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
    }
    actual = germline_variant_calling_workflow.get_input_files("gatk4_hc", "discover")(wildcards)
    assert actual == expected


def test_gatk4_hc_discover_step_part_get_output_files(germline_variant_calling_workflow):
    gvcf = "work/{mapper}.gatk4_hc.{library_name}/out/{mapper}.gatk4_hc.{library_name}.discover.g.vcf.gz"
    expected = {
        "gvcf": gvcf,
        "gvcf_md5": gvcf + ".md5",
        "gvcf_tbi": gvcf + ".tbi",
        "gvcf_tbi_md5": gvcf + ".tbi.md5",
        "output_links": [],
    }
    actual = germline_variant_calling_workflow.get_output_files("gatk4_hc", "discover")
    assert actual == expected


def test_gatk4_hc_genotype_step_part_get_input_files(germline_variant_calling_workflow):
    wildcards = Wildcards(fromdict={"mapper": "bwa", "library_name": "P001-T1-DNA1-WGS1"})
    gvcf = f"work/bwa.gatk4_hc.P001-T1-DNA1-WGS1/out/bwa.gatk4_hc.P001-T1-DNA1-WGS1.discover.g.vcf.gz"
    expected = {
        "reference": "/path/to/ref.fa",
        "gvcf": gvcf,
        "gvcf_md5": gvcf + ".md5",
        "gvcf_tbi": gvcf + ".tbi",
        "gvcf_tbi_md5": gvcf + ".tbi.md5",
    }
    actual = germline_variant_calling_workflow.get_input_files("gatk4_hc", "genotype")(wildcards)
    assert actual == expected


def test_gatk4_hc_genotype_step_part_get_output_files(germline_variant_calling_workflow):
    tpl = "work/{mapper}.gatk4_hc.{library_name}/out/{mapper}.gatk4_hc.{library_name}"
    out = "output/{mapper}.gatk4_hc.{library_name}/out/{mapper}.gatk4_hc.{library_name}"
    log = "output/{mapper}.gatk4_hc.{library_name}/log/{mapper}.gatk4_hc.{library_name}_genotype"
    expected = {
        "gvcf": tpl + ".g.vcf.gz",
        "gvcf_md5": tpl + ".g.vcf.gz.md5",
        "gvcf_tbi": tpl + ".g.vcf.gz.tbi",
        "gvcf_tbi_md5": tpl + ".g.vcf.gz.tbi.md5",
        "vcf": tpl + ".vcf.gz",
        "vcf_md5": tpl + ".vcf.gz.md5",
        "vcf_tbi": tpl + ".vcf.gz.tbi",
        "vcf_tbi_md5": tpl + ".vcf.gz.tbi.md5",
        "output_links": [
            out + ".g.vcf.gz",
            out + ".g.vcf.gz.md5",
            out + ".g.vcf.gz.tbi",
            out + ".g.vcf.gz.tbi.md5",
            out + ".vcf.gz",
            out + ".vcf.gz.md5",
            out + ".vcf.gz.tbi",
            out + ".vcf.gz.tbi.md5",
            log + ".log",
            log + ".log.md5",
            log + ".conda_info.txt",
            log + ".conda_info.txt.md5",
            log + ".conda_list.txt",
            log + ".conda_list.txt.md5",
            log + ".wrapper.py",
            log + ".wrapper.py.md5",
            log + ".environment.yaml",
            log + ".environment.yaml.md5",
        ],
    }
    actual = germline_variant_calling_workflow.get_output_files("gatk4_hc", "genotype")
    assert actual == expected


def test_gatk4_hc_step_part_get_log_file(germline_variant_calling_workflow):
    for action in ("discover", "genotype"):
       base_out = f"work/{{mapper}}.gatk4_hc.{{library_name}}/log/{{mapper}}.gatk4_hc.{{library_name}}_{action}"
       expected = get_expected_log_files_dict(base_out=base_out, extended=True)
       actual = germline_variant_calling_workflow.get_log_file("gatk4_hc", action)
       assert actual == expected


def test_gatk4_hc_step_part_get_args(germline_variant_calling_workflow):
    expected = {
        "window_length": 10000000,
        "num_threads": 16,
        "allow_seq_dict_incompatibility": False,
        "step_key": "variant_calling",
        "caller_key": "gatk4_hc",
        "ignore_chroms": ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"],
    }
    for action in ("discover", "genotype"):
       actual = germline_variant_calling_workflow.get_args("gatk4_hc", action)
       assert actual == expected


def test_gatk4_hc_step_part_get_resources(germline_variant_calling_workflow):
    expected = {
        "threads": 16,
        "memory": "88G",
        "time": "2-00:00:00",
    }
    for action in ("discover", "genotype"):
        for resource, value in expected.items():
            actual = germline_variant_calling_workflow.get_resource("gatk4_hc", action, resource)()
            assert actual == value


# Tests for GermlineCallingWorkflow -------------------------------------------------------


def test_somatic_variant_annotation_workflow(germline_variant_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["gatk4_hc"]
    actual = list(sorted(germline_variant_calling_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{caller}.E002-BS1-TS1-LIB{i}/{dir_}/"
        "{mapper}.{caller}.E002-BS1-TS1-LIB{i}{ext}"
    )
    expected = [
        tpl.format(
            mapper=mapper, caller=caller, i=i, ext=ext, dir_="out"
        )
        for i in (1, 2)
        for ext in (".vcf.gz", ".vcf.gz.md5", ".vcf.gz.tbi", ".vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for caller in ("gatk4_hc",)
   ]
    expected += [
        tpl.format(
            mapper=mapper, caller=caller, i=i, ext=ext, dir_="out"
        )
        for i in (1, 2)
        for ext in (".g.vcf.gz", ".g.vcf.gz.md5", ".g.vcf.gz.tbi", ".g.vcf.gz.tbi.md5")
        for mapper in ("bwa",)
        for caller in ("gatk4_hc",)
    ]
    expected += [
        tpl.format(
            mapper=mapper, caller=caller, i=i, ext="_genotype" + ext, dir_="log"
        )
        for i in (1, 2)
        for ext in (
            ".conda_info.txt",
            ".conda_list.txt",
            ".log",
            ".wrapper.py",
            ".environment.yaml",
            ".conda_info.txt.md5",
            ".conda_list.txt.md5",
            ".log.md5",
            ".wrapper.py.md5",
            ".environment.yaml.md5",
        )
        for mapper in ("bwa",)
        for caller in ("gatk4_hc",)
    ]
    expected = list(sorted(expected))
    actual = list(sorted(germline_variant_calling_workflow.get_result_files()))
    assert expected == actual
