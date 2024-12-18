# -*- coding: utf-8 -*-
"""Tests for the somatic_purity_ploidy_estimate workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.guess_sex import GuessSexWorkflow

from .conftest import patch_module_fs
from .common import get_expected_log_files_dict


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

          guess_sex:
            tools: ['samtools']

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
def guess_sex_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return GuessSexWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return GuessSexWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for SamtoolsStepPart --------------------------------------------------------------------------


def test_samtools_step_part_get_input_files_run(guess_sex_workflow):
    """Tests SamtoolsStepPart._get_input_files_run()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1", "mapper": "bwa"})
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
    }
    actual = guess_sex_workflow.get_input_files("samtools", "run")(wildcards)
    assert actual == expected


def test_samtools_step_part_get_output_files_run(guess_sex_workflow):
    """Tests SamtoolsStepPart._get_output_files_run()"""
    expected = {
        "table": "work/{mapper}.samtools.{library_name}/out/{mapper}.samtools.{library_name}.tsv",
        "decision": "work/{mapper}.samtools.{library_name}/out/{mapper}.samtools.{library_name}.txt",
        "table_md5": "work/{mapper}.samtools.{library_name}/out/{mapper}.samtools.{library_name}.tsv.md5",
        "decision_md5": "work/{mapper}.samtools.{library_name}/out/{mapper}.samtools.{library_name}.txt.md5",
    }
    actual = guess_sex_workflow.get_output_files("samtools", "run")
    assert actual == expected


def test_samtools_step_part_get_args_run(guess_sex_workflow):
    """Tests SamtoolsStepPart._get_args_run()"""
    expected = {
        "min_X_female": 1.5,
        "max_Y_female": 0.33,
        "max_X_male": 1.5,
        "min_Y_male": 0.5,
    }
    actual = guess_sex_workflow.get_args("samtools", "run")
    assert actual == expected


def test_samtools_step_part_get_log_file_run(guess_sex_workflow):
    """Tests SamtoolsStepPart._get_log_file_run()"""
    base_out = "work/{mapper}.samtools.{library_name}/log/{mapper}.samtools.{library_name}"
    expected = get_expected_log_files_dict(base_out=base_out)
    expected["script"] = base_out + ".sh"
    expected["script_md5"] = expected["script"] + ".md5"
    actual = guess_sex_workflow.get_log_file("samtools", "run")
    assert actual == expected


# Tests for GuessSexWorkflow ----------------------------------------------------


def test_guess_sex_workflow(guess_sex_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "samtools"]
    actual = list(sorted(guess_sex_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction for ascat
    log_exts = ("conda_info.txt", "conda_list.txt", "log", "sh")
    mapper = "bwa"
    expected = []
    for lib in ((1, "N1"), (1, "T1"), (2, "N1"), (2, "T1"), (2, "T2")):
        tpl = f"{mapper}.samtools.P00{lib[0]}-{lib[1]}-DNA1-WGS1"
        expected.append(f"output/{tpl}/out/{tpl}.txt")
        expected.append(f"output/{tpl}/out/{tpl}.txt.md5")
        for log_ext in log_exts:
            expected.append(f"output/{tpl}/log/{tpl}.{log_ext}")
            expected.append(f"output/{tpl}/log/{tpl}.{log_ext}.md5")

    expected = set(expected)
    actual = set(guess_sex_workflow.get_result_files())
    assert actual == expected
