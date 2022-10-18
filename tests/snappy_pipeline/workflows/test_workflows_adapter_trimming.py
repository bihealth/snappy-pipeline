# -*- coding: utf-8 -*-
"""Tests for the hla_typing workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.adapter_trimming import AdapterTrimmingWorkflow

from .conftest import patch_module_fs

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


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

        step_config:
           adapter_trimming:
             tools: ["bbduk", "fastp"]
             bbduk:
               adapter_sequences: /path/to/adapter_sequences.fa
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
def adapter_trimming_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return AdapterTrimmingWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Construct the workflow object
    return AdapterTrimmingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for LinkOutFastqStepPart ---------------------------------------------------------------


def test_link_out_fastq_step_part_get_output_files(adapter_trimming_workflow):
    """Tests LinkOutFastqStepPart.get_output_files()"""
    expected = [
        "output/{trimmer}/{library_name}/log/.done",
        "output/{trimmer}/{library_name}/report/.done",
        "output/{trimmer}/{library_name}/out/.done",
    ]
    actual = adapter_trimming_workflow.get_output_files("link_out_fastq", "run")
    assert actual == expected


def test_link_out_fastq_step_part_get_shell_cmd(adapter_trimming_workflow):
    """Tests LinkOutFastqStepPart.get_output_files()"""
    wildcards = Wildcards(fromdict={"trimmer": "bbduk", "library_name": "P001-N1-DNA1-WGS1"})
    expected = textwrap.dedent(
        r"""
        din_=$(dirname work/bbduk.P001-N1-DNA1-WGS1/log/.done) ; dout=$(dirname output/bbduk/P001-N1-DNA1-WGS1/log/.done) ; fns=$(find $din_ -type f -printf '%P\n') ; for fn in $fns ; do     if [[ ! -L $din_/$fn ]] ; then       mkdir -p $(dirname $dout/$fn) ; ln -sr $din_/$fn $dout/$fn   ; fi ; done
        din_=$(dirname work/bbduk.P001-N1-DNA1-WGS1/report/.done) ; dout=$(dirname output/bbduk/P001-N1-DNA1-WGS1/report/.done) ; fns=$(find $din_ -type f -printf '%P\n') ; for fn in $fns ; do     if [[ ! -L $din_/$fn ]] ; then       mkdir -p $(dirname $dout/$fn) ; ln -sr $din_/$fn $dout/$fn   ; fi ; done
        din_=$(dirname work/bbduk.P001-N1-DNA1-WGS1/out/.done) ; dout=$(dirname output/bbduk/P001-N1-DNA1-WGS1/out/.done) ; fns=$(find $din_ -type f -printf '%P\n') ; for fn in $fns ; do     if [[ ! -L $din_/$fn ]] ; then       mkdir -p $(dirname $dout/$fn) ; ln -sr $din_/$fn $dout/$fn   ; fi ; done
        """
    ).strip()
    actual = adapter_trimming_workflow.get_shell_cmd("link_out_fastq", "run", wildcards)
    assert actual == expected


# Tests for BbdukStepPart ----------------------------------------------------------------------


def test_bbduk_step_part_get_input_files(adapter_trimming_workflow):
    """Tests BBdukStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = adapter_trimming_workflow.get_input_files("bbduk", "run")
    assert actual == expected


def test_bbduk_step_part_get_output_files(adapter_trimming_workflow):
    """Tests BbdukStepPart.get_output_files()"""
    expected = {
        "out_done": "work/bbduk.{library_name}/out/.done",
        "report_done": "work/bbduk.{library_name}/report/.done",
        "rejected_done": "work/bbduk.{library_name}/rejected/.done",
    }
    actual = adapter_trimming_workflow.get_output_files("bbduk", "run")
    assert actual == expected


def test_bbduk_step_part_get_args_input(adapter_trimming_workflow):
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1"})
    actual = adapter_trimming_workflow.substep_dispatch("bbduk", "get_args", "run")(wildcards)
    expected = {
        "library_name": "P001-T1-DNA1-WGS1",
        "input": {
            "reads_left": {
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001-T1-DNA1-WGS1_R1.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001-T1-DNA1-WGS1_R1.fastq.gz",
                },
            },
            "reads_right": {
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001-T1-DNA1-WGS1_R2.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001-T1-DNA1-WGS1_R2.fastq.gz",
                },
            },
        },
    }
    assert actual == expected


def test_bbduk_step_part_get_log_file(adapter_trimming_workflow):
    """Tests BbdukStepPart.get_log_file()"""
    expected = {
        "done": "work/bbduk.{library_name}/log/.done",
        "done_md5": "work/bbduk.{library_name}/log/.done.md5",
    }
    key_ext = (
        ("log", ".log"),
        ("conda_info", ".conda_info.txt"),
        ("conda_list", ".conda_list.txt"),
    )
    for key, ext in key_ext:
        expected[key] = "work/bbduk.{library_name}/log/bbduk.{library_name}" + ext
        expected[key + "_md5"] = "work/bbduk.{library_name}/log/bbduk.{library_name}" + ext + ".md5"
    actual = adapter_trimming_workflow.get_log_file("bbduk", "run")
    assert actual == expected


def test_bbduk_step_part_get_resource_usage(adapter_trimming_workflow):
    """Tests BbdukStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 8, "time": "12:00:00", "memory": "24000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = adapter_trimming_workflow.get_resource("bbduk", "run", resource)
        assert actual == expected, msg_error


# Tests for FastpStepPart ----------------------------------------------------------------------


def test_fastp_step_part_get_input_files(adapter_trimming_workflow):
    """Tests FastpStepPart.get_input_files()"""
    expected = {"done": "work/input_links/{library_name}/.done"}
    actual = adapter_trimming_workflow.get_input_files("fastp", "run")
    assert actual == expected


def test_fastp_step_part_get_output_files(adapter_trimming_workflow):
    """Tests FastpStepPart.get_output_files()"""
    expected = {
        "out_done": "work/fastp.{library_name}/out/.done",
        "report_done": "work/fastp.{library_name}/report/.done",
        "rejected_done": "work/fastp.{library_name}/rejected/.done",
    }
    actual = adapter_trimming_workflow.get_output_files("fastp", "run")
    assert actual == expected


def test_fastp_step_part_get_args_input(adapter_trimming_workflow):
    wildcards = Wildcards(fromdict={"library_name": "P001-T1-DNA1-WGS1"})
    actual = adapter_trimming_workflow.substep_dispatch("fastp", "get_args", "run")(wildcards)
    expected = {
        "library_name": "P001-T1-DNA1-WGS1",
        "input": {
            "reads_left": {
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001-T1-DNA1-WGS1_R1.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001-T1-DNA1-WGS1_R1.fastq.gz",
                },
            },
            "reads_right": {
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001-T1-DNA1-WGS1_R2.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001-T1-DNA1-WGS1_R2.fastq.gz",
                },
            },
        },
    }
    assert actual == expected


def test_fastp_step_part_get_log_file(adapter_trimming_workflow):
    """Tests FastpStepPart.get_log_file()"""
    expected = {
        "done": "work/fastp.{library_name}/log/.done",
        "done_md5": "work/fastp.{library_name}/log/.done.md5",
    }
    key_ext = (
        ("log", ".log"),
        ("conda_info", ".conda_info.txt"),
        ("conda_list", ".conda_list.txt"),
    )
    for key, ext in key_ext:
        expected[key] = "work/fastp.{library_name}/log/fastp.{library_name}" + ext
        expected[key + "_md5"] = "work/fastp.{library_name}/log/fastp.{library_name}" + ext + ".md5"
    actual = adapter_trimming_workflow.get_log_file("fastp", "run")
    assert actual == expected


def test_fastp_step_part_get_resource_usage(adapter_trimming_workflow):
    """Tests FastpStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 4, "time": "12:00:00", "memory": "24000M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = adapter_trimming_workflow.get_resource("fastp", "run", resource)
        assert actual == expected, msg_error


# Tests for AdapterTrimmingWorkflow ------------------------------------------------------------


def test_adapter_trimming_workflow_get_results(adapter_trimming_workflow):
    """Tests the ouput of AdapterTrimmingWorkflow.get_result_files()"""
    libraries = (
        "P001-N1-DNA1-WGS1",
        "P001-T1-DNA1-WGS1",
        "P001-T1-RNA1-mRNA_seq1",
        "P002-N1-DNA1-WGS1",
        "P002-T1-DNA1-WGS1",
        "P002-T1-DNA1-WGS2",
        "P002-T2-DNA1-WGS1",
        "P002-T2-RNA1-mRNA_seq1",
    )
    tpl = "output/{tool}/{library}/{sub_dir}/.done"
    expected = []
    for library in libraries:
        for tool in ["bbduk", "fastp"]:
            for sub_dir in ["out", "report", "log"]:
                expected.append(tpl.format(tool=tool, library=library, sub_dir=sub_dir))
    actual = adapter_trimming_workflow.get_result_files()
    assert actual == expected
