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
               adapter_sequences:
                 - /path/to/adapter_sequences.fa
             fastp:
               num_threads: 4
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
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001_T1_DNA1_WGS1_R1.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001_T1_DNA1_WGS1_R1.fastq.gz",
                },
            },
            "reads_right": {
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001_T1_DNA1_WGS1_R2.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001_T1_DNA1_WGS1_R2.fastq.gz",
                },
            },
        },
        "config": {
            "adapter_sequences": ["/path/to/adapter_sequences.fa"],
            "num_threads": 8,
            "interleaved": "auto",
            "qin": "auto",
            "copyundefined": False,
            "nzo": True,
            "qout": "auto",
            "statscolumns": 3,
            "rename": False,
            "refnames": False,
            "trd": False,
            "ordered": False,
            "gcbins": "auto",
            "maxhistlen": 6000,
            "histbefore": True,
            "idbins": 100,
            "k": 21,
            "rcomp": True,
            "maskmiddle": True,
            "minkmerhits": 1,
            "minkmerfraction": 0.0,
            "mincovfraction": 0.0,
            "hammingdistance": 1,
            "qhdist": 0,
            "editdistance": 0,
            "hammingdistance2": 0,
            "qhdist2": 0,
            "editdistance2": 0,
            "forbidn": False,
            "removeifeitherbad": True,
            "trimfailures": False,
            "findbestmatch": False,
            "skipr1": False,
            "skipr2": False,
            "ecco": False,
            "ktrim": "r",
            "kmask": "",
            "maskfullycovered": False,
            "ksplit": False,
            "mink": 11,
            "qtrim": "rl",
            "trimq": 25,
            "minlength": 35,
            "mlf": 0,
            "minavgquality": 0,
            "maqb": 0,
            "minbasequality": 0,
            "maxns": -1,
            "mcb": 0,
            "ottm": False,
            "tp": 0,
            "tbo": False,
            "strictoverlap": True,
            "minoverlap": 14,
            "mininsert": 40,
            "tpe": False,
            "forcetrimleft": 0,
            "forcetrimright": 0,
            "forcetrimright2": 0,
            "forcetrimmod": 5,
            "restrictleft": 0,
            "restrictright": 0,
            "mingc": 0.0,
            "maxgc": 1.0,
            "gcpairs": True,
            "tossjunk": False,
            "swift": False,
            "chastityfilter": False,
            "barcodefilter": "f",
            "barcodes": "",
            "xmin": -1,
            "ymin": -1,
            "xmax": -1,
            "ymax": -1,
            "trimpolya": 0,
            "trimpolygleft": 0,
            "trimpolygright": 8,
            "trimpolyg": 0,
            "filterpolyg": 8,
            "entropy": -1.0,
            "entropywindow": 50,
            "entropyk": 5,
            "minbasefrequency": 0.0,
            "entropytrim": "f",
            "entropymask": "f",
            "entropymark": False,
            "cardinality": False,
            "cardinalityout": False,
            "loglogk": 31,
            "loglogbuckets": 2048,
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
        actual = adapter_trimming_workflow.get_resource("bbduk", "run", resource)()
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
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001_T1_DNA1_WGS1_R1.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001_T1_DNA1_WGS1_R1.fastq.gz",
                },
            },
            "reads_right": {
                "work/input_links/P001-T1-DNA1-WGS1/FCXXXXXX/L001/P001_T1_DNA1_WGS1_R2.fastq.gz": {
                    "relative_path": "FCXXXXXX/L001",
                    "filename": "P001_T1_DNA1_WGS1_R2.fastq.gz",
                },
            },
        },
        "config": {
            "num_threads": 4,
            "trim_front1": 0,
            "trim_tail1": 0,
            "max_len1": 0,
            "trim_front2": 0,
            "trim_tail2": 0,
            "max_len2": 0,
            "dedup": False,
            "dup_calc_accuracy": 0,
            "dont_eval_duplication": True,
            "trim_poly_g": True,
            "poly_g_min_len": 8,
            "trim_poly_x": False,
            "poly_x_min_len": 10,
            "cut_front": False,
            "cut_tail": False,
            "cut_right": False,
            "cut_front_window_size": 4,
            "cut_front_mean_quality": 20,
            "cut_tail_window_size": 4,
            "cut_tail_mean_quality": 20,
            "cut_right_window_size": 4,
            "cut_right_mean_quality": 20,
            "disable_quality_filtering": False,
            "qualified_quality_phred": 15,
            "unqualified_percent_limit": 40,
            "n_base_limit": 5,
            "average_qual": 0,
            "disable_length_filtering": False,
            "length_required": 15,
            "length_limit": 0,
            "low_complexity_filter": False,
            "complexity_threshold": 30,
            "filter_by_index1": "",
            "filter_by_index2": "",
            "filter_by_index_threshold": 0,
            "correction": False,
            "overlap_len_require": 30,
            "overlap_diff_limit": 5,
            "overlap_diff_percent_limit": 20,
            "umi": False,
            "umi_loc": "",
            "umi_len": 0,
            "umi_prefix": "",
            "umi_skip": 0,
            "overrepresentation_analysis": False,
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
        actual = adapter_trimming_workflow.get_resource("fastp", "run", resource)()
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
