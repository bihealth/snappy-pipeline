# -*- coding: utf-8 -*-
"""CUBI Pipeline targeted_seq_cnv_calling step Snakefile"""

import os

from snappy_pipeline import expand_ref
from snappy_pipeline.workflows.targeted_seq_cnv_calling import TargetedSeqCnvCallingWorkflow



# Configuration ===============================================================


configfile: "config.yaml"


# Expand "$ref" JSON pointers in configuration (also works for YAML)
config, lookup_paths, config_paths = expand_ref("config.yaml", config)

# WorkflowImpl Object Setup ===================================================

wf = TargetedSeqCnvCallingWorkflow(
    workflow, config, lookup_paths, config_paths, os.getcwd()
)

# Rules =======================================================================

rule all:
    input:
        wf.substep_getattr("gcnv", "get_cnv_model_result_files"),


rule targeted_seq_cnv_calling_gcnv_preprocess_intervals:
    input:
        unpack(wf.get_input_files("gcnv", "preprocess_intervals")),
    output:
        **wf. get_output_files("gcnv","preprocess_intervals"),

    threads: wf.get_resources("gcnv", "preprocess_intervals", "threads")
    resources:
        time=wf.get_resource("gcnv", "preprocess_intervals", "time"),
        mem=wf.get_resource("gcnv", "preprocess_intervals", "memory"),
        partition=wf.get_resources("gcnv", "preprocess_intervals", "partition"),
    log:
        wf.get_log_file("gcnv", "preprocess_intervals"),
    wrapper:
        wf.wrapper_path("gcnv/preprocess_intervals")


rule targeted_seq_cnv_calling_gcnv_annotate_gc:
    input:
        unpack(wf.get_input_files("gcnv", "annotate_gc")),
    output:
        **wf. get_output_files("gcnv","annotate_gc"),

    threads: wf.get_resources("gcnv", "annotate_gc", "threads")
    resources:
        time=wf.get_resource("gcnv", "annotate_gc", "time"),
        mem=wf.get_resource("gcnv", "annotate_gc", "memory"),
        partition=wf.get_resources("gcnv", "annotate_gc", "partition"),
    log:
        wf.get_log_file("gcnv", "annotate_gc"),
    wrapper:
        wf.wrapper_path("gcnv/annotate_gc")


rule targeted_seq_cnv_calling_gcnv_coverage:
    input:
        unpack(wf.get_input_files("gcnv", "coverage")),
    output:
        **wf. get_output_files("gcnv","coverage"),
    threads: wf.get_resources("gcnv", "coverage", "threads")
    resources:
        time=wf.get_resource("gcnv", "coverage", "time"),
        mem=wf.get_resource("gcnv", "coverage", "memory"),
        partition=wf.get_resources("gcnv", "coverage", "partition"),
    log:
        wf.get_log_file("gcnv", "coverage"),
    wrapper:
        wf.wrapper_path("gcnv/coverage")


rule targeted_seq_cnv_calling_gcnv_filter_intervals:
    input:
        unpack(wf.get_input_files("gcnv", "filter_intervals")),
    output:
        **wf. get_output_files("gcnv","filter_intervals"),
    threads: wf.get_resources("gcnv", "filter_intervals", "threads")
    resources:
        time=wf.get_resource("gcnv", "filter_intervals", "time"),
        mem=wf.get_resource("gcnv", "filter_intervals", "memory"),
        partition=wf.get_resources("gcnv", "filter_intervals", "partition"),
    log:
        wf.get_log_file("gcnv", "filter_intervals"),
    wrapper:
        wf.wrapper_path("gcnv/filter_intervals")


rule targeted_seq_cnv_calling_gcnv_contig_ploidy:
    input:
        unpack(wf.get_input_files("gcnv", "contig_ploidy")),
    output:
        **wf. get_output_files("gcnv","contig_ploidy"),
    threads: wf.get_resources("gcnv", "contig_ploidy", "threads")
    resources:
        time=wf.get_resource("gcnv", "contig_ploidy", "time"),
        mem=wf.get_resource("gcnv", "contig_ploidy", "memory"),
        partition=wf.get_resources("gcnv", "contig_ploidy", "partition"),
    log:
        wf.get_log_file("gcnv", "contig_ploidy"),
    wrapper:
        wf.wrapper_path("gcnv/contig_ploidy")
