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
    workflow, config, cluster_config, lookup_paths, config_paths, os.getcwd()
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
    log:
        wf.get_log_file("gcnv", "preprocess_intervals"),
    wrapper:
        wf.wrapper_path("gcnv/preprocess_intervals")


rule targeted_seq_cnv_calling_gcnv_annotate_gc:
    input:
        unpack(wf.get_input_files("gcnv", "annotate_gc")),
    output:
        **wf. get_output_files("gcnv","annotate_gc"),
    log:
        wf.get_log_file("gcnv", "annotate_gc"),
    wrapper:
        wf.wrapper_path("gcnv/annotate_gc")


rule targeted_seq_cnv_calling_gcnv_coverage:
    input:
        unpack(wf.get_input_files("gcnv", "coverage")),
    output:
        **wf. get_output_files("gcnv","coverage"),
    params:
        args=wf.get_params("gcnv", "coverage"),
    log:
        wf.get_log_file("gcnv", "coverage"),
    wrapper:
        wf.wrapper_path("gcnv/coverage")


rule targeted_seq_cnv_calling_gcnv_filter_intervals:
    input:
        unpack(wf.get_input_files("gcnv", "filter_intervals")),
    output:
        **wf. get_output_files("gcnv","filter_intervals"),
    log:
        wf.get_log_file("gcnv", "filter_intervals"),
    wrapper:
        wf.wrapper_path("gcnv/filter_intervals")


checkpoint targeted_seq_cnv_calling_gcnv_scatter_intervals:
    input:
        unpack(wf.get_input_files("gcnv", "scatter_intervals")),
    output:
        directory(wf.get_output_files("gcnv", "scatter_intervals")),
    log:
        wf.get_log_file("gcnv", "scatter_intervals"),
    wrapper:
        wf.wrapper_path("gcnv/scatter_intervals")


rule targeted_seq_cnv_calling_gcnv_contig_ploidy:
    # TODO: output should be created/linked per-sample
    input:
        unpack(wf.get_input_files("gcnv", "contig_ploidy")),
    output:
        **wf. get_output_files("gcnv","contig_ploidy"),
    log:
        wf.get_log_file("gcnv", "contig_ploidy"),
    wrapper:
        wf.wrapper_path("gcnv/contig_ploidy")
