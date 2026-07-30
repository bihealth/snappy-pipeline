# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline in a directory

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse
import logging
import os
import sys

from snakemake.cli import main as snakemake_main

from .. import __version__
from ..workflows import (
    adapter_trimming,
    cbioportal_export,
    combine_variants,
    create_proteome,
    gene_expression_quantification,
    gene_expression_report,
    helper_gcnv_model_targeted,
    helper_gcnv_model_wgs,
    hla_typing,
    homologous_recombination_deficiency,
    igv_session_generation,
    ngs_data_qc,
    ngs_mapping,
    panel_of_normals,
    repeat_expansion,
    somatic_cnv_checking,
    somatic_gene_fusion_calling,
    somatic_hla_loh_calling,
    somatic_msi_calling,
    somatic_neoepitope_prediction,
    somatic_purity_ploidy_estimate,
    somatic_targeted_seq_cnv_calling,
    somatic_variant_signatures,
    somatic_wgs_cnv_calling,
    sv_calling_targeted,
    sv_calling_wgs,
    tumor_mutational_burden,
    varfish_export,
    variant_annotation,
    variant_calling,
    variant_checking,
    variant_export_external,
    variant_filtration,
    variant_phasing,
    wgs_cnv_export_external,
    wgs_sv_export_external,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Configuration file names
CONFIG_FILES = ("config.yaml", "config.json")

#: Shell to use
SHELL = "/bin/bash"

#: Mapping from step name to module
STEP_TO_MODULE = {
    "adapter_trimming": adapter_trimming,
    "cbioportal_export": cbioportal_export,
    "combine_variants": combine_variants,
    "create_proteome": create_proteome,
    "gene_expression_quantification": gene_expression_quantification,
    "gene_expression_report": gene_expression_report,
    "helper_gcnv_model_targeted": helper_gcnv_model_targeted,
    "helper_gcnv_model_wgs": helper_gcnv_model_wgs,
    "hla_typing": hla_typing,
    "homologous_recombination_deficiency": homologous_recombination_deficiency,
    "igv_session_generation": igv_session_generation,
    "ngs_mapping": ngs_mapping,
    "ngs_data_qc": ngs_data_qc,
    "panel_of_normals": panel_of_normals,
    "repeat_analysis": repeat_expansion,
    "somatic_cnv_checking": somatic_cnv_checking,
    "somatic_gene_fusion_calling": somatic_gene_fusion_calling,
    "somatic_hla_loh_calling": somatic_hla_loh_calling,
    "somatic_msi_calling": somatic_msi_calling,
    "somatic_neoepitope_prediction": somatic_neoepitope_prediction,
    "somatic_purity_ploidy_estimate": somatic_purity_ploidy_estimate,
    "somatic_targeted_seq_cnv_calling": somatic_targeted_seq_cnv_calling,
    "somatic_variant_signatures": somatic_variant_signatures,
    "somatic_wgs_cnv_calling": somatic_wgs_cnv_calling,
    "sv_calling_targeted": sv_calling_targeted,
    "sv_calling_wgs": sv_calling_wgs,
    "tumor_mutational_burden": tumor_mutational_burden,
    "varfish_export": varfish_export,
    "variant_annotation": variant_annotation,
    "variant_calling": variant_calling,
    "variant_checking": variant_checking,
    "variant_export_external": variant_export_external,
    "variant_filtration": variant_filtration,
    "variant_phasing": variant_phasing,
    "wgs_cnv_export_external": wgs_cnv_export_external,
    "wgs_sv_export_external": wgs_sv_export_external,
}


def setup_logging(args):
    """Setup logger."""
    logging.basicConfig(
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
    )
    logger = logging.getLogger("")
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


def run(wrapper_args, snakemake_args):
    """Launch the CUBI Pipeline wrapper for the given arguments"""
    # Point to the master orchestrator Snakefile
    orchestrator_snakefile = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Snakefile"
    )

    snakemake_argv = [
        "--directory",
        wrapper_args.directory,
        "--snakefile",
        orchestrator_snakefile,
    ]

    config_args = ["--config"]
    if wrapper_args.task:
        config_args.append(f"task={wrapper_args.task}")
    if wrapper_args.all_tasks:
        config_args.append("all_tasks=True")
    if wrapper_args.verbose:
        config_args.append("dump_orchestrator=True")

    if len(config_args) > 1:
        snakemake_argv.extend(config_args)

    # Configure profile if snappy pipeline profile is requested
    if wrapper_args.profile_snappy_pipeline:
        profile_path = os.path.join(os.path.dirname(__file__), "tpls", "profile")
        snakemake_argv += ["--profile", profile_path]

    # Append all user-provided snakemake arguments directly
    snakemake_argv += snakemake_args

    logging.info("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    return snakemake_main(snakemake_argv)


def main(argv=None):
    """Main program entry point, starts parsing command line arguments"""
    if argv is None:
        argv = sys.argv[1:]

    # Split arguments at '--' to cleanly separate snappy and snakemake arguments
    try:
        separator_idx = argv.index("--")
        snappy_args_list = argv[:separator_idx]
        snakemake_args = argv[separator_idx + 1 :]
    except ValueError:
        snappy_args_list = argv
        snakemake_args = []

    parser = argparse.ArgumentParser(
        usage="%(prog)s [--version] [-v] [-d directory] [--profile-snappy-pipeline] [--task TASK] [--all-tasks] [--] [snakemake arguments]",
        allow_abbrev=False,
    )

    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verbosity level")
    parser.add_argument(
        "-d", "--directory", default=os.getcwd(), help="Path to directory to run in, default is cwd"
    )
    parser.add_argument(
        "--profile-snappy-pipeline",
        action="store_true",
        help="Uses the profile defined in the snappy pipeline",
    )
    parser.add_argument(
        "--task",
        type=str,
        metavar="TASK",
        default=None,
        help="The specific task name from config.yaml to run",
    )
    parser.add_argument(
        "--all-tasks",
        action="store_true",
        default=False,
        help=(
            "Target all tasks (not just leaf tasks). "
            "By default only leaf tasks (tasks not depended on by any other task) are targeted."
        ),
    )

    # Only parse the arguments meant for snappy
    wrapper_args = parser.parse_args(snappy_args_list)

    # Setup logging
    setup_logging(wrapper_args)

    if not wrapper_args.task:
        if wrapper_args.all_tasks:
            logging.info("No specific --task provided. Targeting all tasks (--all-tasks).")
        else:
            logging.info("No specific --task provided. Targeting leaf tasks only (default).")

    return run(wrapper_args, snakemake_args)


if __name__ == "__main__":
    sys.exit(main())
