# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline in a directory

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse
import datetime
import logging
import os
import subprocess
import sys

import ruamel.yaml as ruamel_yaml
from snakemake import RERUN_TRIGGERS
from snakemake import main as snakemake_main

from .. import __version__
from ..workflows import (
    adapter_trimming,
    cbioportal_export,
    gene_expression_quantification,
    gene_expression_report,
    germline_snvs,
    guess_sex,
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
    somatic_purity_ploidy_estimate,
    somatic_targeted_seq_cnv_calling,
    somatic_variant_annotation,
    somatic_variant_calling,
    somatic_variant_filtration,
    somatic_variant_signatures,
    somatic_variants_for_cnv,
    somatic_wgs_cnv_calling,
    somatic_wgs_sv_calling,
    sv_calling_targeted,
    sv_calling_wgs,
    tumor_mutational_burden,
    varfish_export,
    variant_annotation,
    variant_calling,
    variant_checking,
    variant_denovo_filtration,
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
    "gene_expression_quantification": gene_expression_quantification,
    "gene_expression_report": gene_expression_report,
    "germline_snvs": germline_snvs,
    "guess_sex": guess_sex,
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
    "somatic_purity_ploidy_estimate": somatic_purity_ploidy_estimate,
    "somatic_targeted_seq_cnv_calling": somatic_targeted_seq_cnv_calling,
    "somatic_variant_annotation": somatic_variant_annotation,
    "somatic_variant_calling": somatic_variant_calling,
    "somatic_variant_filtration": somatic_variant_filtration,
    "somatic_variant_signatures": somatic_variant_signatures,
    "somatic_variants_for_cnv": somatic_variants_for_cnv,
    "somatic_wgs_cnv_calling": somatic_wgs_cnv_calling,
    "somatic_wgs_sv_calling": somatic_wgs_sv_calling,
    "sv_calling_targeted": sv_calling_targeted,
    "sv_calling_wgs": sv_calling_wgs,
    "tumor_mutational_burden": tumor_mutational_burden,
    "varfish_export": varfish_export,
    "variant_annotation": variant_annotation,
    "variant_calling": variant_calling,
    "variant_checking": variant_checking,
    "variant_denovo_filtration": variant_denovo_filtration,
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


def binary_available(name):
    retcode = subprocess.call(["which", "mamba"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return retcode == 0


def run(wrapper_args):  # noqa: C901
    """Launch the CUBI Pipeline wrapper for the given arguments"""
    # Setup logging
    setup_logging(wrapper_args)

    mamba_available = binary_available("mamba")

    # Map from step to module
    # Build arguments for wrapped "snakemake" call and parse arguments
    module = STEP_TO_MODULE[wrapper_args.step]
    snakemake_argv = [
        "--directory",
        wrapper_args.directory,
        "--snakefile",
        os.path.join(os.path.dirname(os.path.abspath(module.__file__)), "Snakefile"),
        # Force using job script for now that overrides the TMPDIR from the cluster scheduler
        "--jobscript",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "tpls", "jobscript.sh"),
    ]
    if wrapper_args.keep_going:
        snakemake_argv.append("--keep-going")
    if wrapper_args.rerun_triggers:
        snakemake_argv += ["--rerun-triggers"] + wrapper_args.rerun_triggers
    if wrapper_args.unlock:
        snakemake_argv.append("--unlock")
    if wrapper_args.rerun_incomplete:
        snakemake_argv.append("--rerun-incomplete")
    if wrapper_args.ignore_incomplete:
        snakemake_argv.append("--ignore-incomplete")
    if wrapper_args.touch:
        snakemake_argv.append("--touch")
    if wrapper_args.detailed_summary:
        snakemake_argv.append("-D")
    if wrapper_args.summary:
        snakemake_argv.append("-S")
    if wrapper_args.printshellcmds:
        snakemake_argv.append("-p")
    if wrapper_args.print_compilation:
        snakemake_argv.append("--print-compilation")
    if wrapper_args.verbose:
        snakemake_argv.append("--verbose")
    if wrapper_args.debug:
        snakemake_argv.append("--debug")
    if wrapper_args.dryrun:
        snakemake_argv.append("--dryrun")
    if wrapper_args.reason:
        snakemake_argv.append("--reason")
    if wrapper_args.batch:
        snakemake_argv += ["--batch", wrapper_args.batch]
    if not wrapper_args.snappy_pipeline_use_profile:
        snakemake_argv += ["--cores", str(wrapper_args.cores or 1)]
    if wrapper_args.conda_create_envs_only:
        snakemake_argv.append("--conda-create-envs-only")
        wrapper_args.use_conda = True
    if wrapper_args.use_conda:
        snakemake_argv.append("--use-conda")
        if mamba_available and wrapper_args.use_mamba:
            snakemake_argv += ["--conda-frontend", "mamba"]
        else:
            snakemake_argv += ["--conda-frontend", "conda"]

    # Configure profile if configured
    if wrapper_args.snappy_pipeline_use_profile:
        snakemake_argv += ["--profile", wrapper_args.snappy_pipeline_use_profile]
        snakemake_argv += ["--jobs", str(wrapper_args.snappy_pipeline_jobs)]
        if wrapper_args.restart_times:
            snakemake_argv += ["--restart-times", str(wrapper_args.restart_times)]
        if wrapper_args.max_jobs_per_second:
            snakemake_argv += ["--max-jobs-per-second", str(wrapper_args.max_jobs_per_second)]
        if wrapper_args.max_status_checks_per_second:
            snakemake_argv += [
                "--max-status-checks-per-second",
                str(wrapper_args.max_status_checks_per_second),
            ]
        # create output log directory, use SLURM_JOB_ID env variable for sub directory name if set
        logdir = os.path.abspath(os.path.join(wrapper_args.directory, "slurm_log"))
        if "SLURM_JOB_ID" in os.environ:
            logdir = os.path.join(logdir, os.environ["SLURM_JOB_ID"])
        else:
            logdir = os.path.join(logdir, datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
        os.makedirs(logdir, exist_ok=True)
        logging.info("Creating directory %s", logdir)
    # Add positional arguments
    if wrapper_args.cleanup_metadata:
        snakemake_argv.append("--cleanup-metadata")
    snakemake_argv += wrapper_args.targets
    logging.info("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    return snakemake_main(snakemake_argv)


def main(argv=None):
    """Main program entry point, starts parsing command line arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)

    group = parser.add_argument_group("Snakemake Basics", "Basic arguments from Snakemake")
    group.add_argument("targets", metavar="FILE", nargs="*", help="The files to create")
    group.add_argument(
        "-k",
        "--keep-going",
        action="store_true",
        default=False,
        help="Go on with independent jobs in case of failure",
    )
    group.add_argument(
        "-S",
        "--summary",
        action="store_true",
        default=False,
        help="Print paths that are to be created",
    )
    group.add_argument(
        "-D",
        "--detailed-summary",
        action="store_true",
        default=False,
        help="Print detailed summary",
    )
    group.add_argument(
        "-d", "--directory", default=os.getcwd(), help="Path to directory to run in, default is cwd"
    )
    group.add_argument("--cores", type=int, help="Number of cores to use for local processing")
    group.add_argument(
        "--unlock", action="store_true", default=False, help="Unlock working directory"
    )
    group.add_argument(
        "--rerun-incomplete", action="store_true", default=False, help="Rerun incomplete jobs"
    )
    group.add_argument(
        "--ignore-incomplete", action="store_true", default=False, help="Ignore incomplete jobs"
    )
    group.add_argument(
        "--cleanup-metadata",
        action="store_true",
        default=False,
        help="Cleanup the metadata of given files. Mark files as complete.",
    )
    parser.add_argument(
        "-t",
        "--touch",
        default=False,
        action="store_true",
        help=(
            "Touch output files (mark them up to date without really "
            "changing them) instead of running their commands. This is "
            "used to pretend that the rules were executed, in order to "
            "fool future invocations of snakemake. Fails if a file does "
            "not yet exist."
        ),
    )
    rerun_triggers_default = ["mtime", "params", "input"]
    group.add_argument(
        "--rerun-triggers",
        nargs="+",
        choices=RERUN_TRIGGERS,
        default=rerun_triggers_default,
        help=f"Expose --rerun-triggers from snakemake and set to {rerun_triggers_default} by default",
    )
    group.add_argument(
        "--batch",
        metavar="RULE=BATCH/BATCHES",
        help="Create the given batch for the given rule.  See Snakemake documentation for more info.",
    )
    group = parser.add_argument_group(
        "Snakemake Verbosity / Debugging",
        "Arguments from Snakemake that are useful for debugging, such as increasing verbosity",
    )
    group.add_argument(
        "-p", "--printshellcmds", action="store_true", default=False, help="Print shell commands"
    )
    group.add_argument("--verbose", action="store_true", default=False, help="Enable verbose mode")
    group.add_argument("--debug", action="store_true", default=False, help="Enable debugging")
    group.add_argument(
        "-n",
        "--dryrun",
        action="store_true",
        default=False,
        help="Simulate/run workflow without executing anything",
    )
    group.add_argument(
        "--print-compilation", action="store_true", default=False, help="Print compilation and stop"
    )
    group.add_argument(
        "-r",
        "--reason",
        action="store_true",
        default=False,
        help="Print the reason for each executed rule",
    )

    group = parser.add_argument_group(
        "Cluster Configuration",
        "Arguments for enabling and controlling cluster execution",
    )
    group.add_argument(
        "--snappy-pipeline-use-profile",
        metavar="PROFILE",
        default=None,
        help="Enables running the pipeline in cluster mode using the given profile",
    )
    group.add_argument(
        "--snappy-pipeline-jobs",
        type=int,
        default=100,
        help="Number of cluster jobs to run in parallel",
    )
    group.add_argument(
        "--default-partition", type=str, default="medium", help="Default partition to use, if any"
    )

    group = parser.add_argument_group(
        "Cluster Robustness",
        "Snakemake settings to increase robustness of executing the pipeline "
        "during execution in cluster mode",
    )
    group.add_argument(
        "--restart-times", type=int, default=5, help="Number of times to restart jobs automatically"
    )
    group.add_argument(
        "--max-jobs-per-second",
        type=int,
        default=10,
        help="Maximal number of jobs to launch per second in cluster mode",
    )
    group.add_argument(
        "--max-status-checks-per-second",
        type=int,
        default=10,
        help="Maximal number of status checks to perform per second",
    )

    group = parser.add_argument_group(
        "SNAPPY Pipeline Miscellaneous", "Various options specific to the SNAPPY pipeline"
    )
    group.add_argument(
        "--step",
        type=str,
        metavar="STEP",
        choices=sorted(STEP_TO_MODULE.keys()),
        help="The type of the step to run",
    )
    group.add_argument(
        "--no-use-mamba",
        dest="use_mamba",
        action="store_false",
        default=True,
        help="Disable usage of mamba when available",
    )
    group.add_argument(
        "--no-use-conda",
        dest="use_conda",
        action="store_false",
        default=True,
        help="Disable usage of conda",
    )
    group.add_argument(
        "--conda-create-envs-only",
        dest="conda_create_envs_only",
        action="store_true",
        default=False,
        help="Prepare all conda environments",
    )

    args = parser.parse_args(argv)

    os.environ["SNAPPY_PIPELINE_PARTITION"] = args.default_partition
    if args.snappy_pipeline_use_profile:
        os.environ["SNAPPY_PIPELINE_SNAKEMAKE_PROFILE"] = args.snappy_pipeline_use_profile

    if not args.step:
        for cfg in CONFIG_FILES:
            path = os.path.join(args.directory, cfg)
            if not os.path.exists(path):
                continue
            with open(path, "rt") as f:
                yaml = ruamel_yaml.YAML()
                data = yaml.load(f.read())
            try:
                args.step = data["pipeline_step"]["name"]
                break
            except KeyError:
                logging.info("Could not pick up pipeline step/name from %s", path)
    if not args.step:
        parser.error("the following arguments are required: --step")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
