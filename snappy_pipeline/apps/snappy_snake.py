# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline in a directory

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse
import datetime
import functools
import logging
import os
import subprocess
import sys
import textwrap

import ruamel.yaml as yaml
from snakemake import main as snakemake_main

from .. import __version__
from ..workflows import (
    cbioportal_export,
    gene_expression_quantification,
    gene_expression_report,
    hla_typing,
    igv_session_generation,
    ngs_data_qc,
    ngs_mapping,
    ngs_sanity_checking,
    panel_of_normals,
    roh_calling,
    somatic_gene_fusion_calling,
    somatic_hla_loh_calling,
    somatic_msi_calling,
    somatic_neoepitope_prediction,
    somatic_ngs_sanity_checking,
    somatic_purity_ploidy_estimate,
    somatic_targeted_seq_cnv_calling,
    somatic_variant_annotation,
    somatic_variant_calling,
    somatic_variant_checking,
    somatic_variant_expression,
    somatic_variant_filtration,
    somatic_variant_signatures,
    somatic_wgs_cnv_calling,
    somatic_wgs_sv_calling,
    targeted_seq_cnv_annotation,
    targeted_seq_cnv_calling,
    targeted_seq_cnv_export,
    tcell_crg_report,
    variant_annotation,
    variant_calling,
    variant_checking,
    variant_combination,
    variant_denovo_filtration,
    variant_export,
    variant_filtration,
    variant_phasing,
    wgs_cnv_annotation,
    wgs_cnv_calling,
    wgs_cnv_export,
    wgs_cnv_filtration,
    wgs_mei_annotation,
    wgs_mei_calling,
    wgs_mei_filtration,
    wgs_sv_annotation,
    wgs_sv_calling,
    wgs_sv_export,
    wgs_sv_filtration,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class DrmaaNotAvailable(Exception):
    """Raised when DRMAA is not available"""


#: Configuration file names
CONFIG_FILES = ("config.yaml", "config.json")

#: Shell to use
SHELL = "/bin/bash"

#: Mapping from step name to module
STEP_TO_MODULE = {
    "gene_expression_quantification": gene_expression_quantification,
    "gene_expression_report": gene_expression_report,
    "cbioportal_export": cbioportal_export,
    "hla_typing": hla_typing,
    "igv_session_generation": igv_session_generation,
    "ngs_mapping": ngs_mapping,
    "ngs_data_qc": ngs_data_qc,
    "panel_of_normals": panel_of_normals,
    "ngs_sanity_checking": ngs_sanity_checking,
    "roh_calling": roh_calling,
    "somatic_gene_fusion_calling": somatic_gene_fusion_calling,
    "somatic_hla_loh_calling": somatic_hla_loh_calling,
    "somatic_msi_calling": somatic_msi_calling,
    "somatic_neoepitope_prediction": somatic_neoepitope_prediction,
    "somatic_ngs_sanity_checking": somatic_ngs_sanity_checking,
    "somatic_purity_ploidy_estimate": somatic_purity_ploidy_estimate,
    "somatic_targeted_seq_cnv_calling": somatic_targeted_seq_cnv_calling,
    "somatic_variant_annotation": somatic_variant_annotation,
    "somatic_variant_calling": somatic_variant_calling,
    "somatic_variant_checking": somatic_variant_checking,
    "somatic_variant_expression": somatic_variant_expression,
    "somatic_variant_filtration": somatic_variant_filtration,
    "somatic_variant_signatures": somatic_variant_signatures,
    "somatic_wgs_cnv_calling": somatic_wgs_cnv_calling,
    "somatic_wgs_sv_calling": somatic_wgs_sv_calling,
    "targeted_seq_cnv_annotation": targeted_seq_cnv_annotation,
    "targeted_seq_cnv_calling": targeted_seq_cnv_calling,
    "targeted_seq_cnv_export": targeted_seq_cnv_export,
    "tcell_crg_report": tcell_crg_report,
    "variant_annotation": variant_annotation,
    "variant_calling": variant_calling,
    "variant_checking": variant_checking,
    "variant_combination": variant_combination,
    "variant_denovo_filtration": variant_denovo_filtration,
    "variant_export": variant_export,
    "variant_filtration": variant_filtration,
    "variant_phasing": variant_phasing,
    "wgs_cnv_annotation": wgs_cnv_annotation,
    "wgs_cnv_calling": wgs_cnv_calling,
    "wgs_cnv_export": wgs_cnv_export,
    "wgs_cnv_filtration": wgs_cnv_filtration,
    "wgs_mei_annotation": wgs_mei_annotation,
    "wgs_mei_calling": wgs_mei_calling,
    "wgs_mei_filtration": wgs_mei_filtration,
    "wgs_sv_annotation": wgs_sv_annotation,
    "wgs_sv_calling": wgs_sv_calling,
    "wgs_sv_export": wgs_sv_export,
    "wgs_sv_filtration": wgs_sv_filtration,
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


def run(wrapper_args):
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
    if wrapper_args.unlock:
        snakemake_argv.append("--unlock")
    if wrapper_args.rerun_incomplete:
        snakemake_argv.append("--rerun-incomplete")
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
    if not wrapper_args.snappy_pipeline_use_drmaa:
        snakemake_argv += ["--cores", str(wrapper_args.cores or 1)]
    if wrapper_args.use_conda:
        snakemake_argv.append("--use-conda")
        if mamba_available and wrapper_args.use_mamba:
            snakemake_argv += ["--conda-frontend", "mamba"]

    # Configure DRMAA if configured so
    if wrapper_args.snappy_pipeline_use_drmaa:
        if wrapper_args.restart_times:
            snakemake_argv += ["--restart-times", str(wrapper_args.restart_times)]
        if wrapper_args.max_jobs_per_second:
            snakemake_argv += ["--max-jobs-per-second", str(wrapper_args.max_jobs_per_second)]
        if wrapper_args.max_status_checks_per_second:
            snakemake_argv += [
                "--max-status-checks-per-second",
                str(wrapper_args.max_status_checks_per_second),
            ]
        if not drmaa_available():
            raise DrmaaNotAvailable()
        # create output log directory, use SLURM_JOB_ID env variable for sub directory name if set
        logdir = os.path.abspath(os.path.join(wrapper_args.directory, "slurm_log"))
        if "SLURM_JOB_ID" in os.environ:
            logdir = os.path.join(logdir, os.environ["SLURM_JOB_ID"])
        else:
            logdir = os.path.join(logdir, datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
        os.makedirs(logdir, exist_ok=True)
        logging.info("Creating directory %s", logdir)
        # extend arguments with DRMAA settings
        snakemake_argv += [
            # arguments to DRMAA, ties in cluster configuration
            "--drmaa",
            (
                " --mem={{cluster.mem}} --time={{cluster.time}} "
                "--cpus-per-task={{cluster.ntasks}} --output={logdir}/slurm-%x-%J.log {drmaa_snippet}"
            ).format(logdir=logdir, drmaa_snippet=" ".join(wrapper_args.drmaa_snippets)),
            # number of parallel jobs to launch
            "-j",
            str(wrapper_args.snappy_pipeline_drmaa_jobs),
        ]
    # Add positional arguments
    if wrapper_args.cleanup_metadata:
        snakemake_argv.append("--cleanup-metadata")
    snakemake_argv += wrapper_args.targets
    logging.info("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    return snakemake_main(snakemake_argv)


def drmaa_available():
    """Return whether DRMAA execution is available"""
    try:
        pass
    except (ImportError, RuntimeError):
        return False
    else:
        return True


def run_self_test():
    """Print snappy_pipeline version and enabled/disabled features"""
    yes_no = {True: "yes", False: "no"}
    vals = {"version": __version__, "drmaa_enabled": yes_no[drmaa_available()]}
    print(
        textwrap.dedent(
            r"""
    CUBI Pipeline
    =============

    Version     {version}

    Features
    --------

    DRMAA       {drmaa_enabled}
    """
        )
        .lstrip()
        .format(**vals),
        file=sys.stderr,
    )


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

    group = parser.add_argument_group(
        "Snakemake Verbosity / Debugging",
        "Arguments from Snakemake that are useful for debugging, such as " "increasing verbosity",
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
        "DRMAA/Cluster Configuration",
        "Arguments for enabling and controllingn basic DRMAA execution",
    )
    group.add_argument(
        "--snappy-pipeline-use-drmaa",
        action="store_true",
        default=False,
        help="Enables running the pipeline with DRMA via Snakemake",
    )
    group.add_argument(
        "--snappy-pipeline-drmaa-jobs", type=int, default=100, help="Number of DRMAA jobs to run"
    )
    group.add_argument(
        "--drmaa-snippet",
        dest="drmaa_snippets",
        default=[],
        action="append",
        nargs="+",
        help="Snippet to include to DRMAA submission",
    )

    group = parser.add_argument_group(
        "DRMAA Robustness",
        "Snakemake settings to increase robustness of executing the pipeline "
        "during execution in DRMAA mode",
    )
    group.add_argument(
        "--restart-times", type=int, default=5, help="Number of times to restart jobs automatically"
    )
    group.add_argument(
        "--max-jobs-per-second",
        type=int,
        default=10,
        help="Maximal number of jobs to launch per second in DRMAA mode",
    )
    group.add_argument(
        "--max-status-checks-per-second",
        type=int,
        default=10,
        help="Maximal number of status checks to perform per second",
    )

    group = parser.add_argument_group(
        "CUBI Pipeline Miscellaneous", "Various options specific to the CUBI pipeline"
    )
    group.add_argument(
        "--snappy-pipeline-self-test",
        action="store_true",
        default=False,
        help="Perform self-test and exit",
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

    args = parser.parse_args(argv)

    args.drmaa_snippets = functools.reduce(lambda x, y: x + y, args.drmaa_snippets, [])

    if args.snappy_pipeline_self_test:
        run_self_test()
    else:
        if not args.step:
            for cfg in CONFIG_FILES:
                path = os.path.join(args.directory, cfg)
                if not os.path.exists(path):
                    continue
                with open(path, "rt") as f:
                    data = yaml.round_trip_load(f.read())
                try:
                    args.step = data["pipeline_step"]["name"]
                    break
                except KeyError:
                    logging.info("Could not pick up pipeline step/name from %s", path)
        if not args.step:
            parser.error("the following arguments are required: --step")
        try:
            return run(args)
        except DrmaaNotAvailable:
            parser.error("DRMAA not available but given --snappy-pipeline-use-drmaa")
            raise e


if __name__ == "__main__":
    sys.exit(main())
