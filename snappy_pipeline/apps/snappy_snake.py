# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline in a directory

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse
import logging
import os
import sys

import ruamel.yaml as ruamel_yaml
from snakemake.cli import main as snakemake_main
from snakemake.settings.enums import RerunTrigger

from .. import __version__
from ..workflows import (
    adapter_trimming,
    cbioportal_export,
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
    somatic_purity_ploidy_estimate,
    somatic_targeted_seq_cnv_calling,
    somatic_variant_annotation,
    somatic_variant_calling,
    somatic_variant_filtration,
    somatic_variant_signatures,
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

# snakemake v8 now has an explicit enum for rerun triggers
RERUN_TRIGGERS = RerunTrigger.all()

#: Configuration file names
CONFIG_FILES = ("config.yaml", "config.json")

#: Shell to use
SHELL = "/bin/bash"

#: Mapping from step name to module
STEP_TO_MODULE = {
    "adapter_trimming": adapter_trimming,
    "gene_expression_quantification": gene_expression_quantification,
    "gene_expression_report": gene_expression_report,
    "cbioportal_export": cbioportal_export,
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


def run(wrapper_args, snakemake_args):  # noqa: C901
    """Launch the CUBI Pipeline wrapper for the given arguments"""
    # Build arguments for wrapped "snakemake" call and parse arguments

    # Map from step to module
    module = STEP_TO_MODULE[wrapper_args.step]

    # Fail for forbidden arguments (conflicts between snappy & snakemake)
    # Note that this code is brittle if snakemake allows abbreviations in its parser
    # (which is does, apparently)
    if "-s" in snakemake_args or "--snakefile" in snakemake_args:
        logging.error("User cannot specify snakefile, it is selected by the wrapper")
        return 1

    snakemake_argv = [
        "--directory",
        wrapper_args.directory,
        "--snakefile",
        os.path.join(os.path.dirname(os.path.abspath(module.__file__)), "Snakefile"),
    ]

    # Force conda usage
    if "--use-conda" not in snakemake_args:
        snakemake_argv += ["--use-conda"]
    for sdm in ("--software-deployment-method", "--deployment-method", "--deployment", "--sdm"):
        try:
            i = snakemake_args.index(sdm)
            if snakemake_args[i + 1] != "conda":
                logging.error(
                    "Software deployment method {} not implemented".format(snakemake_args[i + 1])
                )
                return 1
            snakemake_args.pop(i + 1)
            snakemake_args.pop(i)
        except ValueError:
            pass
    snakemake_argv += ["--software-deployment-method", "conda"]
    if "--conda-frontend" not in snakemake_args:
        snakemake_argv += ["--conda-frontend", "conda"]

    # Increase verbosity levels
    if wrapper_args.verbose:
        snakemake_argv += ["--verbose"]

    # Configure profile if snappy pipeline profile is requested
    if wrapper_args.profile_snappy_pipeline:
        if "--profile" in snakemake_args:
            logging.error("--profile-snappy-pipeline & --profile are mutually exclusive")
            return 1
        profile_path = os.path.join(os.path.dirname(__file__), "tpls", "profile")
        snakemake_argv += ["--profile", profile_path]

    # Add cores when missing
    if "--cores" not in snakemake_args:
        snakemake_argv += ["--cores", "1"]

    snakemake_argv += snakemake_args
    logging.info("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    return snakemake_main(snakemake_argv)


def main(argv=None):
    """Main program entry point, starts parsing command line arguments"""
    parser = argparse.ArgumentParser(
        usage="%(prog)s [--version] [-v] [-d directory] [--profile-snappy-pipeline] [snakemake arguments]",
        allow_abbrev=False,
    )

    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verobsity level")
    parser.add_argument(
        "-d", "--directory", default=os.getcwd(), help="Path to directory to run in, default is cwd"
    )
    parser.add_argument(
        "--profile-snappy-pipeline",
        action="store_true",
        help="Uses the profile defined in the snappy pipeline",
    )
    parser.add_argument(
        "--step",
        type=str,
        metavar="STEP",
        choices=sorted(STEP_TO_MODULE.keys()),
        help="The type of the step to run",
    )

    wrapper_args, snakemake_args = parser.parse_known_args(argv)

    # Setup logging
    setup_logging(wrapper_args)

    if not wrapper_args.step:
        for cfg in CONFIG_FILES:
            path = os.path.join(wrapper_args.directory, cfg)
            if not os.path.exists(path):
                continue
            with open(path, "rt") as f:
                yaml = ruamel_yaml.YAML()
                data = yaml.load(f.read())
            try:
                wrapper_args.step = data["pipeline_step"]["name"]
                break
            except KeyError:
                logging.info("Could not pick up pipeline step/name from %s", path)
    if not wrapper_args.step:
        parser.error("the following arguments are required: --step")
    return run(wrapper_args, snakemake_args)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
