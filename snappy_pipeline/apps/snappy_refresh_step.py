# -*- coding: utf-8 -*-
"""Tool for refreshing pipeline_job.sh file for a step.

Global configuration will not be touched.
"""

import argparse
import os
import sys

import ruamel.yaml as yaml

from .snappy_snake import STEP_TO_MODULE
from .impl.logging import LVL_SUCCESS, log
from .impl.fsmanip import (
    assume_path_existing,
    create_directory,
    create_from_tpl,
)
from .. import __version__

#: Allowed steps
STEPS = tuple(sorted(STEP_TO_MODULE))

#: File name for pipeline job file
FILENAME_PIPELINE_JOB_SH = "pipeline_job.sh"

#: Configuration sub directory
CONFIG_SUBDIR = ".snappy_pipeline"

#: Configuration file name
CONFIG_FILENAME = "config.yaml"


class RefreshStepAppException(Exception):
    """Raised in case of problem with refreshing a step."""


class RefreshStepApp:
    """Implementation of ``snappy-refresh-step``."""

    def __init__(self, step, directory, args):
        #: The step name
        self.step = step
        #: Directory to create in
        self.directory = directory
        #: Parsed command line arguments
        self.args = args

    def run(self):
        """Actually perform the step."""
        log("")
        log(
            'Refresh step "{step}" in sub-directory "{directory}" of project dir "{project_dir}"',
            args={"step": self.step, "directory": self.directory, "project_dir": self.directory},
        )

        dest_dir = os.path.join(self.args.project_directory, self.directory)
        if not assume_path_existing(dest_dir):
            return 1  # pragma: nocover

        # Load project-wide configuration and check that the step already exists.
        try:
            config_yaml = self._load_config_yaml()
        except RefreshStepAppException:
            return 1

        # Re-setup the step sub directory.
        self._resetup_step_dir(dest_dir, config_yaml)

        log("all done, have a nice day!", level=LVL_SUCCESS)

    def _load_config_yaml(self):
        """Load configuration."""
        config_filename = os.path.join(self.args.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
        with open(config_filename, "rt") as f:
            return yaml.round_trip_load(f.read())

    def _resetup_step_dir(self, dest_dir, config_yaml):
        """Re-setup the step sub directory."""
        create_directory(dest_dir, exist_ok=True)
        create_directory(os.path.join(dest_dir, "slurm_log"), exist_ok=True)

        create_from_tpl(
            src_path=os.path.join(os.path.dirname(__file__), "tpls", "step_config.yaml"),
            dest_path=os.path.join(dest_dir, CONFIG_FILENAME),
            format_args={"step_name": self.step, "step_version": 1, "config_subdir": CONFIG_SUBDIR},
            message="creating step-wide configuration in {path}",
            message_args={"path": os.path.join(dest_dir, CONFIG_FILENAME)},
        )

        create_from_tpl(
            src_path=os.path.join(os.path.dirname(__file__), "tpls", FILENAME_PIPELINE_JOB_SH),
            dest_path=os.path.join(dest_dir, FILENAME_PIPELINE_JOB_SH),
            format_args={
                "line_m": (
                    "##SBATCH --mail-type ALL" if not self.args.email else "#SBATCH --mail-type ALL"
                ),
                "line_M": (
                    "##SBATCH --mail-user your.name@mdc-berlin.de"
                    if not self.args.email
                    else "##SBATCH --mail-user {}".format(self.args.email)
                ),
                "partition": self.args.partition,
                "step_name": self.step,
            },
            message="creating SGE job shell file in {path}",
            message_args={"path": os.path.join(dest_dir, FILENAME_PIPELINE_JOB_SH)},
        )


def run_refresh_step(step, directory, args):
    """Run ``snappy-refresh-step``."""
    return RefreshStepApp(step, directory, args).run()


def run(args):
    """Program entry point after argument parsing"""
    log("CUBI Pipeline -- refresh_step")
    log("=============================")
    for step, directory in zip(args.steps, args.directories):
        res = run_refresh_step(step, directory, args)
        if res:
            return res


def main(argv=None):
    """Program entry point including command line argument parsing"""
    epilog = 'Valid choices for "--step": {}'.format(", ".join(STEPS))
    parser = argparse.ArgumentParser(epilog=epilog)

    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)

    parser.add_argument(
        "--step",
        type=str,
        metavar="STEP",
        required=True,
        nargs="+",
        dest="steps",
        choices=sorted(STEP_TO_MODULE.keys()),
        default=[],
        action="append",
        help="The type of the step to run",
    )

    parser.add_argument(
        "--project-directory",
        default=os.getcwd(),
        help="Project directory, defaults to current working directory",
    )

    parser.add_argument(
        "--directory",
        type=str,
        nargs="+",
        dest="directories",
        action="append",
        default=[],
        help="Name of directory/ies to create.  Defaults to step name(s)",
    )

    parser.add_argument(
        "--email",
        type=str,
        help=(
            "Email address for pipeline_job.sh file.  You can also set the environment variable "
            "SNAPPY_EMAIL"
        ),
    )

    parser.add_argument("--partition", type=str, help="The partition to run in", default="medium")

    args = parser.parse_args(argv)
    # Flatten ``--step`` and ``-directory`` argument
    args.steps = [item for sublist in args.steps for item in sublist]
    args.directories = [item for sublist in args.directories for item in sublist]

    if not args.directories:
        args.directories = list(args.steps)
    elif len(args.directories) != len(args.steps):
        msg = (
            "Either leave --directory empty or give same number of argments as for --step! "
            "{} vs {} values".format(len(args.directories), len(args.steps))
        )
        raise RuntimeError(msg)

    args.email = args.email or os.environ.get("SNAPPY_EMAIL")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
