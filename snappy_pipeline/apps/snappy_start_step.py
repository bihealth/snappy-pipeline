# -*- coding: utf-8 -*-
"""Tool for starting a pipeline step instance

The tool will update the existing global configuration with sensible default
settings for the given step.
"""

import argparse
import io
import os
import sys

import ruamel.yaml as ruamel_yaml
from ruamel.yaml.comments import CommentedMap

from .. import __version__
from .impl.fsmanip import (
    assume_path_nonexisting,
    backup_file,
    create_directory,
    create_from_tpl,
    update_file,
)
from .impl.logging import LVL_ERROR, LVL_IMPORTANT, LVL_SUCCESS, log
from .impl.yaml_utils import remove_non_required, remove_yaml_comment_lines
from .snappy_snake import STEP_TO_MODULE

#: Allowed steps
STEPS = tuple(sorted(STEP_TO_MODULE))

#: File name for pipeline job file
FILENAME_PIPELINE_JOB_SH = "pipeline_job.sh"

#: Configuration sub directory
CONFIG_SUBDIR = ".snappy_pipeline"

#: Configuration file name
CONFIG_FILENAME = "config.yaml"


class StartStepAppException(Exception):
    """Raised in case of problem with starting a step."""


class StartStepApp:
    """Implementation of ``snappy-start-step``."""

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
            'Starting step "{step}" in sub-directory "{directory}" of project dir "{project_dir}"',
            args={
                "step": self.step,
                "directory": self.directory,
                "project_dir": self.args.project_directory,
            },
        )

        dest_dir = os.path.join(self.args.project_directory, self.directory)
        if not assume_path_nonexisting(dest_dir):
            return 1  # pragma: nocover

        # Load project-wide configuration and check that the step does not exist yet
        # (if configured).
        try:
            config_yaml = self._load_config_yaml()
        except StartStepAppException:
            return 1

        # Setup the step sub directory.
        self._setup_step_dir(dest_dir, config_yaml)

        # Setup the configuration
        if self.args.manage_config:
            self._setup_configuration(config_yaml)

        log(
            "\nDo not forget to fill out the REQUIRED fields in the project configuration file!\n",
            level=LVL_IMPORTANT,
        )
        log("Step {step} created.", args={"step": self.step}, level=LVL_SUCCESS)

    def _load_config_yaml(self):
        """Load configuration."""
        config_filename = os.path.join(self.args.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
        with open(config_filename, "rt") as f:
            yaml = ruamel_yaml.YAML()
            config_yaml = yaml.load(f.read())
        if (
            self.args.manage_config
            and "step_config" in config_yaml
            and self.step in config_yaml["step_config"]
        ):
            log(
                "configuration for step {step} (/step_config/{step}) already present!",
                args={"step": self.step},
                level=LVL_ERROR,
            )
            log(
                (
                    "please comment out in {path}, re-run start_step, and merge configuration "
                    "settings manually"
                ),
                args={"path": config_filename},
                level=LVL_ERROR,
            )
            raise StartStepAppException("Config already exists")
        return config_yaml

    def _setup_step_dir(self, dest_dir, config_yaml):
        """Create and setup the step sub directory."""
        create_directory(dest_dir)
        create_directory(os.path.join(dest_dir, "slurm_log"))

        create_from_tpl(
            src_path=os.path.join(os.path.dirname(__file__), "tpls", "step_config.yaml"),
            dest_path=os.path.join(dest_dir, CONFIG_FILENAME),
            format_args={"step_name": self.step, "step_version": 1, "config_subdir": CONFIG_SUBDIR},
            message="Creating step-wide configuration in {path}",
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
                "conda": self.args.conda,
                "step_name": self.step,
            },
            message="Creating SGE job shell file in {path}",
            message_args={"path": os.path.join(dest_dir, FILENAME_PIPELINE_JOB_SH)},
        )

    def _setup_configuration(self, config_yaml):
        """Setup configuration settings."""
        # Ensure that a "step_config" setting is present and block style is used for it
        if "step_config" not in config_yaml:
            config_yaml["step_config"] = CommentedMap()
        config_yaml["step_config"].fa.set_block_style()

        # Load default configuration, remove comment lines and lines not marked as required;
        # preserve comments
        yaml = ruamel_yaml.YAML()
        default_config_yaml = yaml.load(
            remove_yaml_comment_lines(STEP_TO_MODULE[self.step].DEFAULT_CONFIG)
        )
        only_required = remove_non_required(default_config_yaml)
        if only_required:
            config_yaml["step_config"][self.step] = only_required["step_config"][self.step]

        # Create backup of config.yaml file and overwrite with new string, showing diff
        config_filename = os.path.join(self.args.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
        backup_file(config_filename)
        yaml = ruamel_yaml.YAML()
        buf = io.StringIO()
        yaml.dump(config_yaml, stream=buf)
        buf.seek(0)
        contents = buf.read()
        update_file(
            path=config_filename,
            contents=contents,
            message="Updating project config with required default config for step {step}",
            message_args={"step": self.step},
        )


def run_start_step(step, directory, args):
    """Run ``snappy-start-step``.

    :param step: Name of step.
    :type step: str

    :param directory: Name of directory to store step configurations and results. Usually the same
    name as the step.
    :type directory: str

    :param args: Arguments provided by the user.
    :type args: argparse.Namespace

    :return: Returns the return value of Start Step run call.
    """
    return StartStepApp(step, directory, args).run()


def run(args):
    """Program entry point after argument parsing"""
    log("CUBI Pipeline -- start_step")
    log("===========================")
    for step, directory in zip(args.steps, args.directories):
        res = run_start_step(step, directory, args)
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

    parser.add_argument("--partition", default="medium", help="Partition to submit into")

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

    parser.add_argument(
        "--no-manage-config",
        dest="manage_config",
        default=True,
        action="store_false",
        help=(
            "Do not check config.yaml for existing configuration or change it (IOW: leave it alone)"
        ),
    )

    parser.add_argument(
        "--conda",
        type=str,
        nargs="?",
        default="",
        help="conda environment to load when submitting job",
    )

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
