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
    backup_file,
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
CONFIG_SUBDIR = ""

#: Configuration file name
CONFIG_FILENAME = "config.yaml"


class StartStepAppException(Exception):
    """Raised in case of problem with starting a step."""


class StartStepApp:
    """Implementation of ``snappy-start-step``."""

    def __init__(self, step, directory, args):
        #: The step name
        self.step = step
        #: Task name (maps to directory parameter)
        self.directory = directory
        #: Parsed command line arguments
        self.args = args

    def run(self):
        """Actually perform the step."""
        log("")
        log(
            'Starting step "{step}" with task name "{task_name}" in project dir "{project_dir}"',
            args={
                "step": self.step,
                "task_name": self.directory,
                "project_dir": self.args.project_directory,
            },
        )

        # Load project-wide configuration
        try:
            config_yaml = self._load_config_yaml()
        except StartStepAppException:
            return 1

        # Check if the task name already exists in the tasks list
        if "tasks" in config_yaml and isinstance(config_yaml["tasks"], list):
            for existing_task in config_yaml["tasks"]:
                if isinstance(existing_task, dict) and existing_task.get("name") == self.directory:
                    log(
                        "task with name {task_name} already present in configuration!",
                        args={"task_name": self.directory},
                        level=LVL_ERROR,
                    )
                    return 1

        # Setup the configuration by appending the task
        if self.args.manage_config:
            self._setup_configuration(config_yaml)

        # Ensure pipeline_job.sh is present at project root
        self._ensure_pipeline_job_sh()

        log(
            "\nDo not forget to fill out the REQUIRED fields in the project configuration file!\n",
            level=LVL_IMPORTANT,
        )
        log(
            "Task {task_name} for step {step} created.",
            args={"task_name": self.directory, "step": self.step},
            level=LVL_SUCCESS,
        )

    def _load_config_yaml(self):
        """Load configuration."""
        config_filename = os.path.join(self.args.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
        with open(config_filename, "rt") as f:
            yaml = ruamel_yaml.YAML()
            config_yaml = yaml.load(f.read())
        return config_yaml

    def _ensure_pipeline_job_sh(self):
        dest_path = os.path.join(self.args.project_directory, FILENAME_PIPELINE_JOB_SH)
        if os.path.exists(dest_path):
            return
        create_from_tpl(
            src_path=os.path.join(os.path.dirname(__file__), "tpls", FILENAME_PIPELINE_JOB_SH),
            dest_path=dest_path,
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
                "step_name": os.path.basename(self.args.project_directory),
            },
            message="Creating master job shell file in {path}",
            message_args={"path": dest_path},
        )

    def _setup_configuration(self, config_yaml):
        """Setup configuration settings."""
        if "tasks" not in config_yaml:
            config_yaml["tasks"] = []

        # Load default configuration, remove comment lines and lines not marked as required
        yaml = ruamel_yaml.YAML()
        default_config_yaml = yaml.load(
            remove_yaml_comment_lines(STEP_TO_MODULE[self.step].DEFAULT_CONFIG)
        )
        only_required = remove_non_required(default_config_yaml)

        step_config_block = None
        if only_required and "step_config" in only_required:
            step_name_key = self.step
            if step_name_key in only_required["step_config"]:
                step_config_block = only_required["step_config"][step_name_key]

        if step_config_block is None:
            step_config_block = CommentedMap()

        task_block = CommentedMap()
        task_block["step"] = self.step
        task_block["name"] = self.directory
        task_block["config"] = step_config_block

        config_yaml["tasks"].append(task_block)

        # Create backup of config.yaml file and overwrite with new string
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
