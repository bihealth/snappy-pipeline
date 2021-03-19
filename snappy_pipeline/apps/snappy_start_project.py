# -*- coding: utf-8 -*-
"""Tool for starting a project (a collection of pipeline step instances)

The tool will create the global configuration and a number of steps for you.
"""

import argparse
import datetime
import os
import sys

from .snappy_snake import STEP_TO_MODULE
from .snappy_start_step import run_start_step, CONFIG_FILENAME, CONFIG_SUBDIR
from .impl.logging import LVL_IMPORTANT, LVL_SUCCESS, log
from .impl.fsmanip import assume_path_nonexisting, create_directory, create_from_tpl
from .. import __version__

#: Allowed steps
STEPS = tuple(sorted(STEP_TO_MODULE))

#: README file name
README_FILENAME = "README.md"


def run(args):
    """Perform the actual creation of stuff"""
    log("SNAPPY Pipeline -- start_project")
    log("================================")
    log("")

    # Check if directory already exists - no overwrite
    if not assume_path_nonexisting(args.project_directory):
        return 1

    # Create project directory and subdirectory for configuration files
    paths = (args.project_directory, os.path.join(args.project_directory, CONFIG_SUBDIR))
    for path in paths:
        create_directory(path)

    # Create config file in subdirectory based on template
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", "project_config.yaml"),
        dest_path=os.path.join(args.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME),
        format_args={
            "created_at": datetime.datetime.now().isoformat(),
            "project_name": (args.project_name or os.path.basename(args.project_directory)),
        },
        message="Creating project-wide configuration in {path}",
        message_args={"path": os.path.join(args.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)},
    )

    # Create readme file in subdirectory based on template
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", README_FILENAME),
        dest_path=os.path.join(args.project_directory, README_FILENAME),
        format_args={
            "created_at": datetime.datetime.now().isoformat(),
            "project_name": (args.project_name or os.path.basename(args.project_directory)),
        },
        message="Creating README file in in {path}",
        message_args={"path": os.path.join(args.project_directory, README_FILENAME)},
    )

    # Create additional steps if any was provided
    for step in args.steps:
        run_start_step(step=step, directory=step, args=args)

    log(
        "\nDo not forget to review .snappy_pipeline/config.yaml and to fill out README.md!\n",
        level=LVL_IMPORTANT,
    )
    log("All done, have a nice day!", level=LVL_SUCCESS)


def main(argv=None):
    epilog = 'Valid choices for "--steps": {}'.format(", ".join(STEPS))
    parser = argparse.ArgumentParser(epilog=epilog)

    parser.add_argument(
        "--version",
        action="version",
        version="%%(prog)s %s" % __version__,
        help="Show version and exit",
    )

    parser.add_argument("--partition", default="medium", help="Partition to submit into")

    parser.add_argument(
        "--directory",
        required=True,
        dest="project_directory",
        help="Path to directory to create for the project",
    )

    parser.add_argument(
        "--project-name",
        type=str,
        help="A string to use for the project name, used in the README file",
    )

    parser.add_argument(
        "--steps",
        metavar="STEP",
        nargs="+",
        choices=STEPS,
        dest="steps",
        action="append",
        default=[],
        help=(
            "List of steps to create instances of automatically (otherwise, use snappy-start-step)"
        ),
    )

    parser.add_argument(
        "--no-manage-config",
        dest="manage_config",
        default=True,
        action="store_false",
        help=(
            "Do not check config.yaml for existing configuration or change it (IOW: "
            "leave it alone)"
        ),
    )

    parser.add_argument(
        "--email",
        type=str,
        help=(
            "Email address for pipeline_job.sh file.  You can also set the environment variable "
            "SNAPPY_PIPELINE_EMAIL"
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
    args.steps = [item for sublist in args.steps for item in sublist]
    args.email = args.email or os.environ.get("SNAPPY_PIPELINE_EMAIL")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
