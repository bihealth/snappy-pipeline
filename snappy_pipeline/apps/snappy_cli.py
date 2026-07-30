# -*- coding: utf-8 -*-
"""Unified Click-based CLI for the snappy pipeline.

Replaces the previous multi-script CLI entry points with a single `snappy` command.
"""

import datetime
import io
import logging
import os
import subprocess
import sys

import click
import ruamel.yaml as ruamel_yaml
from ruamel.yaml.comments import CommentedMap
from snakemake.cli import main as snakemake_main

from .. import __version__
from .impl.fsmanip import (
    assume_path_existing,
    assume_path_nonexisting,
    backup_file,
    create_directory,
    create_from_tpl,
    update_file,
)
from .impl.logging import LVL_ERROR, LVL_IMPORTANT, LVL_SUCCESS, log
from .impl.yaml_utils import remove_non_required, remove_yaml_comment_lines
from .snappy_snake import STEP_TO_MODULE


def _get_cwd():
    return os.getcwd()


#: Allowed steps
STEPS = tuple(sorted(STEP_TO_MODULE))


class TaskParamType(click.ParamType):
    name = "task"

    def convert(self, value, param, ctx):
        if isinstance(value, tuple) and len(value) == 2:
            return value
        if "=" in value:
            task_name, step = value.split("=", 1)
        else:
            task_name = step = value

        if step not in STEPS:
            self.fail(
                f"Invalid step type: {step}. Must be one of {', '.join(STEPS)}",
                param,
                ctx,
            )
        return task_name, step


#: README file name
README_FILENAME = "README.md"

#: File name for pipeline job file
FILENAME_PIPELINE_JOB_SH = "pipeline_job.sh"

#: Configuration sub directory
CONFIG_SUBDIR = ""

#: Configuration file name
CONFIG_FILENAME = "config.yaml"


class AddTaskAppException(Exception):
    """Raised in case of problem with adding a task."""


class AddTaskApp:
    """Implementation of task addition logic."""

    def __init__(
        self,
        step,
        task_name,
        project_directory,
        manage_config=True,
        partition="medium",
        email=None,
        conda="",
    ):
        self.step = step
        self.task_name = task_name
        self.project_directory = project_directory
        self.manage_config = manage_config
        self.partition = partition
        self.email = email
        self.conda = conda

    def run(self):
        log("")
        log(
            "Starting task {task_name} (using step {step}) in project dir {project_dir}",
            args={
                "step": self.step,
                "task_name": self.task_name,
                "project_dir": self.project_directory,
            },
        )

        try:
            config_yaml = self._load_config_yaml()
        except AddTaskAppException:
            return 1

        # Check if the task name already exists in the tasks list
        if "tasks" in config_yaml and isinstance(config_yaml["tasks"], list):
            for existing_task in config_yaml["tasks"]:
                if isinstance(existing_task, dict) and existing_task.get("name") == self.task_name:
                    log(
                        "Task with name {task_name} already present in configuration!",
                        args={"task_name": self.task_name},
                        level=LVL_ERROR,
                    )
                    return 1

        # Setup the configuration by appending the task
        if self.manage_config:
            self._setup_configuration(config_yaml)

        # Ensure pipeline_job.sh is present at project root
        self._ensure_pipeline_job_sh()

        log(
            "\nDo not forget to fill out the REQUIRED fields in the project configuration file!\n",
            level=LVL_IMPORTANT,
        )
        log(
            "Task {task_name} (using step {step}) created.",
            args={"task_name": self.task_name, "step": self.step},
            level=LVL_SUCCESS,
        )
        return 0

    def _load_config_yaml(self):
        config_filename = os.path.join(self.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
        if not os.path.exists(config_filename):
            raise AddTaskAppException(f"Configuration file does not exist at {config_filename}")
        with open(config_filename, "rt") as f:
            yaml = ruamel_yaml.YAML()
            config_yaml = yaml.load(f.read())
        return config_yaml

    def _ensure_pipeline_job_sh(self):
        dest_path = os.path.join(self.project_directory, FILENAME_PIPELINE_JOB_SH)
        if os.path.exists(dest_path):
            return
        create_from_tpl(
            src_path=os.path.join(os.path.dirname(__file__), "tpls", FILENAME_PIPELINE_JOB_SH),
            dest_path=dest_path,
            format_args={
                "line_m": "##SBATCH --mail-type ALL"
                if not self.email
                else "#SBATCH --mail-type ALL",
                "line_M": (
                    "##SBATCH --mail-user your.name@mdc-berlin.de"
                    if not self.email
                    else "##SBATCH --mail-user {}".format(self.email)
                ),
                "partition": self.partition,
                "conda": self.conda,
                "project_name": os.path.basename(self.project_directory),
            },
            message="Creating master job shell file in {path}",
            message_args={"path": dest_path},
        )

    def _setup_configuration(self, config_yaml):
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
        task_block["name"] = self.task_name
        task_block["config"] = step_config_block

        config_yaml["tasks"].append(task_block)

        # Create backup of config.yaml file and overwrite with new string
        config_filename = os.path.join(self.project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
        backup_file(config_filename)
        yaml = ruamel_yaml.YAML()
        buf = io.StringIO()
        yaml.dump(config_yaml, stream=buf)
        buf.seek(0)
        contents = buf.read()
        update_file(
            path=config_filename,
            contents=contents,
            message="Updating project config with required default config for task {task_name} (step {step})",
            message_args={"task_name": self.task_name, "step": self.step},
        )


def setup_logging(verbose):
    """Setup logger."""
    logging.basicConfig(
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
    )
    logger = logging.getLogger("")
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


@click.group()
@click.version_option(__version__, "--version", help="Show version and exit")
def main():
    """SNAPPY Nucleic Acid Processing in Python."""
    pass


@main.command()
@click.option(
    "--directory",
    "-d",
    "project_directory",
    required=True,
    type=click.Path(),
    help="Path to directory to create for the project",
)
@click.option(
    "--project-name",
    type=str,
    help="A string to use for the project name, used in the README file",
)
@click.option("--partition", default="medium", help="Partition to submit into")
@click.option(
    "--task",
    "-t",
    "tasks",
    type=TaskParamType(),
    multiple=True,
    metavar="[NAME=]STEP_TYPE",
    help="List of tasks to create automatically, formatted as [name=]step_type",
)
@click.option(
    "--manage-config/--no-manage-config",
    default=True,
    help="Manage the project config.yaml file automatically",
)
@click.option("--email", type=str, help="Email address for pipeline_job.sh file")
@click.option(
    "--conda",
    type=str,
    default="",
    help="conda environment to load when submitting job",
)
def init(project_directory, project_name, partition, tasks, manage_config, email, conda):
    """Initialize a new snappy project directory."""
    log("SNAPPY Pipeline -- start_project")
    log("================================")
    log("")

    # Check if directory already exists - no overwrite
    if not assume_path_nonexisting(project_directory):
        sys.exit(1)

    # Create project directory and subdirectory for configuration files
    paths = [project_directory]
    if CONFIG_SUBDIR:
        paths.append(os.path.join(project_directory, CONFIG_SUBDIR))
    for path in paths:
        create_directory(path)

    # Create config file in subdirectory based on template
    config_dest_path = os.path.join(project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", "project_config.yaml"),
        dest_path=config_dest_path,
        format_args={
            "created_at": datetime.datetime.now().isoformat(),
            "project_name": (project_name or os.path.basename(project_directory)),
        },
        message="Creating project-wide configuration in {path}",
        message_args={"path": config_dest_path},
    )

    # Create readme file in subdirectory based on template
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", README_FILENAME),
        dest_path=os.path.join(project_directory, README_FILENAME),
        format_args={
            "created_at": datetime.datetime.now().isoformat(),
            "project_name": (project_name or os.path.basename(project_directory)),
        },
        message="Creating README file in in {path}",
        message_args={"path": os.path.join(project_directory, README_FILENAME)},
    )

    # Create master job shell file in project directory based on template
    dest_path = os.path.join(project_directory, "pipeline_job.sh")
    email_val = email or os.environ.get("SNAPPY_PIPELINE_EMAIL")
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", "pipeline_job.sh"),
        dest_path=dest_path,
        format_args={
            "line_m": "##SBATCH --mail-type ALL" if not email_val else "#SBATCH --mail-type ALL",
            "line_M": (
                "##SBATCH --mail-user your.name@mdc-berlin.de"
                if not email_val
                else "##SBATCH --mail-user {}".format(email_val)
            ),
            "partition": partition,
            "conda": conda,
            "project_name": (project_name or os.path.basename(project_directory)),
        },
        message="Creating master job shell file in {path}",
        message_args={"path": dest_path},
    )

    # Create additional tasks if any was provided
    for task_name, step in tasks:
        app = AddTaskApp(
            step=step,
            task_name=task_name,
            project_directory=project_directory,
            manage_config=manage_config,
            partition=partition,
            email=email_val,
            conda=conda,
        )
        app.run()

    log(
        "\nDo not forget to review config.yaml and to fill out README.md!\n",
        level=LVL_IMPORTANT,
    )
    log("All done, have a nice day!", level=LVL_SUCCESS)


@main.group()
def task():
    """Manage tasks within a snappy project."""
    pass


@task.command(name="add")
@click.argument(
    "tasks", nargs=-1, required=True, type=TaskParamType(), metavar="[NAME=]STEP_TYPE..."
)
@click.option(
    "--project-directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.option("--partition", default="medium", help="Partition to submit into")
@click.option("--email", type=str, help="Email address for pipeline_job.sh file")
@click.option(
    "--manage-config/--no-manage-config",
    default=True,
    help="Manage the project config.yaml file automatically",
)
@click.option(
    "--conda",
    type=str,
    default="",
    help="conda environment to load when submitting job",
)
def task_add(
    tasks,
    project_directory,
    partition,
    email,
    manage_config,
    conda,
):
    """Add tasks to the project configuration.

    TASKS is a list of tasks to add, formatted as [task_name=]step_type.
    """
    email_val = email or os.environ.get("SNAPPY_EMAIL") or os.environ.get("SNAPPY_PIPELINE_EMAIL")

    for task_name, step in tasks:
        app = AddTaskApp(
            step=step,
            task_name=task_name,
            project_directory=project_directory,
            manage_config=manage_config,
            partition=partition,
            email=email_val,
            conda=conda,
        )
        res = app.run()
        if res:
            sys.exit(res)


@task.command(name="list")
@click.option(
    "--project-directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
def task_list(project_directory):
    """List all defined tasks in the project."""
    config_filename = os.path.join(project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
    if not os.path.exists(config_filename):
        click.echo(f"Error: Configuration file not found at {config_filename}", err=True)
        sys.exit(1)

    with open(config_filename, "rt") as f:
        yaml = ruamel_yaml.YAML()
        config_yaml = yaml.load(f.read())

    tasks = config_yaml.get("tasks", [])
    if not tasks:
        click.echo("No tasks defined in config.yaml")
        return

    click.echo(f"Defined tasks in {config_filename}:")
    for task in tasks:
        if isinstance(task, dict):
            name = task.get("name", "unnamed")
            step = task.get("step", "unknown")
            click.echo(f"  - {name} ({step})")


@main.command(
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.option(
    "-d",
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Path to directory to run in, default is cwd",
)
@click.option(
    "--slurm",
    is_flag=True,
    help="Enable SLURM executor profile (adds executor, job limits, resource defaults)",
)
@click.option(
    "--task",
    "task_name",
    type=str,
    default=None,
    help="The specific task name from config.yaml to run",
)
@click.option(
    "--all-tasks",
    "all_tasks",
    is_flag=True,
    default=False,
    help=(
        "Target all tasks, not just leaf tasks. "
        "By default only leaf tasks (tasks not depended on by any other task) are targeted."
    ),
)
@click.option("-v", "--verbose", is_flag=True, help="Increase verbosity level")
@click.pass_context
def run(ctx, directory, slurm, task_name, all_tasks, verbose):
    """Run snappy pipeline workflows."""
    setup_logging(verbose)
    if task_name:
        logging.info("Targeting single task: %s", task_name)
    elif all_tasks:
        logging.info("Targeting all tasks (--all-tasks flag set).")
    else:
        logging.info("No specific --task provided. Targeting leaf tasks only (default).")

    snakemake_args = list(ctx.args)

    # Point to the master orchestrator Snakefile
    orchestrator_snakefile = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Snakefile"
    )

    snakemake_argv = [
        "--directory",
        directory,
        "--snakefile",
        orchestrator_snakefile,
    ]

    config_args = ["--config"]
    if task_name:
        config_args.append(f"task={task_name}")
    if all_tasks:
        config_args.append("all_tasks=True")
    if verbose:
        config_args.append("dump_orchestrator=True")

    if len(config_args) > 1:
        snakemake_argv.extend(config_args)

    # Always pass the default conda workflow profile
    default_profile = os.path.join(os.path.dirname(__file__), "profile")
    snakemake_argv += ["--workflow-profile", default_profile]

    # Layer the SLURM profile on top when --slurm is requested
    if slurm:
        slurm_profile = os.path.join(os.path.dirname(__file__), "profile-slurm")
        snakemake_argv += ["--workflow-profile", slurm_profile]

    # Append all user-provided snakemake arguments directly
    snakemake_argv += snakemake_args

    logging.info("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    res = snakemake_main(snakemake_argv)
    if res != 0:
        ctx.exit(res)


@main.command()
@click.option(
    "--directory",
    "-d",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.option(
    "--db-path",
    type=click.Path(),
    default=None,
    help="Custom path to the snkmt database. Defaults to <directory>/.snakemake/log/snkmt.sqlite",
)
def watch(directory, db_path):
    """Launch snkmt TUI to monitor a running snakemake workflow.

    Opens the snkmt interactive console for the snkmt database in the
    project directory (or a custom path).  The database is written by
    Snakemake when the profile has ``logger: snkmt`` (the default in
    snappy's workflow profile).
    """
    path = db_path or os.path.join(directory, ".snakemake", "log", "snkmt.sqlite")
    if not os.path.exists(path):
        log(
            "snkmt database not found at {path}.\n\n"
            "Either:\n"
            "  - Start a snakemake run first (snappy run) to generate it\n"
            "  - Pass --db-path with the correct database location",
            {"path": path},
            level=LVL_ERROR,
        )
        sys.exit(1)

    cmd = [sys.executable, "-m", "snkmt", "console", "--db-path", path]
    log(
        "Launching snkmt console for database at {path}",
        {"path": path},
        level=LVL_IMPORTANT,
    )
    res = subprocess.run(cmd)
    if res.returncode != 0:
        sys.exit(res.returncode)


@main.command()
@click.option(
    "--project-directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.option("--partition", default="medium", help="Partition to submit into")
@click.option("--email", type=str, help="Email address for pipeline_job.sh file")
@click.option(
    "--conda",
    type=str,
    default="",
    help="conda environment to load when submitting job",
)
def refresh(project_directory, partition, email, conda):
    """Recreate the master pipeline_job.sh and ensure slurm_log exists."""
    log("CUBI Pipeline -- refresh")
    log("========================")

    if not assume_path_existing(project_directory):
        sys.exit(1)

    config_filename = os.path.join(project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
    if not os.path.exists(config_filename):
        log(
            "Configuration file does not exist at {path}",
            {"path": config_filename},
            level=LVL_ERROR,
        )
        sys.exit(1)

    email_val = email or os.environ.get("SNAPPY_EMAIL") or os.environ.get("SNAPPY_PIPELINE_EMAIL")

    create_directory(os.path.join(project_directory, "slurm_log"), exist_ok=True)

    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", FILENAME_PIPELINE_JOB_SH),
        dest_path=os.path.join(project_directory, FILENAME_PIPELINE_JOB_SH),
        format_args={
            "line_m": "##SBATCH --mail-type ALL" if not email_val else "#SBATCH --mail-type ALL",
            "line_M": (
                "##SBATCH --mail-user your.name@mdc-berlin.de"
                if not email_val
                else "##SBATCH --mail-user {}".format(email_val)
            ),
            "partition": partition,
            "conda": conda,
            "project_name": os.path.basename(project_directory),
        },
        message="creating SGE job shell file in {path}",
        message_args={"path": os.path.join(project_directory, FILENAME_PIPELINE_JOB_SH)},
    )
    log("all done, have a nice day!", level=LVL_SUCCESS)


@main.command()
@click.argument("jobid")
def status(jobid):
    """Check the status of a SLURM job ID."""
    output = str(
        subprocess.check_output(
            "sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid,
            shell=True,
        ).strip()
    )

    running_status = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
    if "COMPLETED" in output:
        click.echo("success")
    elif any(r in output for r in running_status):
        click.echo("running")
    else:
        click.echo("failed")


@main.command(name="pull-sheet")
@click.option(
    "-o",
    "--output",
    type=click.File("wt"),
    default="-",
    help="Destination file, default is stdout.",
)
@click.option("--api-key", required=True, help="API key to use.")
@click.option("--sodar-host", required=True, help="SODAR host to use.")
@click.option("--project-uuid", required=True, help="UUID of project to query")
@click.option(
    "--project-name",
    required=True,
    help="Name of the project for output (machine-readable).",
)
@click.option(
    "--project-title",
    default=".",
    help="Title of the project (human-readable).",
)
@click.option(
    "--project-description",
    default=".",
    help="Description of the project (human-readable).",
)
@click.option("--library-types", help="Library type(s) to use, comma-separated")
def pull_sheet(
    output,
    api_key,
    sodar_host,
    project_uuid,
    project_name,
    project_title,
    project_description,
    library_types,
):
    """Pull biomedsheet sample sheet from SODAR API."""
    library_types_list = library_types.split(",") if library_types else []

    class PullSheetArgs:
        def __init__(
            self,
            output,
            api_key,
            sodar_host,
            project_uuid,
            project_name,
            project_title,
            project_description,
            library_types,
        ):
            self.output = output
            self.api_key = api_key
            self.sodar_host = sodar_host
            self.project_uuid = project_uuid
            self.project_name = project_name
            self.project_title = project_title
            self.project_description = project_description
            self.library_types = library_types

    args = PullSheetArgs(
        output,
        api_key,
        sodar_host,
        project_uuid,
        project_name,
        project_title,
        project_description,
        library_types_list,
    )
    from .snappy_pull_sheet import run as run_pull_sheet

    run_pull_sheet(args)


if __name__ == "__main__":
    main()
