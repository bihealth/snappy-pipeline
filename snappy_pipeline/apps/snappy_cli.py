# -*- coding: utf-8 -*-
"""Unified Click-based CLI for the snappy pipeline.

Replaces the previous multi-script CLI entry points with a single `snappy` command.
"""

import datetime
import enum
import io
import logging
import os
import shutil
import subprocess
import sys

import rich_click as click
import ruamel.yaml as ruamel_yaml
from rich.console import Console
from rich.table import Table
from ruamel.yaml.comments import CommentedMap
from snakemake.cli import main as snakemake_main
from snakemake.common.configfile import load_configfile

from .. import __version__
from ..log_archive import build_archive, latest_snakemake_log
from ..orchestration import load_project, missing_files
from ..workflow_registry import WORKFLOW_REGISTRY
from .impl.fsmanip import (
    backup_file,
    create_directory,
    create_from_tpl,
    update_file,
)
from .impl.logging import LVL_ERROR, LVL_IMPORTANT, LVL_SUCCESS, log
from .impl.yaml_utils import remove_non_required


class ConfigMode(str, enum.Enum):
    """Configuration generation mode for tasks."""

    COMMENTED = "commented"
    FULL = "full"
    MINIMAL = "minimal"


class ShellType(str, enum.Enum):
    """Supported shells for autocompletion."""

    BASH = "bash"
    ZSH = "zsh"
    FISH = "fish"


class JobStatus(str, enum.Enum):
    """SLURM job status strings."""

    SUCCESS = "success"
    RUNNING = "running"
    FAILED = "failed"


# Configure rich-click settings for formatting
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True


def _get_cwd():
    return os.getcwd()


#: Allowed steps
STEPS = tuple(sorted(WORKFLOW_REGISTRY))


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

    def shell_complete(self, ctx, param, incomplete):
        from click.shell_completion import CompletionItem

        if "=" in incomplete:
            prefix, step_part = incomplete.split("=", 1)
            return [
                CompletionItem(f"{prefix}={step}") for step in STEPS if step.startswith(step_part)
            ]
        return [CompletionItem(step) for step in STEPS if step.startswith(incomplete)]


def _complete_task_names(ctx, param, incomplete):
    from click.shell_completion import CompletionItem

    directory = ctx.params.get("directory")
    if callable(directory):
        directory = directory()
    directory = directory or _get_cwd()

    config_filename = os.path.join(directory, CONFIG_SUBDIR, CONFIG_FILENAME)
    task_names = set()
    if os.path.exists(config_filename):
        try:
            with open(config_filename, "rt") as f:
                yaml = ruamel_yaml.YAML()
                config_yaml = yaml.load(f.read())
                if (
                    config_yaml
                    and "tasks" in config_yaml
                    and isinstance(config_yaml["tasks"], list)
                ):
                    for task_item in config_yaml["tasks"]:
                        if isinstance(task_item, dict) and "name" in task_item:
                            task_names.add(task_item["name"])
        except Exception:
            pass
    if not task_names:
        task_names = set(STEPS)

    return [CompletionItem(name) for name in sorted(task_names) if name.startswith(incomplete)]


def _get_registered_steps_info():
    info = []
    for step_name in STEPS:
        cls = WORKFLOW_REGISTRY[step_name]
        doc = cls.__doc__ or ""
        summary = doc.strip().split("\n")[0] if doc else "No description available."
        info.append((step_name, summary))
    return info


def _print_steps_table():
    steps_info = _get_registered_steps_info()
    console = Console()
    table = Table(
        title="Available Workflow Steps in SNAPPY", show_header=True, header_style="bold magenta"
    )
    table.add_column("Step Name", style="cyan", no_wrap=True)
    table.add_column("Description", style="white")
    for step_name, summary in steps_info:
        table.add_row(step_name, summary)
    console.print(table)


#: README file name
README_FILENAME = "README.md"

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
        step: str,
        task_name: str,
        project_directory: str,
        manage_config: bool = True,
        config_mode: ConfigMode | str = ConfigMode.COMMENTED,
    ):
        self.step = step
        self.task_name = task_name
        self.project_directory = project_directory
        self.manage_config = manage_config
        self.config_mode = ConfigMode(config_mode)

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

    def _setup_configuration(self, config_yaml):
        if "tasks" not in config_yaml:
            config_yaml["tasks"] = []

        yaml = ruamel_yaml.YAML()
        step_cls = WORKFLOW_REGISTRY[self.step]

        comment_optional = self.config_mode == ConfigMode.COMMENTED
        step_config_yaml_str = step_cls.config_model_class.default_config_yaml_string(
            comment_optional=comment_optional
        )

        if self.config_mode == ConfigMode.MINIMAL:
            step_config_block = yaml.load(step_config_yaml_str)
            if step_config_block:
                step_config_block = remove_non_required(step_config_block)
        else:
            indented_lines = [
                "    " + line if line.strip() else "" for line in step_config_yaml_str.splitlines()
            ]
            step_config_block = yaml.load("\n".join(indented_lines))

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
@click.argument(
    "directory",
    default=_get_cwd,
    type=click.Path(),
    required=False,
    metavar="[DIRECTORY]",
)
@click.option(
    "--project-name",
    type=str,
    help="A string to use for the project name, used in the README file",
)
@click.option(
    "--task",
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
@click.option(
    "--config-mode",
    type=click.Choice(["commented", "full", "minimal"]),
    default="commented",
    help="Style for task config: 'commented' (defaults commented out), 'full' (all defaults active), 'minimal' (strip default options).",
)
def init(directory, project_name, tasks, manage_config, config_mode):
    """Initialize a new snappy project directory."""
    project_directory = directory() if callable(directory) else directory
    log("SNAPPY Pipeline -- init")
    log("================================")
    log("")

    # Check if config file already exists - do not overwrite existing project
    config_dest_path = os.path.join(project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
    if os.path.exists(config_dest_path):
        log(
            "Project already initialized at {path} ({config} exists)",
            {"path": project_directory, "config": CONFIG_FILENAME},
            level=LVL_ERROR,
        )
        sys.exit(1)

    # Create project directory and subdirectory for configuration files
    paths = [project_directory]
    if CONFIG_SUBDIR:
        paths.append(os.path.join(project_directory, CONFIG_SUBDIR))
    for path in paths:
        create_directory(path, exist_ok=True)

    # Create config file in subdirectory based on template
    config_dest_path = os.path.join(project_directory, CONFIG_SUBDIR, CONFIG_FILENAME)
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", "project_config.yaml"),
        dest_path=config_dest_path,
        format_args={
            "created_at": datetime.datetime.now().isoformat(),
            "project_name": (project_name or os.path.basename(os.path.abspath(project_directory))),
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
            "project_name": (project_name or os.path.basename(os.path.abspath(project_directory))),
        },
        message="Creating README file in in {path}",
        message_args={"path": os.path.join(project_directory, README_FILENAME)},
    )

    # Create additional tasks if any was provided
    for task_name, step in tasks:
        app = AddTaskApp(
            step=step,
            task_name=task_name,
            project_directory=project_directory,
            manage_config=manage_config,
            config_mode=config_mode,
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
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.option(
    "--manage-config/--no-manage-config",
    default=True,
    help="Manage the project config.yaml file automatically",
)
@click.option(
    "--config-mode",
    type=click.Choice(["commented", "full", "minimal"]),
    default="commented",
    help="Style for task config: 'commented' (defaults commented out), 'full' (all defaults active), 'minimal' (strip default options).",
)
def task_add(
    tasks,
    directory,
    manage_config,
    config_mode,
):
    """Add tasks to the project configuration.

    TASKS is a list of tasks to add, formatted as [task_name=]step_type.
    """
    project_directory = directory() if callable(directory) else directory

    for task_name, step in tasks:
        app = AddTaskApp(
            step=step,
            task_name=task_name,
            project_directory=project_directory,
            manage_config=manage_config,
            config_mode=config_mode,
        )
        res = app.run()
        if res:
            sys.exit(res)


@task.command(name="list")
@click.option(
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.option(
    "--available",
    is_flag=True,
    default=False,
    help="List all available registered workflow steps in SNAPPY instead of configured tasks.",
)
def task_list(directory, available):
    """List defined tasks in the project (or available steps with --available)."""
    if available:
        _print_steps_table()
        return

    project_directory = directory() if callable(directory) else directory
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
    for task_item in tasks:
        if isinstance(task_item, dict):
            name = task_item.get("name", "unnamed")
            step = task_item.get("step", "unknown")
            click.echo(f"  - {name} ({step})")


@task.command(name="steps")
def task_steps():
    """List all available workflow steps registered in SNAPPY."""
    _print_steps_table()


@main.command(
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.option(
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Path to directory to run in, default is current working directory",
)
@click.option(
    "--task",
    "task_name",
    type=str,
    default=None,
    shell_complete=_complete_task_names,
    help="The specific task name from config.yaml to run",
)
@click.option(
    "--frozen",
    is_flag=True,
    default=False,
    help=(
        "With --task: run only that task's rules. Outputs of upstream tasks must already exist; "
        "they are used as input files and never rebuilt."
    ),
)
@click.option("--verbose", is_flag=True, help="Increase verbosity level")
@click.pass_context
def run(ctx, directory, task_name, frozen, verbose):
    """Run snappy pipeline workflows.

    The executor and its settings, such as SLURM, come from a Snakemake profile of the system or
    the user (``--profile`` after ``--``, or ``SNAKEMAKE_PROFILE``)."""
    if frozen and not task_name:
        raise click.UsageError("--frozen requires --task")
    directory_path = directory() if callable(directory) else directory
    setup_logging(verbose)
    if task_name:
        logging.info("Targeting single task: %s", task_name)
    else:
        logging.info("No specific --task provided. Targeting all tasks (default).")

    snakemake_args = list(ctx.args)

    config_args = []
    if task_name:
        config_args.append(f"task={task_name}")
    if frozen:
        config_args.append("frozen=True")

    snakemake_argv = _snakemake_argv(directory_path, tuple(config_args))

    # Append all user-provided snakemake arguments directly
    snakemake_argv += snakemake_args

    logging.info("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    res = snakemake_main(snakemake_argv)
    if res != 0:
        ctx.exit(res)


@main.command()
@click.option(
    "--directory",
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
    directory_path = directory() if callable(directory) else directory
    path = db_path or _snkmt_db_path(directory_path)
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

    cmd = _snkmt_console_cmd(path)
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
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.argument("snakemake_log", required=False, type=click.Path(exists=True, dir_okay=False))
def logs(directory, snakemake_log):
    """Archive the logs of a run in logs/<run start>.tar.gz.

    snappy run does this when it ends. Use this command for a run whose Snakemake process was
    killed (walltime, scancel), so that its hooks did not run. Without SNAKEMAKE_LOG, the
    newest log in <directory>/.snakemake/log is used.
    """
    directory_path = directory() if callable(directory) else directory
    snakemake_log = snakemake_log or latest_snakemake_log(directory_path)
    if snakemake_log is None:
        log("No Snakemake log found below {path}", {"path": directory_path}, level=LVL_ERROR)
        sys.exit(1)
    archive = build_archive(directory_path, snakemake_log)
    log("Logs archived in {path}", {"path": str(archive)}, level=LVL_IMPORTANT)


@main.command(name="check-files")
@click.option(
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
def check_files(directory):
    """List the files that config.yaml names but that do not exist.

    These are the files in tool sections, the external_data paths, the sample sheets and search
    paths, the .fai and .dict of a genome and the index of a VCF or of the mappability BED file.
    Snakemake checks the inputs of the jobs it runs anyway; this command checks all of them, for
    example after moving a project to another machine. It exits with 1 if a file is missing.
    """
    directory_path = directory() if callable(directory) else directory
    config_path = os.path.abspath(os.path.join(directory_path, CONFIG_SUBDIR, CONFIG_FILENAME))
    project = load_project(load_configfile(config_path), os.path.dirname(config_path))
    missing = missing_files(project)
    for name, paths in missing.items():
        for path in paths:
            click.echo(f"{name}\t{path}")
    if missing:
        sys.exit(1)
    log("All files of {path} exist", {"path": config_path}, level=LVL_SUCCESS)


def _snakemake_argv(directory_path: str, config: tuple = ()) -> list[str]:
    """Return the Snakemake arguments that snappy run and snappy report share."""
    apps = os.path.dirname(os.path.abspath(__file__))
    argv = [
        "--directory",
        directory_path,
        # The project directory of relative paths is that of the config file, see the Snakefile
        "--configfile",
        os.path.abspath(os.path.join(directory_path, CONFIG_SUBDIR, CONFIG_FILENAME)),
        # The orchestrator Snakefile, which loads the tasks of config.yaml
        "--snakefile",
        os.path.join(os.path.dirname(apps), "Snakefile"),
    ]
    if config:
        argv += ["--config", *config]
    # Conda, the logger and the database persistence backend; the executor comes from a profile
    # of the system or the user
    argv += ["--workflow-profile", os.path.join(apps, "profile")]
    # Without an explicit path, the snkmt logger writes to one database per user, shared by all
    # projects. The path must be absolute because the logger resolves it against the process
    # working directory, not against --directory.
    return argv + ["--logger-snkmt-db", _snkmt_db_path(directory_path)]


@main.command(
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@click.option(
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)
@click.option(
    "--output",
    default="report.zip",
    show_default=True,
    help="Report file below the project directory. A .zip report also holds HTML directories "
    "such as the FastQC reports.",
)
@click.pass_context
def report(ctx, directory, output):
    """Build the Snakemake report of the QC outputs of the project's tasks.

    Arguments after -- go to Snakemake, as with snappy run.
    """
    directory_path = directory() if callable(directory) else directory
    target = os.path.join(os.path.abspath(directory_path), output)
    argv = _snakemake_argv(directory_path) + ["--report", target] + list(ctx.args)
    logging.info("Executing snakemake %s", " ".join(map(repr, argv)))
    res = snakemake_main(argv)
    if res != 0:
        ctx.exit(res)


def _snkmt_db_path(directory_path: str) -> str:
    """Return the absolute path of the project's snkmt database."""
    return os.path.join(os.path.abspath(directory_path), ".snakemake", "log", "snkmt.sqlite")


def _snkmt_console_cmd(db_path: str) -> list[str]:
    """Return the command line to launch the snkmt console for ``db_path``.

    ``python -m snkmt.cli`` silently does nothing because the package has no
    ``__main__`` guard; the ``snkmt`` console script is the proper entry point
    (``snkmt = snkmt.cli:main``).  Prefer it on PATH, and fall back to calling
    ``main()`` directly with the current interpreter.
    """
    snkmt_bin = shutil.which("snkmt")
    if snkmt_bin:
        return [snkmt_bin, "console", "--db-path", db_path]
    return [
        sys.executable,
        "-c",
        "from snkmt.cli import main; main()",
        "console",
        "--db-path",
        db_path,
    ]


@main.command(name="pull-sheet")
@click.option(
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


@main.command()
@click.argument("shell", type=click.Choice(["bash", "zsh", "fish"]), required=False)
def completion(shell):
    """Generate shell autocompletion script.

    Usage:
        eval "$(snappy completion zsh)"   # for zsh
        eval "$(snappy completion bash)"  # for bash
        eval "$(snappy completion fish)"  # for fish
    """
    if not shell:
        shell_path = os.environ.get("SHELL", "")
        if "zsh" in shell_path:
            shell = "zsh"
        elif "fish" in shell_path:
            shell = "fish"
        else:
            shell = "bash"

    from click.shell_completion import get_completion_class

    comp_cls = get_completion_class(shell)
    if comp_cls is None:
        click.echo(f"Error: Unsupported shell '{shell}'", err=True)
        sys.exit(1)

    comp = comp_cls(
        cli=main,
        ctx_args={},
        prog_name="snappy",
        complete_var="_SNAPPY_COMPLETE",
    )
    click.echo(comp.source())


if __name__ == "__main__":
    main()
