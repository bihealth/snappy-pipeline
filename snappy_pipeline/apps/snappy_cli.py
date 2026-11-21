# -*- coding: utf-8 -*-
"""Unified Click-based CLI for the snappy pipeline.

Replaces the previous multi-script CLI entry points with a single `snappy` command.
"""

import datetime
import importlib
import io
import json
import logging
import os
import shutil
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import rich_click as click
from click.core import ParameterSource
import ruamel.yaml as ruamel_yaml
from rich.console import Console
from rich.table import Table
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.util import load_yaml_guess_indent
from snakemake.cli import main as snakemake_main
from snakemake.common.configfile import load_configfile

from .. import __version__
from ..base import ConfigError
from ..log_archive import build_archive, latest_snakemake_log
from ..models import ConfigMode
from ..orchestration import create_task_instances, load_project, missing_files
from ..sheet_summary import missing_reads
from ..workflow_registry import WORKFLOW_REGISTRY
from .impl.fsmanip import (
    create_directory,
    create_from_tpl,
    update_file,
)
from .impl.logging import LVL_ERROR, LVL_IMPORTANT, LVL_INFO, LVL_SUCCESS, log

logger = logging.getLogger(__name__)

# Configure rich-click settings for formatting
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True


def _get_cwd():
    return os.getcwd()


#: The project directory of the commands that work on a project
directory_option = click.option(
    "--directory",
    type=click.Path(),
    default=_get_cwd,
    help="Project directory, defaults to current working directory",
)

#: The style of the task configurations that init and task add write
config_mode_option = click.option(
    "--config-mode",
    type=click.Choice([mode.value for mode in ConfigMode]),
    default=ConfigMode.COMMENTED.value,
    help="Style of the task configuration: 'commented' (every field with its description, the optional ones commented out), 'full' (every field active), 'minimal' (the required fields only, without descriptions).",
)


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

    directory = ctx.params.get("directory") or _get_cwd()

    task_names = set()
    try:
        with open(_project_config(directory), "rt") as f:
            yaml = ruamel_yaml.YAML()
            config_yaml = yaml.load(f.read())
            if config_yaml and "tasks" in config_yaml and isinstance(config_yaml["tasks"], list):
                for task_item in config_yaml["tasks"]:
                    if isinstance(task_item, dict) and "name" in task_item:
                        task_names.add(task_item["name"])
    except Exception:  # also without a config file
        pass
    if not task_names:
        task_names = set(STEPS)

    return [CompletionItem(name) for name in sorted(task_names) if name.startswith(incomplete)]


def _get_registered_steps_info():
    """Return (step name, title) of every step; the title is the first line of the step's module
    docstring, as on its docs page, else of its class docstring."""
    info = []
    for step_name in STEPS:
        cls = WORKFLOW_REGISTRY[step_name]
        doc = importlib.import_module(cls.__module__).__doc__ or cls.__doc__ or ""
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
GITIGNORE_FILENAME = ".gitignore"
SHEET_FILENAME = "samplesheet.tsv"
#: The sample sheet formats of ``snappy init --sheet-type`` and the data set type of each
SHEET_TYPES = {"generic": "generic", "germline": "germline_variants", "cancer": "matched_cancer"}

#: Configuration file name
CONFIG_FILENAME = "config.yaml"


class AddTaskApp:
    """Implementation of task addition logic."""

    def __init__(
        self,
        step: str,
        task_name: str,
        project_directory: str,
        config_mode: ConfigMode | str = ConfigMode.COMMENTED,
    ):
        self.step = step
        self.task_name = task_name
        self.project_directory = project_directory
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

        config_yaml = self._load_config_yaml()

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
        """Load ``config.yaml``, keeping its comments and remembering its indentation of
        sequences (``- name`` or ``  - name``), so that the file is written back in its style."""
        self.config_filename = _project_config(self.project_directory)
        with open(self.config_filename, "rt") as f:
            text = f.read()
        _, sequence, offset = load_yaml_guess_indent(text)
        # A file without sequences (``tasks: []``) gets the style of the docs
        self.sequence_indent = sequence if sequence is not None else 4
        self.sequence_offset = offset if offset is not None else 2
        return self._yaml().load(text)

    def _yaml(self):
        """A YAML round-tripper that keeps the file's quotes, line breaks and sequence style."""
        yaml = ruamel_yaml.YAML()
        yaml.preserve_quotes = True
        yaml.width = 4096
        yaml.indent(
            mapping=self.sequence_indent - self.sequence_offset,
            sequence=self.sequence_indent,
            offset=self.sequence_offset,
        )
        return yaml

    def _setup_configuration(self, config_yaml):
        if "tasks" not in config_yaml:
            config_yaml["tasks"] = []

        yaml = self._yaml()
        step_cls = WORKFLOW_REGISTRY[self.step]
        step_config_yaml_str = step_cls.config_model_class.default_config_yaml_string(
            self.config_mode
        )
        # ruamel keeps the column of a comment, so the block is loaded at the column where it
        # ends up: below ``config:`` of a task entry
        indent = " " * (self.sequence_indent + self.sequence_indent - self.sequence_offset)
        step_config_block = yaml.load(textwrap.indent(step_config_yaml_str, indent))

        if step_config_block is None:
            step_config_block = CommentedMap()

        task_block = CommentedMap()
        task_block["name"] = self.task_name
        task_block["step"] = self.step
        task_block["config"] = step_config_block

        config_yaml["tasks"].append(task_block)

        buf = io.StringIO()
        yaml.dump(config_yaml, stream=buf)
        buf.seek(0)
        contents = buf.read()
        update_file(
            path=self.config_filename,
            contents=contents,
            message="Updating project config with required default config for task {task_name} (step {step})",
            message_args={"task_name": self.task_name, "step": self.step},
        )


def setup_logging(verbose):
    """Print the records of the ``snappy_pipeline`` loggers to stderr.

    The root logger stays untouched: a handler on it would print every record of Snakemake a
    second time and the setup messages of the libraries that its plugins use (alembic).
    """
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("[snappy] %(levelname)s: %(message)s"))
    package_logger = logging.getLogger("snappy_pipeline")
    package_logger.addHandler(handler)
    package_logger.propagate = False
    package_logger.setLevel(logging.DEBUG if verbose else logging.INFO)


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
    help="The name of the project, in README.md, config.yaml and the sample sheet; default: the directory name",
)
@click.option(
    "--external-id",
    type=str,
    default="",
    help="The identifier of the project in an external data management system, such as the SODAR project UUID, which tracks where the raw data lives",
)
@click.option(
    "--description",
    type=str,
    default="",
    help="A description of the project, in README.md, config.yaml and the sample sheet",
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
    "--sheet-type",
    type=click.Choice(SHEET_TYPES),
    default="generic",
    show_default=True,
    help="The biomedsheets format of the sample sheet samplesheet.tsv, which is created if it does not exist and named in config.yaml",
)
@click.option(
    "--track-sheet/--no-track-sheet",
    default=False,
    help="Track samplesheet.tsv in git. Only if it holds no personal data, or the remote is approved for such data (such as the internal GitLab).",
)
@config_mode_option
@click.option(
    "--git/--no-git",
    default=True,
    help="Track the configuration in a git repository: .gitignore, config.yaml and README.md, nothing else unless --track-sheet.",
)
@click.option(
    "--remote",
    type=str,
    metavar="URL",
    help="Add URL as the remote 'origin' of the new repository; nothing is pushed.",
)
@click.option(
    "--interactive/--no-interactive",
    default=None,
    help="Ask for what no option gives; the default on a terminal.",
)
@click.pass_context
def init(
    ctx,
    directory,
    project_name,
    external_id,
    description,
    tasks,
    sheet_type,
    track_sheet,
    config_mode,
    git,
    remote,
    interactive,
):
    """Initialize a new snappy project directory.

    Creates config.yaml, README.md and the sample sheet samplesheet.tsv, and puts the
    configuration under git. On a terminal, asks for what no option gives.
    """
    project_directory = directory
    default_name = os.path.basename(os.path.abspath(project_directory))
    if interactive is None:
        interactive = sys.stdin.isatty() and sys.stdout.isatty()
    if interactive:

        def asked(name: str) -> bool:
            return ctx.get_parameter_source(name) == ParameterSource.DEFAULT

        if asked("project_name"):
            project_name = click.prompt("Project name", default=project_name or default_name)
        if asked("external_id"):
            external_id = click.prompt(
                "External id of the project, such as the SODAR project UUID",
                default="",
                show_default=False,
            )
        if asked("description"):
            description = click.prompt("Description", default="", show_default=False)
        if asked("sheet_type"):
            sheet_type = click.prompt(
                "Sample sheet type", type=click.Choice(SHEET_TYPES), default=sheet_type
            )
        if asked("tasks") and click.confirm("Add tasks now?", default=False):
            tasks = _prompt_tasks()
        if git:
            if asked("track_sheet"):
                track_sheet = click.confirm(
                    "Track samplesheet.tsv in git? Only if it holds no personal data, or the "
                    "remote is approved for such data (such as the internal GitLab)",
                    default=False,
                )
            if asked("remote"):
                remote = click.prompt(
                    "Git remote URL (origin; nothing is pushed)", default="", show_default=False
                )
            if track_sheet and remote:
                track_sheet = click.confirm(
                    f"Is {remote} approved for the data in the sample sheet?", default=False
                )
                if not track_sheet:
                    log("samplesheet.tsv stays untracked", level=LVL_IMPORTANT)
    project_name = project_name or default_name

    log("SNAPPY Pipeline -- init")
    log("================================")
    log("")

    # Check if config file already exists - do not overwrite existing project
    config_dest_path = os.path.join(project_directory, CONFIG_FILENAME)
    if os.path.exists(config_dest_path):
        log(
            "Project already initialized at {path} ({config} exists)",
            {"path": project_directory, "config": CONFIG_FILENAME},
            level=LVL_ERROR,
        )
        sys.exit(1)

    create_directory(project_directory, exist_ok=True)

    tpls = os.path.join(os.path.dirname(__file__), "tpls")
    format_args = {
        "created_at": datetime.datetime.now().isoformat(),
        "project_name": project_name,
        "external_id": external_id,
        "description": description,
        # YAML double-quoted scalars
        "project_name_yaml": json.dumps(project_name, ensure_ascii=False),
        "external_id_yaml": json.dumps(external_id, ensure_ascii=False),
        "description_yaml": json.dumps(description, ensure_ascii=False),
        "sheet_file": SHEET_FILENAME,
        "sheet_type": SHEET_TYPES[sheet_type],
    }

    # Create config file in subdirectory based on template; the minimal one has a short header
    # that points to the docs instead of explaining the sections
    template = "project_config.yaml"
    if ConfigMode(config_mode) is ConfigMode.MINIMAL:
        template = "project_config_minimal.yaml"
    create_from_tpl(
        src_path=os.path.join(tpls, template),
        dest_path=config_dest_path,
        format_args=format_args,
        message="Creating project-wide configuration in {path}",
        message_args={"path": config_dest_path},
    )

    # Create readme file in subdirectory based on template
    create_from_tpl(
        src_path=os.path.join(tpls, README_FILENAME),
        dest_path=os.path.join(project_directory, README_FILENAME),
        format_args={
            **format_args,
            "description": description or "TODO: describe your project here",
        },
        message="Creating README file in {path}",
        message_args={"path": os.path.join(project_directory, README_FILENAME)},
    )

    # The sample sheet, unless the project directory has one
    sheet_path = os.path.join(project_directory, SHEET_FILENAME)
    if os.path.exists(sheet_path):
        log("Keeping the existing sample sheet {path}", {"path": sheet_path}, level=LVL_INFO)
    else:
        create_from_tpl(
            src_path=os.path.join(tpls, f"samplesheet_{sheet_type}.tsv"),
            dest_path=sheet_path,
            format_args=format_args,
            message="Creating the {sheet_type} sample sheet {path}",
            message_args={"sheet_type": sheet_type, "path": sheet_path},
        )

    # Add the tasks; one that fails (a duplicate name) fails init once the others are added
    failed = 0
    for task_name, step in tasks:
        app = AddTaskApp(
            step=step,
            task_name=task_name,
            project_directory=project_directory,
            config_mode=config_mode,
        )
        failed = failed or app.run()

    if git:
        init_git(project_directory, remote or None, track_sheet)
    if failed:
        sys.exit(failed)

    log(
        "\nDo not forget to review config.yaml and to fill out README.md and the sample sheet!\n",
        level=LVL_IMPORTANT,
    )
    log("All done, have a nice day!", level=LVL_SUCCESS)


def _prompt_tasks() -> tuple:
    """Ask for tasks, ``[name=]step_type`` each, until an empty line."""
    click.echo("Step types: `snappy task list --available`")
    tasks = []
    while True:
        value = click.prompt(
            "Task ([name=]step_type, empty line to finish)", default="", show_default=False
        )
        if not value:
            return tuple(tasks)
        try:
            tasks.append(TaskParamType().convert(value, None, None))
        except click.BadParameter as e:
            click.echo(e.format_message())


def init_git(project_directory: str, remote: str | None, track_sheet: bool = False) -> None:
    """Put the configuration under version control: the first commit of a new repository with
    .gitignore, config.yaml and README.md (and the sample sheet with ``track_sheet``), or these
    files staged in the repository that already contains the directory (no nested repository,
    no remote). A remote is added without pushing. Without git or without a commit identity
    the files stay staged and the gap is reported; init does not fail."""
    gitignore = os.path.join(project_directory, GITIGNORE_FILENAME)
    create_from_tpl(
        src_path=os.path.join(os.path.dirname(__file__), "tpls", "gitignore"),
        dest_path=gitignore,
        format_args={},
        message="Creating {path}",
        message_args={"path": gitignore},
    )
    files = [GITIGNORE_FILENAME, CONFIG_FILENAME, README_FILENAME]
    if track_sheet:
        with open(gitignore, "at") as f:
            f.write(f"!/{SHEET_FILENAME}\n")
        files.append(SHEET_FILENAME)
    git = shutil.which("git")
    if git is None:
        log("git is not installed: the project is not under version control", level=LVL_IMPORTANT)
        return

    def run(*args, check=True):
        return subprocess.run(
            [git, "-C", project_directory, *args], capture_output=True, text=True, check=check
        )

    try:
        inside = run("rev-parse", "--is-inside-work-tree", check=False).returncode == 0
        if inside:
            log(
                "{path} is inside a git repository: staging the configuration there",
                {"path": project_directory},
                level=LVL_INFO,
            )
        else:
            run("init", "-q")
        run("add", "--", *files)
        if inside:
            if remote:
                log("The remote is not added to an existing repository", level=LVL_IMPORTANT)
            return
        commit = run("commit", "-q", "-m", "Initialize the snappy project", check=False)
        if commit.returncode != 0:
            log(
                "The configuration is staged, the commit failed: {error}",
                {"error": commit.stderr.strip()},
                level=LVL_IMPORTANT,
            )
        if remote:
            run("remote", "add", "origin", remote)
    except subprocess.CalledProcessError as e:
        log(
            "git failed ({command}): {error}",
            {"command": " ".join(e.cmd[3:]), "error": e.stderr.strip()},
            level=LVL_ERROR,
        )


@main.group()
def task():
    """Manage tasks within a snappy project."""
    pass


@task.command(name="add")
@click.argument(
    "tasks", nargs=-1, required=True, type=TaskParamType(), metavar="[NAME=]STEP_TYPE..."
)
@directory_option
@config_mode_option
def task_add(tasks, directory, config_mode):
    """Add tasks to the project configuration.

    TASKS is a list of tasks to add, formatted as [task_name=]step_type.
    """
    for task_name, step in tasks:
        app = AddTaskApp(
            step=step,
            task_name=task_name,
            project_directory=directory,
            config_mode=config_mode,
        )
        res = app.run()
        if res:
            sys.exit(res)


@task.command(name="list")
@directory_option
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

    config_filename = _project_config(directory)

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


@main.command(
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    )
)
@directory_option
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
@click.option(
    "--verbose", is_flag=True, help="Log the Snakemake command line and other debug messages"
)
@click.pass_context
def run(ctx, directory, task_name, frozen, verbose):
    """Run snappy pipeline workflows.

    The executor and its settings, such as SLURM, come from a Snakemake profile of the system or
    the user (``--profile`` after ``--``, or ``SNAKEMAKE_PROFILE``)."""
    if frozen and not task_name:
        raise click.UsageError("--frozen requires --task")
    setup_logging(verbose)

    snakemake_args = list(ctx.args)

    config_args = []
    if task_name:
        config_args.append(f"task={task_name}")
    if frozen:
        config_args.append("frozen=True")

    snakemake_argv = _snakemake_argv(directory, tuple(config_args))

    # Append all user-provided snakemake arguments directly
    snakemake_argv += snakemake_args

    logger.debug("Executing snakemake %s", " ".join(map(repr, snakemake_argv)))
    res = snakemake_main(snakemake_argv)
    if res != 0:
        ctx.exit(res)


@main.command()
@directory_option
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
    path = db_path or _snkmt_db_path(directory)
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
@directory_option
@click.argument("snakemake_log", required=False, type=click.Path(exists=True, dir_okay=False))
def logs(directory, snakemake_log):
    """Archive the logs of a run in logs/<run start>.tar.gz.

    snappy run does this when it ends. Use this command for a run whose Snakemake process was
    killed (walltime, scancel), so that its hooks did not run. Without SNAKEMAKE_LOG, the
    newest log in <directory>/.snakemake/log is used.
    """
    snakemake_log = snakemake_log or latest_snakemake_log(directory)
    if snakemake_log is None:
        log("No Snakemake log found below {path}", {"path": directory}, level=LVL_ERROR)
        sys.exit(1)
    archive = build_archive(directory, snakemake_log)
    log("Logs archived in {path}", {"path": str(archive)}, level=LVL_IMPORTANT)


@main.command(name="check-files")
@directory_option
def check_files(directory):
    """List the files that config.yaml names but that do not exist, and the libraries whose
    reads are not found.

    The files are those in tool sections, the external_data paths, the sample sheets and search
    paths, the .fai and .dict of a genome and the index of a VCF or of the mappability BED file.
    The reads are checked for every library of the tasks that read FASTQ files, below the search
    paths of its data set or external_data task. Snakemake checks the inputs of the jobs it runs
    anyway; this command checks all of them, for example after moving a project to another
    machine. It exits with 1 if a file is missing.
    """
    config_path = _project_config(directory)
    try:
        project = load_project(load_configfile(config_path), os.path.dirname(config_path))
        missing = missing_files(project)
        for name, paths in missing.items():
            for path in paths:
                click.echo(f"{name}\t{path}")
        # The task objects load the sample sheets and find the reads as the Snakefile does
        problems = missing_reads(create_task_instances(None, project))
    except ConfigError as e:
        log(str(e), level=LVL_ERROR)
        sys.exit(1)
    for problem in problems:
        click.echo("\t".join(problem))
    if missing or problems:
        sys.exit(1)
    log("All files of {path} exist", {"path": config_path}, level=LVL_SUCCESS)


def _project_config(directory_path: str) -> str:
    """Return the absolute path of the project's config file.
    Without one, fail with the usage and a hint instead of a traceback from Snakemake."""
    path = os.path.abspath(os.path.join(directory_path, CONFIG_FILENAME))
    if not os.path.isfile(path):
        raise click.UsageError(
            f"No {CONFIG_FILENAME} in {os.path.dirname(path)}. Change into the project directory, "
            "pass --directory <project>, or create a project with snappy init."
        )
    return path


def _snakemake_argv(directory_path: str, config: tuple = ()) -> list[str]:
    """Return the Snakemake arguments that snappy run and snappy report share."""
    apps = os.path.dirname(os.path.abspath(__file__))
    argv = [
        "--directory",
        directory_path,
        # The project directory of relative paths is that of the config file, see the Snakefile
        "--configfile",
        _project_config(directory_path),
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
@directory_option
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
    target = os.path.join(os.path.abspath(directory), output)
    argv = _snakemake_argv(directory) + ["--report", target] + list(ctx.args)
    setup_logging(verbose=False)
    logger.debug("Executing snakemake %s", " ".join(map(repr, argv)))
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
    from .snappy_pull_sheet import run as run_pull_sheet

    run_pull_sheet(
        SimpleNamespace(
            output=output,
            api_key=api_key,
            sodar_host=sodar_host,
            project_uuid=project_uuid,
            project_name=project_name,
            project_title=project_title,
            project_description=project_description,
            library_types=library_types.split(",") if library_types else [],
        )
    )


@main.command()
@click.argument("shell", type=click.Choice(["bash", "zsh", "fish"]), required=False)
def completion(shell):
    """Print the shell completion script.

    Load it in every shell:
        eval "$(snappy completion zsh)"   # zsh
        eval "$(snappy completion bash)"  # bash
        eval "$(snappy completion fish)"  # fish

    Or write it once to a path the shell loads on demand, without the eval at every start:
        snappy completion zsh > ~/.zfunc/_snappy   # fpath+=(~/.zfunc) before compinit
        snappy completion bash > ~/.local/share/bash-completion/completions/snappy
        snappy completion fish > ~/.config/fish/completions/snappy.fish
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
