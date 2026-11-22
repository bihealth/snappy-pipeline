# -*- coding: utf-8 -*-
"""Unified Click-based CLI for the snappy pipeline.

Replaces the previous multi-script CLI entry points with a single `snappy` command.
"""

import contextlib
import datetime
import importlib
import io
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import textwrap
from types import SimpleNamespace
from typing import Any

import rich_click as click
from click.core import ParameterSource
import ruamel.yaml as ruamel_yaml
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.util import load_yaml_guess_indent
from snakemake.cli import main as snakemake_main
from snakemake.common.configfile import load_configfile

from .. import __version__
from ..base import ConfigError
from ..log_archive import SNKMT_DB, build_archive, latest_snakemake_log
from ..models import ConfigMode
from ..orchestration import CONFIG_FILE, create_task_instances, load_project, missing_files
from ..sheet_summary import missing_reads
from ..workflow_registry import WORKFLOW_REGISTRY
from .impl.fsmanip import (
    create_directory,
    create_from_tpl,
    update_file,
)
from .impl.logging import LVL_ERROR, LVL_IMPORTANT, LVL_INFO, LVL_SUCCESS, log
from .impl.output import Column, listing_options, print_listing
from .impl.prompts import ask_confirm, ask_select, ask_text
from .impl.task_chain import (
    EXTERNAL,
    Key,
    candidates,
    entries_of,
    entry_for,
    external_data_settings,
    providers,
    required_keys,
    step_tools,
)

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


#: The columns of the listings (:func:`print_listing`)
STEP_COLUMNS = [Column("step", "Step", no_wrap=True), Column("description", "Description")]
TASK_COLUMNS = [
    Column("name", "Task", no_wrap=True),
    Column("step", "Step"),
    Column("tool", "Tool"),
    Column("depends_on", "Depends on"),
]
CHECK_COLUMNS = [
    Column("kind", "Kind"),
    Column("task", "Task"),
    Column("item", "Missing"),
    Column("problem", "Problem"),
]
RULE_COLUMNS = [
    Column("rule", "Rule", no_wrap=True),
    Column("task", "Task"),
    Column("step", "Step"),
    Column("part", "Part"),
    Column("action", "Action"),
]


def _load_config_yaml(config_filename: str):
    with open(config_filename, "rt") as f:
        return ruamel_yaml.YAML().load(f.read())


def _task_rows_of(config_yaml) -> list[dict[str, Any]]:
    """The rows of the task listing from the loaded ``config.yaml``."""
    tasks = (config_yaml or {}).get("tasks") or []
    return [_task_row(item) for item in tasks if isinstance(item, dict)]


def _task_row(item: dict) -> dict[str, Any]:
    """A task entry of ``config.yaml`` as a row: name, step, tool and ``depends_on``."""
    config = item.get("config") or {}
    return {
        "name": item.get("name", ""),
        "step": item.get("step", ""),
        "tool": config.get("tool", ""),
        "depends_on": dict(config.get("depends_on") or {}),
    }


#: README file name
README_FILENAME = "README.md"
GITIGNORE_FILENAME = ".gitignore"
SHEET_FILENAME = "samplesheet.tsv"
#: The sample sheet formats of ``snappy init --sheet-type`` and the data set type of each
SHEET_TYPES = {"generic": "generic", "germline": "germline_variants", "cancer": "matched_cancer"}
#: The sheet types as prompt choices: value, title, description
SHEET_TYPE_CHOICES = [
    ("generic", "generic", "Libraries without relations between donors"),
    ("germline", "germline", "Pedigrees: index, parents and siblings, with sex and affection"),
    ("cancer", "cancer", "The tumor and normal libraries of each donor"),
]

#: Configuration file name
CONFIG_FILENAME = CONFIG_FILE


class AddTaskApp:
    """Implementation of task addition logic.

    ``settings`` are values set in the generated block, as the questions of a terminal give
    them: ``tool``, the ``depends_on`` keys, ``produces`` and the file keys of ``external_data``.
    """

    def __init__(
        self,
        step: str,
        task_name: str,
        project_directory: str,
        config_mode: ConfigMode | str = ConfigMode.COMMENTED,
        settings: dict[str, Any] | None = None,
    ):
        self.step = step
        self.task_name = task_name
        self.project_directory = project_directory
        self.config_mode = ConfigMode(config_mode)
        self.settings = settings or {}

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
        for key, value in self.settings.items():
            # A key of the block keeps its comments; a commented-out key is added anew
            if isinstance(value, dict) and isinstance(step_config_block.get(key), CommentedMap):
                for sub_key, sub_value in value.items():
                    step_config_block[key][sub_key] = sub_value
            else:
                step_config_block[key] = value

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
    add_tasks = False  # ask for the tasks once config.yaml exists
    if interactive:

        def asked(name: str) -> bool:
            return ctx.get_parameter_source(name) == ParameterSource.DEFAULT

        if asked("project_name"):
            project_name = ask_text("Project name", default=project_name or default_name)
        if asked("external_id"):
            external_id = ask_text("External id of the project, such as the SODAR project UUID")
        if asked("description"):
            description = ask_text("Description")
        if asked("sheet_type"):
            sheet_type = ask_select("Sample sheet type", SHEET_TYPE_CHOICES, default=sheet_type)
        if git:
            if asked("track_sheet"):
                track_sheet = ask_confirm(
                    "Track samplesheet.tsv in git? Only if it holds no personal data, or the "
                    "remote is approved for such data (such as the internal GitLab)",
                    default=False,
                )
            if asked("remote"):
                remote = ask_text("Git remote URL (origin; nothing is pushed)")
            if track_sheet and remote:
                track_sheet = ask_confirm(
                    f"Is {remote} approved for the data in the sample sheet?", default=False
                )
                if not track_sheet:
                    log("samplesheet.tsv stays untracked", level=LVL_IMPORTANT)
        add_tasks = asked("tasks") and ask_confirm("Add tasks now?", default=False)
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
    if add_tasks:
        _prompt_tasks(project_directory, config_mode)

    if git:
        init_git(project_directory, remote or None, track_sheet)
    if failed:
        sys.exit(failed)

    log(
        "\nDo not forget to review config.yaml and to fill out README.md and the sample sheet!\n",
        level=LVL_IMPORTANT,
    )
    log("All done, have a nice day!", level=LVL_SUCCESS)


#: The answers of a key question besides the tasks that can provide the key
ADD_PROVIDER = "+"
LEAVE_EMPTY = ""


class _TaskQuestions:
    """The task questions of a terminal: one task after the other, each written into
    ``config.yaml`` at once, with the tasks of the project as the answers to its ``depends_on``
    keys, or a new task added on the spot for a key nothing provides yet."""

    def __init__(self, project_directory: str, config_mode: ConfigMode | str):
        self.project_directory = project_directory
        self.config_mode = config_mode
        self.config_filename = _project_config(project_directory)
        self.entries = entries_of(_load_config_yaml(self.config_filename))
        self.summaries = dict(_get_registered_steps_info())
        self.added: list[str] = []

    def run(self) -> list[str]:
        """Ask for tasks until no more are wanted; print the task table; return the added names."""
        while True:
            self.task(ask_select("Step", self._step_choices(STEPS)))
            if not ask_confirm("Add another task?", default=True):
                break
        print_listing(_task_rows_of(_load_config_yaml(self.config_filename)), TASK_COLUMNS, "Tasks")
        return self.added

    def _step_choices(self, steps) -> list[tuple[str, str, str | None]]:
        return [(step, step, self.summaries.get(step)) for step in steps]

    def task(self, step: str, provides: Key | None = None) -> str:
        """Ask for one task of ``step`` and write it; ``provides`` is the key of the task that
        asked for it, which names the task by default and sets the data of ``external_data``."""
        taken = {entry.name for entry in self.entries}
        default = provides.name if provides else step
        name = ""
        while not name or name in taken:
            name = ask_text("Task name", default="" if default in taken else default).strip()
            if name in taken:
                log("A task named {name} exists", {"name": name}, level=LVL_ERROR)
        settings: dict[str, Any] = {}
        tools = step_tools(step)
        tool = ask_select("Tool", [(tool, tool, None) for tool in tools]) if tools else None
        if tool is not None:
            settings["tool"] = tool
        if step == EXTERNAL and provides is not None:
            settings |= external_data_settings(provides)
        depends_on: dict[str, str] = {}
        for key in required_keys(step, tool):
            choices = [(name_, name_, None) for name_ in candidates(key, self.entries)]
            choices.append((ADD_PROVIDER, "add a task that provides it now", None))
            choices.append((LEAVE_EMPTY, "leave empty, fill in config.yaml later", None))
            message = f"{key.name}: {key.description}" if key.description else key.name
            answer = ask_select(message, choices)
            if answer == ADD_PROVIDER:
                provider = ask_select("Step of the task", self._step_choices(providers(key)))
                answer = self.task(provider, provides=key)
            if answer:
                depends_on[key.name] = answer
        if depends_on:
            settings["depends_on"] = depends_on
        app = AddTaskApp(
            step=step,
            task_name=name,
            project_directory=self.project_directory,
            config_mode=self.config_mode,
            settings=settings,
        )
        if app.run():
            sys.exit(1)
        self.entries.append(entry_for(name, step, tool, settings))
        self.added.append(name)
        return name


def _prompt_tasks(project_directory: str, config_mode: ConfigMode | str) -> list[str]:
    """Ask for tasks on a terminal, one by one (see :class:`_TaskQuestions`)."""
    return _TaskQuestions(project_directory, config_mode).run()


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
    "tasks", nargs=-1, required=False, type=TaskParamType(), metavar="[NAME=]STEP_TYPE..."
)
@directory_option
@config_mode_option
@click.option(
    "--interactive/--no-interactive",
    default=None,
    help="Without tasks, ask for them one by one; the default on a terminal.",
)
def task_add(tasks, directory, config_mode, interactive):
    """Add tasks to the project configuration, or ask for them one by one on a terminal.

    TASKS is a list of tasks to add, formatted as [task_name=]step_type. Without tasks, each
    one is asked for: step, name, tool and, per required depends_on key, the task that provides
    it, among those of the project or added on the spot.
    """
    if not tasks:
        if interactive is None:
            interactive = sys.stdin.isatty() and sys.stdout.isatty()
        if not interactive:
            raise click.UsageError(
                "Give the tasks as [NAME=]STEP_TYPE, or run on a terminal to be asked for them"
            )
        _prompt_tasks(directory, config_mode)
        return
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
@listing_options
def task_list(directory, available, as_json, plain):
    """List the tasks of the project: name, step, tool and the tasks it depends on
    (with --available, the steps a task can use instead).
    """
    if available:
        rows = [
            {"step": name, "description": summary} for name, summary in _get_registered_steps_info()
        ]
        print_listing(rows, STEP_COLUMNS, "Steps", as_json=as_json, plain=plain)
        return
    config_filename = _project_config(directory)
    rows = _task_rows_of(_load_config_yaml(config_filename))
    print_listing(
        rows,
        TASK_COLUMNS,
        "Tasks",
        as_json=as_json,
        plain=plain,
        empty=f"No tasks in {config_filename}",
    )


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
@listing_options
def check_files(directory, as_json, plain):
    """List the files that config.yaml names but that do not exist, and the libraries whose
    reads are not found.

    The files are those in tool sections, the external_data paths, the sample sheets and search
    paths, the .fai and .dict of a genome and the index of a VCF or of the mappability BED file.
    The reads are checked for every library of the tasks that read FASTQ files, below the search
    paths of its data set or external_data task. Snakemake checks the inputs of the jobs it runs
    anyway; this command checks all of them, for example after moving a project to another
    machine. Each row has the kind (file or reads), the task, the file or the library, and the
    problem. It exits with 1 if a file is missing.
    """
    config_path = _project_config(directory)
    try:
        project = load_project(load_configfile(config_path), os.path.dirname(config_path))
        missing = missing_files(project)
        # The task objects load the sample sheets and find the reads as the Snakefile does
        problems = missing_reads(create_task_instances(None, project))
    except ConfigError as e:
        log(str(e), level=LVL_ERROR)
        sys.exit(1)
    rows = [
        {"kind": "file", "task": name, "item": path, "problem": "does not exist"}
        for name, paths in missing.items()
        for path in paths
    ]
    rows += [
        {"kind": "reads", "task": task, "item": library, "problem": problem}
        for task, library, problem in problems
    ]
    print_listing(rows, CHECK_COLUMNS, "Missing files and reads", as_json=as_json, plain=plain)
    if rows:
        sys.exit(1)
    log("All files of {path} exist", {"path": config_path}, level=LVL_SUCCESS)


def _rule_names(directory_path: str) -> list[str]:
    """The rule names of the project from Snakemake's rule listing, which prints to stdout like
    the library summary of the Snakefile; only the lines of one word (with a docstring in
    parentheses) are rule names."""
    argv = _snakemake_argv(directory_path, profile=False) + ["--list"]
    listing = io.StringIO()
    with contextlib.redirect_stdout(listing):
        try:
            snakemake_main(argv)  # ends with sys.exit, also after a listing
        except SystemExit as e:
            if e.code:
                sys.exit(e.code)
    names = []
    for line in listing.getvalue().splitlines():
        if match := re.fullmatch(r"(\w+)(?: \(.*\))?", line):
            names.append(match.group(1))
    return names


def _rule_rows(names: list[str], instances: dict[str, Any]) -> list[dict[str, str]]:
    """Rule name -> task, step, part and action: a rule is ``<task>_<step>_<part>_<action>`` for
    the jobs of a step part or ``<task>_results``; a rule outside the pattern (``all``) keeps the
    fields empty."""
    rows = []
    for name in names:
        row = {"rule": name, "task": "", "step": "", "part": "", "action": ""}
        task = max((t for t in instances if name.startswith(t + "_")), key=len, default=None)
        if task is not None:
            step = instances[task]
            row["task"], row["step"] = task, step.name
            rest = name.removeprefix(f"{task}_{step.name}_")
            for part in sorted(step.sub_steps.values(), key=lambda p: -len(p.name)):
                action = rest.removeprefix(part.name + "_")
                if action != rest and action in part.actions:
                    row["part"], row["action"] = part.name, action
                    break
        rows.append(row)
    return rows


@main.command()
@directory_option
@click.option(
    "--task",
    "task_names",
    multiple=True,
    shell_complete=_complete_task_names,
    help="Only the rules of this task (repeatable).",
)
@click.option(
    "--part",
    "parts",
    multiple=True,
    help="Only the rules of this step part, such as calls_report or qc_summary (repeatable).",
)
@click.option(
    "--names", is_flag=True, help="Only the rule names, one per line, for snakemake --forcerun."
)
@listing_options
def rules(directory, task_names, parts, names, as_json, plain):
    """List the rules of the project as Snakemake names them: <task>_<step>_<part>_<action> for
    the jobs of a step part, <task>_results for the job that links the results of a task.

    Rerun the report jobs of every task:
        snappy run -- --forcerun $(snappy rules --part calls_report --part qc_summary --names)
    """
    if names and (as_json or plain):
        raise click.UsageError("--names excludes --json and --plain")
    config_path = _project_config(directory)
    try:
        project = load_project(load_configfile(config_path), os.path.dirname(config_path))
        instances = create_task_instances(None, project)
    except ConfigError as e:
        log(str(e), level=LVL_ERROR)
        sys.exit(1)
    rows = _rule_rows(_rule_names(directory), instances)
    if task_names:
        rows = [row for row in rows if row["task"] in task_names]
    if parts:
        rows = [row for row in rows if row["part"] in parts]
    if names:
        for row in rows:
            click.echo(row["rule"])
        return
    print_listing(rows, RULE_COLUMNS, "Rules", as_json=as_json, plain=plain, empty="No rules")


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


def _snakemake_argv(directory_path: str, config: tuple = (), profile: bool = True) -> list[str]:
    """Return the Snakemake arguments that snappy run and snappy report share; without
    ``profile``, those of a listing (no workflow profile, so no logger database)."""
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
    if not profile:
        return argv
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
    return os.path.join(os.path.abspath(directory_path), SNKMT_DB)


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
