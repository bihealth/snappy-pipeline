# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# Register a ``tsv`` lexer alias so ``.. code-block:: tsv`` sample-sheet blocks render.
from pygments.lexers._mapping import LEXERS
from pygments.lexers.special import TextLexer
import pygments.lexers as _pygments_lexers


class TsvLexer(TextLexer):
    name = "TSV"
    aliases = ["tsv"]
    filenames = ["*.tsv"]
    mimetypes = ["text/tab-separated-values"]


LEXERS["TsvLexer"] = (
    "pygments.lexers.special",
    "TSV",
    ("tsv",),
    ("*.tsv",),
    ("text/tab-separated-values",),
)
_pygments_lexers._lexer_cache["TSV"] = TsvLexer

project = "snappy-pipeline"
copyright = "2015-2024, CUBI, Berlin Institute of Health"
author = "CUBI, Berlin Institute of Health"

# Get the project root dir, which is the parent dir of this
import os, sys

cwd = os.getcwd()
project_root = os.path.dirname(cwd)

# Insert the project root dir as the first element in the PYTHONPATH.
# This lets us ensure that the source package is imported, and that its
# version is used.
sys.path.insert(0, project_root)

import snappy_pipeline

# The short X.Y version.
version = snappy_pipeline.__version__
# The full version, including alpha/beta/rc tags.
release = snappy_pipeline.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx.ext.graphviz",
]
graphviz_output_format = "svg"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"
myst_enable_extensions = ["deflist"]


man_pages = [("index", "snappy_pipeline", "SNAPPY Pipeline Documentation", ["Manuel Holtgrewe"], 1)]

# -- Options for LaTeX output ------------------------------------------

latex_documents = [
    (
        "index",
        "snappy_pipeline.tex",
        "SNAPPY Pipeline Documentation",
        "Manuel Holtgrewe",
        "manual",
    )
]

texinfo_documents = [
    (
        "index",
        "snappy_pipeline",
        "SNAPPY Pipeline Documentation",
        "Manuel Holtgrewe",
        "snappy_pipeline",
        "One line description of project.",
        "Miscellaneous",
    )
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
pygments_style = "sphinx"
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["graphviz.css"]
# Output file base name for HTML help builder.
htmlhelp_basename = "snappy_pipelinedoc"
html_last_updated_fmt = "%b %d, %Y"


# -- Step pages ----------------------------------------------------------------
# One page per registered step, generated into step/ (ignored by git). The title and the text
# come from the step's module docstring, whose first line is the title. The inputs table comes
# from the step's ``depends_on`` model, the default configuration from its config model.

import importlib
import re
import textwrap

from snappy_pipeline import step_overview
from snappy_pipeline.models import RequiredBy, reader_markers
from snappy_pipeline.orchestration import ON_HOLD_STEPS
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import TOOL, field_contract

STEP_DIR = "step"

INPUTS_INTRO = """\
Each key of ``depends_on`` names an upstream task.
That task must produce the listed data type with all listed tags (``a|b``: one of them, ``-a``: not ``a``).
A key is an error for a task that none of the readers in "Read by" applies to."""


def _inputs_section(model_cls) -> str:
    field = model_cls.model_fields.get("depends_on")
    if field is None:
        return "Inputs\n------\n\nThis step reads no other task.\n"
    rows = []
    for key, info in field.annotation.model_fields.items():
        signature = field_contract(info, TOOL).signature
        markers = reader_markers(info)
        read_by = "; ".join(str(marker) for marker in markers) or "all"
        required_by = [str(marker) for marker in markers if isinstance(marker, RequiredBy)]
        if info.is_required():
            required = "yes"
        elif required_by:
            required = "for " + "; ".join(required_by)
        else:
            required = "no"
        description = " ".join((info.description or "").split())
        row = f"   * - ``{key}``\n     - ``{signature}``\n     - {read_by}\n     - {required}"
        rows.append(f"{row}\n     - {description}".rstrip())
    header = "   * - Key\n     - Requires\n     - Read by\n     - Required\n     - Description"
    table = "\n".join([".. list-table::", "   :header-rows: 1", "", header, *rows])
    return f"Inputs\n------\n\n{INPUTS_INTRO}\n\n{table}\n"


def _step_page(name: str, workflow_cls) -> str:
    doc = importlib.import_module(workflow_cls.__module__).__doc__ or name
    title, _, body = textwrap.dedent(doc.strip()).partition("\n")
    body = textwrap.dedent(body).strip()
    # The inputs go before the docstring's first section
    first_section = re.search(r"^[^\n]+\n-{3,}$", body, re.MULTILINE)
    split = first_section.start() if first_section else len(body)
    config = textwrap.indent(workflow_cls.default_config_yaml().strip(), "   ")
    return (
        "\n\n".join(
            part
            for part in (
                f".. _step_{name}:",
                f"{title}\n{'=' * len(title)}",
                body[:split].strip(),
                _inputs_section(workflow_cls.config_model_class).strip(),
                body[split:].strip(),
                f"Default configuration\n---------------------\n\n.. code-block:: yaml\n\n{config}",
            )
            if part
        )
        + "\n"
    )


# -- Step overview -------------------------------------------------------------
# The data contracts of all steps as a graph and a table (step_overview.rst, ignored by git)

OVERVIEW_PAGE = "step_overview.rst"

OVERVIEW_INTRO = """\
Every step declares what it produces and what it requires from the tasks named in its ``depends_on`` (:ref:`dev_contracts`).
The graph shows the data types that flow between steps.
Inputs from resources (reference data, regions, databases, indexes and existing files, from ``external_data``, ``reference_download`` and ``reference_index`` tasks) are listed in each box, optional ones in parentheses.
The table lists every ``depends_on`` key with the data type and tags it requires; the page of a step has the details."""


def _overview_page() -> str:
    contracts = step_overview.step_contracts()
    graph = textwrap.indent(step_overview.dot(contracts), "   ")
    return (
        ".. _step_overview:\n\nStep overview\n=============\n\n"
        f"{OVERVIEW_INTRO}\n\n"
        ".. graphviz::\n   :alt: The steps and the data types they read and write\n\n"
        f"{graph}\n{step_overview.rst_table(contracts)}\n"
    )


# Write changed pages only, so that Sphinx rereads only those, and remove pages of old steps
pages = {
    os.path.join(STEP_DIR, f"{name}.rst"): _step_page(name, workflow_cls)
    for name, workflow_cls in WORKFLOW_REGISTRY.items()
    if name not in ON_HOLD_STEPS
}
os.makedirs(STEP_DIR, exist_ok=True)
for entry in os.scandir(STEP_DIR):
    if entry.path not in pages:
        os.remove(entry.path)
pages[OVERVIEW_PAGE] = _overview_page()
for path, page in pages.items():
    if not os.path.exists(path) or open(path).read() != page:
        with open(path, "wt") as f:
            f.write(page)


# -- Command table -------------------------------------------------------------
# The ``snappy-commands`` directive renders the commands of ``snappy`` with the short help of each,
# so that the docs list what the CLI has (usage.md)

from docutils import nodes
from docutils.parsers.rst import Directive

from snappy_pipeline.apps import snappy_cli


def _commands(group, prefix="snappy"):
    """Yield (command line, short help) of the commands below ``group`` in definition order."""
    for name, command in group.commands.items():
        if hasattr(command, "commands"):
            yield from _commands(command, f"{prefix} {name}")
        elif not command.hidden:
            yield f"{prefix} {name}", command.get_short_help_str(limit=1000)


def _help(text):
    """The help text with its options as literals, which Sphinx's smart quotes leave alone."""
    parts = re.split(r"(--[\w-]+)", text)
    return [nodes.literal(text=p) if i % 2 else nodes.Text(p) for i, p in enumerate(parts) if p]


def _row(*cells):
    row = nodes.row()
    for cell in cells:
        row += nodes.entry("", nodes.paragraph("", "", *cell))
    return row


class SnappyCommands(Directive):
    def run(self):
        # Sphinx reads the page again when the CLI changes
        self.state.document.settings.env.note_dependency(snappy_cli.__file__)
        tgroup = nodes.tgroup(cols=2)
        tgroup += [nodes.colspec(colwidth=1), nodes.colspec(colwidth=3)]
        tgroup += nodes.thead("", _row([nodes.Text("Command")], [nodes.Text("Purpose")]))
        tgroup += nodes.tbody(
            "", *(_row([nodes.literal(text=c)], _help(h)) for c, h in _commands(snappy_cli.main))
        )
        return [nodes.table("", tgroup)]


def setup(app):
    app.add_directive("snappy-commands", SnappyCommands)
