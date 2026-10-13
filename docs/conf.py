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
]

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

from snappy_pipeline.orchestration import ON_HOLD_STEPS
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import DataSignature

STEP_DIR = "step"

INPUTS_INTRO = """\
Each key of ``depends_on`` names an upstream task.
That task must produce the listed data type with all listed tags (``a|b``: one of them, ``-a``: not ``a``)."""


def _inputs_section(model_cls) -> str:
    field = model_cls.model_fields.get("depends_on")
    if field is None:
        return "Inputs\n------\n\nThis step reads no other task.\n"
    rows = []
    for key, info in field.annotation.model_fields.items():
        signature = next(m for m in info.metadata if isinstance(m, DataSignature))
        required = "yes" if info.is_required() else "no"
        description = " ".join((info.description or "").split())
        row = f"   * - ``{key}``\n     - ``{signature}``\n     - {required}\n     - {description}"
        rows.append(row.rstrip())
    header = "   * - Key\n     - Requires\n     - Required\n     - Description"
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
    return "\n\n".join(
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
    ) + "\n"


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
for path, page in pages.items():
    if not os.path.exists(path) or open(path).read() != page:
        with open(path, "wt") as f:
            f.write(page)
