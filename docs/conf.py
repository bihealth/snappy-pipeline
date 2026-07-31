# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

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
    "sphinx_mdinclude",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "step/DEFAULT_CONFIG_*.rst"]
source_suffix = ".rst"
master_doc = "index"


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


import textwrap

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

for name, workflow_cls in sorted(WORKFLOW_REGISTRY.items()):
    cfg = workflow_cls.default_config_yaml()
    if not cfg:
        continue  # swallow
    with open("step/DEFAULT_CONFIG_{}.rst".format(name), "wt") as outf:
        outf.write("::\n\n{}\n\n".format(textwrap.indent(cfg.strip(), "    ")))
