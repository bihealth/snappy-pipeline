# -*- coding: utf-8 -*-
"""Helper code for file system (directory/file) manipulation

The routines in this package provides helpers for creating files from templates, for example,
and creating directories and files while printing log messages.
"""

import os

from .logging import LVL_INFO, log


def create_directory(path, msg_lvl=LVL_INFO, exist_ok=False):
    """Create directory ``path``

    Switch off messaging by setting msg_lvl to ``None``
    """
    if msg_lvl:
        log("Creating directory {path}", {"path": path}, level=msg_lvl)
    os.makedirs(path, exist_ok=exist_ok)


def create_from_tpl(src_path, dest_path, format_args, message, message_args, msg_lvl=LVL_INFO):
    """Create file at ``dest_path`` from template at ``src_path``

    Put in format arguments from format_args (properly escape ``%``).

    Disable printing by setting ``msg_lvl`` to ``None``
    """
    if msg_lvl:
        log(message, message_args, level=msg_lvl)
    # Read in the template and fill in values
    with open(src_path, "rt") as f:
        contents = f.read()
    formatted = contents % format_args
    with open(dest_path, "wt") as f:
        f.write(formatted)


def update_file(path, contents):
    """Write ``contents`` to the file at ``path``; the caller tells what changed."""
    with open(path, "wt") as f:
        f.write(contents)
