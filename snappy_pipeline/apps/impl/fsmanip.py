# -*- coding: utf-8 -*-
"""Helper code for file system (directory/file) manipulation

The routines in this package provides helpers for creating files from templates, for example,
and creating directories and files while printing log messages.
"""

import datetime
import difflib
import os
import shutil

from .logging import log, LVL_INFO, LVL_ERROR


def file_mtime(path):
    """Return mtime of file at ``path``"""
    t = datetime.datetime.fromtimestamp(os.stat(path).st_mtime)
    return t.isoformat()


def assume_path_nonexisting(path, msg_lvl=LVL_ERROR):
    """Check whether path does exist or not, print message if it does

    Return ``True`` if it does not exist and ``False`` otherwise.

    Switch off messaging by setting msg_lvl to ``None``
    """
    if os.path.exists(path):
        if msg_lvl:
            log("directory {path} already exists", {"path": path}, level=msg_lvl)
        return False
    else:
        return True


def assume_path_existing(path, msg_lvl=LVL_ERROR):
    """Check whether path does exist or not, print message if it does not.

    Return ``True`` if it does exist and ``False`` otherwise.

    Switch off messaging by setting msg_lvl to ``None``
    """
    if not os.path.exists(path):
        if msg_lvl:
            log("directory {path} does not exist yet", {"path": path}, level=msg_lvl)
        return False
    else:
        return True


def create_directory(path, msg_lvl=LVL_INFO, exist_ok=False):
    """Create directory ``path``

    Switch off messaging by setting msg_lvl to ``None``
    """
    if msg_lvl:
        log("creating directory {path}", {"path": path}, level=msg_lvl)
    os.makedirs(path, exist_ok=exist_ok)


def create_from_tpl(
    src_path, dest_path, format_args, message, message_args, msg_lvl=LVL_INFO, diff_context=3
):
    """Create file at ``dest_path`` from template at ``src_path``

    Put in format arguments from format_args (properly escape ``%``).

    Disable printing by setting ``msg_lvl`` to ``None``
    """
    if msg_lvl:
        log("creating project-wide configuration in {path}", message_args, level=msg_lvl)
    # Read in the template and fill in values
    with open(src_path, "rt") as f:
        contents = f.read()
    print("----\n%s\n----\n%s\n----\n" % (contents, format_args))
    formatted = contents % format_args
    # Create diff output and print
    if msg_lvl:
        diff = "".join(
            difflib.unified_diff(
                [],
                formatted.splitlines(True),
                "/dev/null",
                dest_path,
                file_mtime("/dev/null"),
                datetime.datetime.now().isoformat(),
                n=diff_context,
            )
        )
    log("applying the following change: \n\n{diff}", {"diff": diff}, level=LVL_INFO)
    with open(dest_path, "wt") as f:
        f.write(formatted)


def backup_file(path):
    """Create backup of file at the given path"""
    dest_path = path + ".bak"
    if os.path.exists(dest_path):
        i = 2
        dest_path = path + ".bak{}".format(i)
        while os.path.exists(dest_path):
            i += 1  # count up
            dest_path = path + ".bak{}".format(i)
    log("Creating backup {src} => {dest}", {"src": path, "dest": dest_path}, level=LVL_INFO)
    shutil.copyfile(path, dest_path)


def update_file(path, contents, message, message_args, msg_lvl=LVL_INFO, diff_context=3):
    """Update file at path with the given contents.

    Message is printed similar as in create_from_tpl.

    Disable printing by setting ``msg_lvl`` to ``None``
    """
    if msg_lvl:
        log(message, message_args, level=msg_lvl)
    # Read in the template and fill in values
    with open(path, "rt") as f:
        old_contents = f.read()
    slash = "/" if not os.path.isabs(path) else ""
    # Create diff output and print
    if msg_lvl:
        diff = "".join(
            difflib.unified_diff(
                old_contents.splitlines(True),
                contents.splitlines(True),
                "a" + slash + path,
                "b" + slash + path,
                file_mtime(path),
                datetime.datetime.now().isoformat(),
                n=diff_context,
            )
        )
    log("applying the following change: \n\n{diff}", {"diff": diff}, level=LVL_INFO)
    with open(path, "wt") as f:
        f.write(contents)
