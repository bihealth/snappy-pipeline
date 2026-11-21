# -*- coding: utf-8 -*-
"""Helper code for logging"""

import sys

import click

#: Message level: ERROR
LVL_ERROR = "ERROR"

#: Message level: INFO
LVL_INFO = "INFO"

#: Message level: IMPORTANT
LVL_IMPORTANT = "IMPORTANT"

#: Message level: SUCCESS
LVL_SUCCESS = "SUCCESS"


def log(msg, args=None, level=None, file=sys.stderr):
    """Print log message for given levels of importance

    For LVL_ERROR, LVL_INFO, LVL_SUCCESS, the message will be prefixed with a colored keyword
    identifying the level.  For IMPORTANT, the message itself will be colored.
    Colors are left out if ``file`` is not a terminal.
    """
    args = args or {}
    if level == LVL_IMPORTANT:
        click.echo(click.style(msg.format(**args), fg="yellow"), file=file)
    else:
        if level == LVL_ERROR:
            prefix = click.style("ERROR: ", fg="red", bold=True)
        elif level == LVL_INFO:
            prefix = click.style("INFO: ", fg="yellow", bold=True)
        elif level == LVL_SUCCESS:
            prefix = click.style("SUCCESS: ", fg="green", bold=True)
        else:
            prefix = ""
        click.echo(prefix + msg.format(**args), file=file)
