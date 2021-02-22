# -*- coding: utf-8 -*-
"""Helper code for logging"""

import sys

from termcolor import colored


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
    """
    args = args or {}
    if level == LVL_IMPORTANT:
        print(colored(msg.format(**args), "yellow"), file=file)
    else:
        if level == LVL_ERROR:
            prefix = colored("ERROR: ", "red", attrs=["bold"])
        elif level == LVL_INFO:
            prefix = colored("INFO: ", "yellow", attrs=["bold"])
        elif level == LVL_SUCCESS:
            prefix = colored("SUCCESS: ", "green", attrs=["bold"])
        else:
            prefix = ""
        print(prefix, msg.format(**args), sep="", file=file)
