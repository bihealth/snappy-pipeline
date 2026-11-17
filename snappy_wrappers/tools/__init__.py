# -*- coding: utf-8 -*-

# TODO: snappy_wrappers.tools is the wrong place for this

import contextlib
import sys


def open_file(path, mode: str):
    """Open ``path`` as a context manager; ``-`` or ``None`` is standard input or output, which
    stays open. The tools take paths on the command line and open them after parsing, as
    ``argparse.FileType`` did before Python 3.14 deprecated it."""
    if path in (None, "-"):
        return contextlib.nullcontext(sys.stdin if "r" in mode else sys.stdout)
    return open(path, mode)
