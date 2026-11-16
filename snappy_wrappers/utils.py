# -*- coding: utf-8 -*-
"""Utility code for snappy_wrappers"""

import os

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def basename(fname, strip_suffix=""):
    """Unix `basename`-like function"""
    fname = os.path.basename(fname)
    if strip_suffix and fname.endswith(strip_suffix):
        return fname[: -len(strip_suffix)]
    else:
        return fname


def listify(gen):
    """Decorator that converts a generator into a function which returns a list

    Use it in the case where a generator is easier to write but you want
    to enforce returning a list::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield i
    """

    def patched(*args, **kwargs):
        """Wrapper function"""
        return list(gen(*args, **kwargs))

    return patched


def dictify(gen):
    """Decorator that converts a generator into a function which returns a dict

    Use it in the case where a generator is easier to write but you want
    to enforce returning a dict::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield 'key{}'.format(i), i
    """

    def patched(*args, **kwargs):
        """Wrapper function"""
        return dict(gen(*args, **kwargs))

    return patched


#: Options of ``mehari annotate seqvars`` that take their values comma-separated, in one occurrence
MEHARI_COMMA_SEPARATED = ("pick_order", "restrict_transcripts")


def mehari_options(config):
    """Return the command line options of ``mehari annotate seqvars`` for the fields of the
    ``Mehari`` config section (``snappy_pipeline.models.annotation``): snake_case fields become
    kebab-case options, booleans flags, ``plugins`` one ``--plugin`` per entry, the transcript
    lists comma-separated; unset fields (``None``, empty) leave mehari's default."""
    options = []
    for key, value in config.items():
        if value is None or value == "" or value == []:
            continue
        option = "--" + key.replace("_", "-")
        if isinstance(value, bool):
            if value:
                options.append(option)
        elif key in MEHARI_COMMA_SEPARATED:
            options.append("{} {}".format(option, ",".join(value)))
        elif key == "plugins":
            options.extend("--plugin {}".format(item) for item in value)
        elif isinstance(value, list):
            options.extend("{} {}".format(option, item) for item in value)
        else:
            options.append("{} {}".format(option, value))
    return options
