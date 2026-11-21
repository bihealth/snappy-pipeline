# -*- coding: utf-8 -*-
"""Utility code for snappy_wrappers"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

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
