"""Helpers of the fgumi wrappers of ``umi_mapping`` (``fgumi_grouping``, ``fgumi_consensus``):
the two wrappers share them, so their shell texts live in a module each, and the params that one
wrapper reads stay apart from the other's.

Compression levels are structural: 0 inside pipes, 1 for the files of the job's temp directory
and for the grouped BAM, the configured level for the final BAM.
"""

import re


def fill(template: str, **values) -> str:
    """``template`` with the placeholders of ``values`` filled in as literal text: ``ShellWrapper``
    formats the result once more, so their braces are doubled, while the other placeholders and
    the doubled braces of the template stay for it."""
    return re.sub(
        r"\{(" + "|".join(map(re.escape, values)) + r")\}",
        lambda match: str(values[match.group(1)]).replace("{", "{{").replace("}", "}}"),
        template,
    )


def flag(value) -> str:
    """A boolean as fgumi reads it."""
    return "true" if value else "false"


def index_prefix(snakemake) -> str:
    """The prefix of the bwa-mem3 index, whose files are inputs with the .amb file first."""
    return snakemake.input.index[0].removesuffix(".amb")
