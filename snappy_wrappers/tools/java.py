"""The Java heap of a wrapper's JVMs."""

# The wrappers run with the Python of their environment, so this stays valid Python 3.7.

from math import ceil


def java_heap_mb(snakemake, shards=1) -> int:
    """75% of the job's memory (``resources.mem_mb``) as Java heap in MB.

    ``shards`` JVMs that run at the same time get an equal part each.
    """
    return ceil(snakemake.resources.mem_mb * 0.75 / shards)
