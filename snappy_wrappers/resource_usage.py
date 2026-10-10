# -*- coding: utf-8 -*-
"""Resource usage definition"""

import re
import typing
from dataclasses import dataclass


@dataclass
class ResourceUsage:
    """Resource usage specification to be used in ``BaseStepPart.default_resource_usage`` and
    ``BaseStepPart.resource_usage.values()``; as well as in the parallel wrappers classes.
    """

    threads: int
    runtime: str
    mem: str
    partition: typing.Optional[str] = None
    tmpdir: typing.Optional[str] = None
    #: Memory added per GB of input files, e.g. ``"2GB"``
    mem_per_gb_input: typing.Optional[str] = None
    #: Runtime added per GB of input files, e.g. ``"30m"``
    runtime_per_gb_input: typing.Optional[str] = None


_MEM_UNITS = {"": 1, "M": 1, "MB": 1, "G": 1024, "GB": 1024, "T": 1024**2, "TB": 1024**2}
_RUNTIME_UNITS = {"": 1, "m": 1, "h": 60, "d": 24 * 60}


def mem_mb(value: str) -> float:
    """Return a memory string such as ``"4GB"``, ``"16G"`` or ``"8000MB"`` in MB."""
    number, unit = re.fullmatch(r"\s*([0-9.]+)\s*([A-Za-z]*)\s*", str(value)).groups()
    return float(number) * _MEM_UNITS[unit.upper()]


def runtime_minutes(value: str) -> float:
    """Return a runtime string such as ``"4h"``, ``"2d"`` or ``"30m"`` in minutes."""
    number, unit = re.fullmatch(r"\s*([0-9.]+)\s*([a-z]*)\s*", str(value)).groups()
    return float(number) * _RUNTIME_UNITS[unit]
