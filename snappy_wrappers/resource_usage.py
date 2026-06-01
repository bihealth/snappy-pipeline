# -*- coding: utf-8 -*-
"""Resource usage definition"""

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
