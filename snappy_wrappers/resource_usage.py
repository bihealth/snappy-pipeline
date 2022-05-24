# -*- coding: utf-8 -*-
"""Resource usage definition"""

import typing

import attr


@attr.s(frozen=True, auto_attribs=True)
class ResourceUsage:
    """Resource usage specification to be used in ``BaseStepPart.default_resource_usage`` and
    ``BaseStepPart.resource_usage.values()``; as well as in the parallel wrappers classes.
    """

    threads: int
    time: str
    memory: str
    partition: typing.Optional[str] = None
