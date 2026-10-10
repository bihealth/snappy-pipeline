# -*- coding: utf-8 -*-
"""Tests for memory and runtime requests that grow with the input and with retries (F7)."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflow_model import Resources
from snappy_pipeline.workflows.abstract import BaseStepPart
from snappy_wrappers.resource_usage import ResourceUsage, mem_mb, runtime_minutes


@pytest.mark.parametrize(
    "value, expected", [("4GB", 4096), ("16G", 16384), ("8000MB", 8000), ("1.5 GB", 1536)]
)
def test_mem_mb(value, expected):
    assert mem_mb(value) == expected


@pytest.mark.parametrize("value, expected", [("4h", 240), ("2d", 2880), ("30m", 30), ("90", 90)])
def test_runtime_minutes(value, expected):
    assert runtime_minutes(value) == expected


def _resource(name, usage, attempt=1, input_mb=0, **resources):
    part = object.__new__(BaseStepPart)
    part.actions = ("run",)
    part.w_config = SimpleNamespace(resources=Resources(**resources))
    part.get_resource_usage = lambda action, **kwargs: usage
    return part.get_resource("run", name)(input=SimpleNamespace(size_mb=input_mb), attempt=attempt)


USAGE = ResourceUsage(threads=1, runtime="4h", mem="8GB")


def test_first_attempt_requests_the_declared_values():
    assert _resource("mem", USAGE) == "8GB"
    assert _resource("runtime", USAGE) == "4h"


def test_retries_grow_by_the_retry_factor():
    assert _resource("mem", USAGE, attempt=2) == "12288MB"
    assert _resource("runtime", USAGE, attempt=3) == "540m"
    assert _resource("mem", USAGE, attempt=2, retry_factor=1) == "8GB"


def test_input_term_adds_per_gb_of_input():
    usage = ResourceUsage(threads=1, runtime="1h", mem="2GB", mem_per_gb_input="1GB")
    assert _resource("mem", usage, input_mb=3 * 1024) == "5120MB"
