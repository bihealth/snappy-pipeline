# -*- coding: utf-8 -*-
"""Tests of the wrapper rules that the shared include generates (``BaseStep.wrapper_rules``)."""

from types import SimpleNamespace

from snappy_pipeline.workflows.abstract import BaseStep, RuleOptions
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType


class _Step(BaseStep):
    name = "example"
    produces = [DataSignature(DataType.QC)]


def _step(parts):
    step = object.__new__(_Step)
    step.step_name = "example"
    step.sub_steps = parts
    return step


def test_one_rule_per_wrapper_action_named_after_the_part_and_action():
    tool = SimpleNamespace(
        name="tool", actions=("prepare", "run", "write"), wrappers={"run": "x", "prepare": "y"}
    )
    tool.rule_options = {}
    helper = SimpleNamespace(name="helper", actions=("run",), wrappers={}, rule_options={})
    rules = _step({"helper": helper, "tool": tool}).wrapper_rules()

    assert [(r.name, r.part, r.action) for r in rules] == [
        ("tool_prepare", "tool", "prepare"),
        ("tool_run", "tool", "run"),
    ]
    assert rules[0].options == RuleOptions()


def test_explicit_rules_are_left_to_the_snakefile_and_cached_rules_are_separate():
    part = SimpleNamespace(
        name="tool",
        actions=("index", "run", "report"),
        wrappers={"index": "i", "run": "r", "report": "p"},
        rule_options={
            "index": RuleOptions(cache=True, retries=2),
            "report": RuleOptions(explicit=True),
            "run": RuleOptions(local=True, constraints={"x": "[^/]+"}),
        },
    )
    step = _step({"tool": part})

    plain = step.wrapper_rules()
    cached = step.wrapper_rules(cached=True)

    assert [r.action for r in plain] == ["run"]
    assert plain[0].options.local and plain[0].options.constraints == {"x": "[^/]+"}
    assert [(r.action, r.options.retries) for r in cached] == [("index", 2)]
