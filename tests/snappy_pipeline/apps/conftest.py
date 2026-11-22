# -*- coding: utf-8 -*-
"""Fixtures of the CLI tests: the prompts, scripted or driven by key presses."""

import pytest
from prompt_toolkit.input import create_pipe_input
from prompt_toolkit.output import DummyOutput

from snappy_pipeline.apps import snappy_cli
from snappy_pipeline.apps.impl import prompts


@pytest.fixture
def keys(monkeypatch):
    """Feed key presses to the real prompts: a pipe input in place of the terminal.
    Text and a filter end with Enter (``\\r``); a prefilled default is cleared first
    (``\\x15``); a confirmation takes one key, ``y`` or ``n``."""
    with create_pipe_input() as pipe:
        monkeypatch.setattr(prompts, "PROMPT_KWARGS", {"input": pipe, "output": DummyOutput()})
        yield pipe.send_text


class Script:
    """Scripted answers for the prompts of the CLI, in the order of the questions;
    ``asked`` records the questions. An answer of a selection must be one of its values."""

    def __init__(self):
        self.pending: list = []
        self.asked: list[str] = []

    def __call__(self, *answers):
        self.pending.extend(answers)

    def answer(self, message, *args, **kwargs):
        self.asked.append(message)
        assert self.pending, f"no answer scripted for {message!r}"
        answer = self.pending.pop(0)
        if args and isinstance(args[0], list):  # the choices of ask_select
            values = [value for value, _, _ in args[0]]
            assert answer in values, f"{answer!r} is no choice of {message!r}: {values}"
        return answer


@pytest.fixture
def answers(monkeypatch):
    """Script the prompts of the CLI (see :class:`Script`)."""
    script = Script()
    for name in ("ask_text", "ask_confirm", "ask_select"):
        monkeypatch.setattr(snappy_cli, name, script.answer)
    return script
