# -*- coding: utf-8 -*-
"""Tests of the prompts through a pipe input."""

from snappy_pipeline.apps.impl.prompts import ask_confirm, ask_select, ask_text

ENTER = "\r"
#: Clears a prefilled default (kill to the start of the line)
CLEAR = "\x15"

CHOICES = [("a", "alpha", "the first"), ("b", "beta", "the second")]


def test_text_select_and_confirm_through_a_pipe(keys):
    keys("My project" + ENTER)
    assert ask_text("Name") == "My project"
    keys(CLEAR + "mine" + ENTER)
    assert ask_text("Name", default="x") == "mine"
    keys(ENTER)
    assert ask_text("Name", default="x") == "x"
    keys("bet" + ENTER)  # a part of a title filters the list
    assert ask_select("Pick", CHOICES) == "b"
    keys(ENTER)
    assert ask_select("Pick", CHOICES, default="b") == "b"
    keys("n")
    assert ask_confirm("Sure?", default=True) is False
