# -*- coding: utf-8 -*-
"""Tests for the BBDuk references of ``adapter_trimming``."""

import pydantic
import pytest

from snappy_pipeline.workflows.adapter_trimming.model import Bbduk


def test_bbduk_takes_builtin_references_without_files():
    config = Bbduk(builtin_references=["adapters", "phix"])

    assert config.adapter_sequences == []
    assert config.builtin_references == ["adapters", "phix"]


def test_bbduk_rejects_unknown_builtin_references():
    with pytest.raises(pydantic.ValidationError):
        Bbduk(builtin_references=["adaptors"])


def test_bbduk_needs_a_reference():
    with pytest.raises(pydantic.ValidationError, match="adapter_sequences or builtin_references"):
        Bbduk()
