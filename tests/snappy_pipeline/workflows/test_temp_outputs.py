# -*- coding: utf-8 -*-
"""Tests for marking scatter fragments and intermediates ``temp()``."""

from types import SimpleNamespace

import pytest
from snakemake.io import is_flagged

from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart


class _Part(BaseStepPart):
    name = "caller"
    actions = ("run",)
    temp_outputs = {"run": ("chunk",)}
    intermediate_outputs = {"run": ("raw",)}
    bam_outputs = {"run": ("bam",)}

    def get_output_files(self, action):
        return {
            "chunk": "work/1-of-2.vcf",
            "raw": "work/raw.vcf",
            "vcf": "work/out.vcf",
            "bam": "work/out.bam",
        }


@pytest.mark.parametrize(
    "delete_intermediates, temp",
    [
        ([], ["chunk"]),
        (["intermediates"], ["chunk", "raw"]),
        (["bam"], ["chunk", "bam"]),
        (["intermediates", "bam"], ["chunk", "raw", "bam"]),
    ],
)
def test_fragments_are_temp_and_intermediates_on_request(delete_intermediates, temp):
    step = object.__new__(BaseStep)
    step.config = SimpleNamespace(delete_intermediates=delete_intermediates)
    step.task_name = "calling"
    part = object.__new__(_Part)
    step._get_sub_step = lambda name: part

    outputs = step.get_output_files("caller", "run")

    assert [key for key, path in outputs.items() if is_flagged(path, "temp")] == temp
