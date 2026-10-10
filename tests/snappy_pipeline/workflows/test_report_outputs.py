# -*- coding: utf-8 -*-
"""Tests for outputs flagged for the Snakemake report (plans.md F5)."""

import os

from click.testing import CliRunner
from snakemake.io import get_flag_value, is_flagged

from snappy_pipeline.apps import snappy_cli
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ReportOutput


class _PlotPart(BaseStepPart):
    name = "plot"
    actions = ("run",)
    report_outputs = {"run": {"plot": ReportOutput("Plots", caption="plots.rst")}}

    def get_output_files(self, action):
        return {"plot": "work/{library_name}/report/{library_name}.png", "table": "work/x.tsv"}


def _outputs():
    step = object.__new__(BaseStep)
    step.task_name = "qc"
    step._get_sub_step = lambda name: object.__new__(_PlotPart)
    return step.get_output_files("plot", "run")


def test_flagged_outputs_go_into_the_task_category_labelled_by_their_wildcards():
    outputs = _outputs()

    assert is_flagged(outputs["plot"], "report") and not is_flagged(outputs["table"], "report")
    flag = get_flag_value(outputs["plot"], "report")
    assert (flag.category, flag.subcategory) == ("qc", "Plots")
    assert flag.labels == {"library_name": "{library_name}"}
    assert flag.caption == os.path.join(os.path.dirname(__file__), "report", "plots.rst")


def test_snappy_report_writes_report_zip_in_the_project(mocker, tmp_path):
    snakemake_main = mocker.patch("snappy_pipeline.apps.snappy_cli.snakemake_main", return_value=0)
    result = CliRunner().invoke(snappy_cli.main, ["report", "--directory", str(tmp_path)])

    assert result.exit_code == 0, result.output
    argv = snakemake_main.call_args.args[0]
    assert argv[argv.index("--report") + 1] == str(tmp_path / "report.zip")


def _subclasses(cls):
    for sub in cls.__subclasses__():
        yield sub
        yield from _subclasses(sub)


def test_declared_captions_exist():
    import inspect

    import snappy_pipeline.workflow_registry  # noqa: F401 (imports every step module)

    missing = []
    for part in _subclasses(BaseStepPart):
        specs = part.__dict__.get("report_outputs")
        # A property depends on the config; test parts are not shipped
        if not isinstance(specs, dict) or not part.__module__.startswith("snappy_pipeline."):
            continue
        for outputs in specs.values():
            for spec in outputs.values():
                path = os.path.join(os.path.dirname(inspect.getfile(part)), "report", spec.caption)
                if spec.caption and not os.path.isfile(path):
                    missing.append(path)
    assert not missing
