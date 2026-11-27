"""The metro-map rule graph of a run for the report."""

import re
from pathlib import Path
from types import SimpleNamespace

from ruamel.yaml import YAML

from snappy_pipeline.rulegraph import fit_labels, rule_graph_dot, rule_label, write_rulegraph

#: Snakemake's ``--rulegraph`` output of a small run
DOT = """digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans, fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
\t0[label = "mapping: bwa run", color = "0.00 0.6 0.85", style="rounded"];
\t1[label = "calling: freebayes run", color = "0.17 0.6 0.85", style="rounded"];
\t2[label = "annotation: mehari run with a label that runs past the drawing", color = "0.33 0.6 0.85", style="rounded"];
\t0 -> 1
\t1 -> 2
}
"""


class Job:
    def __init__(self, rule: str) -> None:
        self.rule = SimpleNamespace(name=rule)


def test_rule_label_names_the_task_and_the_rest_of_the_rule():
    assert rule_label("mapping__bwa_run") == "mapping: bwa run"
    assert rule_label("mapping_qc__qc_summary_run") == "mapping_qc: qc summary run"
    assert rule_label("mapping__results") == "mapping: results"
    assert rule_label("all") == "all"


def test_rule_graph_dot_lists_the_rules_and_their_edges():
    bwa_a, bwa_b, calls, all_ = Job("bwa"), Job("bwa"), Job("calls__freebayes_run"), Job("all")
    dag = SimpleNamespace(
        jobs=[all_, calls, bwa_a, bwa_b],
        dependencies={all_: [calls], calls: [bwa_a, bwa_b], bwa_a: [], bwa_b: []},
    )
    assert rule_graph_dot(dag).splitlines() == [
        "digraph snakemake_dag {",
        '\t0[label = "all", style="rounded"];',
        '\t1[label = "bwa", style="rounded"];',
        '\t2[label = "calls: freebayes run", style="rounded"];',
        "\t2 -> 0",
        "\t1 -> 2",
        "}",
    ]
    assert rule_graph_dot(dag, skip=["all"]).splitlines() == [
        "digraph snakemake_dag {",
        '\t0[label = "bwa", style="rounded"];',
        '\t1[label = "calls: freebayes run", style="rounded"];',
        "\t0 -> 1",
        "}",
    ]


def test_write_rulegraph_draws_the_rules_and_fits_the_longest_label(tmp_path):
    svg = tmp_path / "report" / "rulegraph.svg"
    write_rulegraph(DOT, str(svg))
    text = svg.read_text()
    assert "<svg" in text
    assert "mapping: bwa run" in text
    assert "calling: freebayes run" in text
    left, _, width, _ = (float(v) for v in re.search(r'viewBox="([^"]*)"', text).group(1).split())
    x = max(float(v) for v in re.findall(r'<text[^>]*x="([0-9.]+)"', text))
    assert left + width >= x + 9 * len(
        "annotation: mehari run with a label that runs past the drawing"
    )


def test_fit_labels_keeps_a_drawing_that_is_wide_enough():
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0,0,500,100">'
        '<text x="10" y="10">short</text></svg>'
    )
    assert 'viewBox="0.0 0.0 500.0 100.0"' in fit_labels(svg)
    assert fit_labels('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0,0,1,1"/>').endswith("/>")


def test_the_rulegraph_job_draws_the_whole_run(tmp_path):
    """In all-tasks mode the main Snakefile draws the rule graph of the full DAG for the report.

    The DAG holds the finished jobs too, so a restarted run draws the whole workflow.
    """
    from snakemake.api import SnakemakeApi
    from snakemake.io import is_flagged
    from snakemake.settings.types import (
        ConfigSettings,
        DAGSettings,
        OutputSettings,
        ResourceSettings,
    )

    import snappy_pipeline

    root = Path(snappy_pipeline.__file__).resolve().parent.parent
    pipeline = root / ".tests" / "test-workflow" / "pipelines" / "snappy-cancer_wes"
    task_names = [t["name"] for t in YAML(typ="safe").load(pipeline / "config.yaml")["tasks"]]
    with SnakemakeApi(OutputSettings()) as api:
        workflow_api = api.workflow(
            resource_settings=ResourceSettings(cores=1),
            config_settings=ConfigSettings(),
            snakefile=Path(snappy_pipeline.__file__).parent / "Snakefile",
            workdir=pipeline,
        )
        dag_api = workflow_api.dag(dag_settings=DAGSettings(forceall=True))
        # No public API returns the DAG's jobs; this mirrors tests/scripts/dump_dag.py
        workflow = dag_api.workflow_api._workflow
        workflow._prepare_dag(forceall=True, ignore_incomplete=True, lock_warn_only=True)
        workflow._build_dag()
        jobs = [job for job in workflow.dag.jobs if job.rule.name == "rulegraph"]
        assert len(jobs) == 1
        assert is_flagged(jobs[0].output[0], "report")
        skip = ["all", "rulegraph", *(f"{name}__results" for name in task_names)]
        dot = rule_graph_dot(workflow.dag, skip=skip)

    labels = re.findall(r'label = "([^"]+)"', dot)
    assert len(labels) > 20
    assert "mapping: bwa run" in labels
    assert "report_final: datavzrd run" in labels
    assert not [label for label in labels if label.endswith("results")]
    svg = tmp_path / "rulegraph.svg"
    write_rulegraph(dot, str(svg))
    assert "mapping: bwa run" in svg.read_text()
