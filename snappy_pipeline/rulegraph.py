"""The rule graph of a run as a metro map for the report, drawn by snakevision.

The start hook of ``snappy run`` draws it in the main process, the only place that holds the
run's DAG (a job's ``run:`` block executes in a spawned Snakemake process with a DAG of that job
alone), and the report's global description embeds the SVG on its Workflow panel.
"""

import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Iterable

from snappy_pipeline.rule_names import split_task

#: The SVG below the project directory
RULEGRAPH_SVG = "report/rulegraph.svg"
#: The config key through which the main Snakefile tells the report description the SVG's path
RULEGRAPH_CONFIG_KEY = "__rulegraph__"

#: The width of a label character in the drawing's units (the labels have the viewer's default
#: font size, 16 units)
CHAR_WIDTH = 9.0
SVG_NS = "http://www.w3.org/2000/svg"


def rule_label(rule: str) -> str:
    """``<task>: <part> <action>`` for a task's rule ``<task>__<part>_<action>``, else ``rule``."""
    split = split_task(rule)
    if split is None:
        return rule
    task, rest = split
    return f"{task}: {rest.replace('_', ' ')}"


def rule_graph_dot(dag, skip: Iterable[str] = ()) -> str:
    """The rule graph of Snakemake's ``dag`` as DOT text, as ``snakemake --rulegraph`` prints it.

    One node per rule with a job, labelled by ``rule_label``, and one edge from a rule to each
    rule that reads its outputs. The rules in ``skip`` are left out.
    """
    skipped = set(skip)
    graph: dict[str, set[str]] = defaultdict(set)
    for job in dag.jobs:
        if job.rule.name in skipped:
            continue
        graph[job.rule.name].update(
            dep.rule.name for dep in dag.dependencies[job] if dep.rule.name not in skipped
        )
    ids = {name: i for i, name in enumerate(sorted(graph))}
    labels = {name: rule_label(name) for name in ids}
    nodes = [f'\t{i}[label = "{labels[name]}", style="rounded"];' for name, i in ids.items()]
    edges = [f"\t{ids[dep]} -> {ids[name]}" for name in ids for dep in sorted(graph[name])]
    return "digraph snakemake_dag {\n" + "\n".join(nodes + edges) + "\n}\n"


def write_rulegraph(dot: str, path: str) -> None:
    """Draw the DOT rule graph ``dot`` as the SVG ``path``, wide enough for its longest label."""
    # Imported here: snakevision pulls in networkx, which the CLI does not need otherwise
    from snakevision.snakevision import SnakeVision

    SnakeVision(dot.splitlines()).write(path)
    with open(path, encoding="utf-8") as f:
        svg = f.read()
    with open(path, "w", encoding="utf-8") as f:
        f.write(fit_labels(svg))


def draw_rulegraph(dag, skip: Iterable[str] = ()) -> bool:
    """Draw the rule graph of the run's ``dag`` into ``RULEGRAPH_SVG`` and return whether it did.

    Only a run that targets the whole workflow draws it, that is, one whose target jobs include
    the final target ``all``: a run with file targets or ``--until`` builds a partial DAG, which
    would draw a part of the workflow.
    """
    if not any(job.rule.name == "all" for job in dag.targetjobs):
        return False
    write_rulegraph(rule_graph_dot(dag, skip=skip), RULEGRAPH_SVG)
    return True


def fit_labels(svg: str) -> str:
    """Widen the viewBox of the SVG ``svg`` to the end of its longest label and size it to that.

    dagviz sizes the drawing by its nodes only, so a long label would be cut off, and it sets
    the width and height to 100 %, which leaves an ``<img>`` without a size.
    """
    ET.register_namespace("", SVG_NS)
    root = ET.fromstring(svg)
    ends = [
        float(text.get("x", "0")) + CHAR_WIDTH * len(text.text or "")
        for text in root.iter(f"{{{SVG_NS}}}text")
    ]
    if not ends:
        return svg
    left, top, width, height = (float(v) for v in root.get("viewBox", "").replace(",", " ").split())
    width = max(width, max(ends) - left)
    root.set("viewBox", f"{left} {top} {width} {height}")
    # The drawing's own size instead of 100 %, so that an <img> in the report has a size to show
    root.set("width", f"{width:g}")
    root.set("height", f"{height:g}")
    return ET.tostring(root, encoding="unicode")
