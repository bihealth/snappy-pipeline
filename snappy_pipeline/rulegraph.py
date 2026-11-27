"""The rule graph of a run as a metro map for the report, drawn by snakevision."""

import os
import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Iterable

from snakevision.snakevision import SnakeVision

from snappy_pipeline.rule_names import split_task

#: The SVG below the project directory, a ``report()`` output of the main Snakefile
RULEGRAPH_SVG = "report/rulegraph.svg"
#: Its caption in the report
RULEGRAPH_CAPTION = os.path.join(os.path.dirname(__file__), "report", "rulegraph.rst")

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
    SnakeVision(dot.splitlines()).write(path)
    with open(path, encoding="utf-8") as f:
        svg = f.read()
    with open(path, "w", encoding="utf-8") as f:
        f.write(fit_labels(svg))


def fit_labels(svg: str) -> str:
    """Widen the viewBox of the SVG ``svg`` to the end of its longest label.

    dagviz sizes the drawing by its nodes only, so a long label would be cut off.
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
    root.set("viewBox", f"{left} {top} {max(width, max(ends) - left)} {height}")
    return ET.tostring(root, encoding="unicode")
