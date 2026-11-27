"""The rule graph of a run as a metro map for the report, drawn by snakevision."""

import os
import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Iterable, Mapping

from snakevision.snakevision import SnakeVision

#: The SVG below the project directory, a ``report()`` output of the main Snakefile
RULEGRAPH_SVG = "report/rulegraph.svg"
#: Its caption in the report
RULEGRAPH_CAPTION = os.path.join(os.path.dirname(__file__), "report", "rulegraph.rst")

#: The width of a label character in the drawing's units (the labels have the viewer's default
#: font size, 16 units)
CHAR_WIDTH = 9.0
SVG_NS = "http://www.w3.org/2000/svg"


def rule_label(rule: str, task_steps: Mapping[str, str]) -> str:
    """``<task>: <part> <action>`` for a task's rule ``<task>_<step>_<part>_<action>``, else ``rule``.

    ``task_steps`` maps the task names to their steps.
    """
    # A task name may be the prefix of another one's (filtration, filtration_dkfz)
    for task in sorted(task_steps, key=len, reverse=True):
        prefix = f"{task}_"
        if rule.startswith(prefix):
            rest = rule.removeprefix(prefix).removeprefix(f"{task_steps[task]}_")
            return f"{task}: {rest.replace('_', ' ')}"
    return rule


def rule_graph_dot(
    dag, task_steps: Mapping[str, str] | None = None, skip: Iterable[str] = ()
) -> str:
    """The rule graph of Snakemake's ``dag`` as DOT text, as ``snakemake --rulegraph`` prints it.

    One node per rule with a job, labelled by ``rule_label`` when ``task_steps`` is given, and
    one edge from a rule to each rule that reads its outputs. The rules in ``skip`` are left out.
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
    labels = {name: rule_label(name, task_steps) if task_steps else name for name in ids}
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
