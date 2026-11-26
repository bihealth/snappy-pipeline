"""The data contracts of all steps, what each produces and what it requires, as the table and the
graph of the docs' step overview page (``docs/conf.py``) and as a Graphviz graph for a slide
(``pixi run -e docs step-overview``).

A requirement is part of the data flow when a step other than the resource steps produces
something that satisfies it; the graph draws it as an arrow from the data type's hub. Every
other requirement is a resource (reference data, regions, databases, indexes, existing files)
and is listed in the step's box."""

import argparse
import textwrap
from dataclasses import dataclass
from html import escape

from snappy_pipeline.orchestration import ON_HOLD_STEPS
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract.protocol import (
    TOOL,
    DataSignature,
    DataType,
    field_contract,
)

#: Steps that provide reference data and existing files, not results of other steps
RESOURCE_STEPS = frozenset({"external_data", "reference_download", "reference_index"})

#: One hue per data type that flows between steps, in a fixed order
TYPE_COLORS = {
    DataType.RAW: "#2a78d6",
    DataType.ALIGNMENTS: "#eb6834",
    DataType.VARIANTS: "#1baf7a",
    DataType.EXPRESSION: "#eda100",
    DataType.TABULAR: "#e87ba4",
    DataType.MODELS: "#008300",
    DataType.EXPORTS: "#6250d6",
    DataType.QC: "#e34948",
}
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
HAIRLINE = "#c3c2b7"
GRID = "#e1e0d9"


@dataclass(frozen=True)
class Requirement:
    key: str
    signature: DataSignature
    required: bool

    @property
    def text(self) -> str:
        """The signature, in parentheses when the key is optional."""
        return str(self.signature) if self.required else f"({self.signature})"


@dataclass(frozen=True)
class StepContract:
    name: str
    produces: tuple[DataSignature, ...]
    requirements: tuple[Requirement, ...]


def step_contracts() -> list[StepContract]:
    """The contracts of the registered steps that are not on hold."""
    contracts = []
    for name, workflow_cls in WORKFLOW_REGISTRY.items():
        if name in ON_HOLD_STEPS:
            continue
        field = workflow_cls.config_model_class.model_fields.get("depends_on")
        requirements = []
        if field is not None:
            for key, info in field.annotation.model_fields.items():
                signature = field_contract(info, TOOL).signature
                requirements.append(Requirement(key, signature, info.is_required()))
        contracts.append(StepContract(name, tuple(workflow_cls.produces), tuple(requirements)))
    return contracts


def producers(requirement: Requirement, consumer: str, contracts: list[StepContract]) -> list[str]:
    """The steps other than ``consumer`` and the resource steps that produce what satisfies
    ``requirement``; the tool placeholder matches any tool."""
    signature = DataSignature(requirement.signature.type, requirement.signature.tags - {TOOL})
    return [
        contract.name
        for contract in contracts
        if contract.name != consumer
        and contract.name not in RESOURCE_STEPS
        and any(produced.satisfies(signature) for produced in contract.produces)
    ]


def _small(text: str, color: str = INK_SECONDARY, italic: bool = False) -> str:
    html = f'<FONT POINT-SIZE="7" COLOR="{color}">{escape(text)}</FONT>'
    return f"<I>{html}</I>" if italic else html


def _node_label(contract: StepContract, needs: list[str]) -> str:
    lines = [f"<B>{contract.name}</B>"]
    if contract.name == "external_data":
        lines.append(_small("any type and tags, as the task declares"))
    for signature in contract.produces:
        color = TYPE_COLORS.get(signature.type, INK_SECONDARY)
        tags = ", ".join(sorted(signature.tags))
        text = signature.type.value + (f" {tags}" if tags else "")
        lines.append(f'<FONT POINT-SIZE="7" COLOR="{color}">●</FONT> ' + _small(text))
    if needs:
        for line in textwrap.wrap("needs " + "; ".join(needs), 56):
            lines.append(_small(line, italic=True))
    return "<" + '<BR ALIGN="LEFT"/>'.join(lines) + '<BR ALIGN="LEFT"/>>'


def dot(contracts: list[StepContract]) -> str:
    """The Graphviz graph: the steps, one hub per data type, arrows for reading and writing."""
    lines = [
        "digraph steps {",
        "    rankdir=TB; newrank=true; concentrate=true; nodesep=0.12; ranksep=0.45;",
        '    bgcolor="transparent"; fontname="Helvetica"; labelloc=b;',
        "    label=<"
        + _small(
            "Arrows: a step reads a data type (in), writes it (out) or both; "
            "needs: inputs from resources, optional ones in parentheses",
            INK_MUTED,
        )
        + ">;",
        '    node [shape=box, style="rounded,filled", fillcolor="#ffffff", '
        f'color="{HAIRLINE}", fontname="Helvetica", fontsize=9, fontcolor="{INK}", '
        'margin="0.08,0.03"];',
        "    edge [arrowsize=0.5, penwidth=1.0];",
    ]
    hubs = set()
    edges = {}
    for contract in contracts:
        needs = []
        read = set()
        for requirement in contract.requirements:
            if contract.name in RESOURCE_STEPS or not producers(
                requirement, contract.name, contracts
            ):
                needs.append(requirement.text)
            else:
                read.add(requirement.signature.type)
        written = (
            set()
            if contract.name in RESOURCE_STEPS
            else {signature.type for signature in contract.produces}
        )
        for data_type in read | written:
            hubs.add(data_type)
            edges[(contract.name, data_type)] = (data_type in read, data_type in written)
        lines.append(f'    "{contract.name}" [label={_node_label(contract, needs)}];')
    lines.append(
        "    subgraph cluster_resources { label=<"
        + _small("resources", INK_MUTED)
        + f'>; color="{GRID}"; style="rounded"; labelloc=t; '
        + " ".join(f'"{name}";' for name in sorted(RESOURCE_STEPS))
        + " }"
    )
    # The FASTQ files of the libraries come from an external_data task's files tables
    edges[("external_data", DataType.RAW)] = (False, True)
    hubs.add(DataType.RAW)
    for data_type in sorted(hubs):
        color = TYPE_COLORS[data_type]
        ink = INK if data_type == DataType.EXPRESSION else "#ffffff"
        lines.append(
            f'    "{data_type.value}" [style="rounded,filled", fillcolor="{color}", '
            f'color="{color}", fontcolor="{ink}", fontsize=10, label=<<B>{data_type.value}</B>>];'
        )
    for (name, data_type), (reads, writes) in sorted(edges.items()):
        color = TYPE_COLORS[data_type] + "b3"
        if reads and writes:
            lines.append(f'    "{data_type.value}" -> "{name}" [color="{color}", dir=both];')
        elif reads:
            lines.append(f'    "{data_type.value}" -> "{name}" [color="{color}"];')
        else:
            lines.append(f'    "{name}" -> "{data_type.value}" [color="{color}"];')
    lines.append("}")
    return "\n".join(lines) + "\n"


def rst_table(contracts: list[StepContract]) -> str:
    """The table of the docs: one row per step with what it produces and every ``depends_on``
    key with the data it requires."""
    rows = []
    for contract in contracts:
        if contract.name == "external_data":
            produces = "as the task declares (``produces: {type, tags}``)"
        else:
            produces = "\n       ".join(f"- ``{signature}``" for signature in contract.produces)
        depends_on = (
            "\n       ".join(
                f"- ``{requirement.key}``: ``{requirement.signature}``"
                + ("" if requirement.required else " (optional)")
                for requirement in contract.requirements
            )
            or "none"
        )
        rows.append(
            f"   * - :ref:`{contract.name} <step_{contract.name}>`\n"
            f"     - {produces}\n"
            f"     - {depends_on}"
        )
    header = "   * - Step\n     - Produces\n     - Depends on"
    return "\n".join(
        [".. list-table::", "   :header-rows: 1", "   :widths: 22 30 48", "", header, *rows]
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Print the step overview graph (Graphviz)")
    parser.add_argument(
        "--rst", action="store_true", help="print the docs table instead of the graph"
    )
    args = parser.parse_args(argv)
    contracts = step_contracts()
    print(rst_table(contracts) if args.rst else dot(contracts), end="")


if __name__ == "__main__":
    main()
