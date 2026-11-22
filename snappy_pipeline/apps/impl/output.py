# -*- coding: utf-8 -*-
"""The output of the listing commands, the same for every subcommand:
a table on a terminal, tab-separated lines on a pipe or with ``--plain``, JSON with ``--json``.
The data goes to stdout; messages go to stderr (:mod:`.logging`).
"""

import json
from dataclasses import dataclass
from typing import Any

import click
from rich import box
from rich.console import Console
from rich.table import Table


@dataclass(frozen=True)
class Column:
    """One column of a listing: the key of the rows and the header of the table."""

    key: str
    header: str
    no_wrap: bool = False


def listing_options(f):
    """``--json`` and ``--plain`` of a listing command."""
    f = click.option(
        "--json", "as_json", is_flag=True, help="JSON: an array with one object per row."
    )(f)
    return click.option(
        "--plain",
        is_flag=True,
        help="Tab-separated lines without header and frame, the output on a pipe.",
    )(f)


def text(value: Any) -> str:
    """A cell: ``None`` empty, a list joined by commas, a mapping as ``key: value`` pairs."""
    if value is None:
        return ""
    if isinstance(value, dict):
        return ", ".join(f"{key}: {text(item)}" for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return ", ".join(text(item) for item in value)
    return str(value)


def print_listing(
    rows: list[dict[str, Any]],
    columns: list[Column],
    title: str,
    *,
    as_json: bool = False,
    plain: bool = False,
    empty: str | None = None,
) -> None:
    """Print ``rows`` as the listing of a command.
    JSON keeps every key of the rows; the table and the lines show ``columns``.
    Without rows, a terminal shows ``empty`` (nothing if ``None``); a pipe gets no line, JSON ``[]``.
    """
    if as_json and plain:
        raise click.UsageError("--json and --plain exclude each other")
    if as_json:
        click.echo(json.dumps(rows, indent=2))
        return
    console = Console()
    if plain or not console.is_terminal:
        for row in rows:
            click.echo("\t".join(text(row.get(column.key)) for column in columns))
        return
    if not rows:
        if empty:
            console.print(empty, style="dim")
        return
    table = Table(
        title=title,
        title_justify="left",
        title_style="bold",
        box=box.SIMPLE_HEAD,
        header_style="bold",
        pad_edge=False,
    )
    for index, column in enumerate(columns):
        table.add_column(
            column.header, style="cyan" if index == 0 else None, no_wrap=column.no_wrap
        )
    for row in rows:
        table.add_row(*(text(row.get(column.key)) for column in columns))
    console.print(table)
