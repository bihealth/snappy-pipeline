# -*- coding: utf-8 -*-
"""The questions of the interactive commands, on questionary, in one style with the listings.

Every prompt takes :data:`PROMPT_KWARGS`, which tests set to a pipe input and a dummy output.
Ctrl-C raises ``KeyboardInterrupt``, which click turns into an abort.
"""

from typing import Any

import questionary
from questionary import Choice, Style

#: Bold questions; answers, cursor and highlight in cyan, like the first column of the tables
STYLE = Style(
    [
        ("qmark", "bold"),
        ("question", "bold"),
        ("answer", "fg:cyan"),
        ("pointer", "fg:cyan bold"),
        ("highlighted", "fg:cyan bold"),
        ("selected", "fg:cyan"),
        ("instruction", "fg:#888888"),
        ("disabled", "fg:#888888 italic"),
    ]
)

#: Keyword arguments of every prompt; tests give ``input`` (a pipe) and ``output``
PROMPT_KWARGS: dict[str, Any] = {}


def ask_text(message: str, default: str = "") -> str:
    """A line of text, with ``default`` prefilled."""
    return questionary.text(message, default=default, style=STYLE, **PROMPT_KWARGS).unsafe_ask()


def ask_confirm(message: str, default: bool = False) -> bool:
    """Yes or no by one key press, Enter for ``default``."""
    return questionary.confirm(message, default=default, style=STYLE, **PROMPT_KWARGS).unsafe_ask()


def ask_select(
    message: str, choices: list[tuple[Any, str, str | None]], default: Any = None
) -> Any:
    """One of ``choices``, ``(value, title, description)`` each, by cursor or by typing a part of
    a title; the description of the highlighted one is shown below the list."""
    items = [
        Choice(title, value=value, description=description) for value, title, description in choices
    ]
    return questionary.select(
        message,
        choices=items,
        default=default,
        style=STYLE,
        use_search_filter=True,
        use_jk_keys=False,
        **PROMPT_KWARGS,
    ).unsafe_ask()
