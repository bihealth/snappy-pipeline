#!/usr/bin/env python3
"""Check or fix the layout of the workflow modules and their step and step part classes.

A module lists its imports, constants and functions, then the step parts in the order the step
class registers them (a base part before its subclasses), then the step class.

A class lists its attributes first: those of its base classes in their order, root first, then
its own; a property that overrides a base attribute counts as that attribute. Then come the
methods: the overrides of base methods in the order of the bases, the action methods by action
(``actions`` order) and within an action inputs, outputs, log file, params, then the class's own
methods in their order.

Usage: class_layout.py [--write] [path ...]; without paths, all workflow modules.
"""

from __future__ import annotations

import argparse
import ast
import importlib
import inspect
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKFLOWS = ROOT / "snappy_pipeline" / "workflows"

#: The per-action methods, in the order of the rule sections they feed
_ACTION_KINDS = ("input_files", "output_files", "log_file", "params")
_ACTION_METHOD = re.compile(r"^_get_(input_files|output_files|log_file|params)(?:_(.+))?$")
_ATTRIBUTE_DECORATORS = {"property", "cached_property", "functools.cached_property"}


def module_paths() -> list[Path]:
    """The workflow modules: the step modules, the shared modules and the base classes."""
    return sorted(
        path
        for path in WORKFLOWS.rglob("*.py")
        if "model" not in path.name and path.name != "protocol.py"
    )


@dataclass
class Member:
    """A statement of a class body or a module with its source text."""

    node: ast.stmt
    name: str
    kind: str  # "attr" or "method"
    comments: list[str]
    body: list[str]
    index: int

    @property
    def text(self) -> str:
        return "\n".join([*self.comments, *self.body])

    @property
    def is_big(self) -> bool:
        return bool(self.comments) or len(self.body) > 1


def _first_line(node: ast.stmt) -> int:
    return min([node.lineno, *(d.lineno for d in getattr(node, "decorator_list", []))])


def _segments(
    lines: list[str], nodes: list[ast.stmt], start: int
) -> list[tuple[list[str], list[str]]]:
    """The comment lines above each node (from the previous node on) and its own lines."""
    segments = []
    previous_end = start
    for node in nodes:
        first = _first_line(node)
        comments = [
            line for line in lines[previous_end : first - 1] if line.strip().startswith("#")
        ]
        segments.append((comments, lines[first - 1 : node.end_lineno]))
        previous_end = node.end_lineno
    return segments


def _docstring(body: list[ast.stmt]) -> ast.stmt | None:
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        return body[0]
    return None


def _classify(node: ast.stmt, base_attributes: set[str]) -> tuple[str, str]:
    """The (name, kind) of a class body statement."""
    if isinstance(node, (ast.Assign, ast.AnnAssign)):
        target = node.target if isinstance(node, ast.AnnAssign) else node.targets[0]
        name = ast.unparse(target)
        return name, "method" if _ACTION_METHOD.match(name) else "attr"
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        decorators = {ast.unparse(d) for d in node.decorator_list}
        if decorators & _ATTRIBUTE_DECORATORS and node.name in base_attributes:
            return node.name, "attr"
        return node.name, "method"
    return ast.unparse(node)[:40], "attr"


def _declared(cls: type) -> list[tuple[str, str]]:
    """The (name, kind) of the members a class declares, in source order."""
    body = ast.parse(inspect.getsource(cls)).body[0].body
    docstring = _docstring(body)
    return [_classify(node, set()) for node in body if node is not docstring]


def _base_orders(cls: type) -> tuple[dict[str, int], dict[str, int]]:
    """The attributes and methods of the base classes of ``cls``, root first, each ranked."""
    attributes: dict[str, int] = {}
    methods: dict[str, int] = {}
    for base in reversed(cls.__mro__[1:]):
        if not base.__module__.startswith("snappy_pipeline"):
            continue
        for name, kind in _declared(base):
            ranks = attributes if kind == "attr" else methods
            ranks.setdefault(name, len(ranks))
    return attributes, methods


def _method_key(member: Member, methods: dict[str, int], actions: tuple[str, ...]):
    if member.name in methods:
        return (0, methods[member.name], 0)
    if match := _ACTION_METHOD.match(member.name):
        kind, action = match.groups()
        if action is None:  # the generic method of all actions, placed like the public one
            return (0, methods.get(f"get_{kind}", len(methods)), _ACTION_KINDS.index(kind))
        if action in actions:
            return (1, actions.index(action), _ACTION_KINDS.index(kind))
        return (1, len(actions) + 1, member.index)  # the actions of subclasses keep their order
    return (2, member.index, 0)


def reorder_class(lines: list[str], node: ast.ClassDef, cls: type) -> list[str]:
    """The lines of the class with its body in the conventional order."""
    attributes, methods = _base_orders(cls)
    body = node.body
    docstring = _docstring(body)
    items = [n for n in body if n is not docstring]
    header_end = docstring.end_lineno if docstring else body[0].lineno - 1
    members = [
        Member(n, *_classify(n, set(attributes)), comments, text, i)
        for i, (n, (comments, text)) in enumerate(zip(items, _segments(lines, items, header_end)))
    ]
    attrs = sorted(
        (m for m in members if m.kind == "attr"),
        key=lambda m: (0, attributes[m.name]) if m.name in attributes else (1, m.index),
    )
    actions = tuple(getattr(cls, "actions", ()) or ())
    defs = sorted(
        (m for m in members if m.kind == "method"),
        key=lambda m: _method_key(m, methods, actions),
    )
    out = lines[_first_line(node) - 1 : header_end]
    previous = None
    for member in _resolved([*attrs, *defs]):
        if member.kind == "method" or (previous is not None and (previous.is_big or member.is_big)):
            out.append("")
        out += [*member.comments, *member.body]
        previous = member
    return out


def _resolved(members: list[Member]) -> list[Member]:
    """The members with those that a class-level assignment names moved before it."""
    by_name = {m.name: m for m in members}
    result: list[Member] = []
    emitted: set[int] = set()

    def emit(member: Member) -> None:
        if member.index in emitted:
            return
        emitted.add(member.index)
        if isinstance(member.node, (ast.Assign, ast.AnnAssign)):
            for name in _mentions(member.node, set(by_name) - {member.name}):
                emit(by_name[name])
        result.append(member)

    for member in members:
        emit(member)
    return result


def _module_of(path: Path):
    names = list(path.resolve().relative_to(ROOT).with_suffix("").parts)
    if names[-1] == "__init__":
        names.pop()
    return importlib.import_module(".".join(names))


def _mentions(node: ast.AST, names: set[str]) -> list[str]:
    """The names of ``names`` that ``node`` mentions, by first appearance in the source."""
    found = sorted(
        (n for n in ast.walk(node) if isinstance(n, ast.Name) and n.id in names),
        key=lambda n: (n.lineno, n.col_offset),
    )
    return list(dict.fromkeys(n.id for n in found))


def _class_order(classes: dict[str, type], step: str | None, module) -> list[str]:
    """The step parts in registration order, each after its bases, then the step class."""
    if step is None:
        return list(classes)
    nodes = {
        n.name: n for n in ast.parse(inspect.getsource(module)).body if isinstance(n, ast.ClassDef)
    }
    mentioned = _mentions(nodes[step], set(classes) - {step})
    unmentioned = [name for name in classes if name not in mentioned and name != step]
    order: list[str] = []

    def place(name: str) -> None:
        """Place the class after its bases and the classes its body (not its methods) names."""
        if name in order:
            return
        order.append(name)
        bases = [b.__name__ for b in reversed(classes[name].__mro__[1:]) if b.__name__ in classes]
        statements = [n for n in nodes[name].body if not isinstance(n, ast.FunctionDef)]
        for dependency in [*bases, *_mentions(ast.Module(statements, []), set(classes))]:
            if dependency != name and dependency not in order:
                place(dependency)
                order.remove(dependency)
                order.insert(order.index(name), dependency)

    for name in [*mentioned, *unmentioned]:
        place(name)
    return [*order, step]


def reorder_module(path: Path) -> str:
    """The module with its statements and classes in the conventional order."""
    from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart

    lines = path.read_text(encoding="utf-8").splitlines()
    tree = ast.parse("\n".join(lines))
    module = _module_of(path)
    items = tree.body
    first_class = next((i for i, n in enumerate(items) if isinstance(n, ast.ClassDef)), None)
    if first_class is None:
        return "\n".join(lines) + "\n"
    # Everything above the first class stays as it is
    preamble_end = items[first_class - 1].end_lineno if first_class else 0
    segments = _segments(lines, items[first_class:], preamble_end)
    classes = {n.name: getattr(module, n.name) for n in items if isinstance(n, ast.ClassDef)}
    steps = [name for name, cls in classes.items() if issubclass(cls, BaseStep)]
    if len(steps) > 1:
        raise SystemExit(f"{path}: several step classes: {', '.join(steps)}")
    order = _class_order(classes, steps[0] if steps else None, module)
    moved, after, class_texts = [], [], {}
    for node, (comments, body) in zip(items[first_class:], segments):
        if isinstance(node, ast.ClassDef):
            cls = classes[node.name]
            if issubclass(cls, (BaseStep, BaseStepPart)):
                body = reorder_class(lines, node, cls)
            class_texts[node.name] = "\n".join([*comments, *body])
        elif _mentions(node, set(classes)):
            after.append("\n".join([*comments, *body]))
        else:
            moved.append("\n".join([*comments, *body]))
    parts = ["\n".join(lines[:preamble_end])] if preamble_end else []
    parts += moved
    parts += [class_texts[name] for name in order if name not in steps]
    parts += after
    parts += [class_texts[name] for name in steps]
    return "\n\n\n".join(parts) + "\n"


def _shape(source: str) -> tuple[list[str], list[str]]:
    """The statements of a module regardless of their order (class bodies included) and its
    comment lines."""
    shape = []
    for node in ast.parse(source).body:
        if isinstance(node, ast.ClassDef):
            members = sorted(ast.dump(n) for n in node.body)
            shape.append(f"class {node.name}({ast.dump(ast.Tuple(node.bases))}): {members}")
        else:
            shape.append(ast.dump(node))
    comments = sorted(line.strip() for line in source.splitlines() if line.strip().startswith("#"))
    return sorted(shape), comments


def check(paths: list[Path], write: bool = False) -> list[Path]:
    """Return the modules whose layout differs from the convention; ``write`` fixes them."""
    changed = []
    for path in paths:
        old = path.read_text(encoding="utf-8")
        new = reorder_module(path)
        if _shape(old) != _shape(new):
            raise SystemExit(f"{path}: the reordered module has other statements or comments")
        if ast.dump(ast.parse(old)) != ast.dump(ast.parse(new)):
            changed.append((path, new))
    # Written after all modules are read: the base classes' source must match the imported code
    if write:
        for path, new in changed:
            path.write_text(new, encoding="utf-8")
    return [path for path, _ in changed]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true", help="fix the layout")
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args()
    changed = check(args.paths or module_paths(), write=args.write)
    for path in changed:
        print(f"{'fixed' if args.write else 'out of order'}: {path.relative_to(ROOT)}")
    sys.exit(0 if args.write or not changed else 1)


if __name__ == "__main__":
    main()
