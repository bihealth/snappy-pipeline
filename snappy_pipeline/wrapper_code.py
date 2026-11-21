"""The code behind a wrapper rule, as a rerun trigger.

Snakemake reruns a job when the rule's code changes, and for a ``wrapper:`` rule that is the wrapper script alone.
The scripts of snappy are thin and import their logic from ``snappy_wrappers`` (the ``ShellWrapper`` base, the tools of the reports), which Snakemake does not see.
So every params function of a step gets a ``code`` entry: a hash of the files that the rule's wrappers import from ``snappy_wrappers``, found by reading the imports, without importing anything (the tools need the wrapper's environment).
The rules of a step part and action are read from the step's Snakefile.
"""

from __future__ import annotations

import ast
import functools
import hashlib
import importlib.util
import os
import re
from pathlib import Path

#: The ``snappy_wrappers`` package, located without importing its modules
PACKAGE = "snappy_wrappers"

_RULES = re.compile(r"^\s*rule\s+\w+\s*:", re.M)
_PART_ACTION = re.compile(r"""wf\.get_(?:params|input_files)\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)""")
_WRAPPER = re.compile(r"""wf\.wrapper_path\(\s*"([^"]+)"\s*\)""")


@functools.cache
def package_root() -> Path:
    spec = importlib.util.find_spec(PACKAGE)
    return Path(spec.origin).parent


@functools.cache
def rule_wrappers(snakefile: str) -> dict[tuple[str, str], tuple[str, ...]]:
    """(step part, action) -> the wrappers of the rules of a Snakefile that take its params or inputs."""
    text = Path(snakefile).read_text(encoding="utf-8")
    starts = [match.start() for match in _RULES.finditer(text)] + [len(text)]
    found: dict[tuple[str, str], set[str]] = {}
    for start, end in zip(starts, starts[1:]):
        block = text[start:end]
        wrappers = _WRAPPER.findall(block)
        if not wrappers:
            continue
        for part, action in set(_PART_ACTION.findall(block)):
            found.setdefault((part, action), set()).update(wrappers)
    return {key: tuple(sorted(value)) for key, value in found.items()}


def _module_file(module: str) -> Path | None:
    """The file of a module of the package, a package's ``__init__.py`` included."""
    parts = module.split(".")[1:]
    if not parts:  # ``import snappy_wrappers``
        return package_root() / "__init__.py"
    relative = Path(*parts)
    for candidate in (
        package_root() / relative.with_suffix(".py"),
        package_root() / relative / "__init__.py",
    ):
        if candidate.is_file():
            return candidate
    return None


def _imports(path: Path) -> set[str]:
    """The modules of the package that ``path`` imports, relative imports resolved."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    relative = path.relative_to(package_root())
    parts = relative.with_suffix("").parts
    package = parts[:-1] if parts[-1] != "__init__" else parts[:-1]
    own = [PACKAGE, *package]
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] == PACKAGE:
                    modules.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = own[: len(own) - node.level + 1]
                module = ".".join([*base, node.module] if node.module else base)
            elif node.module and node.module.split(".")[0] == PACKAGE:
                module = node.module
            else:
                continue
            modules.add(module)
            for alias in node.names:  # ``from . import x`` names modules
                modules.add(f"{module}.{alias.name}")
    return modules


@functools.cache
def imported_files(wrapper: str) -> tuple[str, ...]:
    """The files of the package that a wrapper's scripts import, transitively, without the
    scripts themselves (Snakemake hashes them), as paths relative to the package."""
    directory = package_root() / "wrappers" / wrapper
    pending = [path for path in sorted(directory.glob("*.py"))]
    own = set(pending)
    seen: set[Path] = set()
    while pending:
        path = pending.pop()
        for module in _imports(path):
            file = _module_file(module)
            if file is not None and file not in seen and file not in own:
                seen.add(file)
                pending.append(file)
    return tuple(sorted(str(path.relative_to(package_root())) for path in seen))


@functools.cache
def wrapper_code(wrappers: tuple[str, ...]) -> str:
    """A hash of the imported files of ``wrappers``: the rerun trigger of their rules."""
    digest = hashlib.sha256()
    files = sorted({file for wrapper in wrappers for file in imported_files(wrapper)})
    for file in files:
        digest.update(file.encode())
        digest.update(b"\\0")
        digest.update((package_root() / file).read_bytes())
        digest.update(b"\\0")
    return digest.hexdigest()[:16]


def code_of(snakefile: str, part: str, action: str) -> str | None:
    """The ``code`` param of the rules of ``part`` and ``action`` in ``snakefile``: the hash of
    what their wrappers import, ``None`` for rules without a wrapper."""
    wrappers = rule_wrappers(snakefile).get((part, action))
    if not wrappers:
        return None
    return wrapper_code(wrappers)


def snakefile_of(step_module_file: str) -> str:
    return os.path.join(os.path.dirname(step_module_file), "Snakefile")
