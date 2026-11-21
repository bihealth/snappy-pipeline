"""The code of a wrapper, as the rerun trigger of its rules.

Snakemake reruns a job when the rule's code changes, but for a ``wrapper:`` rule it records no
code at all: its code trigger covers shell commands and ``run:`` blocks, and its script mtime
check covers ``script:`` and ``notebook:`` rules.
So every rule that runs a wrapper gets a ``code`` param (``BaseStep.get_params``): a hash of the
wrapper's own files, of the files they import from ``snappy_wrappers`` (the ``ShellWrapper``
base, the tools of the reports) and of the scripts of the package they run by path or as a
console script.
Python files are hashed without docstrings and comments, so an edit to those reruns nothing.
The files are found by reading the imports, without importing anything (the tools need the
wrapper's environment).
"""

from __future__ import annotations

import ast
import functools
import hashlib
import importlib.metadata
import importlib.util
import re
from pathlib import Path

#: The ``snappy_wrappers`` package, located without importing its modules
PACKAGE = "snappy_wrappers"

#: Files of a wrapper directory that are no code: the environment with its pins and post-deploy
#: script (Snakemake's software trigger)
_NOT_CODE = re.compile(r"^(environment\..*|meta\.yaml|__init__\.py)$")

#: A script of the package that a wrapper runs by path, such as ``tools/gcnv_merge_vcfs.py``
_SCRIPT_PATH = re.compile(r"tools/[\w/]+\.py")


@functools.cache
def package_root() -> Path:
    spec = importlib.util.find_spec(PACKAGE)
    return Path(spec.origin).parent


def _is_wrapper(directory: Path) -> bool:
    return (directory / "wrapper.py").is_file()


def _code_files(directory: Path) -> list[Path]:
    """The code files of a wrapper directory: its own files and those of its helper
    subdirectories (not of the wrappers below it)."""
    files = []
    for path in sorted(directory.iterdir()):
        if path.is_dir():
            if path.name != "__pycache__" and not _is_wrapper(path):
                files += _code_files(path)
        elif not _NOT_CODE.match(path.name):
            files.append(path)
    return files


@functools.cache
def wrapper_files(wrapper: str) -> tuple[Path, ...]:
    """The files of a wrapper: its directory (``wrappers/<wrapper>``) and the plain files that
    its parent directories share between their wrappers (``cbioportal/helper_functions.R``)."""
    directory = package_root() / "wrappers" / wrapper
    if not directory.is_dir():
        raise ValueError(f"no wrapper {wrapper!r} below {package_root() / 'wrappers'}")
    files = _code_files(directory)
    for parent in directory.parents:
        if parent == package_root() / "wrappers":
            break
        files += [
            path
            for path in sorted(parent.iterdir())
            if path.is_file() and path.name != "wrapper.py" and not _NOT_CODE.match(path.name)
        ]
    return tuple(files)


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


@functools.cache
def _console_scripts() -> dict[str, Path]:
    """The console scripts of the package (``snappy-vcf_filter_to_info``) and their modules"""
    scripts = {}
    for entry in importlib.metadata.entry_points(group="console_scripts"):
        module = entry.value.split(":")[0]
        if module.split(".")[0] == PACKAGE and (file := _module_file(module)):
            scripts[entry.name] = file
    return scripts


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


def _referenced(path: Path) -> set[Path]:
    """The files of the package that a Python file imports, runs by path or as a console script"""
    text = path.read_text(encoding="utf-8")
    files = {file for module in _imports(path) if (file := _module_file(module))}
    files |= {package_root() / match for match in _SCRIPT_PATH.findall(text)}
    files |= {file for name, file in _console_scripts().items() if name in text}
    return {file for file in files if file.is_file()}


@functools.cache
def code_files(wrapper: str) -> tuple[str, ...]:
    """The files behind a wrapper: its own files (``wrapper_files``) and, transitively, the files
    of the package they reference, as paths relative to the package."""
    pending = list(wrapper_files(wrapper))
    seen: set[Path] = set(pending)
    while pending:
        path = pending.pop()
        if path.suffix != ".py":
            continue
        for file in _referenced(path) - seen:
            seen.add(file)
            pending.append(file)
    return tuple(sorted(str(path.relative_to(package_root())) for path in seen))


def _normalized(path: Path) -> bytes:
    """The content of a code file; a Python file without its docstrings and comments"""
    if path.suffix != ".py":
        return path.read_bytes()
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if ast.get_docstring(node, clean=False) is not None:
                node.body = node.body[1:] or [ast.Pass()]
    return ast.unparse(tree).encode()


@functools.cache
def wrapper_code(wrapper: str) -> str:
    """The hash of the code behind a wrapper: the rerun trigger of its rules."""
    digest = hashlib.sha256()
    for file in code_files(wrapper):
        digest.update(file.encode())
        digest.update(b"\0")
        digest.update(_normalized(package_root() / file))
        digest.update(b"\0")
    return digest.hexdigest()[:16]


def code_of(wrapper: str | None) -> str | None:
    """The ``code`` param of a rule that runs ``wrapper``, ``None`` for a rule without one."""
    return wrapper_code(wrapper) if wrapper else None
