"""Check or format the shell scripts and the shell templates of the wrappers with shfmt.

A wrapper's template is the literal string it passes to ``ShellWrapper(snakemake).run()`` or
``shell()``, directly, through a variable with literal assignments or with ``.replace()`` calls. Its Snakemake
placeholders (``{snakemake.input.bam}``) become words and its doubled braces single ones, so that
shfmt parses it as Bash; formatting writes it back the other way. A template that is no raw string
is only checked.
"""

import argparse
import ast
import string
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).parents[2]

#: shfmt options: Bash, 4 spaces, pipes and && at the start of continued lines, indented case
#: patterns, a space after redirections (the options that change the existing wrappers least)
SHFMT = [
    "shfmt",
    "--language-dialect",
    "bash",
    "--indent",
    "4",
    "--binary-next-line",
    "--case-indent",
    "--space-redirects",
]


@dataclass
class Template:
    path: Path
    node: ast.Constant
    #: The template as Bash, and the placeholder words with their original text
    script: str
    placeholders: dict[str, str]

    @property
    def raw(self) -> bool:
        return self.source_segment().startswith(("r", "R"))

    def source_segment(self) -> str:
        lines = self.path.read_text().splitlines(keepends=True)
        start = sum(map(len, lines[: self.node.lineno - 1])) + self.node.col_offset
        end = sum(map(len, lines[: self.node.end_lineno - 1])) + self.node.end_col_offset
        return "".join(lines)[start:end]


def to_bash(template: str) -> tuple[str, dict[str, str]]:
    """The template as Bash: each placeholder a word, doubled braces single."""
    parts, placeholders = [], {}
    for literal, field, spec, conversion in string.Formatter().parse(template):
        parts.append(literal)
        if field is not None:
            word = f"__SMK{len(placeholders)}__"
            original = "{" + field + (f"!{conversion}" if conversion else "")
            placeholders[word] = original + (f":{spec}" if spec else "") + "}"
            parts.append(word)
    return "".join(parts), placeholders


def from_bash(script: str, placeholders: dict[str, str]) -> str:
    """The template of a Bash script from ``to_bash``."""
    template = script.replace("{", "{{").replace("}", "}}")
    for word, original in placeholders.items():
        template = template.replace(word, original)
    return template


def _templates(path: Path) -> tuple[list[Template], list[str]]:
    """The templates of a wrapper, and the shell calls whose template is not literal."""
    tree = ast.parse(path.read_text())
    assignments: dict[str, list[ast.expr]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    assignments.setdefault(target.id, []).append(node.value)
    found, skipped = [], []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and node.args):
            continue
        func = node.func
        is_shell = isinstance(func, ast.Name) and func.id == "shell"
        is_run = (
            isinstance(func, ast.Attribute)
            and func.attr == "run"
            and isinstance(func.value, ast.Call)
            and getattr(func.value.func, "id", None) == "ShellWrapper"
        )
        if not (is_shell or is_run):
            continue
        # A variable with literal assignments (if/else), and a literal with .replace() calls
        args = assignments.get(node.args[0].id, []) if isinstance(node.args[0], ast.Name) else []
        for arg in args or [node.args[0]]:
            while isinstance(arg, ast.Call) and getattr(arg.func, "attr", None) == "replace":
                arg = arg.func.value
            if not (isinstance(arg, ast.Constant) and isinstance(arg.value, str)):
                skipped.append(f"{path.relative_to(ROOT)}:{node.lineno}: template is not a literal")
                continue
            try:
                found.append(Template(path, arg, *to_bash(arg.value)))
            except ValueError as e:
                skipped.append(f"{path.relative_to(ROOT)}:{node.lineno}: no format template ({e})")
    return found, skipped


def shfmt(script: str, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(SHFMT + list(args), input=script, capture_output=True, text=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Format instead of checking")
    parser.add_argument("paths", nargs="*", type=Path, help="Files (default: the repository's)")
    args = parser.parse_args()

    paths = [p.resolve() for p in args.paths] or [
        ROOT / p
        for p in subprocess.run(
            ["git", "ls-files", "*.sh", "snappy_wrappers/wrappers/**/*.py"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.split()
    ]
    problems = []
    for path in paths:
        name = path.relative_to(ROOT)
        if path.suffix == ".sh":
            result = subprocess.run(SHFMT + ["--write" if args.write else "--diff", str(path)])
            if result.returncode:
                problems.append(f"{name}: shfmt")
            continue
        if path.suffix != ".py" or not path.is_relative_to(ROOT / "snappy_wrappers" / "wrappers"):
            continue
        templates, skipped = _templates(path)
        problems += skipped
        replacements = []
        for template in templates:
            result = shfmt(template.script)
            if result.returncode:
                line = template.node.lineno
                problems.append(f"{name}:{line}: {result.stderr.strip()} (line in the template)")
                continue
            formatted = result.stdout
            if template.node.value.startswith("\n"):
                formatted = "\n" + formatted
            if not template.node.value.endswith("\n"):
                formatted = formatted.rstrip("\n")
            if formatted == template.script:
                continue
            if not args.write or not template.raw:
                problems.append(f"{name}:{template.node.lineno}: not formatted")
                continue
            new = from_bash(formatted, template.placeholders)
            segment = template.source_segment()
            quote = segment[1:4] if segment[1:4] in ('"""', "'''") else segment[1]
            replacements.append((template, f"r{quote}{new}{quote}"))
        if replacements:
            text = path.read_text()
            for template, literal in replacements:
                text = text.replace(template.source_segment(), literal, 1)
            path.write_text(text)
    for problem in problems:
        print(problem, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
