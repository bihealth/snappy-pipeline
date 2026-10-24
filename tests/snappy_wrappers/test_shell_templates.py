"""The shell templates of the wrappers become Bash for shfmt and back without change."""

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).parents[2] / "tests/scripts/shell_templates.py"


def _module():
    spec = importlib.util.spec_from_file_location("shell_templates", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TEMPLATE = r"""
run_vep() {{
    local out="$1"
    awk '{{print $1}}' {snakemake.input.vcf} > "$out"
}}
for ((i = 0; i < ${{#reads[@]}}; i++)); do
    echo {args[name]!r} {threads:>4} ${{reads[$i]}}
done
"""


def test_templates_become_bash_and_back():
    module = _module()

    script, placeholders = module.to_bash(TEMPLATE)

    assert "run_vep() {\n" in script and "${#reads[@]}" in script
    assert "{snakemake" not in script
    assert set(placeholders.values()) == {
        "{snakemake.input.vcf}",
        "{args[name]!r}",
        "{threads:>4}",
    }
    assert module.from_bash(script, placeholders) == TEMPLATE
