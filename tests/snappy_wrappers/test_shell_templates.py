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


def test_the_shell_text_of_the_tool_modules_is_checked():
    module = _module()
    root = SCRIPT.parents[2]
    found = {
        str(path): [t.node.lineno for t in module._module_templates(root / path)[0]]
        for path in module.SHELL_MODULES
    }
    assert {path: len(lines) for path, lines in found.items()} == {
        "snappy_wrappers/tools/bwa_mem.py": 3,  # READ_FUNCTIONS, MAPPING_PIPELINE, TEMPLATE
        "snappy_wrappers/tools/fgumi_consensus.py": 3,  # CLIP, NO_CLIP, CONSENSUS
        "snappy_wrappers/tools/fgumi_grouping.py": 4,  # EXTRACT, CORRECT, NO_CORRECT, GROUP
        "snappy_wrappers/tools/shards.py": 3,  # INTERVALS, RUN_SHARDS, GATHER
    }
