# -*- coding: utf-8 -*-
"""Tests of the code hash of a rule's wrappers: the rules of a Snakefile, the imported files
and the hash as a rerun trigger."""

from pathlib import Path
from types import SimpleNamespace

from snappy_pipeline import wrapper_code
from snappy_pipeline.wrapper_code import code_of, imported_files, rule_wrappers, snakefile_of
from snappy_pipeline.workflows import variant_filtration

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = str(REPO / "snappy_pipeline" / "workflows" / "variant_filtration" / "Snakefile")


def test_the_rules_of_a_snakefile_name_the_wrappers_of_each_part_and_action():
    rules = rule_wrappers(SNAKEFILE)

    assert rules[("calls_report", "run")] == ("calls_report",)
    assert set(rules[("filter", "run")]) >= {
        "bcftools/filter",
        "vembrane/tag",
        "varlociraptor/control_fdr",
    }
    assert ("filter", "write_panel") not in rules or rules[("filter", "write_panel")]


def test_imported_files_follow_the_imports_of_the_package_transitively():
    files = imported_files("calls_report")

    assert "snappy_wrapper.py" in files  # the base class of every wrapper
    assert "tools/calls_report/pedigree.py" in files  # via ``from . import pedigree``
    assert "tools/qc_summary/__init__.py" in files  # the framework
    assert "tools/vcf_events.py" in files  # via igv.py
    assert not any(file.startswith("wrappers/") for file in files)  # the scripts themselves
    assert "tools/vcf_events.py" in imported_files("varlociraptor/control_fdr")
    assert imported_files("bcftools/concat_variants") == ("snappy_wrapper.py",)


def test_the_code_hash_changes_with_a_tool_and_is_part_of_the_params(monkeypatch, tmp_path):
    before = code_of(SNAKEFILE, "calls_report", "run")

    assert before and len(before) == 16
    assert code_of(SNAKEFILE, "calls_report", "run") == before
    assert code_of(SNAKEFILE, "no_such_part", "run") is None
    # Another content of an imported file gives another hash
    root = tmp_path / "snappy_wrappers"
    for file in (
        "__init__.py",
        "snappy_wrapper.py",
        "tools/__init__.py",
        "tools/calls_report/__init__.py",
    ):
        (root / file).parent.mkdir(parents=True, exist_ok=True)
        (root / file).write_text((wrapper_code.package_root() / file).read_text(encoding="utf-8"))
    (root / "wrappers" / "calls_report").mkdir(parents=True)
    (root / "wrappers" / "calls_report" / "wrapper.py").write_text(
        "from snappy_wrappers.tools.calls_report import summarize\n"
    )
    (root / "tools" / "calls_report" / "__init__.py").write_text("x = 1\n")
    monkeypatch.setattr(wrapper_code, "package_root", lambda: root)
    wrapper_code.imported_files.cache_clear()
    wrapper_code.wrapper_code.cache_clear()
    changed = wrapper_code.wrapper_code(("calls_report",))
    (root / "tools" / "calls_report" / "__init__.py").write_text("x = 2\n")
    wrapper_code.imported_files.cache_clear()
    wrapper_code.wrapper_code.cache_clear()
    assert wrapper_code.wrapper_code(("calls_report",)) != changed
    wrapper_code.imported_files.cache_clear()
    wrapper_code.wrapper_code.cache_clear()


def test_the_steps_params_carry_the_code():
    assert snakefile_of(variant_filtration.__file__).endswith("variant_filtration/Snakefile")
    part = SimpleNamespace(
        get_params=lambda action: lambda wildcards: {"title": "t"},
        with_extra_arguments=lambda action, params: params,
    )
    step = object.__new__(variant_filtration.VariantFiltrationWorkflow)
    step._get_sub_step = lambda name: part

    params = step.get_params("calls_report", "run")(None)

    assert params["title"] == "t" and len(params["code"]) == 16
