# -*- coding: utf-8 -*-
"""Tests of the code hash of a wrapper: its files, the files they reference and the hash as the
rerun trigger of the step's params."""

from pathlib import Path
from types import SimpleNamespace

from snappy_pipeline import wrapper_code
from snappy_pipeline.wrapper_code import code_files, code_of, wrapper_files
from snappy_pipeline.workflows import variant_filtration

ROOT = wrapper_code.package_root()


def _relative(paths) -> set[str]:
    return {str(Path(path).relative_to(ROOT)) for path in paths}


def test_the_files_of_a_wrapper_are_its_directory_and_the_shared_files_of_its_parents():
    assert _relative(wrapper_files("calls_report")) == {"wrappers/calls_report/wrapper.py"}
    # The helper files of a wrapper directory belong to the wrapper, sub-wrappers do not
    assert "wrappers/vcf2maf/vcf_to_table/parser.py" in _relative(
        wrapper_files("vcf2maf/vcf_to_table")
    )
    assert "wrappers/control_freec/transform/wrapper.py" not in _relative(
        wrapper_files("control_freec")
    )
    # A parent directory's plain files are shared by its wrappers, its own wrapper.py is not
    assert _relative(wrapper_files("cbioportal/merge_tables")) == {
        "wrappers/cbioportal/merge_tables/wrapper.py",
        "wrappers/cbioportal/merge_tables/script.R",
        "wrappers/cbioportal/helper_functions.R",
    }
    assert "wrappers/bcftools/wrapper.py" not in _relative(wrapper_files("bcftools/filter"))
    # The environment is Snakemake's software trigger, not code
    assert not any("environment" in file for file in _relative(wrapper_files("bwa")))


def test_the_code_files_follow_imports_scripts_and_console_scripts():
    files = code_files("calls_report")

    assert "snappy_wrapper.py" in files  # the base class of every wrapper
    assert "tools/calls_report/pedigree.py" in files  # via ``from . import pedigree``
    assert "tools/qc_summary/__init__.py" in files  # the framework
    assert "tools/vcf_events.py" in files  # via igv.py
    assert "tools/vcf_events.py" in code_files("varlociraptor/control_fdr")
    assert code_files("bcftools/concat_variants") == (
        "snappy_wrapper.py",
        "wrappers/bcftools/concat_variants/wrapper.py",
    )
    # A script run by path, and the console scripts of the package
    assert "tools/gcnv_merge_vcfs.py" in code_files("gcnv/merge_multikit_families")
    assert {"tools/vcf_filter_to_info.py", "tools/vcf_filter_from_info.py"} <= set(
        code_files("gatk_read_backed_phasing")
    )


def _fake_package(root: Path, wrapper_source: str, tool_source: str) -> None:
    for file in ("__init__.py", "snappy_wrapper.py", "tools/__init__.py"):
        (root / file).parent.mkdir(parents=True, exist_ok=True)
        (root / file).write_text((ROOT / file).read_text(encoding="utf-8"))
    (root / "wrappers" / "example").mkdir(parents=True)
    (root / "wrappers" / "example" / "wrapper.py").write_text(wrapper_source)
    (root / "tools" / "example.py").write_text(tool_source)


def _clear_caches() -> None:
    for function in (
        wrapper_code.wrapper_files,
        wrapper_code.code_files,
        wrapper_code.wrapper_code,
    ):
        function.cache_clear()


def test_the_hash_follows_the_code_but_not_the_docstrings_and_comments(monkeypatch, tmp_path):
    root = tmp_path / "snappy_wrappers"
    _fake_package(root, "from snappy_wrappers.tools.example import run\n", "x = 1\n")
    monkeypatch.setattr(wrapper_code, "package_root", lambda: root)
    _clear_caches()
    before = wrapper_code.wrapper_code("example")
    assert len(before) == 16

    # A comment or docstring changes nothing
    (root / "tools" / "example.py").write_text('"""The tool."""\n\nx = 1  # one\n')
    _clear_caches()
    assert wrapper_code.wrapper_code("example") == before
    # Another value in an imported tool changes the hash
    (root / "tools" / "example.py").write_text("x = 2\n")
    _clear_caches()
    changed = wrapper_code.wrapper_code("example")
    assert changed != before
    # So does the wrapper script itself
    (root / "wrappers" / "example" / "wrapper.py").write_text(
        "from snappy_wrappers.tools.example import run\nrun()\n"
    )
    _clear_caches()
    assert wrapper_code.wrapper_code("example") not in (before, changed)
    _clear_caches()


def test_the_steps_params_carry_the_code_of_the_actions_wrapper():
    assert code_of(None) is None
    part = SimpleNamespace(
        wrappers={"run": "calls_report"},
        get_params=lambda action: lambda wildcards: {"title": "t"},
        with_extra_arguments=lambda action, params: params,
    )
    step = object.__new__(variant_filtration.VariantFiltrationWorkflow)
    step._get_sub_step = lambda name: part

    params = step.get_params("calls_report", "run")(None)

    assert params["title"] == "t" and params["code"] == code_of("calls_report")
    part.wrappers = {}
    assert step.get_params("calls_report", "run")(None) == {"title": "t"}
