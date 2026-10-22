# -*- coding: utf-8 -*-
"""Checks for the coding conventions in ``docs/dev_conventions.md``."""

from __future__ import annotations

import ast
import collections
import enum
import importlib
import json
import os
import re
from collections.abc import Iterator
from pathlib import Path

from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

REPO = Path(__file__).resolve().parents[2]
WORKFLOWS = REPO / "snappy_pipeline" / "workflows"
SNAPSHOTS = Path(__file__).resolve().parent / "snapshots"

#: Path convention violations that already exist. The test fails for new ones and
#: for fixed ones, so the list only shrinks.
PATH_BASELINE = SNAPSHOTS / "path_convention_violations.json"

#: Input functions that a Snakefile calls with Snakemake's ``checkpoints`` object.
CHECKPOINT_INPUT_FUNCTIONS = {
    ("GcnvModelStepPart", "_get_input_files_post_germline_calls"),
    ("PostGermlineCallsMixin", "_get_input_files_post_germline_calls"),
}


def _class_members() -> Iterator[tuple[Path, ast.ClassDef, ast.stmt]]:
    for path in sorted(WORKFLOWS.rglob("*.py")):
        for cls in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(cls, ast.ClassDef):
                for node in cls.body:
                    yield path, cls, node


def _location(path: Path, node: ast.stmt) -> str:
    return f"{path.relative_to(REPO)}:{node.lineno}"


def test_no_get_args():
    found = [
        _location(path, node)
        for path, _, node in _class_members()
        if isinstance(node, ast.FunctionDef) and re.match(r"_?get_args", node.name)
    ]
    assert not found, f"get_args is called get_params: {found}"


def test_input_and_params_functions_take_wildcards():
    allowed = (["self", "wildcards"], ["self", "wildcards", "input"])
    found = []
    for path, cls, node in _class_members():
        if not isinstance(node, ast.FunctionDef):
            continue
        if not re.fullmatch(r"_get_(input_files|params)_\w+", node.name):
            continue
        if (cls.name, node.name) in CHECKPOINT_INPUT_FUNCTIONS:
            continue
        args = node.args
        if [a.arg for a in args.args] not in allowed or args.vararg or args.kwarg:
            found.append(f"{_location(path, node)} {cls.name}.{node.name}")
    assert not found, f"use (self, wildcards) or (self, wildcards, input): {found}"


def test_actions_are_tuples():
    found = [
        f"{_location(path, node)} {cls.name}"
        for path, cls, node in _class_members()
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "actions" for t in node.targets)
        and not isinstance(node.value, ast.Tuple)
    ]
    assert not found, f"actions must be a tuple: {found}"


def test_snakefiles_pass_input_and_params_functions():
    found = []
    rule_files = sorted([*WORKFLOWS.glob("*/Snakefile"), *WORKFLOWS.glob("*/*.rules")])
    for path in rule_files:
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            unpacked = re.search(r"\*\*\(?wf\.get_(input_files|params)\(", line)
            stray_params = "wf.get_params(" in line and "args=wf.get_params(" not in line
            if unpacked or stray_params:
                found.append(f"{path.relative_to(REPO)}:{lineno}: {line.strip()}")
    assert not found, f"use unpack(wf.get_input_files(...)) and args=wf.get_params(...): {found}"


# Path conventions ---------------------------------------------------------------------------------


def _tool_names_by_step() -> dict[str, set[str]]:
    """Return the lower-case values of each step's ``Tool`` enums."""
    result = {}
    for step in WORKFLOW_REGISTRY:
        module = importlib.import_module(f"snappy_pipeline.workflows.{step}.model")
        result[step] = {
            str(member.value).lower()
            for name, obj in vars(module).items()
            if isinstance(obj, type) and issubclass(obj, enum.Enum) and "Tool" in name
            for member in obj
        }
    return result


def _path_violations() -> dict[str, list[str]]:
    """Return step -> violations of the path conventions in all DAG snapshots.

    Paths are normalized back to templates by replacing wildcard values with ``{name}``.
    """
    tools = _tool_names_by_step()
    steps = sorted(WORKFLOW_REGISTRY, key=len, reverse=True)
    found: dict[str, set[str]] = {}
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            if job["rule"] == "all" or job["rule"].endswith("_publish"):
                continue
            step = next(s for s in steps if job["rule"].startswith(s + "_"))
            values = sorted(job["wildcards"].items(), key=lambda kv: (-len(str(kv[1])), kv[0]))
            for path in job["output"] + job["log"]:
                template = re.sub(r"^tasks/[^/]+/", "", path)
                for name, value in values:
                    template = template.replace(str(value), "{%s}" % name)
                tokens = set(re.split(r"[/._\-]", re.sub(r"\{\w+\}", "", template)))
                if tokens & tools[step]:
                    found.setdefault(step, set()).add(f"tool name: {template}")
                if re.search(r"\}_", template):
                    found.setdefault(step, set()).add(f"underscore after wildcard: {template}")
    return {step: sorted(violations) for step, violations in sorted(found.items())}


def test_path_conventions_only_improve():
    """No tool names in paths and ``{entity}.{suffix}`` naming."""
    current = _path_violations()
    if os.environ.get("SNAPPY_UPDATE_SNAPSHOTS"):
        PATH_BASELINE.write_text(json.dumps(current, indent=1) + "\n", encoding="utf-8")
        return

    baseline = json.loads(PATH_BASELINE.read_text(encoding="utf-8"))
    flat_current = {(step, v) for step, vs in current.items() for v in vs}
    flat_baseline = {(step, v) for step, vs in baseline.items() for v in vs}
    new = sorted(flat_current - flat_baseline)
    fixed = sorted(flat_baseline - flat_current)
    assert not new, f"new path convention violations: {new}"
    assert not fixed, (
        f"{len(fixed)} violations are fixed; shrink the baseline with SNAPPY_UPDATE_SNAPSHOTS=1"
    )


def test_step_snakefiles_fetch_their_workflow_object():
    found = []
    for path in sorted(WORKFLOWS.glob("*/Snakefile")):
        text = path.read_text(encoding="utf-8")
        if path.parent.name == "external_data":
            continue  # configuration carrier without rules
        if 'wf = task_instance(config["__task_name__"])' not in text or re.search(
            r"Workflow\(", text
        ):
            found.append(str(path.relative_to(REPO)))
    assert not found, f"use wf = task_instance(config['__task_name__']) only: {found}"


#: Allowed ``depends_on`` keys: named after the data, with role prefixes where a step needs two
#: inputs of one kind.
DEPENDS_ON_KEYS = {
    "reads",
    "kits",
    "regions",
    "excluded_regions",
    "mappable_regions",
    "par_regions",
    "microsatellites",
    "snp_loci",
    "gene_regions",
    "common_variants",
    "germline_resource",
    "database",
    "alignments",
    "variants",
    "somatic_variants",
    "germline_variants",
    "combined_variants",
    "annotated_variants",
    "phased_variants",
    "copy_number",
    "structural_variants",
    "fusions",
    "proteome",
    "expression",
    "hla_types",
    "strandedness",
    "panel_of_normals",
    "access",
    "reference",
    "features",
    "dbsnp",
    "exons",
    "index",
    "gcnv_model",
}


def test_depends_on_keys_use_the_vocabulary():
    found = []
    for step, workflow_cls in WORKFLOW_REGISTRY.items():
        depends_on = workflow_cls.config_model_class.model_fields.get("depends_on")
        if depends_on is not None:
            keys = set(depends_on.annotation.model_fields) - DEPENDS_ON_KEYS
            found += [f"{step}.{key}" for key in sorted(keys)]
    assert not found, f"depends_on keys outside the vocabulary in docs/dev_conventions.md: {found}"


def test_tool_is_required():
    found = [
        step
        for step, workflow_cls in WORKFLOW_REGISTRY.items()
        if (tool := workflow_cls.config_model_class.model_fields.get("tool")) is not None
        and not tool.is_required()
    ]
    assert not found, f"tool must be set explicitly, without a default: {found}"


def test_upstream_paths_come_from_contracts():
    """Consumers ask providers for named outputs instead of building their paths."""
    found = [
        f"{_location(path, node)}"
        for path in sorted(WORKFLOWS.rglob("*.py"))
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "upstream"
    ]
    assert not found, f"use get_upstream_paths(field, ...) instead of upstream(field): {found}"


#: Params a wrapper reads only in a case that the snapshot jobs do not take.
#: Keys are ``(wrapper, name)``; ``name`` is a key of ``params.args`` or, with the
#: ``params.`` prefix, a top-level param.
WRAPPER_PARAM_EXCEPTIONS = {
    ("cnvkit/report", "breaks"): "read only for a breaks output",
    ("cnvkit/report", "genemetrics"): "read only for a genemetrics output",
    ("cnvkit/report", "segmetrics"): "read only for a segmetrics output",
    ("mbcs", "barcode_config"): "read only with mbcs.use_barcodes",
    ("vembrane/tag", "expression"): "read only in filter mode",
}

_ARGS_ALIAS = re.compile(
    r"""^[ \t]*(\w+)[ \t]*(?::[ \t]*[\w\[\], |]+)?=[ \t]*(?:dict\(\s*)?"""
    r"""(?:snakemake\.params\.args\b|getattr\(\s*snakemake\.params\s*,\s*["']args["'])""",
    re.M,
)
_PARAMS_ATTR = re.compile(r"""\bsnakemake\.params\.(?!args\b)(\w+)\b(?!\s*\()""")
_PARAMS_ITEM = re.compile(r"""\bsnakemake\.params\[\s*["'](\w+)["']\s*\]""")
_PARAMS_GETATTR = re.compile(r"""\bgetattr\(\s*snakemake\.params\s*,\s*["'](\w+)["']\s*\)""")
_PARAMS_GUARD = re.compile(
    r"""\bhasattr\(\s*snakemake\.params\s*,\s*["'](\w+)["']"""
    r"""|\bgetattr\(\s*snakemake\.params\s*,\s*["'](\w+)["']\s*,"""
    r"""|\bsnakemake\.params\.get\(\s*["'](\w+)["']"""
    r"""|["'](\w+)["']\s+(?:not\s+)?in\s+snakemake\.params\b"""
)
_ARGS_GUARD = (
    r"""["'](\w+)["']\s+(?:not\s+)?in\s+(?:snakemake\.params\.)?{name}\b"""
    r"""|(?:snakemake\.params\.)?{name}\.get\(\s*["'](\w+)["']"""
)


def _wrapper_reads(wrapper: str) -> tuple[set[str], set[str]]:
    """Return the (args keys, top-level params) that a wrapper reads without guard or default.

    Reads are found textually in the wrapper's Python files:

    - ``snakemake.params.<name>``, ``snakemake.params["<name>"]`` and
      ``getattr(snakemake.params, "<name>")`` for top-level params;
    - ``args[...]`` and ``{args[...]}`` in command templates, for ``params.args`` keys, with the
      same for every variable assigned ``snakemake.params.args`` (also via ``getattr`` or
      ``dict(...)``).

    Reads with a default or a guard (``.get``, ``in``, ``hasattr``, three-argument ``getattr``) are
    optional. A name guarded anywhere in the wrapper's files is optional in all of them.
    """
    args_keys, top, guarded_args, guarded_top = set(), set(), set(), set()
    for path in (REPO / "snappy_wrappers" / "wrappers" / wrapper).glob("*.py"):
        text = path.read_text(encoding="utf-8")
        for name in {"args"} | set(_ARGS_ALIAS.findall(text)):
            literal = re.compile(
                rf"""(?:^|[^\w.]|snakemake\.params\.){name}\[\s*["'](\w+)["']\s*\]""", re.M
            )
            template = re.compile(rf"\{{(?:snakemake\.params\.)?{name}\[(\w+)\]")
            guard = re.compile(_ARGS_GUARD.replace("{name}", name))
            args_keys |= set(literal.findall(text)) | set(template.findall(text))
            guarded_args |= {a or b for a, b in guard.findall(text)}
        top |= (
            set(_PARAMS_ATTR.findall(text))
            | set(_PARAMS_ITEM.findall(text))
            | set(_PARAMS_GETATTR.findall(text))
        )
        guarded_top |= {next(g for g in groups if g) for groups in _PARAMS_GUARD.findall(text)}
    return args_keys - guarded_args, top - guarded_top


def test_wrappers_get_the_params_they_read():
    """Every param a wrapper reads is in the params of the snapshot jobs that run it.

    A read counts when it is unconditional: ``snakemake.params.<name>``, ``snakemake.params[...]``,
    ``getattr(snakemake.params, "<name>")`` and ``args[...]`` (``args`` being ``params.args``).
    Reads with a default or a guard (``.get(...)``, ``in``, ``hasattr``, three-argument ``getattr``)
    may be absent. A name guarded anywhere in the wrapper's files counts as optional.
    ``WRAPPER_PARAM_EXCEPTIONS`` lists reads that are conditional in another way.
    """
    found, used_exceptions = set(), set()
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            wrapper = (job.get("wrapper") or "").partition("snappy_wrappers/wrappers/")[2]
            params = job.get("params")
            args = (params or {}).get("args", {})
            if not wrapper or not isinstance(params, dict) or not isinstance(args, dict):
                continue  # params not known when dumping
            args_keys, top = _wrapper_reads(wrapper)
            missing = [(k, f"args[{k}]") for k in sorted(args_keys - set(args))]
            missing += [(f"params.{k}", f"params.{k}") for k in sorted(top - set(params))]
            for key, label in missing:
                if (wrapper, key) in WRAPPER_PARAM_EXCEPTIONS:
                    used_exceptions.add((wrapper, key))
                else:
                    found.add(f"{wrapper} reads {label}, missing in {job['rule']}")
    assert not found, f"wrappers read params that their rules do not pass: {sorted(found)}"
    stale = sorted(set(WRAPPER_PARAM_EXCEPTIONS) - used_exceptions)
    assert not stale, f"remove these entries from WRAPPER_PARAM_EXCEPTIONS: {stale}"


#: Wrappers that read their params in a way ``test_wrappers_read_the_params_they_get`` does not see
WRAPPERS_READING_PARAMS_DYNAMICALLY = {
    "cbioportal/case_lists": "iterates over its params (cbioportal is postponed)",
    "cbioportal/meta_files": "iterates over its params (cbioportal is postponed)",
    "mehari/bam_qc": "its keys are library names (varfish export, left for later)",
    "varfish_annotator/annotate_svs": "varfish export, left for later",
}


def test_wrappers_read_the_params_they_get():
    """Every param of a wrapper job occurs in the wrapper's Python files.

    A param that no wrapper reads is dead, or a config option that never reaches the tool, while
    Snakemake reruns the job when it changes. ``extra_arguments`` is read by the base class.
    """
    found, used = set(), set()
    for _, job in _snapshot_jobs():
        wrapper = (job.get("wrapper") or "").partition("snappy_wrappers/wrappers/")[2]
        params = job.get("params")
        if not wrapper or not isinstance(params, dict):
            continue
        if wrapper in WRAPPERS_READING_PARAMS_DYNAMICALLY:
            used.add(wrapper)
            continue
        paths = (REPO / "snappy_wrappers" / "wrappers" / wrapper).glob("*.py")
        text = "\n".join(path.read_text(encoding="utf-8") for path in paths)
        args = params.get("args", {})
        keys = (
            {key for key in args if key != "extra_arguments"} if isinstance(args, dict) else set()
        )
        keys |= {f"params.{key}" for key in params if key != "args"}
        for key in keys:
            if not re.search(rf"\b{re.escape(key.removeprefix('params.'))}\b", text):
                found.add(f"{wrapper} does not read {key}")
    assert not found, f"rules pass params that their wrappers do not read: {sorted(found)}"
    stale = sorted(set(WRAPPERS_READING_PARAMS_DYNAMICALLY) - used)
    assert not stale, f"remove these entries from WRAPPERS_READING_PARAMS_DYNAMICALLY: {stale}"


#: Params that hold file paths. The test fails for new ones and for fixed ones, so
#: the list only shrinks.
PARAMS_FILES_BASELINE = SNAPSHOTS / "params_file_violations.json"


def _flat_params(value, key=""):
    if isinstance(value, dict):
        for name, item in value.items():
            yield from _flat_params(item, f"{key}.{name}")
    elif isinstance(value, list):
        for item in value:
            yield from _flat_params(item, f"{key}[]")
    else:
        yield key, value


def _params_file_violations() -> list[str]:
    """Return ``rule key`` for each param of a wrapper job whose value is a file path."""
    found = set()
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            params = job.get("params") or {}
            if not job.get("wrapper") or not isinstance(params.get("args", {}), dict):
                continue
            inputs = set(job["input"])
            for key, value in _flat_params(params):
                if isinstance(value, str) and (
                    value in inputs or value.startswith(("<repo>/", "<project>/", "/", "tasks/"))
                ):
                    found.add(f"{job['rule']} {key}")
    return sorted(found)


def test_files_are_inputs_not_params():
    """Wrappers get the files they read as inputs, so the DAG and rerun triggers see them."""
    current = _params_file_violations()
    if os.environ.get("SNAPPY_UPDATE_SNAPSHOTS"):
        PARAMS_FILES_BASELINE.write_text(json.dumps(current, indent=1) + "\n", encoding="utf-8")
        return

    baseline = json.loads(PARAMS_FILES_BASELINE.read_text(encoding="utf-8"))
    new = sorted(set(current) - set(baseline))
    fixed = sorted(set(baseline) - set(current))
    assert not new, f"params that hold file paths; pass these files as inputs: {new}"
    assert not fixed, (
        f"{len(fixed)} file params are fixed; shrink the baseline with SNAPPY_UPDATE_SNAPSHOTS=1"
    )


#: Wrappers that do not run through a base class yet, with the reason
WRAPPERS_WITHOUT_BASE_CLASS = {"pvactools/create_proteome": "on hold with create_proteome"}


def test_wrappers_use_a_base_class():
    """Wrappers run through ShellWrapper, RWrapper or PythonWrapper (log, conda info, md5)."""
    wrappers = REPO / "snappy_wrappers" / "wrappers"
    found, base = [], re.compile(r"\b(ShellWrapper|RWrapper|PythonWrapper)\(")
    for path in sorted(wrappers.rglob("wrapper.py")):
        name = str(path.parent.relative_to(wrappers))
        if not base.search(path.read_text(encoding="utf-8")):
            found.append(name)
    assert sorted(set(found) - set(WRAPPERS_WITHOUT_BASE_CLASS)) == []
    stale = sorted(set(WRAPPERS_WITHOUT_BASE_CLASS) - set(found))
    assert not stale, f"remove these entries from WRAPPERS_WITHOUT_BASE_CLASS: {stale}"


def _rule_blocks(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (rule name, rule text) for the rules of a Snakefile."""
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"^(\s*)(?:rule|checkpoint) (\w+):", lines[i])
        if not m:
            i += 1
            continue
        indent, j = len(m.group(1)), i + 1
        while j < len(lines) and (
            not lines[j].strip() or len(lines[j]) - len(lines[j].lstrip()) > indent
        ):
            j += 1
        yield m.group(2), "\n".join(lines[i:j])
        i = j


def test_wrapper_rules_write_a_benchmark():
    """Benchmarks record wall time and memory per job to calibrate resources.

    Snakemake does not allow benchmarks for rules that are cached between workflows.
    """
    snakefiles = sorted(WORKFLOWS.glob("*/Snakefile")) + sorted(WORKFLOWS.glob("*/*.rules"))
    found = [
        f"{path.relative_to(WORKFLOWS)}: {name}"
        for path in snakefiles
        for name, text in _rule_blocks(path)
        if "wrapper:" in text and "benchmark:" not in text and "cache:" not in text
    ]
    assert not found, f"wrapper rules without benchmark: wf.get_benchmark_file(...): {found}"


#: Wrappers that run without the ``extra_arguments`` of a tool, by reason; the list only shrinks
WRAPPERS_WITHOUT_EXTRA_ARGUMENTS = {
    wrapper: reason
    for reason, wrappers in {
        "Python code, no command line": (
            "cbioportal/case_lists",
            "cbioportal/clinical_data",
            "cbioportal/meta_files",
            "pvactools/features",
            "pvactools/hla_types",
            "reference_download/ensembl",
        ),
        "an R script without command line arguments": (
            "ascat/run_ascat",
            "cbioportal/generate_cna",
            "cbioportal/merge_tables",
            "cnvkit/postprocess",
            "control_freec/transform",
            "rnaqc/dupradar",
            "scarHRD/run",
            "sequenza/run",
            "signatures/deconstruct_sigs",
        ),
        "installs software": ("scarHRD/install", "sequenza/install", "singularity"),
        "a helper (scatter, gather, merge, conversion, sidecars), not a tool's main command": (
            "arcashla/prepare_reference",
            "ascat/build_baf",
            "ascat/build_cnv",
            "bcftools/TMB",
            "bcftools/concat",
            "delly2/germline/merge_calls",
            "delly2/germline/merge_genotypes",
            "gcnv/joint_germline_cnv_segmentation",
            "gcnv/merge_multikit_families",
            "melt/merge_vcf",
            "mutect2/gather",
            "picard/prepare",
            "popdel/concat_calls",
            "popdel/reorder_vcf",
            "pvactools/normalize",
            "reference_download/sidecars",
            "reference_index/common",
            "rnaqc/stats",
            "sequenza/gcreference",
            "split_regions",
        ),
        "a step of several small tools without a main command (somatic_cnv_checking)": (
            "bcftools/heterozygous_variants",
            "bcftools/pileups",
            "somatic_cnv_checking",
            "vcfpy/add_bed",
        ),
        "a secondary action; the step part's extra_arguments_action gets the arguments": (
            "cnvetti/wgs/segment",
            "cnvetti/wgs/tumor_normal_ratio",
            "cnvkit/antitarget",
            "cnvkit/call",
            "cnvkit/coverage",
            "cnvkit/export",
            "cnvkit/fix",
            "cnvkit/plot",
            "cnvkit/report",
            "cnvkit/target",
            "delly2/germline/genotype",
            "gatk4_hc/combine_gvcfs",
            "gatk4_hc/genotype",
            "gcnv/annotate_gc",
            "gcnv/annotate_gc_wgs",
            "gcnv/coverage",
            "gcnv/filter_intervals",
            "gcnv/post_germline_calls",
            "gcnv/post_germline_calls_case_mode",
            "gcnv/preprocess_intervals",
            "gcnv/preprocess_intervals_wgs",
            "gcnv/scatter_intervals",
            "melt/genotype",
            "melt/group_analysis",
            "melt/make_vcf",
            "melt/preprocess",
            "popdel/profile",
            "purecn/coverage",
            "purecn/prepare",
            "pvactools/junction",
            "scramble/cluster",
        ),
        "the tool section depends on the prediction wildcard": (
            "pvactools/netchop",
            "pvactools/netstab",
        ),
        "puts the mapper section's extra_arguments on bwa mem itself": ("mbcs",),
        "varfish export, left for later": (
            "bcftools/gvcf_to_vcf",
            "bcftools/merge_vcf",
            "mehari/annotate_seqvars",
            "mehari/annotate_strucvars",
            "mehari/bam_qc",
            "varfish_annotator/annotate",
            "varfish_annotator/annotate_svs",
        ),
    }.items()
    for wrapper in wrappers
}


def test_tool_wrappers_get_extra_arguments():
    """Every wrapper runs with the ``extra_arguments`` of its tool section, unless
    ``WRAPPERS_WITHOUT_EXTRA_ARGUMENTS`` lists it with a reason."""
    with_arguments, without = set(), set()
    for _, job in _snapshot_jobs():
        if job["wrapper"]:
            wrapper = job["wrapper"].split("/snappy_wrappers/wrappers/", 1)[1]
            args = job["params"].get("args")
            has = isinstance(args, dict) and "extra_arguments" in args
            (with_arguments if has else without).add(wrapper)
    without -= with_arguments
    listed = set(WRAPPERS_WITHOUT_EXTRA_ARGUMENTS)
    assert not without - listed, f"wrappers without extra_arguments: {sorted(without - listed)}"
    assert not listed - without, f"listed, but they get extra_arguments: {sorted(listed - without)}"


def test_wrappers_use_the_extra_arguments_they_get():
    """A rule whose params carry the tool's ``extra_arguments`` runs a wrapper that puts them on
    the command line, as ``{extra_arguments}`` (filled in by the wrapper base classes)."""
    found = set()
    for _, job in _snapshot_jobs():
        args = job["params"].get("args")
        if not job["wrapper"] or not isinstance(args, dict) or "extra_arguments" not in args:
            continue
        wrapper = job["wrapper"].split("/snappy_wrappers/wrappers/", 1)[1]
        code = "".join(
            p.read_text(encoding="utf-8")
            for p in (REPO / "snappy_wrappers" / "wrappers" / wrapper).glob("*.py")
        )
        if "{extra_arguments}" not in code:
            found.add(wrapper)
    assert not found, f"wrappers that ignore their extra_arguments: {sorted(found)}"


def test_wrappers_leave_extra_arguments_to_the_base_class():
    """An f-string that formats ``{extra_arguments}`` needs the name in the wrapper; without it,
    ``{{extra_arguments}}`` leaves the placeholder to the wrapper base classes."""
    found = []
    for path in sorted((REPO / "snappy_wrappers" / "wrappers").rglob("wrapper.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        assigned = {
            target.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Assign)
            for target in node.targets
            if isinstance(target, ast.Name)
        }
        formatted = {
            value.value.id
            for node in ast.walk(tree)
            if isinstance(node, ast.JoinedStr)
            for value in node.values
            if isinstance(value, ast.FormattedValue) and isinstance(value.value, ast.Name)
        }
        if "extra_arguments" in formatted - assigned:
            found.append(str(path.parent.relative_to(REPO / "snappy_wrappers" / "wrappers")))
    assert not found, f"f-strings that format an undefined extra_arguments: {found}"


def test_logs_and_benchmarks_are_unique_per_job():
    """Two jobs of one DAG never share a log or benchmark file."""
    found = []
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        jobs = json.loads(snapshot.read_text(encoding="utf-8"))
        paths = [p for job in jobs for p in job["log"] + [job.get("benchmark")] if p]
        found += [f"{snapshot.stem}: {p}" for p, n in collections.Counter(paths).items() if n > 1]
    assert not found, f"files shared by several jobs: {found}"


#: Rules of steps on hold; their own wrappers still write conda files and md5 files
ON_HOLD_RULES = ("create_proteome",)

_CONDA_OR_COPY = re.compile(r"\.(conda_info\.txt|conda_list\.txt|wrapper\.py|environment\.yaml)$")


def _snapshot_jobs() -> Iterator[tuple[str, dict]]:
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            if not job["rule"].startswith(ON_HOLD_RULES):
                yield snapshot.stem, job


def _snapshot_dags() -> Iterator[tuple[str, list[dict], set[str]]]:
    """Yield (snapshot, jobs, published files) per DAG snapshot, without the steps on hold.

    The published files are the inputs of the ``<task>_publish`` jobs: ``work/`` files that the
    publish job links into ``output/``.
    """
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        jobs = [
            job
            for job in json.loads(snapshot.read_text(encoding="utf-8"))
            if not job["rule"].startswith(ON_HOLD_RULES)
        ]
        published = {
            path for job in jobs if job["rule"].endswith("_publish") for path in job["input"]
        }
        yield snapshot.stem, jobs, published


def test_only_the_publish_job_writes_into_output():
    """``work/`` is the contract between tasks; ``output/`` is publication only."""
    found = {
        f"{name}: {job['rule']}: {path}"
        for name, jobs, _ in _snapshot_dags()
        for job in jobs
        if not job["rule"].endswith("_publish")
        for path in job["output"]
        if "/output/" in path
    }
    assert not found, f"jobs writing into output/: {sorted(found)[:20]}"


def test_no_job_reads_from_output():
    """Consumers read the ``work/`` paths that ``get_output_paths`` returns."""
    found = {
        f"{name}: {job['rule']}: {path}"
        for name, jobs, _ in _snapshot_dags()
        for job in jobs
        if job["rule"] != "all"
        for path in job["input"]
        if "/output/" in path and not path.startswith(("tasks/create_proteome/",))
    }
    assert not found, f"jobs reading from output/: {sorted(found)[:20]}"


def test_logs_hold_only_logs():
    """A rule's log is its log (plus a script copy); Snakemake records the rest."""
    found = {
        f"{name}: {path}"
        for name, job in _snapshot_jobs()
        for path in job["log"]
        if path.endswith(".md5") or _CONDA_OR_COPY.search(path)
    }
    assert not found, f"md5, conda or copied files declared as logs: {sorted(found)[:20]}"


def test_logs_are_not_published():
    """Logs go into the run's log archive, not into output/."""
    found = {
        f"{name}: {path}"
        for name, _, published in _snapshot_dags()
        for path in published
        if "/log/" in path
        or path.removesuffix(".md5").endswith(".log")
        or _CONDA_OR_COPY.search(path)
    }
    assert not found, f"published logs: {sorted(found)[:20]}"


def test_md5_files_only_for_published_files():
    """SODAR needs an md5 per published file; nothing reads the md5 of an intermediate file.

    The wrapper base classes write ``<file>.md5`` when it is declared next to ``<file>``. Only
    tasks with a publish job in the snapshot are checked.
    """
    found = set()
    for name, jobs, published in _snapshot_dags():
        targets = {_task_of(path) for path in published}
        for job in jobs:
            for path in job["output"]:
                if not path.endswith(".md5") or _task_of(path) not in targets:
                    continue
                payload = path.removesuffix(".md5")
                if payload not in job["output"]:
                    found.add(f"{name}: {path} (its file is no output of the job)")
                elif payload not in published:
                    found.add(f"{name}: {path} (its file is not published)")
    assert not found, f"md5 files of unpublished files: {sorted(found)[:20]} ({len(found)} total)"


def test_wrappers_read_no_md5_or_conda_files():
    found = {
        f"{name}: {job['rule']}: {path}"
        for name, job in _snapshot_jobs()
        if job.get("wrapper")
        for path in job["input"]
        if path.endswith(".md5") or _CONDA_OR_COPY.search(path)
    }
    assert not found, f"wrapper jobs reading md5 or conda files: {sorted(found)[:20]}"


#: Rules that build project-wide resources, which later samples need again: never temp
_PROJECT_RESOURCES = re.compile(
    r"_(install|prepare_reference|gcreference|write_pedigree_run|write_trio_pedigree_run|"
    r"write_panel)$"
)

#: Intermediates that stay, with the reason
NOT_TEMP = {
    r"\.scatter_intervals$": "checkpoint output that is globbed while the DAG is built",
    r"\.shards\.txt$": "checkpoint output that is read while the DAG is built",
    r"^tasks/gcnv_model\w*/work/[^/]+/out/[^/.]+\.interval_list$": "the callers read it",
    r"gatk4_hc_gvcf_combine_gvcfs\.[^/]+/out/[^/]+\.g\.vcf\.gz(\.tbi)?$": (
        "the published {library}.g.vcf.gz of the genotype rule links to it"
    ),
}


def _task_of(path: str) -> str:
    return path.split("/")[1]


def _intermediate_jobs() -> Iterator[tuple[str, dict, list[dict], set[str], dict]]:
    """Yield (snapshot, job, jobs, published files, readers) for jobs of the published tasks."""
    for name, jobs, published in _snapshot_dags():
        targets = {_task_of(path) for path in published}
        readers = collections.defaultdict(set)
        for job in jobs:
            for path in job["input"]:
                readers[path].add(job["rule"])
        rule_task = {job["rule"]: _task_of(job["output"][0]) for job in jobs if job["output"]}
        readers = {
            path: {rule_task.get(rule) for rule in rules if rule != "all"}
            for path, rules in readers.items()
        }
        for job in jobs:
            if job["output"] and _task_of(job["output"][0]) in targets:
                yield name, job, jobs, published, readers


def test_intermediates_are_temp():
    """Files that only later rules of the same task read are temp()."""
    found = set()
    for name, job, _, published, readers in _intermediate_jobs():
        if _PROJECT_RESOURCES.search(job["rule"]):
            continue
        for path in job["output"]:
            tasks = readers.get(path)
            if (
                path in job.get("temp", ())
                or "/work/" not in path
                or path.endswith(".md5")
                or path in published
                or tasks != {_task_of(path)}
                or any(re.search(pattern, path) for pattern in NOT_TEMP)
            ):
                continue
            found.add(f"{name}: {path}")
    assert not found, f"intermediates that are not temp: {sorted(found)[:20]} ({len(found)})"


def test_temp_files_are_intermediates():
    """A temp file is unpublished and has a declared reader; its index is temp, too."""
    found = set()
    for name, job, _, published, readers in _intermediate_jobs():
        temp = set(job.get("temp", ()))
        for path in temp:
            if path in published:
                found.add(f"{name}: {path} is published")
            if not readers.get(path):
                found.add(f"{name}: {path} has no reader, Snakemake deletes it right away")
            elif readers[path] != {_task_of(path)}:
                found.add(f"{name}: {path} is read by another task")
            for ext in (".tbi", ".csi", ".bai"):
                if path + ext in job["output"] and path + ext not in temp:
                    found.add(f"{name}: the index {path}{ext} is not temp")
    assert not found, f"temp files that are no intermediates: {sorted(found)[:20]} ({len(found)})"


#: Steps that publish without md5 files: steps whose results are not uploaded to SODAR
#: (reference indices, gCNV models, cBioPortal files)
UNCHECKED_PUBLISHING = (
    "reference_index",
    "gcnv_model",
    "cbioportal_export",
)


def test_published_files_come_with_md5_and_index():
    """SODAR needs an md5 for every published file; bgzipped VCF and BED files need their index.

    ``.done`` markers of directories have none.
    """
    found = set()
    for name, _, published in _snapshot_dags():
        for path in published:
            if path.endswith((".md5", ".done")) or path.startswith(
                tuple(f"tasks/{step}" for step in UNCHECKED_PUBLISHING)
            ):
                continue
            companions = [path + ".md5"]
            if re.search(r"\.(vcf|bed)\.gz$", path):
                companions += [path + ".tbi", path + ".tbi.md5"]
            found.update(f"{name}: {c}" for c in companions if c not in published)
    assert not found, f"published files without md5 or index: {sorted(found)[:20]} ({len(found)})"


#: Wrappers on hold, not pinned yet
UNPINNED_WRAPPERS = {"pvactools/create_proteome"}

#: Environments that conda cannot solve as written (unavailable versions, a defunct channel,
#: no dependencies); each needs new versions before it can be pinned
UNSOLVABLE_ENVS = {
    "ascat/run_ascat",
    "cbioportal/case_lists",
    "cbioportal/clinical_data",
    "cbioportal/meta_files",
    "gatk_phase_by_transmission",
    "gatk_read_backed_phasing",
    "star_fusion",
}

_EXACT_DEPENDENCY = re.compile(r"^(?:[\w-]+::)?([\w.-]+)==([\w.]+)$")


def test_wrapper_envs_are_pinned():
    """Snakemake builds a wrapper environment from its pin file, so versions never drift.

    ``snakedeploy pin-conda-envs`` writes ``environment.linux-64.pin.txt``; a symlinked
    ``environment.yaml`` needs a symlink to its target's pin file.
    """
    import yaml

    wrappers = REPO / "snappy_wrappers" / "wrappers"
    found = []
    for env in sorted(wrappers.rglob("environment.yaml")):
        name = str(env.parent.relative_to(wrappers))
        if name in UNPINNED_WRAPPERS or name in UNSOLVABLE_ENVS:
            continue
        pin = env.with_name("environment.linux-64.pin.txt")
        if not pin.exists() or pin.resolve() != env.resolve().with_name(pin.name):
            found.append(f"{name}: no pin file of its environment next to it")
            continue
        pinned = pin.read_text(encoding="utf-8")
        for dependency in yaml.safe_load(env.read_text(encoding="utf-8"))["dependencies"]:
            if not isinstance(dependency, str):
                continue
            package = re.split(r"[=<>!~ ]", dependency.split("::")[-1], maxsplit=1)[0]
            # <channel>/<subdir>/<name>-<version>-<build>.conda; names may contain "-"
            if not re.search(rf"/{re.escape(package)}-[^/-]+-[^/-]+\.(conda|tar\.bz2)", pinned):
                found.append(f"{name}: the pin file does not hold {package}")
                continue
            m = _EXACT_DEPENDENCY.match(dependency.replace(" ", ""))
            # conda reads ==4.0 as 4.0.0
            if m and not re.search(rf"/{re.escape(m[1])}-{re.escape(m[2])}(\.[\w.]+)?-", pinned):
                found.append(f"{name}: the pin file does not hold {dependency}")
    assert not found, f"wrapper environments to pin again: {found}"
