# -*- coding: utf-8 -*-
"""Checks for the coding conventions in ``docs/dev_conventions.rst``."""

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

#: Path convention violations that already exist (plans.md C3). The test fails for new ones and
#: for fixed ones, so the list only shrinks.
PATH_BASELINE = SNAPSHOTS / "path_convention_violations.json"

#: Input functions that a Snakefile calls with Snakemake's ``checkpoints`` object.
CHECKPOINT_INPUT_FUNCTIONS = {
    ("BuildGcnvTargetSeqModelStepPart", "_get_input_files_post_germline_calls"),
    ("BuildGcnvWgsModelStepPart", "_get_input_files_post_germline_calls"),
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
            if job["rule"] == "all":
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
    """No tool names in paths and ``{entity}.{suffix}`` naming (plans.md C3)."""
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
    "expression",
    "hla_types",
    "strandedness",
    "panel_of_normals",
    "reference",
    "features",
    "dbsnp",
    "index",
}


def test_depends_on_keys_use_the_vocabulary():
    found = []
    for step, workflow_cls in WORKFLOW_REGISTRY.items():
        depends_on = workflow_cls.config_model_class.model_fields.get("depends_on")
        if depends_on is not None:
            keys = set(depends_on.annotation.model_fields) - DEPENDS_ON_KEYS
            found += [f"{step}.{key}" for key in sorted(keys)]
    assert not found, f"depends_on keys outside the vocabulary in docs/dev_conventions.rst: {found}"


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


#: Keys a wrapper reads from ``params.args`` only in a case that the snapshot jobs do not take.
WRAPPER_PARAM_EXCEPTIONS = {
    ("cnvkit/report", "breaks"): "read only for a breaks output",
    ("cnvkit/report", "genemetrics"): "read only for a genemetrics output",
    ("cnvkit/report", "segmetrics"): "read only for a segmetrics output",
    ("mbcs", "barcode_config"): "read only with mbcs.use_barcodes",
    ("vembrane/tag", "expression"): "read only in filter mode",
}

_ARGS_LITERAL = re.compile(r"""(?:^|[^\w.]|snakemake\.params\.)args\[\s*["'](\w+)["']\s*\]""", re.M)
_ARGS_TEMPLATE = re.compile(r"\{(?:snakemake\.params\.)?args\[(\w+)\]")
_ARGS_GUARD = re.compile(
    r"""["'](\w+)["']\s+(?:not\s+)?in\s+(?:snakemake\.params\.)?args\b"""
    r"""|(?:snakemake\.params\.)?args\.get\(\s*["'](\w+)["']"""
)


def _wrapper_args_keys(wrapper: str) -> set[str]:
    """Return the ``params.args`` keys that a wrapper reads without a guard or default."""
    keys, guarded = set(), set()
    for path in (REPO / "snappy_wrappers" / "wrappers" / wrapper).glob("*.py"):
        text = path.read_text(encoding="utf-8")
        keys |= set(_ARGS_LITERAL.findall(text)) | set(_ARGS_TEMPLATE.findall(text))
        guarded |= {quoted or got for quoted, got in _ARGS_GUARD.findall(text)}
    return keys - guarded


def test_wrappers_get_the_params_they_read():
    """Every ``args[...]`` key a wrapper reads is in the params of the snapshot jobs that run it."""
    found, used_exceptions = set(), set()
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            wrapper = (job.get("wrapper") or "").partition("snappy_wrappers/wrappers/")[2]
            args = (job.get("params") or {}).get("args", {})
            if not wrapper or not isinstance(args, dict):  # params not known when dumping
                continue
            for key in sorted(_wrapper_args_keys(wrapper) - set(args)):
                if (wrapper, key) in WRAPPER_PARAM_EXCEPTIONS:
                    used_exceptions.add((wrapper, key))
                else:
                    found.add(f"{wrapper} reads args[{key}], missing in {job['rule']}")
    assert not found, f"wrappers read params that their rules do not pass: {sorted(found)}"
    stale = sorted(set(WRAPPER_PARAM_EXCEPTIONS) - used_exceptions)
    assert not stale, f"remove these entries from WRAPPER_PARAM_EXCEPTIONS: {stale}"


#: Params that hold file paths (plans.md F11). The test fails for new ones and for fixed ones, so
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
WRAPPERS_WITHOUT_BASE_CLASS = {
    name: "on hold with somatic_neoepitope_prediction"
    for name in (
        "pvactools/combine",
        "pvactools/create_proteome",
        "pvactools/junction",
        "pvactools/netchop",
        "pvactools/normalize",
        "pvactools/phasing",
        "pvactools/pileup",
        "pvactools/pvacfuse",
        "pvactools/pvacseq",
        "pvactools/pvacsplice",
    )
}


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


#: Wrapper rules without a benchmark: their wrappers do not exist (legacy)
RULES_WITHOUT_BENCHMARK = {
    "gene_expression_report_aggregate_counts",
    "gene_expression_report_compute_ranks",
    "gene_expression_report_compute_signatures",
    "gene_expression_report_plot_expression_distribution",
}


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
    """Benchmarks record wall time and memory per job to calibrate resources (plans.md F7).

    Snakemake does not allow benchmarks for rules that are cached between workflows.
    """
    snakefiles = sorted(WORKFLOWS.glob("*/Snakefile")) + sorted(WORKFLOWS.glob("*/*.rules"))
    found = [
        f"{path.relative_to(WORKFLOWS)}: {name}"
        for path in snakefiles
        for name, text in _rule_blocks(path)
        if "wrapper:" in text
        and "benchmark:" not in text
        and "cache:" not in text
        and name not in RULES_WITHOUT_BENCHMARK
    ]
    assert not found, f"wrapper rules without benchmark: wf.get_benchmark_file(...): {found}"


def test_logs_and_benchmarks_are_unique_per_job():
    """Two jobs of one DAG never share a log or benchmark file."""
    found = []
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        jobs = json.loads(snapshot.read_text(encoding="utf-8"))
        paths = [p for job in jobs for p in job["log"] + [job.get("benchmark")] if p]
        found += [f"{snapshot.stem}: {p}" for p, n in collections.Counter(paths).items() if n > 1]
    assert not found, f"files shared by several jobs: {found}"


#: Rules of steps on hold (plans.md); their own wrappers still write conda files and md5 files
ON_HOLD_RULES = ("somatic_neoepitope_prediction", "create_proteome")

_CONDA_OR_COPY = re.compile(r"\.(conda_info\.txt|conda_list\.txt|wrapper\.py|environment\.yaml)$")


def _snapshot_jobs() -> Iterator[tuple[str, dict]]:
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        for job in json.loads(snapshot.read_text(encoding="utf-8")):
            if not job["rule"].startswith(ON_HOLD_RULES):
                yield snapshot.stem, job


def test_logs_hold_only_logs():
    """A rule's log is its log (plus a script copy); Snakemake records the rest (plans.md W2)."""
    found = {
        f"{name}: {path}"
        for name, job in _snapshot_jobs()
        for path in job["log"]
        if path.endswith(".md5") or _CONDA_OR_COPY.search(path)
    }
    assert not found, f"md5, conda or copied files declared as logs: {sorted(found)[:20]}"


def test_logs_are_not_published():
    """Logs go into the run's log archive, not into output/ (plans.md F3, W2)."""
    found = {
        f"{name}: {path}"
        for name, job in _snapshot_jobs()
        for path in job["output"]
        if "/output/" in path
        and (
            "/log/" in path
            or path.removesuffix(".md5").endswith(".log")
            or _CONDA_OR_COPY.search(path)
        )
    }
    assert not found, f"logs linked into output/: {sorted(found)[:20]}"


def test_md5_files_only_for_published_files():
    """SODAR needs an md5 per published file; nothing reads the md5 of an intermediate file.

    The wrapper base classes write ``<file>.md5`` when it is declared next to ``<file>``. Only
    the target tasks of a snapshot link their files into output/, so only they are checked. A
    link may sit in another directory than its file, so file names are compared. The orchestrator
    does not target ``.tbi`` files, so an index counts as published with the file it indexes.
    """
    found = set()
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        jobs = [
            job
            for job in json.loads(snapshot.read_text(encoding="utf-8"))
            if not job["rule"].startswith(ON_HOLD_RULES)
        ]
        targets = {
            path.split("/")[1] for job in jobs if job["rule"] == "all" for path in job["input"]
        }
        published = {
            (path.split("/")[1], path.rsplit("/", 1)[1])
            for job in jobs
            for path in job["output"]
            if "/output/" in path
        }
        for job in jobs:
            for path in job["output"]:
                if not path.endswith(".md5") or "/work/" not in path:
                    continue
                if path.split("/")[1] not in targets:
                    continue
                payload = path.removesuffix(".md5")
                if payload not in job["output"]:
                    found.add(f"{snapshot.stem}: {path} (its file is no output of the job)")
                elif (path.split("/")[1], payload.rsplit("/", 1)[1].removesuffix(".tbi")) not in (
                    published
                ):
                    found.add(f"{snapshot.stem}: {path} (its file is not published)")
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


#: Steps that publish without md5 files: a legacy step whose wrappers do not exist, and the
#: reference indices, which are not uploaded to SODAR
UNCHECKED_PUBLISHING = ("gene_expression_report", "reference_index")


def test_published_files_come_with_md5_and_index():
    """SODAR needs an md5 for every published file; bgzipped VCF and BED files need their index.

    ``PublishStepPart`` links a file together with these companions. ``.done`` markers of
    directories have none.
    """
    found = set()
    for snapshot in sorted((SNAPSHOTS / "dag").glob("*.json")):
        jobs = [
            job
            for job in json.loads(snapshot.read_text(encoding="utf-8"))
            if not job["rule"].startswith(ON_HOLD_RULES + UNCHECKED_PUBLISHING)
        ]
        targets = {
            path.split("/")[1] for job in jobs if job["rule"] == "all" for path in job["input"]
        }
        published = {
            path
            for job in jobs
            for path in job["output"]
            if "/output/" in path and path.split("/")[1] in targets
        }
        for path in published:
            if path.endswith((".md5", ".done")):
                continue
            companions = [path + ".md5"]
            if re.search(r"\.(vcf|bed)\.gz$", path):
                companions += [path + ".tbi", path + ".tbi.md5"]
            found.update(f"{snapshot.stem}: {c}" for c in companions if c not in published)
    assert not found, f"published files without md5 or index: {sorted(found)[:20]} ({len(found)})"


#: Wrapper directories on hold (plans.md), not pinned yet
UNPINNED_WRAPPERS = {"pvactools"}

#: Environments that conda cannot solve as written (unavailable versions, a defunct channel,
#: no dependencies); each needs new versions before it can be pinned
UNSOLVABLE_ENVS = {
    "ascat/run_ascat",
    "cbioportal/case_lists",
    "cbioportal/clinical_data",
    "cbioportal/meta_files",
    "gatk_phase_by_transmission",
    "gatk_read_backed_phasing",
    "hera/quant",
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
        if env.relative_to(wrappers).parts[0] in UNPINNED_WRAPPERS or name in UNSOLVABLE_ENVS:
            continue
        pin = env.with_name("environment.linux-64.pin.txt")
        if not pin.exists() or pin.resolve() != env.resolve().with_name(pin.name):
            found.append(f"{name}: no pin file of its environment next to it")
            continue
        pinned = pin.read_text(encoding="utf-8")
        for dependency in yaml.safe_load(env.read_text(encoding="utf-8"))["dependencies"]:
            m = isinstance(dependency, str) and _EXACT_DEPENDENCY.match(dependency.replace(" ", ""))
            # conda reads ==4.0 as 4.0.0
            if m and not re.search(rf"/{re.escape(m[1])}-{re.escape(m[2])}(\.[\w.]+)?-", pinned):
                found.append(f"{name}: the pin file does not hold {dependency}")
    assert not found, f"wrapper environments to pin again: {found}"
