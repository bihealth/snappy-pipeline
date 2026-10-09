from __future__ import annotations

import difflib
import json
import os
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from enum import Enum
from pathlib import Path
from typing import Any

import pytest
import yaml

from tests.scripts.generate_task_configs import build_all_tasks, load_yaml


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _get_task_names() -> list[str]:
    root = _repo_root()
    base_config_path = root / "tests/snappy_pipeline/fixtures/base_config.yaml"
    base_config = load_yaml(base_config_path)
    # Redirect stdout to suppress print logs during test collection
    import io
    import sys as sys_orig

    f = io.StringIO()
    sys_orig.stdout = f
    try:
        tasks = build_all_tasks(base_config, base_config_path)
    finally:
        sys_orig.stdout = sys_orig.__stdout__
    return [t["name"] for t in tasks]


TASK_NAMES = _get_task_names()


RNA_TASKS = {
    "ngs_mapping_star",
    "somatic_gene_fusion_calling_arriba",
    "somatic_gene_fusion_calling_fusioncatcher",
    "somatic_gene_fusion_calling_jaffa",
    "somatic_gene_fusion_calling_defuse",
    "somatic_gene_fusion_calling_hera",
    "somatic_gene_fusion_calling_pizzly",
    "somatic_gene_fusion_calling_star_fusion",
    "gene_expression_quantification_strandedness",
    "gene_expression_quantification_featurecounts",
    "gene_expression_quantification_dupradar",
    "gene_expression_quantification_duplication",
    "gene_expression_quantification_rnaseqc",
    "gene_expression_quantification_salmon",
    "gene_expression_quantification_stats",
    "gene_expression_report",
    "somatic_neoepitope_prediction_pvacseq",
    "somatic_neoepitope_prediction_pvacfuse",
    "somatic_neoepitope_prediction_pvacsplice",
    "create_proteome",
    "hla_typing_arcashla",
}

#: Tasks that need a germline trio: germline callers and their pedigree-based consumers.
GERMLINE_TASKS = {
    "variant_calling_bcftools_call",
    "variant_calling_gatk3_hc",
    "variant_calling_gatk3_ug",
    "variant_calling_gatk4_hc_gvcf",
    "variant_calling_gatk4_hc_joint",
    "variant_phasing",
    "igv_session_generation",
}


#: Closures whose DAG cannot be built yet, with the reason. Strict xfail, so a fix shows up.
KNOWN_BROKEN = {
    "igv_session_generation": "still builds tool-prefixed upstream paths (plans.md K3, C3)",
}

_GERMLINE_CALLER_RESULTS = (
    "germline callers declare only work/ outputs, and get_result_files keeps only output/ paths "
    "(plans.md W1)"
)

#: Closures that target no files, with the reason. Every other closure must target at least one
#: file, otherwise its snapshot checks nothing.
EXPECTED_EMPTY = {
    "external_data": "provides existing files; has no rules",
    "external_vcf": "provides existing files; has no rules",
    "external_cnv": "provides existing files; has no rules",
    "external_sv": "provides existing files; has no rules",
    "ngs_mapping_minimap2": "minimap2 maps only long-read libraries; there is no long-read fixture",
    "variant_calling_bcftools_call": _GERMLINE_CALLER_RESULTS,
    "variant_calling_gatk3_hc": _GERMLINE_CALLER_RESULTS,
    "variant_calling_gatk3_ug": _GERMLINE_CALLER_RESULTS,
    "variant_calling_gatk4_hc_gvcf": _GERMLINE_CALLER_RESULTS,
    "variant_calling_gatk4_hc_joint": _GERMLINE_CALLER_RESULTS,
}


def _fixture_dir() -> Path:
    return _repo_root() / "tests" / "snappy_pipeline" / "fixtures"


def _task_sample_sheet(task_name: str) -> Path:
    if task_name in RNA_TASKS:
        return _fixture_dir() / "samplesheet_rna.tsv"
    if task_name in GERMLINE_TASKS:
        return _fixture_dir() / "samplesheet_germline.tsv"
    return _fixture_dir() / "samplesheet.tsv"


def _task_raw_folders(task_name: str) -> tuple[str, ...]:
    if task_name in GERMLINE_TASKS:
        return (
            "case001subregion-N1-DNA1-WES1",
            "case001subregionFather-N1-DNA1-WES1",
            "case001subregionMother-N1-DNA1-WES1",
        )
    if task_name in RNA_TASKS:
        return (
            "case001subregion-N1-DNA1-WES1",
            "case001subregion-T1-DNA1-WES1",
            "case001subregion-T1-RNA1-mRNA_seq1",
        )
    return ("case001subregion-N1-DNA1-WES1", "case001subregion-T1-DNA1-WES1")


#: Committed DAG snapshots, one JSON file per generated task.
SNAPSHOT_DIR = Path(__file__).resolve().parent / "snapshots" / "dag"


def _run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    # The default slurm partition ends up in job resources; keep snapshots independent of the
    # caller's environment.
    env["SNAPPY_PIPELINE_PARTITION"] = "medium"
    if "PYTHONPATH" not in env:
        env["PYTHONPATH"] = str(_repo_root())
    else:
        env["PYTHONPATH"] = str(_repo_root()) + os.pathsep + env["PYTHONPATH"]
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False, env=env)


def _tail(text: str, n: int = 40) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def to_plain_obj(obj: Any) -> Any:
    """Convert ruamel/pydantic/path-like values to plain YAML-safe builtins."""
    if obj is None:
        return None
    if isinstance(obj, bool):
        return bool(obj)
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return float(obj)
    if isinstance(obj, str):
        return str(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, Enum):
        return to_plain_obj(obj.value)
    if isinstance(obj, Mapping):
        return {str(k): to_plain_obj(v) for k, v in obj.items()}
    if isinstance(obj, Sequence) and not isinstance(obj, (str, bytes, bytearray)):
        return [to_plain_obj(v) for v in obj]
    if isinstance(obj, set):
        return [to_plain_obj(v) for v in sorted(obj, key=lambda x: str(x))]
    return str(obj)


@pytest.fixture(scope="session")
def generated_task_config(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    root = _repo_root()
    out_dir = tmp_path_factory.mktemp("generated-task-configs")

    gen = _run(
        [sys.executable, "tests/scripts/generate_task_configs.py", "--out-dir", str(out_dir)],
        cwd=root,
    )
    gen_output = (gen.stdout or "") + "\n" + (gen.stderr or "")
    if gen.returncode != 0:
        raise AssertionError(
            f"generate_task_configs.py failed\nstdout/stderr tail:\n{_tail(gen_output)}"
        )

    return {
        "root": root,
        "out_dir": out_dir,
        "config_path": out_dir / "all-workflows" / "config.yaml",
        "generate_output": gen_output,
    }


def dependency_closure(task_name: str, tasks_by_name: dict[str, dict[str, Any]]) -> set[str]:
    seen = set()
    stack = [task_name]
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        task = tasks_by_name.get(current)
        if not task:
            continue
        task_config = task.get("config", {}) if isinstance(task, dict) else {}
        depends_on = task_config.get("depends_on", {}) if isinstance(task_config, dict) else {}
        if not isinstance(depends_on, dict):
            continue
        for dep_task_name in depends_on.values():
            if isinstance(dep_task_name, str) and dep_task_name and dep_task_name not in seen:
                stack.append(dep_task_name)
    return seen


def _write_closure_project(
    task_name: str, generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    """Write the config of ``task_name`` and its upstream tasks, with dummy FASTQs, to tmp_path."""
    config_path = generated_task_config["config_path"]

    # Load config.yaml
    config = load_yaml(config_path)
    tasks = config.get("tasks", [])
    tasks_by_name = {t["name"]: t for t in tasks if isinstance(t, dict) and "name" in t}

    # Calculate closure
    closure = dependency_closure(task_name, tasks_by_name)
    tasks_subset = [t for t in tasks if isinstance(t, dict) and t.get("name") in closure]

    # Construct closure config
    closure_config = {
        "static_data_config": config.get("static_data_config", {}),
        "tasks": tasks_subset,
        "data_sets": config.get("data_sets", {}),
    }

    # Redirect search_paths to a temporary raw directory inside the tmp_path
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    if "data_sets" in closure_config:
        for ds_name, ds_config in closure_config["data_sets"].items():
            if isinstance(ds_config, dict):
                ds_config["file"] = str(_task_sample_sheet(task_name))
                if task_name in GERMLINE_TASKS:
                    ds_config["type"] = "germline_variants"
                ds_config["search_paths"] = [str(raw_dir)]
    for task in closure_config["tasks"]:
        if task["step"] == "external_data":
            task["config"]["search_paths"] = [str(raw_dir)]

    # Touch dummy FASTQ files, and VCFs for the external_data tasks
    for folder in _task_raw_folders(task_name):
        f_dir = raw_dir / folder
        f_dir.mkdir(parents=True, exist_ok=True)
        for suffix in (".R1.fastq.gz", ".R2.fastq.gz", ".vcf.gz", ".vcf.gz.tbi"):
            (f_dir / f"{folder}{suffix}").touch()
        (raw_dir / f"{folder}.R1.fastq.gz").touch()
        (raw_dir / f"{folder}.R2.fastq.gz").touch()

    # Write config.yaml directly in tmp_path (no .snappy_pipeline subfolder!)
    closure_config_path = tmp_path / "config.yaml"
    closure_config_plain = to_plain_obj(closure_config)
    with closure_config_path.open("wt", encoding="utf-8") as f:
        yaml.safe_dump(closure_config_plain, f, sort_keys=False)

    # Ensure we emitted parseable YAML before invoking snappy/snakemake.
    try:
        reloaded = yaml.safe_load(closure_config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        raise AssertionError(
            f"Generated closure config is invalid YAML for {task_name}: {e}"
        ) from e
    assert isinstance(reloaded, dict), f"Generated closure config is not a mapping for {task_name}"


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize(
    "task_name",
    [
        pytest.param(name, marks=pytest.mark.xfail(reason=KNOWN_BROKEN[name], strict=True))
        if name in KNOWN_BROKEN
        else name
        for name in TASK_NAMES
    ],
    ids=TASK_NAMES,
)
def test_generated_config_task_closure_passes(
    task_name: str, generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    root = generated_task_config["root"]
    _write_closure_project(task_name, generated_task_config, tmp_path)

    # Build the DAG that `snappy run --task <task_name>` would build. A broken DAG fails here,
    # like a dry-run would; the job list (paths, params, resources, wrappers) must match the
    # committed snapshot.
    dump_path = tmp_path / "dag.json"
    cmd = [
        sys.executable,
        "tests/scripts/dump_dag.py",
        "--directory",
        str(tmp_path),
        "--task",
        task_name,
        "--output",
        str(dump_path),
    ]
    dump = _run(cmd, cwd=root)
    dump_output = (dump.stdout or "") + "\n" + (dump.stderr or "")

    assert dump.returncode == 0, (
        f"building the DAG failed for {task_name}\nstdout/stderr excerpt:\n{dump_output}"
    )
    actual = dump_path.read_text(encoding="utf-8")
    targets = next(job["input"] for job in json.loads(actual) if job["rule"] == "all")
    if task_name in EXPECTED_EMPTY:
        assert not targets, f"{task_name} now targets files; remove it from EXPECTED_EMPTY"
    else:
        assert targets, f"{task_name} targets no files, so its snapshot would check nothing"
    for path, placeholder in ((str(tmp_path), "<project>"), (str(root), "<repo>")):
        actual = actual.replace(path, placeholder)
    _check_snapshot(task_name, actual)


@pytest.mark.integration
def test_frozen_run_uses_existing_upstream_outputs(
    generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    """``snappy run --task X --frozen`` loads only X's rules; missing upstream files fail fast."""
    _write_closure_project("variant_annotation_vep", generated_task_config, tmp_path)
    cmd = [sys.executable, "tests/scripts/dump_dag.py", "--directory", str(tmp_path)]
    cmd += ["--task", "variant_annotation_vep", "--output", str(tmp_path / "dag.json"), "--frozen"]
    dump = _run(cmd, cwd=generated_task_config["root"])
    output = (dump.stdout or "") + (dump.stderr or "")

    assert dump.returncode != 0
    assert "MissingInputException" in output
    assert "tasks/variant_calling_gatk4_hc_gvcf/output/" in output


def _mapping_job_with_reads_from(
    source_task: dict[str, Any], generated_task_config: dict[str, Any], tmp_path: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return the first BWA job and all jobs when ngs_mapping_bwa reads from ``source_task``."""
    _write_closure_project("ngs_mapping_bwa", generated_task_config, tmp_path)
    config = yaml.safe_load((tmp_path / "config.yaml").read_text(encoding="utf-8"))
    config["tasks"] = [t for t in config["tasks"] if t["name"] != source_task["name"]]
    config["tasks"].insert(0, source_task)
    mapping = next(t for t in config["tasks"] if t["name"] == "ngs_mapping_bwa")
    mapping["config"]["depends_on"]["reads"] = source_task["name"]
    (tmp_path / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    dump_path = tmp_path / "dag.json"
    cmd = [sys.executable, "tests/scripts/dump_dag.py", "--directory", str(tmp_path)]
    cmd += ["--task", "ngs_mapping_bwa", "--output", str(dump_path)]
    dump = _run(cmd, cwd=generated_task_config["root"])
    assert dump.returncode == 0, _tail((dump.stdout or "") + (dump.stderr or ""))
    jobs = json.loads(dump_path.read_text(encoding="utf-8"))
    return next(j for j in jobs if j["rule"].endswith("_bwa_run")), jobs


@pytest.mark.integration
def test_mapping_reads_trimmed_fastqs(generated_task_config: dict[str, Any], tmp_path: Path):
    trimming = {
        "step": "adapter_trimming",
        "name": "trimming",
        "config": {"depends_on": {"reads": "data_sets"}, "tool": "fastp", "fastp": {}},
    }
    job, jobs = _mapping_job_with_reads_from(trimming, generated_task_config, tmp_path)
    lib = job["wildcards"]["library_name"]

    assert job["input"] == [f"tasks/trimming/output/{lib}/out/.done"]
    assert job["params"]["args"]["input"]["reads_left"] == [
        f"tasks/trimming/output/{lib}/out/{lib}.R1.fastq.gz"
    ]
    assert any(j["rule"].startswith("trimming_adapter_trimming_fastp") for j in jobs)


@pytest.mark.integration
def test_mapping_reads_from_an_external_data_task(
    generated_task_config: dict[str, Any], tmp_path: Path
):
    external_reads = {
        "step": "external_data",
        "name": "external_reads",
        "config": {
            "produces": {"type": "raw"},
            "search_paths": [str(tmp_path / "raw")],
            "search_patterns": [
                {
                    "left": r"(?P<readgroup>.+)\.R1\.fastq\.gz",
                    "right": r"(?P<readgroup>.+)\.R2\.fastq\.gz",
                }
            ],
        },
    }
    job, _ = _mapping_job_with_reads_from(external_reads, generated_task_config, tmp_path)
    lib = job["wildcards"]["library_name"]

    assert job["params"]["args"]["input"]["reads_left"] == [
        str(tmp_path / "raw" / lib / f"{lib}.R1.fastq.gz")
    ]
    assert str(tmp_path / "raw" / lib / f"{lib}.R1.fastq.gz") in job["input"]


def _check_snapshot(task_name: str, actual: str) -> None:
    """Compare a normalized DAG dump with its snapshot; rewrite it if SNAPPY_UPDATE_SNAPSHOTS is set."""
    snapshot_path = SNAPSHOT_DIR / f"{task_name}.json"
    if os.environ.get("SNAPPY_UPDATE_SNAPSHOTS"):
        snapshot_path.parent.mkdir(parents=True, exist_ok=True)
        snapshot_path.write_text(actual, encoding="utf-8")
        return

    assert snapshot_path.exists(), (
        f"No DAG snapshot for {task_name}; create it with SNAPPY_UPDATE_SNAPSHOTS=1"
    )
    expected = snapshot_path.read_text(encoding="utf-8")
    if actual != expected:
        diff = difflib.unified_diff(
            expected.splitlines(keepends=True),
            actual.splitlines(keepends=True),
            f"snapshots/dag/{task_name}.json",
            "actual",
            n=2,
        )
        diff_head = "".join(list(diff)[:150])
        raise AssertionError(
            f"DAG of {task_name} differs from its snapshot. If the change is intended, "
            f"update it with SNAPPY_UPDATE_SNAPSHOTS=1.\n{diff_head}"
        )


@pytest.mark.integration
def test_generated_config_audit_regression_guard(generated_task_config: dict[str, Any]) -> None:
    gen_output = generated_task_config["generate_output"]
    unresolved = re.findall(r"unresolved validation errors, kept best-effort config", gen_output)

    cfg_path = generated_task_config["config_path"]
    cfg: dict[str, Any] = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}

    auto_count = 0

    def walk(v: Any) -> None:
        nonlocal auto_count
        if isinstance(v, dict):
            for vv in v.values():
                walk(vv)
        elif isinstance(v, list):
            for vv in v:
                walk(vv)
        elif isinstance(v, str) and v == "AUTO":
            auto_count += 1

    walk(cfg)

    # Guard against silent regressions: keep generated-config quality from getting worse.
    assert len(unresolved) <= int(os.environ.get("SNAPPY_MAX_UNRESOLVED_CONFIGS", "35"))
    assert auto_count <= int(os.environ.get("SNAPPY_MAX_AUTO_PLACEHOLDERS", "70"))
