from __future__ import annotations

import difflib
import json
import os
import re
import subprocess
import sys
from collections.abc import Callable, Mapping, Sequence
from enum import Enum
from pathlib import Path
from typing import Any

import pytest
import yaml

from snappy_pipeline.models import reader_markers
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from tests.scripts.generate_task_configs import build_all_tasks, build_config, dump_yaml, load_yaml


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


BASE_CONFIG_PATH = _repo_root() / "tests/snappy_pipeline/fixtures/base_config.yaml"
BASE_CONFIG = load_yaml(BASE_CONFIG_PATH)

#: The generated tasks, built once per process: the collection needs their names, the session
#: fixture writes them.
TASKS = build_all_tasks(BASE_CONFIG, BASE_CONFIG_PATH)
TASK_NAMES = [task["name"] for task in TASKS]


RNA_TASKS = {
    "ngs_mapping_star",
    "somatic_gene_fusion_calling_arriba",
    "somatic_gene_fusion_calling_star_fusion",
    "gene_expression_quantification_strandedness",
    "gene_expression_quantification_featurecounts",
    "gene_expression_quantification_dupradar",
    "gene_expression_quantification_duplication",
    "gene_expression_quantification_rnaseqc",
    "gene_expression_quantification_salmon",
    "gene_expression_quantification_stats",
    "somatic_neoepitope_prediction_pvactools",
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
    "variant_calling_freebayes_germline",
    "variant_calling_varlociraptor_germline",
    "variant_calling_varlociraptor_sv",
    "variant_filtration_varlociraptor_germline",
    "variant_report_datavzrd",
    "variant_report_igv",
    "variant_phasing_gatk_phase_by_transmission",
    "variant_phasing_gatk_read_backed_phasing",
    "igv_session_generation",
}


#: Closures whose DAG cannot be built yet, with the reason. Strict xfail, so a fix shows up.
KNOWN_BROKEN: dict[str, str] = {}

#: Closures that target no files, with the reason. Every other closure must target at least one
#: file, otherwise its snapshot checks nothing.
EXPECTED_EMPTY = {
    "external_data": "provides existing files; has no rules",
    "reads": "provides existing files; has no rules",
    "genome": "provides existing files; has no rules",
    "features": "provides existing files; has no rules",
    "dbsnp": "provides existing files; has no rules",
    "exons": "provides existing files; has no rules",
    "kits": "provides existing files; has no rules",
    "regions": "provides existing files; has no rules",
    "tagged_regions": "provides existing files; has no rules",
    "mappable_regions": "provides existing files; has no rules",
    "population_variants": "provides existing files; has no rules",
    "umis": "provides existing files; has no rules",
    "salmon_index": "provides existing files; has no rules",
    "star_index": "provides existing files; has no rules",
    "proteome": "provides existing files; has no rules",
    "vep_cache": "provides existing files; has no rules",
    "mehari_transcripts": "provides existing files; has no rules",
    "mehari_db": "provides existing files; has no rules",
    "arriba_database": "provides existing files; has no rules",
    "star_fusion_database": "provides existing files; has no rules",
    "control_freec_database": "provides existing files; has no rules",
    "expansionhunter_database": "provides existing files; has no rules",
    "external_vcf": "provides existing files; has no rules",
    "external_cnv": "provides existing files; has no rules",
    "external_sv": "provides existing files; has no rules",
    "ngs_mapping_minimap2": "minimap2 maps only long-read libraries; there is no long-read fixture",
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


def _dump(project: Path, task: str, frozen: bool = False) -> list[dict[str, Any]]:
    """Return the jobs of the DAG that ``snappy run --task <task>`` builds in ``project``."""
    output = project / "dag.json"
    cmd = [sys.executable, "tests/scripts/dump_dag.py", "--directory", str(project)]
    cmd += ["--task", task, "--output", str(output)] + (["--frozen"] if frozen else [])
    dump = _run(cmd, cwd=_repo_root())
    assert dump.returncode == 0, f"building the DAG of {task} failed\n{dump.stdout}\n{dump.stderr}"
    return json.loads(output.read_text(encoding="utf-8"))


def _edit_config(project: Path, edit: Callable[[dict[str, Any]], None]) -> dict[str, Any]:
    """Apply ``edit`` to the project's ``config.yaml`` and return the edited config."""
    path = project / "config.yaml"
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    edit(config)
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config


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
    config_path = tmp_path_factory.mktemp("generated-task-configs") / "config.yaml"
    dump_yaml(config_path, build_config(BASE_CONFIG, TASKS, BASE_CONFIG_PATH))
    return {"root": _repo_root(), "config_path": config_path}


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
        for value in depends_on.values():
            for dep_task_name in value if isinstance(value, list) else [value]:
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
        "project": config["project"],
        "sample_sheets": [str(_task_sample_sheet(task_name))],
        "tasks": tasks_subset,
    }

    # The files tables of the external_data tasks list dummy files below tmp_path/raw
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    for task in closure_config["tasks"]:
        if task["step"] == "external_data" and "library_files" in task["config"]:
            reads = task["config"]["produces"]["type"] == "reads"
            task["config"]["library_files"] = [
                str(tmp_path / ("reads.tsv" if reads else "vcfs.tsv"))
            ]
        # The snapshots record every declared intermediate as temp, not only scatter fragments
        task["config"]["delete_intermediates"] = ["intermediates"]

    # Touch dummy FASTQ files and VCFs, and list them in the files tables
    _write_files_tables(tmp_path, _task_raw_folders(task_name))

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


def _write_files_tables(project_dir: Path, folders: tuple[str, ...], append: bool = False) -> None:
    """Touch the dummy FASTQs and VCFs of the libraries ``folders`` below ``project_dir/raw`` and
    list them in the files tables ``reads.tsv`` and ``vcfs.tsv`` of the project."""
    raw_dir = project_dir / "raw"
    reads, vcfs = [], []
    for folder in folders:
        f_dir = raw_dir / folder
        f_dir.mkdir(parents=True, exist_ok=True)
        for suffix in (".R1.fastq.gz", ".R2.fastq.gz", ".vcf.gz", ".vcf.gz.tbi"):
            (f_dir / f"{folder}{suffix}").touch()
        reads.append(
            f"{folder}\t{folder}\traw/{folder}/{folder}.R1.fastq.gz\traw/{folder}/{folder}.R2.fastq.gz"
        )
        vcfs.append(f"{folder}\traw/{folder}/{folder}.vcf.gz\traw/{folder}/{folder}.vcf.gz.tbi")
    mode = "a" if append else "w"
    with (project_dir / "reads.tsv").open(mode) as f:
        f.write("" if append else "library\treadgroup\tleft\tright\n")
        f.write("".join(line + "\n" for line in reads))
    with (project_dir / "vcfs.tsv").open(mode) as f:
        f.write("" if append else "library\tvcf\tvcf_tbi\n")
        f.write("".join(line + "\n" for line in vcfs))


def _upstream_roots(upstream: str, tasks: dict[str, Any], project_dir: Path) -> list[str]:
    """Return the path prefixes of the files that a task reads from ``upstream``."""
    if tasks[upstream]["step"] != "external_data":
        return [f"tasks/{upstream}/"]
    config = tasks[upstream]["config"]
    if config.get("library_files"):
        return [str(project_dir / "raw")]
    paths = []
    for value in config.get("files", {}).values():
        paths += value if isinstance(value, list) else [value]
    for kit in config.get("kits", {}).values():
        paths += kit.values()
    return [os.path.normpath(os.path.join(project_dir, path)) for path in paths]


#: Required ``depends_on`` fields that a task reads where the DAG of ``snappy run --task`` does not
#: show it, by task and field, with the reason
READ_OUTSIDE_THE_DAG = {
    ("ngs_mapping_bwa_mem2", "reference"): "load_project() checks that the index is of this genome",
    ("ngs_mapping_bwa_mem3", "reference"): "load_project() checks that the index is of this genome",
    ("sv_calling_targeted_gcnv", "alignments"): "read after the checkpoint of the gCNV shards",
    ("sv_calling_wgs_gcnv", "alignments"): "read after the checkpoint of the gCNV shards",
}


def _unread_dependencies(task_name: str, project_dir: Path, jobs: list[dict]) -> list[str]:
    """Return the required or marked (``ReadBy``, ``RequiredBy``) ``depends_on`` fields of
    ``task_name`` that no job of the task reads; unmarked optional fields are skipped."""
    config = yaml.safe_load((project_dir / "config.yaml").read_text(encoding="utf-8"))
    tasks = {task["name"]: task for task in config["tasks"]}
    own = [job for job in jobs if job["rule"].startswith(f"{task_name}__")]
    inputs = [path for job in own for path in job["input"]]
    model = WORKFLOW_REGISTRY[tasks[task_name]["step"]].config_model_class
    if "depends_on" not in model.model_fields:
        return []
    fields = model.model_fields["depends_on"].annotation.model_fields
    unread = []
    for field, value in tasks[task_name]["config"].get("depends_on", {}).items():
        checked = fields[field].is_required() or reader_markers(fields[field])
        if not checked or (task_name, field) in READ_OUTSIDE_THE_DAG:
            continue
        for upstream in value if isinstance(value, list) else [value]:
            roots = _upstream_roots(upstream, tasks, project_dir)
            if not any(path.startswith(root) for path in inputs for root in roots):
                unread.append(f"{field}: {upstream}")
    return unread


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
    jobs = _dump(tmp_path, task_name)
    unread = [] if task_name in EXPECTED_EMPTY else _unread_dependencies(task_name, tmp_path, jobs)
    assert not unread, f"{task_name} sets depends_on fields that none of its jobs reads: {unread}"
    targets = next(job["input"] for job in jobs if job["rule"] == "all")
    if task_name in EXPECTED_EMPTY:
        assert not targets, f"{task_name} now targets files; remove it from EXPECTED_EMPTY"
    else:
        assert targets, f"{task_name} targets no files, so its snapshot would check nothing"
    actual = json.dumps(jobs, indent=1, sort_keys=True) + "\n"
    for path, placeholder in ((str(tmp_path), "<project>"), (str(root), "<repo>")):
        actual = actual.replace(path, placeholder)
    # The code hash of a rule's wrappers changes with every edit of the tools they import
    actual = re.sub(r'"code": "[0-9a-f]{16}"', '"code": "<hash>"', actual)
    _check_snapshot(task_name, actual)


@pytest.mark.integration
def test_frozen_run_uses_existing_upstream_outputs(
    generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    """``snappy run --task X --frozen`` loads only X's rules; missing upstream files fail fast."""
    _write_closure_project("variant_annotation_vep", generated_task_config, tmp_path)
    with pytest.raises(AssertionError, match="MissingInputException") as failure:
        _dump(tmp_path, "variant_annotation_vep", frozen=True)
    assert "tasks/variant_calling_gatk4_hc_gvcf/work/" in str(failure.value)


@pytest.mark.integration
def test_job_group_groups_the_jobs_of_its_tasks(
    generated_task_config: dict[str, Any], tmp_path: Path
) -> None:
    """The non-local jobs of the tasks with ``job_group`` join that group; local rules stay local,
    and the jobs of other tasks have no group."""
    _write_closure_project("variant_filtration_bcftools", generated_task_config, tmp_path)
    grouped = {"variant_filtration_bcftools"}

    def join_group(config):
        tasks = {task["name"]: task for task in config["tasks"]}
        grouped.add(tasks["variant_filtration_bcftools"]["config"]["depends_on"]["variants"])
        for name in grouped:
            tasks[name]["config"]["job_group"] = "tail"

    names = [task["name"] for task in _edit_config(tmp_path, join_group)["tasks"]]
    jobs = _dump(tmp_path, "variant_filtration_bcftools")

    def task_of(job):
        return next((name for name in names if job["rule"].startswith(f"{name}_")), None)

    in_group = [job for job in jobs if task_of(job) in grouped and not job.get("local")]
    assert in_group and all(job.get("group") == "tail" for job in in_group)
    assert not [job for job in jobs if job.get("local") and "group" in job]
    assert not [job for job in jobs if task_of(job) not in grouped and "group" in job]


def _job_with_reads_from(
    task_name: str,
    rule_suffix: str,
    source_task: dict[str, Any],
    generated_task_config: dict[str, Any],
    tmp_path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return the first job of the rule ending in ``rule_suffix`` and all jobs when the task
    ``task_name`` reads from ``source_task``."""
    _write_closure_project(task_name, generated_task_config, tmp_path)

    def read_from_source(config):
        config["tasks"] = [t for t in config["tasks"] if t["name"] != source_task["name"]]
        config["tasks"].insert(0, source_task)
        consumer = next(t for t in config["tasks"] if t["name"] == task_name)
        consumer["config"]["depends_on"]["reads"] = source_task["name"]

    _edit_config(tmp_path, read_from_source)
    jobs = _dump(tmp_path, task_name)
    return next(j for j in jobs if j["rule"].endswith(rule_suffix)), jobs


#: An adapter trimming task that reads the FASTQs of the reads task
TRIMMING = {
    "step": "adapter_trimming",
    "name": "trimming",
    "config": {"depends_on": {"reads": "reads"}, "tool": "fastp", "fastp": {}},
}


@pytest.mark.integration
def test_mapping_reads_trimmed_fastqs(generated_task_config: dict[str, Any], tmp_path: Path):
    job, jobs = _job_with_reads_from(
        "ngs_mapping_bwa", "_bwa_run", TRIMMING, generated_task_config, tmp_path
    )
    lib = job["wildcards"]["library_name"]

    # Only the directory is a job output, so the FASTQs are params, not inputs.
    assert f"tasks/trimming/work/{lib}/out" in job["input"]
    assert not [path for path in job["input"] if path.endswith(".fastq.gz")]
    assert job["params"]["args"]["input"]["reads_left"] == [
        f"tasks/trimming/work/{lib}/out/{lib}.R1.fastq.gz"
    ]
    assert any(j["rule"].startswith("trimming__fastp") for j in jobs)


@pytest.mark.integration
def test_mapping_bam_is_temp_when_the_task_deletes_it(
    generated_task_config: dict[str, Any], tmp_path: Path
):
    """With ``delete_intermediates: [bam]`` the BAM is temp, has no .md5 and is no result file."""
    _write_closure_project("ngs_mapping_bwa", generated_task_config, tmp_path)

    def delete_bam(config):
        mapping = next(t for t in config["tasks"] if t["name"] == "ngs_mapping_bwa")
        mapping["config"]["delete_intermediates"] = ["intermediates", "bam"]

    _edit_config(tmp_path, delete_bam)
    jobs = _dump(tmp_path, "ngs_mapping_bwa")

    job = next(j for j in jobs if j["rule"].endswith("_bwa_run"))
    lib = job["wildcards"]["library_name"]
    bam = f"tasks/ngs_mapping_bwa/work/{lib}/out/{lib}.bam"
    assert {bam, bam + ".bai"} <= set(job["temp"])
    assert not [path for path in job["output"] if path.endswith((".bam.md5", ".bai.md5"))]
    results = next(j for j in jobs if j["rule"] == "ngs_mapping_bwa__results")
    assert not [path for path in results["input"] if path.endswith((".bam", ".bai"))]
    assert any(path.endswith(".bamstats.txt") for path in results["input"])


@pytest.mark.integration
@pytest.mark.parametrize("tool", ["arriba", "star_fusion"])
def test_fusion_calling_reads_trimmed_fastqs(
    tool: str, generated_task_config: dict[str, Any], tmp_path: Path
):
    task_name = f"somatic_gene_fusion_calling_{tool}"
    job, _ = _job_with_reads_from(
        task_name, f"_{tool}_run", TRIMMING, generated_task_config, tmp_path
    )
    lib = job["wildcards"]["library_name"]

    # The trimmed FASTQs are no job outputs, so they are params; the input is their directory.
    assert f"tasks/trimming/work/{lib}/out" in job["input"]
    assert not [path for path in job["input"] if path.endswith(".fastq.gz")]
    assert job["params"]["args"]["input"]["reads_left"] == [
        f"tasks/trimming/work/{lib}/out/{lib}.R1.fastq.gz"
    ]


@pytest.mark.integration
def test_mapping_reads_from_an_external_data_task(
    generated_task_config: dict[str, Any], tmp_path: Path
):
    external_reads = {
        "step": "external_data",
        "name": "external_reads",
        "config": {
            "produces": {"type": "reads"},
            "library_files": [str(tmp_path / "reads.tsv")],
        },
    }
    job, _ = _job_with_reads_from(
        "ngs_mapping_bwa", "_bwa_run", external_reads, generated_task_config, tmp_path
    )
    lib = job["wildcards"]["library_name"]

    assert "input" not in job["params"]["args"]
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
def test_somatic_trio_calls_the_relapse_against_the_primary_tumor(
    generated_task_config: dict[str, Any], tmp_path: Path
):
    # A sample sheet column (timepoint) selects the tumors; a relationship pairs the relapse with
    # the primary tumor instead of the normal
    _write_closure_project("variant_calling_mutect2", generated_task_config, tmp_path)
    _write_files_tables(tmp_path, ("case001subregion-T2-DNA1-WES1",), append=True)
    primary_tumor = "role == 'tumor' and extraction_type == 'dna' and timepoint == 'primary'"
    tasks = {
        "calling_primary": {"library_selection": "timepoint == 'primary'"},
        "calling_relapse": {
            "library_selection": "timepoint == 'relapse'",
            "relationships": {"matched_normal_lib": {"via": "donor_name", "target": primary_tumor}},
        },
    }

    def add_tasks(config):
        config["sample_sheets"] = [str(_fixture_dir() / "samplesheet_somatic_trio.tsv")]
        calling = next(t for t in config["tasks"] if t["name"] == "variant_calling_mutect2")
        for name, settings in tasks.items():
            config["tasks"].append(
                {**calling, "name": name, "config": {**calling["config"], **settings}}
            )

    _edit_config(tmp_path, add_tasks)
    bams = {}
    for name in tasks:
        for job in _dump(tmp_path, name):
            if job["rule"] == f"{name}__mutect2_run":
                libraries = [
                    Path(p).name.removesuffix(".bam") for p in job["input"] if p.endswith(".bam")
                ]
                bams.setdefault(name, set()).add(tuple(libraries))

    # The tumor first, then the library it is compared with
    assert bams == {
        "calling_primary": {("case001subregion-T1-DNA1-WES1", "case001subregion-N1-DNA1-WES1")},
        "calling_relapse": {("case001subregion-T2-DNA1-WES1", "case001subregion-T1-DNA1-WES1")},
    }
