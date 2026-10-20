# -*- coding: utf-8 -*-
"""Tests for loading a project and selecting targets in ``snappy_pipeline.orchestration``."""

import logging
from typing import Annotated

import pydantic
import pytest

from snappy_pipeline import orchestration
from snappy_pipeline.orchestration import Project, load_project, select_target_tasks
from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflow_model import ConfigModel
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import (
    TOOL,
    DataSignature,
    DataType,
    ExpectedReference,
)

WORK_DIR = "/projects/p1"


#: The reference genome task that every test project has
GENOME = (
    "external_data",
    "genome",
    {
        "produces": {"type": "raw", "tags": ["reference", "dna"]},
        "files": {"fasta": "/refs/genome.fa"},
    },
)


GENES = ("external_data", "genes", {"produces": {"type": "raw", "tags": ["features"]}})
GENES[2]["files"] = {"gtf": "/refs/genes.gtf"}
DBSNP = ("external_data", "dbsnp", {"produces": {"type": "variants", "tags": ["dbsnp"]}})
DBSNP[2]["files"] = {"vcf": "/refs/dbsnp.vcf.gz"}


#: The BWA index that ``_mapping`` tasks read
BWA_INDEX = ("external_data", "bwa_index", {"produces": {"type": "index", "tags": ["bwa", "dna"]}})
BWA_INDEX[2]["files"] = {"index": "/refs/genome"}

#: The VEP cache that ``_annotation`` tasks read
VEP_CACHE = ("external_data", "vep_cache", {"produces": {"type": "database", "tags": ["vep"]}})
VEP_CACHE[2]["files"] = {"cache": "/refs/vep"}

#: The enrichment kits that ``_tmb`` tasks read
KITS = (
    "external_data",
    "kits",
    {
        "produces": {"type": "regions", "tags": ["kits"]},
        "kits": {"MedExome": {"targets": "/refs/t.bed", "baits": "/refs/b.bed"}},
    },
)


def _config(*tasks):
    return {
        "tasks": [
            {"step": step, "name": name, "config": config}
            for step, name, config in (GENOME, BWA_INDEX, *tasks, VEP_CACHE)
        ],
        "data_sets": {},
    }


def _mapping(name="mapping", depends_on=None, **config):
    depends_on = {"reads": "data_sets", "reference": "genome", "index": "bwa_index"} | (
        depends_on or {}
    )
    config = {"tool": "bwa", "bwa": {}, **config}
    return ("ngs_mapping", name, {"depends_on": depends_on, **config})


def _calling(name="calling", mapping="mapping", tool="mutect2", **depends_on):
    depends_on = {"alignments": mapping, "reference": "genome", **depends_on}
    config = {"depends_on": depends_on, "tool": tool}
    config[tool] = {"contamination": {}} if tool == "mutect2" else {}
    return ("variant_calling", name, config)


def _annotation(name="annotation", variants="calling"):
    return (
        "variant_annotation",
        name,
        {
            "depends_on": {"variants": variants, "reference": "genome", "database": "vep_cache"},
            "tool": "vep",
            "vep": {},
        },
    )


def _filtration(name="filtration", variants="annotation", **depends_on):
    config = {"depends_on": {"variants": variants, **depends_on}, "tool": "bcftools"}
    config["bcftools"] = {"include": "QUAL > 10"}
    return ("variant_filtration", name, config)


def _tmb(name="tmb", variants="filtration"):
    return ("tumor_mutational_burden", name, {"depends_on": {"variants": variants, "kits": "kits"}})


# load_project ------------------------------------------------------------------------------------


def test_load_project_resolves_dependencies_and_orders_tasks():
    # Listed out of order on purpose.
    project = load_project(_config(_filtration(), _annotation(), _calling(), _mapping()), WORK_DIR)

    assert [task.name for task in project.tasks] == [
        "genome",
        "bwa_index",
        "mapping",
        "calling",
        "vep_cache",
        "annotation",
        "filtration",
    ]
    assert project.dependencies["filtration"] == {"variants": "annotation"}
    assert project.dependencies["mapping"] == {"reference": "genome", "index": "bwa_index"}
    assert project.task_configs["calling"].tool == "mutect2"
    assert project.project_dir == WORK_DIR
    assert project.config_paths == (f"{WORK_DIR}/config.yaml",)


def test_load_project_rejects_unknown_step():
    with pytest.raises(ValueError, match="unknown step 'variant_caling'"):
        load_project(_config(("variant_caling", "calling", {})), WORK_DIR)


def test_load_project_rejects_duplicate_task_names():
    with pytest.raises(ValueError, match="duplicates: mapping"):
        load_project(_config(_mapping(), _mapping()), WORK_DIR)


def test_load_project_rejects_dependency_on_missing_task():
    with pytest.raises(
        ValueError, match="depends_on.alignments is 'mapping2', which is not a task"
    ):
        load_project(_config(_mapping(), _calling(mapping="mapping2")), WORK_DIR)


def test_load_project_rejects_dependency_cycles():
    tasks = (_annotation(variants="filtration"), _filtration(variants="annotation"))
    with pytest.raises(ValueError, match="Dependency cycle between tasks"):
        load_project(_config(*tasks), WORK_DIR)


def test_load_project_uses_defaults_for_keys_without_value(caplog):
    with caplog.at_level(logging.INFO, logger="snappy_pipeline.orchestration"):
        mapping = _mapping(bwa={"mask_duplicates": None})
        project = load_project(_config(mapping), WORK_DIR)

    assert project.task_configs["mapping"].bwa.mask_duplicates is True
    assert "tasks[2].config.bwa.mask_duplicates has no value" in caplog.text


def test_load_project_validates_step_configs():
    with pytest.raises(pydantic.ValidationError, match="contamination"):
        load_project(
            _config(("variant_calling", "calling", {"tool": "mutect2", "mutect2": {}})), WORK_DIR
        )


def test_load_project_passes_variant_tags_through_annotation_and_filtration():
    tasks = (KITS, _mapping(), _calling(), _annotation(), _filtration(), _tmb())
    project = load_project(_config(*tasks), WORK_DIR)

    somatic = frozenset({"somatic", "snv", "indel"})
    assert project.signatures["mapping"] == (
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
    )
    assert project.signatures["filtration"] == (
        DataSignature(DataType.VARIANTS, somatic | {"annotated", "filtered"}),
    )


def test_load_project_rejects_germline_variants_for_tmb():
    calling = _calling(tool="gatk4_hc_gvcf", dbsnp="dbsnp")
    tasks = (KITS, DBSNP, _mapping(), calling, _annotation(), _filtration(), _tmb())
    with pytest.raises(
        ValueError,
        match=r"Task 'tmb': depends_on.variants requires variants \[somatic\], "
        r"but task 'filtration' produces variants \[annotated, filtered, germline, indel, snv\]",
    ):
        load_project(_config(*tasks), WORK_DIR)


STAR_INDEX = (
    "external_data",
    "star_index",
    {"produces": {"type": "index", "tags": ["star", "rna"]}},
)
STAR_INDEX[2]["files"] = {"index": "/refs/star"}
EXONS = ("external_data", "exons", {"produces": {"type": "regions", "tags": ["exons"]}})
EXONS[2]["files"] = {"bed": "/refs/exons.bed"}


def test_load_project_rejects_rna_alignments_for_variant_calling():
    depends_on = {"features": "genes", "index": "star_index", "exons": "exons"}
    star = _mapping("star", depends_on=depends_on, tool="star", star={})
    with pytest.raises(
        ValueError,
        match=r"Task 'calling': depends_on.alignments requires alignments \[dna\], "
        r"but task 'star' produces alignments \[rna\]",
    ):
        load_project(_config(GENES, STAR_INDEX, EXONS, star, _calling(mapping="star")), WORK_DIR)


def test_load_project_rejects_dna_alignments_for_expression_quantification():
    strandedness = {"tool": "strandedness", "strandedness": {}}
    expression = ("gene_expression_quantification", "expression", strandedness)
    expression[2]["depends_on"] = {"alignments": "mapping", "exons": "exons"}
    with pytest.raises(
        ValueError,
        match=r"Task 'expression': depends_on.alignments requires alignments \[rna\], "
        r"but task 'mapping' produces alignments \[dna\]",
    ):
        load_project(_config(EXONS, _mapping(), expression), WORK_DIR)


def test_load_project_requires_depends_on_keys():
    # No default task name such as "ngs_mapping" fills in a missing key.
    calling = ("variant_calling", "calling", {"tool": "mutect2", "mutect2": {"contamination": {}}})
    with pytest.raises(pydantic.ValidationError, match="depends_on\n  Field required"):
        load_project(_config(_mapping(), calling), WORK_DIR)
    mapping = _mapping()
    del mapping[2]["depends_on"]["reads"]
    with pytest.raises(pydantic.ValidationError, match="depends_on.reads\n  Field required"):
        load_project(_config(mapping), WORK_DIR)


def test_load_project_requires_tool():
    calling = ("variant_calling", "calling", {"depends_on": {"alignments": "mapping"}})
    with pytest.raises(pydantic.ValidationError, match="tool\n  Field required"):
        load_project(_config(_mapping(), calling), WORK_DIR)


def test_load_project_accepts_reads_from_data_sets():
    project = load_project(_config(_mapping(depends_on={"reads": "data_sets"})), WORK_DIR)
    assert project.dependencies["mapping"] == {"reference": "genome", "index": "bwa_index"}


def test_load_project_reserves_data_sets_for_reads():
    with pytest.raises(
        ValueError, match="depends_on.alignments is 'data_sets', which is not a task"
    ):
        load_project(_config(_mapping(), _calling(mapping="data_sets")), WORK_DIR)
    with pytest.raises(ValueError, match="Task name 'data_sets' is reserved"):
        load_project(_config(_mapping(name="data_sets")), WORK_DIR)


def _external(name, data_type, tags, **config):
    return ("external_data", name, {"produces": {"type": data_type, "tags": tags}, **config})


BAM_PATTERN = {"bam": r".+\.bam", "bai": r".+\.bam\.bai"}
#: A reads pattern with only a left file
LEFT = {"left": r"(?P<readgroup>.+)_1\.fq\.gz"}


def test_external_alignments_feed_variant_calling():
    bams = _external(
        "bams", "alignments", ["dna"], search_paths=["/data"], search_patterns=[BAM_PATTERN]
    )
    project = load_project(_config(bams, _calling(mapping="bams")), WORK_DIR)
    assert project.signatures["bams"] == (DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),)


def test_external_data_is_checked_against_requirements_and_schemas():
    vcfs = _external("vcfs", "variants", ["germline"], files={"vcf": "/data/x.vcf.gz"})
    with pytest.raises(ValueError, match=r"depends_on.variants requires variants \[somatic\]"):
        load_project(_config(KITS, vcfs, _tmb(variants="vcfs")), WORK_DIR)

    bams = _external(
        "bams", "alignments", ["dna"], search_paths=["/data"], search_patterns=[{"bam": r".+\.bam"}]
    )
    with pytest.raises(
        ValueError, match="depends_on.alignments needs bai, which task 'bams' does not provide"
    ):
        load_project(_config(bams, _calling(mapping="bams")), WORK_DIR)


def test_regions_filter_reads_any_bed_file_of_regions():
    vcfs = _external(
        "vcfs", "variants", ["germline"], files={"vcf": "/x.vcf.gz", "vcf_tbi": "/x.tbi"}
    )
    config = {
        "depends_on": {"variants": "vcfs"},
        "tool": "regions",
        "regions": {"action": "exclude"},
    }
    filtration = ("variant_filtration", "filtration", config)
    with pytest.raises(pydantic.ValidationError, match="depends_on.regions is required"):
        load_project(_config(vcfs, EXONS, KITS, filtration), WORK_DIR)

    config["depends_on"]["regions"] = "exons"
    load_project(_config(vcfs, EXONS, KITS, filtration), WORK_DIR)

    config["depends_on"]["regions"] = "genome"
    with pytest.raises(ValueError, match=r"depends_on.regions requires regions, but task 'genome'"):
        load_project(_config(vcfs, EXONS, KITS, filtration), WORK_DIR)

    config["depends_on"]["regions"] = "kits"
    with pytest.raises(
        ValueError, match="depends_on.regions needs bed, which task 'kits' does not"
    ):
        load_project(_config(vcfs, EXONS, KITS, filtration), WORK_DIR)


@pytest.mark.parametrize(
    "config, error",
    [
        ({}, "Set one of files"),
        (
            {"files": {"x": "/a"}, "search_paths": ["/b"], "search_patterns": [LEFT]},
            "Set one of files",
        ),
        ({"files": {"x": "/a"}, "kits": {"k": {"x": "/a"}}}, "Set one of files"),
        ({"search_paths": ["/b"]}, "search_paths need search_patterns"),
        ({"search_paths": ["/b"], "search_patterns": [{"left": r".+\.fastq\.gz"}]}, "readgroup"),
        ({"kits": {"Kit V1": {"targets": "/a"}}}, "Kit names may only contain"),
        ({"kits": {"a": {"targets": "/a"}, "b": {"baits": "/b"}}}, "same keys"),
    ],
)
def test_external_data_config_errors(config, error):
    with pytest.raises(pydantic.ValidationError, match=error):
        load_project(_config(_external("reads", "raw", [], **config)), WORK_DIR)


# create_task_instances ---------------------------------------------------------------------------


def test_create_task_instances_builds_each_task_once(monkeypatch):
    project = load_project(_config(_calling(), _mapping()), WORK_DIR)
    built = []

    class Stub:
        def __init__(self, workflow, project, task_name):
            built.append(task_name)

    for step in ("ngs_mapping", "variant_calling"):
        monkeypatch.setitem(WORKFLOW_REGISTRY, step, Stub)
    instances = orchestration.create_task_instances(object(), project)

    assert built == ["mapping", "calling"]
    assert orchestration.task_instance("calling") is instances["calling"]


# select_target_tasks -----------------------------------------------------------------------------


def _project(**dependencies):
    """Return a project with the given task -> {field: upstream} dependencies."""
    tasks = [{"step": "dummy", "name": name, "config": {}} for name in dependencies]
    model = ConfigModel(tasks=tasks, data_sets={})
    return Project(
        config={},
        model=model,
        project_dir=WORK_DIR,
        config_paths=(),
        tasks=tuple(model.tasks),
        task_configs={},
        dependencies=dependencies,
        signatures={},
    )


#: mapping -> calling -> annotation -> filtration, plus a QC task on the mapping.
CHAIN = _project(
    mapping={},
    calling={"alignments": "mapping"},
    annotation={"variants": "calling"},
    filtration={"variants": "annotation", "alignments": "mapping"},
    qc={"alignments": "mapping"},
)


def test_default_targets_every_task_in_config_order():
    assert select_target_tasks(CHAIN) == ["mapping", "calling", "annotation", "filtration", "qc"]


def test_single_task():
    assert select_target_tasks(CHAIN, target_task="calling") == ["calling"]


def test_unknown_single_task_raises():
    with pytest.raises(ValueError, match="Unknown task 'calls'"):
        select_target_tasks(CHAIN, target_task="calls")


def test_load_project_rejects_upstream_without_the_output_keys_a_consumer_needs():
    cnv = ("somatic_wgs_cnv_calling", "cnv", {"tool": "cnvetti", "cnvetti": {}})
    cnv[2]["depends_on"] = {"alignments": "mapping", "reference": "genome"}
    check = ("somatic_cnv_checking", "cnv_check", {})
    check[2]["depends_on"] = {"alignments": "mapping", "copy_number": "cnv", "reference": "genome"}
    with pytest.raises(
        ValueError,
        match=r"Task 'cnv_check': depends_on.copy_number needs dnacopy_seg, which task 'cnv' "
        r"does not provide; it provides bcf, bcf_csi",
    ):
        load_project(_config(_mapping(), cnv, check), WORK_DIR)


def test_external_data_files_may_name_lists():
    dbs = _external("dbs", "database", ["mehari"], files={"transcript_db": ["/a.bin", "/b.bin"]})
    project = load_project(_config(dbs), WORK_DIR)
    config = project.task_configs["dbs"]
    paths = WORKFLOW_REGISTRY["external_data"].get_output_paths(config)
    assert paths == {"transcript_db": ["/a.bin", "/b.bin"]}


def test_requirement_with_tool_placeholder_names_the_task_tool():
    # ngs_mapping requires the index of its tool (mbcs: of its mapping_tool)
    bwa_mem2 = _mapping(tool="bwa_mem2", bwa_mem2={})
    with pytest.raises(
        ValueError,
        match=r"depends_on.index requires index \[bwa_mem2, dna\|rna\], but task 'bwa_index'",
    ):
        load_project(_config(bwa_mem2), WORK_DIR)

    mbcs = _mapping(tool="mbcs", mbcs={"mapping_tool": "bwa", "recalibrate": False})
    load_project(_config(mbcs), WORK_DIR)


def test_tool_placeholder_needs_a_tool():
    required = DataSignature(DataType.DATABASE, frozenset({TOOL}))
    assert required.for_tool("vep") == DataSignature(DataType.DATABASE, frozenset({"vep"}))
    with pytest.raises(ValueError, match="names the tool, but the task has none"):
        required.for_tool(None)


@pytest.mark.parametrize(
    "fasta, sequence_dict",
    [
        ("/refs/genome.fa", "/refs/genome.dict"),
        ("/refs/genome.fasta", "/refs/genome.dict"),
        ("/refs/genome.fna.gz", "/refs/genome.dict"),
        ("/refs/hs37d5.fa.gz", "/refs/hs37d5.dict"),
    ],
)
def test_reference_sidecars_follow_the_samtools_and_gatk_names(fasta, sequence_dict):
    reference = ExpectedReference(fasta=fasta)
    assert reference.fai == fasta + ".fai"
    assert reference.sequence_dict == sequence_dict


def test_config_paths_are_resolved_against_the_project_directory_with_pathvars():
    relative = _external("relative", "regions", [], files={"bed": "beds/a.bed"})
    named = _external("named", "regions", [], files={"bed": "<refs>/b.bed"})
    default = _external("default", "regions", [], files={"bed": "<resources>/c.bed"})
    config = _config(relative, named, default) | {"pathvars": {"refs": "/data/refs"}}
    data_set = {"file": "sheets/sheet.tsv", "search_paths": ["raw"], "search_patterns": [LEFT]}
    config["data_sets"] = {"first": data_set}
    project = load_project(config, WORK_DIR)

    assert project.task_configs["relative"].files == {"bed": f"{WORK_DIR}/beds/a.bed"}
    assert project.task_configs["named"].files == {"bed": "/data/refs/b.bed"}
    # Snakemake's default pathvars apply, too
    assert project.task_configs["default"].files == {"bed": f"{WORK_DIR}/resources/c.bed"}
    data_set = project.model.data_sets["first"]
    assert data_set.file == f"{WORK_DIR}/sheets/sheet.tsv"
    assert data_set.search_paths == [f"{WORK_DIR}/raw"]
    assert project.file_paths["relative"] == (f"{WORK_DIR}/beds/a.bed",)


def test_undefined_pathvar_fails_when_the_project_loads():
    named = _external("named", "regions", [], files={"bed": "<refs>/b.bed"})
    with pytest.raises(pydantic.ValidationError, match="Undefined pathvar 'refs'"):
        load_project(_config(named), WORK_DIR)


def test_missing_files_lists_paths_and_genome_sidecars(tmp_path):
    (tmp_path / "genome.fa").touch()
    (tmp_path / "genome.fa.fai").touch()
    (tmp_path / "index.amb").touch()
    genome = (
        "external_data",
        "genome",
        {"produces": {"type": "raw", "tags": ["reference", "dna"]}},
    )
    genome[2]["files"] = {"fasta": "genome.fa"}
    bwa_index = (
        "external_data",
        "bwa_index",
        {"produces": {"type": "index", "tags": ["bwa", "dna"]}},
    )
    bwa_index[2]["files"] = {"index": "index"}  # a prefix: index.amb, ...
    tasks = (genome, bwa_index, _mapping())
    config = {"tasks": [{"step": s, "name": n, "config": c} for s, n, c in tasks], "data_sets": {}}
    project = load_project(config, str(tmp_path))

    assert orchestration.missing_files(project) == {"genome": [str(tmp_path / "genome.dict")]}


class _HlaListDependsOn(SnappyModel):
    hla_types: Annotated[list[str], DataSignature(DataType.TABULAR, frozenset({"hla"}))]


class _HlaListConfig(SnappyStepModel):
    depends_on: _HlaListDependsOn


class _HlaListStep(BaseStep):
    """A step whose ``depends_on.hla_types`` takes several tasks."""

    name = "hla_list"
    produces = [DataSignature(DataType.TABULAR, frozenset({"neoepitope"}))]
    config_model_class = _HlaListConfig


def _hla(name, tags=("hla",)):
    produces = {"type": "tabular", "tags": list(tags)}
    return ("external_data", name, {"produces": produces, "files": {"txt": f"/refs/{name}.txt"}})


@pytest.fixture
def hla_list_step(monkeypatch):
    monkeypatch.setitem(WORKFLOW_REGISTRY, "hla_list", _HlaListStep)


def test_a_list_field_depends_on_each_of_its_tasks(hla_list_step):
    tasks = (
        ("hla_list", "consumer", {"depends_on": {"hla_types": ["hla_rna", "hla_dna"]}}),
        _hla("hla_rna"),
        _hla("hla_dna"),
    )
    project = load_project(_config(*tasks), WORK_DIR)
    assert project.dependencies["consumer"] == {"hla_types": ("hla_rna", "hla_dna")}
    order = [task.name for task in project.tasks]
    assert order.index("consumer") > max(order.index("hla_rna"), order.index("hla_dna"))


@pytest.mark.parametrize(
    "hla_types, error",
    [
        (["hla_dna", "hla_x"], "depends_on.hla_types contains 'hla_x', which is not a task"),
        (["hla_dna", "hla_dna"], "depends_on.hla_types names a task twice"),
        (["hla_dna", "other"], "requires tabular \\[hla\\], but task 'other' produces"),
    ],
)
def test_each_task_of_a_list_field_is_checked(hla_list_step, hla_types, error):
    tasks = (
        ("hla_list", "consumer", {"depends_on": {"hla_types": hla_types}}),
        _hla("hla_dna"),
        _hla("other", tags=("other",)),
    )
    with pytest.raises(ValueError, match=error):
        load_project(_config(*tasks), WORK_DIR)
