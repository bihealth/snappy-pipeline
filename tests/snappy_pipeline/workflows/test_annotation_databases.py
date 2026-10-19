# -*- coding: utf-8 -*-
"""Tests for the ``depends_on.database`` of variant annotation and VarFish export."""

import pydantic
import pytest

from snappy_pipeline.workflows.variant_annotation import MehariStepPart, VepStepPart
from snappy_pipeline.workflows.variant_annotation.model import VariantAnnotation
from snappy_pipeline.workflows.varfish_export.model import VarfishExport

from ..test_orchestration import (
    WORK_DIR,
    _annotation,
    _calling,
    _config,
    _external,
    _mapping,
    load_project,
)
from .test_contracts import _step


def _annotation_inputs(part_cls, project):
    step = _step(project, "annotation")
    wildcards = type("Wildcards", (), {"library_name": "L1"})
    return part_cls(step)._get_input_files_run(wildcards)


def _mehari_annotation():
    annotation = _annotation()
    annotation[2].update(tool="mehari", mehari={})
    annotation[2]["depends_on"]["database"] = "mehari_db"
    return annotation


def test_vep_needs_a_database():
    depends_on = {"variants": "calling", "reference": "genome"}
    with pytest.raises(pydantic.ValidationError, match="database"):
        VariantAnnotation(tool="vep", vep={}, depends_on=depends_on)
    VariantAnnotation(tool="vep", vep={}, depends_on=depends_on | {"database": "vep_cache"})


def test_vep_cache_is_an_input_of_the_job():
    project = load_project(_config(_mapping(), _calling(), _annotation()), WORK_DIR)
    assert _annotation_inputs(VepStepPart, project)["cache"] == "/refs/vep"


def test_mehari_database_may_list_transcript_databases():
    databases = _external(
        "mehari_db",
        "database",
        ["mehari"],
        files={
            "transcript_db": ["/db/a.bin.zst", "/db/b.bin.zst"],
            "frequency_db": "/db/frequencies",
        },
    )
    project = load_project(
        _config(_mapping(), _calling(), databases, _mehari_annotation()), WORK_DIR
    )
    inputs = _annotation_inputs(MehariStepPart, project)
    assert inputs["transcripts"] == ["/db/a.bin.zst", "/db/b.bin.zst"]
    assert inputs["frequencies"] == "/db/frequencies"
    assert "clinvar" not in inputs


def test_mehari_database_needs_transcripts():
    databases = _external("mehari_db", "database", ["mehari"], files={"clinvar_db": "/db/clinvar"})
    with pytest.raises(ValueError, match="transcript_db"):
        load_project(_config(_mapping(), _calling(), databases, _mehari_annotation()), WORK_DIR)


def test_varfish_export_needs_all_mehari_keys():
    depends_on = {
        "alignments": "mapping",
        "variants": "calling",
        "reference": "genome",
        "database": "mehari_db",
    }
    VarfishExport(depends_on=depends_on)
    with pytest.raises(pydantic.ValidationError, match="database"):
        VarfishExport(depends_on={k: v for k, v in depends_on.items() if k != "database"})

    keys = {
        "hgnc_tsv": "/db/hgnc.tsv",
        "transcript_db": "/db/txs.bin.zst",
        "clinvar_db": "/db/clinvar",
        "frequency_db": "/db/frequencies",
    }
    incomplete = _external("mehari_db", "database", ["mehari"], files={"hgnc_tsv": "/db/hgnc.tsv"})
    export = (
        "varfish_export",
        "export",
        {"depends_on": depends_on},
    )
    tasks = (_mapping(), _calling(tool="bcftools_call"), export)
    with pytest.raises(ValueError, match="clinvar_db"):
        load_project(_config(*tasks, incomplete), WORK_DIR)
    complete = _external("mehari_db", "database", ["mehari"], files=keys)
    project = load_project(_config(*tasks, complete), WORK_DIR)
    assert _step(project, "export").get_upstream_paths("database").model_dump() == keys
