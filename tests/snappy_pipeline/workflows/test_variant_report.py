# -*- coding: utf-8 -*-
"""Tests of the ``variant_report`` step: the calls of the VCFs of an upstream task as tables."""

from types import SimpleNamespace

import pytest
from pydantic import ValidationError
from snakemake.io import is_flagged

from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_report import DatavzrdStepPart, VariantReportWorkflow
from snappy_pipeline.workflows.variant_report.model import VariantReport


def _upstream_paths(field, library_name):
    if field == "alignments":
        return SimpleNamespace(
            bam=f"tasks/mapping/work/{library_name}/out/{library_name}.bam",
            bai=f"tasks/mapping/work/{library_name}/out/{library_name}.bam.bai",
        )
    assert field == "variants"
    return SimpleNamespace(
        vcf=f"tasks/fdr/work/{library_name}/out/{library_name}.vcf.gz",
        vcf_tbi=f"tasks/fdr/work/{library_name}/out/{library_name}.vcf.gz.tbi",
    )


def _part(tags=frozenset({"germline", "filtered"}), depends_on=None, **settings):
    part = object.__new__(DatavzrdStepPart)
    depends_on = {"variants": "fdr", **(depends_on or {})}
    part.config = VariantReport(depends_on=depends_on, tool="datavzrd", **settings)
    part.parent = SimpleNamespace(
        task_name="report",
        config=part.config,
        get_upstream_paths=_upstream_paths,
        upstream_signature=lambda field: DataSignature(DataType.VARIANTS, frozenset(tags)),
        effective_group_by="library",
    )
    return part


def test_the_report_reads_the_vcf_of_the_entity_from_the_upstream_task():
    part = _part(page_size=10)
    wildcards = SimpleNamespace(library_name="P001-N1-DNA1-WES1")

    assert part._get_input_files_run(wildcards) == {
        "vcf": "tasks/fdr/work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz",
        "vcf_tbi": "tasks/fdr/work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz.tbi",
    }
    assert part._get_params_run(wildcards) == {
        "title": "report: calls of P001-N1-DNA1-WES1",
        "page_size": 10,
        "max_in_memory_rows": 50000,
        "linkouts": [
            "gnomad",
            "clinvar",
            "dbsnp",
            "ensembl",
            "genecards",
            "omim",
            "sequence_ontology",
            "varlociraptor_inspect",
        ],
        "assembly": None,
        "igv": None,
        "pedigree": [],
    }


def test_the_steps_extra_arguments_reach_the_datavzrd_command():
    part = _part(extra_arguments=["--debug"])

    assert part.tool_section() is part.config
    params = part.with_extra_arguments("run", part._get_params_run)
    assert params(SimpleNamespace(library_name="P001-N1-DNA1-WES1"))["extra_arguments"] == [
        "--debug"
    ]


def test_the_pedigree_of_a_related_entity_names_every_donor_with_parents_and_status():
    import pandas as pd

    part = _part()
    part.parent.effective_group_by = "related"
    part.parent.pedigree_donors = lambda entity, include_background=False: pd.DataFrame(
        [
            {
                "donor_name": "idx",
                "father_name": "fa",
                "mother_name": "mo",
                "sex": "female",
                "library_name": "idx-N1-DNA1-WES1",
            },
            {
                "donor_name": "fa",
                "father_name": "0",
                "mother_name": "0",
                "sex": "male",
                "library_name": "fa-N1-DNA1-WES1",
            },
            {
                "donor_name": "mo",
                "father_name": "0",
                "mother_name": "gm",
                "sex": "female",
                "library_name": "mo-N1-DNA1-WES1",
            },
            {
                "donor_name": "gm",
                "father_name": "0",
                "mother_name": "0",
                "sex": "female",
                "library_name": "",
            },
        ]
    )
    part.parent.build_library_dataframe = lambda: pd.DataFrame(
        [
            {"library_name": "idx-N1-DNA1-WES1", "disease_state": 2},
            {"library_name": "fa-N1-DNA1-WES1", "disease_state": 1},
            {"library_name": "mo-N1-DNA1-WES1", "disease_state": 0},
        ]
    )

    members = part.pedigree("FAM_idx")

    assert members == [
        {
            "name": "idx-N1-DNA1-WES1",
            "label": "idx",
            "father": "fa-N1-DNA1-WES1",
            "mother": "mo-N1-DNA1-WES1",
            "sex": "female",
            "affected": True,
            "sample": "idx-N1-DNA1-WES1",
            "event_suffix": "idx",
        },
        {
            "name": "fa-N1-DNA1-WES1",
            "label": "fa",
            "father": "0",
            "mother": "0",
            "sex": "male",
            "affected": False,
            "sample": "fa-N1-DNA1-WES1",
            "event_suffix": "fa",
        },
        {
            "name": "mo-N1-DNA1-WES1",
            "label": "mo",
            "father": "0",
            "mother": "gm",
            "sex": "female",
            "affected": None,
            "sample": "mo-N1-DNA1-WES1",
            "event_suffix": "mo",
        },
        {
            "name": "gm",
            "label": "gm",
            "father": "0",
            "mother": "0",
            "sex": "female",
            "affected": None,
            "sample": None,
            "event_suffix": None,
        },
    ]
    assert part._get_params_run(SimpleNamespace(library_name="FAM_idx"))["pedigree"] == members
    assert _part().pedigree("P001-N1-DNA1-WES1") == []  # per library: no pedigree


def test_the_assembly_is_the_tag_of_the_upstream_variants():
    tagged = _part(tags={"germline", "annotated", "mehari", "GRCh37"}, linkouts=["gnomad"])

    assert tagged.assembly() == "GRCh37"
    params = tagged._get_params_run(SimpleNamespace(library_name="FAM_P1"))
    assert (params["assembly"], params["linkouts"]) == ("GRCh37", ["gnomad"])
    assert _part().assembly() is None


def test_the_report_writes_the_calls_table_and_the_site():
    outputs = _part()._get_output_files_run()

    assert outputs["tsv"] == "work/{library_name}/report/calls/{library_name}.calls.tsv"
    assert outputs["tsv_md5"] == outputs["tsv"] + ".md5"
    assert is_flagged(outputs["html"], "directory")
    assert str(outputs["html"]) == "work/{library_name}/report/calls/html"
    assert DatavzrdStepPart.result_keys == ("tsv", "tsv_md5")
    spec = DatavzrdStepPart.report_outputs["run"]["html"]
    assert (spec.subcategory, spec.htmlindex) == ("calls", "index.html")
    assert _part()._get_log_file_run() == {
        "log": "work/{library_name}/log/{library_name}.calls_report.log"
    }
    assert VariantReportWorkflow.get_output_paths(None, library_name="FAM_P1") == {
        "tsv": "work/FAM_P1/report/calls/FAM_P1.calls.tsv"
    }


def _part_with_igv(matched_normal=None, **igv):
    part = _part(depends_on={"alignments": "mapping", "reference": "genome"}, igv=igv)
    rows = [{"library_name": "P001-N1-DNA1-WES1", "matched_normal_lib": ""}]
    if matched_normal is not None:
        rows = [{"library_name": "P001-T1-DNA1-WES1", "matched_normal_lib": matched_normal}]
    normals = {row["library_name"]: row["matched_normal_lib"] for row in rows}
    part.parent.matched_normal = lambda library: normals.get(library, "")
    part.parent.member_libraries = lambda entity: [rows[0]["library_name"]]
    part.parent.reference_inputs = lambda: {
        "reference": "genome.fa",
        "reference_fai": "genome.fa.fai",
        "reference_dict": "genome.dict",
    }
    return part


def test_igv_views_read_the_bams_of_the_members_and_their_matched_normals():
    part = _part_with_igv(matched_normal="P001-N1-DNA1-WES1", top=5, views=["somatic_tumor_high"])
    wildcards = SimpleNamespace(library_name="P001-T1-DNA1-WES1")

    inputs = part._get_input_files_run(wildcards)
    params = part._get_params_run(wildcards)

    assert part.igv_libraries("P001-T1-DNA1-WES1") == ["P001-T1-DNA1-WES1", "P001-N1-DNA1-WES1"]
    assert inputs["vcf"] == "tasks/fdr/work/P001-T1-DNA1-WES1/out/P001-T1-DNA1-WES1.vcf.gz"
    assert inputs["bam"] == [
        "tasks/mapping/work/P001-T1-DNA1-WES1/out/P001-T1-DNA1-WES1.bam",
        "tasks/mapping/work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.bam",
    ]
    assert inputs["bai"][1].endswith("P001-N1-DNA1-WES1.bam.bai")
    assert inputs["reference"] == "genome.fa" and inputs["reference_dict"] == "genome.dict"
    assert params["igv"] == {
        "top": 5,
        "views": ["somatic_tumor_high"],
        "flanking": 300,
        "max_span": 10000,
        "max_size_mb": 50,
        "extra_arguments": [],
        "libraries": ["P001-T1-DNA1-WES1", "P001-N1-DNA1-WES1"],
    }


def test_without_igv_the_report_reads_only_the_vcf():
    part = _part()
    wildcards = SimpleNamespace(library_name="P001-N1-DNA1-WES1")

    assert set(part._get_input_files_run(wildcards)) == {"vcf", "vcf_tbi"}
    assert part._get_params_run(wildcards)["igv"] is None


def test_igv_views_need_the_alignments_and_the_genome():
    settings = {"depends_on": {"variants": "fdr"}, "tool": "datavzrd", "igv": {}}
    with pytest.raises(
        ValidationError, match="igv needs depends_on.alignments and depends_on.reference"
    ):
        VariantReport(**settings)
    settings["depends_on"] |= {"alignments": "mapping", "reference": "genome"}
    assert VariantReport(**settings).igv.top == 20
    # Without IGV views, nothing reads the BAMs
    del settings["igv"]
    with pytest.raises(ValidationError, match="depends_on.alignments is read only by igv"):
        VariantReport(**settings)
