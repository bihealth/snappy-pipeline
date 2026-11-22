# -*- coding: utf-8 -*-
"""Tests of the shared calls report step part of the VCF steps."""

from types import SimpleNamespace

from snakemake.io import is_flagged

from snappy_pipeline.models.calls_report import CallsReport
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.calls_report import CallsReportStepPart


def _paths(config, library_name):
    return {
        "vcf": f"work/{library_name}/out/{library_name}.vcf.gz",
        "vcf_tbi": f"work/{library_name}/out/{library_name}.vcf.gz.tbi",
    }


def _part(tags=frozenset({"germline", "filtered"}), **settings):
    part = object.__new__(CallsReportStepPart)
    signatures = {"fdr": (DataSignature(DataType.VARIANTS, frozenset(tags)),)}
    part.parent = SimpleNamespace(
        task_name="fdr",
        config=SimpleNamespace(calls_report=CallsReport(**settings)),
        get_output_paths=_paths,
        project=SimpleNamespace(signatures=signatures),
        effective_group_by="library",
    )
    return part


def test_the_report_reads_the_tasks_own_vcf_of_the_entity():
    part = _part(page_size=10)
    wildcards = SimpleNamespace(library_name="P001-N1-DNA1-WES1")

    assert part._get_input_files_run(wildcards) == {
        "vcf": "work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz",
        "vcf_tbi": "work/P001-N1-DNA1-WES1/out/P001-N1-DNA1-WES1.vcf.gz.tbi",
    }
    assert part._get_params_run(wildcards) == {
        "title": "fdr: calls of P001-N1-DNA1-WES1",
        "page_size": 10,
        "max_in_memory_rows": 50000,
        "linkouts": [
            "gnomad",
            "clinvar",
            "dbsnp",
            "ensembl",
            "genecards",
            "omim",
            "varlociraptor_inspect",
        ],
        "assembly": None,
        "igv": None,
        "pedigree": [],
    }


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


def test_the_assembly_is_the_tag_of_the_tasks_own_variants():
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
    assert CallsReportStepPart.result_keys == ("tsv", "tsv_md5")
    spec = CallsReportStepPart.report_outputs["run"]["html"]
    assert (spec.subcategory, spec.htmlindex) == ("calls", "index.html")
    assert _part()._get_log_file_run() == {
        "log": "work/{library_name}/log/{library_name}.calls_report.log"
    }


def _part_with_igv(matched_normal=None, **igv):
    part = _part(igv=igv)
    rows = [{"library_name": "P001-N1-DNA1-WES1", "matched_normal_lib": ""}]
    if matched_normal is not None:
        rows = [{"library_name": "P001-T1-DNA1-WES1", "matched_normal_lib": matched_normal}]
    normals = {row["library_name"]: row["matched_normal_lib"] for row in rows}
    part.parent.matched_normal = lambda library: normals.get(library, "")
    part.parent.member_libraries = lambda entity: [rows[0]["library_name"]]
    part.parent.get_upstream_paths = lambda field, library_name: SimpleNamespace(
        bam=f"tasks/mapping/work/{library_name}/out/{library_name}.bam",
        bai=f"tasks/mapping/work/{library_name}/out/{library_name}.bam.bai",
    )
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
    import pytest
    from pydantic import ValidationError

    from snappy_pipeline.workflows.variant_filtration.model import VariantFiltration

    settings = {
        "depends_on": {"variants": "calling"},
        "tool": "varlociraptor",
        "varlociraptor": {"events": ["present"]},
        "calls_report": {"igv": {}},
    }
    with pytest.raises(
        ValidationError,
        match="calls_report.igv needs depends_on.alignments and depends_on.reference",
    ):
        VariantFiltration(**settings)
    settings["depends_on"] |= {"alignments": "mapping", "reference": "genome"}
    assert VariantFiltration(**settings).calls_report.igv.top == 20
