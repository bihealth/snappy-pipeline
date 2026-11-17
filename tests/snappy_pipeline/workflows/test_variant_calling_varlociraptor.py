"""Tests for varlociraptor calling in ``variant_calling``: the scenario presets, the libraries
observed, the settings passed to the wrappers and the contract of the tool."""

from types import SimpleNamespace

import pandas as pd
import pytest
import yaml
from pydantic import ValidationError
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.workflows.abstract import MATCHED_NORMAL_REQUIREMENT
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_calling import (
    VariantCallingWorkflow,
    VarlociraptorStepPart,
)
from snappy_pipeline.workflows.variant_calling.model import VariantCalling, Varlociraptor
from snappy_pipeline.workflows.variant_calling.scenario import (
    germline_scenario,
    render_template,
    somatic_scenario,
    to_yaml,
)

SPECIES = Varlociraptor(scenario="germline").species.model_dump()


def _row(name, sex="unknown", father="0", mother="0", **extra):
    return {
        "library_name": name,
        "sex": sex,
        "father_name": father,
        "mother_name": mother,
        "matched_normal_lib": "",
        "related_group": name.split("-")[0],
        **extra,
    }


#: A trio, the child's RNA library, a tumor with its normal (purity 0.8) and a relapse, and a
#: tumor without a normal
ROWS = {
    "child-DNA1": _row("child-DNA1", "male", "father-DNA1", "mother-DNA1"),
    "father-DNA1": _row("father-DNA1", "male"),
    "mother-DNA1": _row("mother-DNA1", "female"),
    "child-RNA1": _row("child-RNA1", "male", "father-DNA1", "mother-DNA1"),
    "P1-T1-DNA1": _row("P1-T1-DNA1", matched_normal_lib="P1-N1-DNA1", purity=0.8),
    "P1-T2-DNA1": _row("P1-T2-DNA1", matched_normal_lib="P1-N1-DNA1", purity=0.5),
    "P1-N1-DNA1": _row("P1-N1-DNA1"),
    "P2-T1-DNA1": _row("P2-T1-DNA1", purity=""),
}
TRIO = [ROWS["child-DNA1"], ROWS["father-DNA1"], ROWS["mother-DNA1"]]
#: The tumor libraries a somatic task takes: one per donor, or with the relapse of P1
TUMORS = pd.DataFrame({"library_name": ["P1-T1-DNA1", "P2-T1-DNA1"], "related_group": ["P1", "P2"]})
WITH_RELAPSE = pd.DataFrame(
    {
        "library_name": ["P1-T1-DNA1", "P1-T2-DNA1", "P2-T1-DNA1"],
        "related_group": ["P1", "P1", "P2"],
    }
)


# --- the presets -------------------------------------------------------------------------


def test_the_germline_preset_of_a_trio_has_mendelian_inheritance_and_a_denovo_event():
    scenario = germline_scenario(TRIO, SPECIES)

    assert scenario["species"]["germline-mutation-rate"] == 0.001
    assert scenario["samples"]["child-DNA1"] == {
        "sex": "male",
        "inheritance": {"mendelian": {"from": ["mother-DNA1", "father-DNA1"]}},
    }
    assert scenario["samples"]["mother-DNA1"] == {"sex": "female"}
    assert scenario["events"] == {
        "denovo": "child-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father-DNA1:0.0",
        "present": (
            "(child-DNA1:]0.0,1.0] | father-DNA1:]0.0,1.0] | mother-DNA1:]0.0,1.0]) & !$denovo"
        ),
    }


def test_the_germline_preset_of_a_quad_ors_the_denovo_events_of_both_children():
    sibling = _row("sibling-DNA1", "female", "father-DNA1", "mother-DNA1")

    scenario = germline_scenario([*TRIO, sibling], SPECIES)

    assert scenario["events"]["denovo"] == (
        "(child-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father-DNA1:0.0)"
        " | (sibling-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father-DNA1:0.0)"
    )


def test_the_germline_preset_without_a_child_has_only_present_and_unknown_sex_is_diploid():
    pair = [_row("a-DNA1", "unknown"), _row("b-DNA1", "female")]

    scenario = germline_scenario(pair, SPECIES)

    assert scenario["samples"] == {"a-DNA1": {"ploidy": 2}, "b-DNA1": {"sex": "female"}}
    assert scenario["events"] == {
        "denovo": False,
        "present": "a-DNA1:]0.0,1.0] | b-DNA1:]0.0,1.0]",
    }


def test_the_germline_preset_models_several_generations():
    grandparents = [_row("gf-DNA1", "male"), _row("gm-DNA1", "female")]
    mother = _row("mother-DNA1", "female", "gf-DNA1", "gm-DNA1")

    members = [ROWS["child-DNA1"], ROWS["father-DNA1"], mother, *grandparents]

    scenario = germline_scenario(members, SPECIES)

    assert scenario["samples"]["mother-DNA1"]["inheritance"] == {
        "mendelian": {"from": ["gm-DNA1", "gf-DNA1"]}
    }
    assert scenario["samples"]["child-DNA1"]["inheritance"]["mendelian"]["from"] == [
        "mother-DNA1",
        "father-DNA1",
    ]
    assert scenario["events"]["denovo"].count(" | ") == 1


@pytest.mark.parametrize(
    "member, message",
    [
        (_row("duo-DNA1", "male", "0", "mother-DNA1"), "parents mother-DNA1 in the sample sheet"),
        (
            _row("child-DNA1", "male", "father-DNA1", "mother-DNA1"),
            "parents father-DNA1 and mother-DNA1 in the sample sheet",
        ),
    ],
)
def test_a_member_with_a_parent_outside_the_task_is_rejected(member, message):
    with pytest.raises(ValueError, match=message):
        germline_scenario([member, ROWS["mother-DNA1"]], SPECIES)


def test_the_somatic_preset_models_the_purity_as_contamination_by_the_normal():
    scenario = somatic_scenario("P1-T1-DNA1", "P1-N1-DNA1", 0.8)

    assert scenario["samples"]["P1-T1-DNA1"]["contamination"] == {
        "by": "P1-N1-DNA1",
        "fraction": 0.2,
    }
    assert scenario["events"]["somatic_tumor_high"] == "P1-T1-DNA1:[0.3,1.0] & P1-N1-DNA1:0.0"
    assert scenario["events"]["germline_het"] == "P1-N1-DNA1:0.5"


def test_ffpe_libraries_get_the_fragment_end_bias_model():
    scenario = somatic_scenario("T", "N", None, frozenset({"T"}))
    assert scenario["samples"]["T"]["fragment_end_bias"] is True
    assert "fragment_end_bias" not in scenario["samples"]["N"]

    ffpe_child = {**ROWS["child-DNA1"], "is_ffpe": "Y"}
    scenario = germline_scenario([ffpe_child, ROWS["father-DNA1"], ROWS["mother-DNA1"]], SPECIES)
    assert scenario["samples"]["child-DNA1"]["fragment_end_bias"] is True
    assert "fragment_end_bias" not in scenario["samples"]["father-DNA1"]


def test_the_step_part_reads_the_ffpe_status_of_the_sheet():
    ROWS["P2-T1-DNA1"]["is_ffpe"] = "Y"
    part = _part(scenario="somatic")

    scenario = yaml.safe_load(part.scenario("P2-T1-DNA1"))

    assert scenario["samples"]["P2-T1-DNA1"]["fragment_end_bias"] is True


def test_the_somatic_preset_of_a_tumor_alone_keeps_an_unobserved_normal_and_no_contamination():
    scenario = somatic_scenario("P2-T1-DNA1", None, None)

    assert list(scenario["samples"]) == ["P2-T1-DNA1", "normal"]
    assert "contamination" not in scenario["samples"]["P2-T1-DNA1"]
    assert scenario["events"]["somatic_normal"] == "normal:]0.0,0.5["


def test_the_yaml_quotes_the_formulas_and_keeps_the_order():
    text = to_yaml(somatic_scenario("T", "N", 0.75))

    assert yaml.safe_load(text)["samples"]["T"]["contamination"]["fraction"] == 0.25
    assert text.startswith("samples:\n  T:\n")


def test_a_template_is_rendered_with_the_samples_and_the_settings():
    template = (
        "samples:\n"
        "  ?for sample in samples:\n"
        "    ?sample['library_name']:\n"
        "      universe: '[0.0,1.0]'\n"
        "events:\n"
        '  present: ?f"{library_name}:]0.0,1.0]"\n'
        "  purity: ?settings['purity']\n"
    )

    scenario = render_template(template, [ROWS["P1-T1-DNA1"]], "P1-T1-DNA1", {"purity": 0.8})

    assert scenario["samples"] == {"P1-T1-DNA1": {"universe": "[0.0,1.0]"}}
    assert scenario["events"] == {"present": "P1-T1-DNA1:]0.0,1.0]", "purity": 0.8}


def test_a_template_without_samples_or_events_is_rejected():
    with pytest.raises(ValueError, match="needs samples and events"):
        render_template("samples: {}\n", [], "x", {})


# --- the step part -------------------------------------------------------------------------


def _upstream(field, **kwargs):
    name = kwargs.get("library_name", "")
    if field == "candidates":
        return SimpleNamespace(vcf=f"candidates/{name}.vcf.gz", vcf_tbi=f"candidates/{name}.tbi")
    if field == "reference":
        return SimpleNamespace(fai="ref.fa.fai")
    return SimpleNamespace(bam=f"{name}.bam", bai=f"{name}.bam.bai")


def _part(group_by="library", libraries=TUMORS, **settings):
    part = object.__new__(VarlociraptorStepPart)
    part.config = SimpleNamespace(varlociraptor=Varlociraptor(**settings))
    part.parent = SimpleNamespace(
        task_name="calling",
        effective_group_by=group_by,
        libraries=libraries,
        primary_dna_libraries=lambda: pd.DataFrame(
            {"library_name": ["child-DNA1", "father-DNA1", "mother-DNA1"]}
        ),
        get_related_libraries=lambda name: ["father-DNA1", "child-DNA1", "mother-DNA1"],
        library_row=lambda name: ROWS[name],
        get_upstream_paths=_upstream,
        reference_inputs=lambda: {"reference": "ref.fa"},
        pedigree_donors=lambda index: pd.DataFrame(),
    )
    part._scenarios = {}
    return part


def test_the_germline_scenario_observes_the_pedigree_or_the_library_alone():
    assert _part("related", scenario="germline").members("child-DNA1") == [
        "child-DNA1",
        "father-DNA1",
        "mother-DNA1",
    ]
    assert _part("library", scenario="germline").members("P1-T1-DNA1") == ["P1-T1-DNA1"]


@pytest.mark.parametrize(
    "mode, expected",
    [
        ("automatic", ["P1-T1-DNA1", "P1-N1-DNA1"]),
        ("paired", ["P1-T1-DNA1", "P1-N1-DNA1"]),
        ("tumor_only", ["P1-T1-DNA1"]),
    ],
)
def test_the_somatic_scenario_observes_the_normal_by_the_tumor_normal_mode(mode, expected):
    part = _part(scenario="somatic", tumor_normal_mode=mode)

    assert part.members("P1-T1-DNA1") == expected
    assert part.members("P2-T1-DNA1") == ["P2-T1-DNA1"]


def test_the_somatic_preset_rejects_a_donor_with_several_tumors_in_the_task():
    part = _part(scenario="somatic", libraries=WITH_RELAPSE)

    with pytest.raises(ValueError, match="tumors P1-T1-DNA1, P1-T2-DNA1 in the task"):
        part.scenario("P1-T1-DNA1")
    assert "P2-T1-DNA1" in part.scenario("P2-T1-DNA1")


def test_a_somatic_template_with_related_grouping_observes_all_tumors_of_the_donor(tmp_path):
    template = tmp_path / "relapse.yaml"
    template.write_text(
        "samples:\n"
        "  ?for sample in samples:\n"
        "    ?sample['library_name']:\n"
        "      universe: '[0.0,1.0]'\n"
        "events:\n"
        "  present: ?' | '.join(f\"{s['library_name']}:]0.0,1.0]\" for s in samples)\n"
    )
    part = _part("related", WITH_RELAPSE, scenario="somatic", template=str(template))

    assert part.members("P1-T1-DNA1") == ["P1-T1-DNA1", "P1-T2-DNA1", "P1-N1-DNA1"]
    scenario = yaml.safe_load(part.scenario("P1-T1-DNA1"))
    assert list(scenario["samples"]) == ["P1-T1-DNA1", "P1-T2-DNA1", "P1-N1-DNA1"]


def test_the_purity_comes_from_the_sheet_column_and_is_checked():
    part = _part(scenario="somatic")

    assert part.purity("P1-T1-DNA1") == 0.8
    assert part.purity("P2-T1-DNA1") is None
    assert part.purity("P1-N1-DNA1") is None
    ROWS["bad"] = _row("bad", purity="high")
    with pytest.raises(ConfigError, match="purity of bad is 'high'"):
        part.purity("bad")


def test_the_call_params_hold_the_sample_names_and_the_rendered_scenario():
    part = _part(scenario="somatic")
    wildcards = Wildcards(fromdict={"library_name": "P1-T1-DNA1", "scatteritem": "1-of-2"})

    params = part._get_params_call(wildcards)
    inputs = part._get_input_files_call(wildcards)

    assert params["obs_names"] == ["P1-T1-DNA1", "P1-N1-DNA1"]
    scenario = yaml.safe_load(params["scenario"])
    assert scenario["samples"]["P1-T1-DNA1"]["contamination"]["fraction"] == 0.2
    assert inputs["observations"] == [
        "work/P1-T1-DNA1/out/P1-T1-DNA1/observations/P1-T1-DNA1.1-of-2.bcf",
        "work/P1-T1-DNA1/out/P1-T1-DNA1/observations/P1-N1-DNA1.1-of-2.bcf",
    ]
    assert part.scenario("P1-T1-DNA1") is params["scenario"]  # rendered once


def test_the_germline_scenario_of_the_step_part_takes_the_pedigree_rows():
    part = _part("related", scenario="germline")

    scenario = yaml.safe_load(part.scenario("child-DNA1"))

    assert set(scenario["samples"]) == {"child-DNA1", "father-DNA1", "mother-DNA1"}
    assert "denovo" in scenario["events"]


def test_a_template_file_replaces_the_preset(tmp_path):
    template = tmp_path / "scenario.yaml"
    template.write_text(
        "samples:\n  ?library_name:\n    universe: '[0.0,1.0]'\n"
        'events:\n  present: ?f"{library_name}:]0.0,1.0]"\n'
    )
    part = _part(scenario="somatic", tumor_normal_mode="tumor_only", template=str(template))

    scenario = yaml.safe_load(part.scenario("P2-T1-DNA1"))

    assert scenario == {
        "samples": {"P2-T1-DNA1": {"universe": "[0.0,1.0]"}},
        "events": {"present": "P2-T1-DNA1:]0.0,1.0]"},
    }


def test_a_template_must_name_every_observed_library(tmp_path):
    template = tmp_path / "scenario.yaml"
    template.write_text("samples:\n  tumor: {}\nevents:\n  present: 'tumor:0.5'\n")
    part = _part(scenario="somatic", template=str(template))

    with pytest.raises(ValueError, match="lacks the samples P1-N1-DNA1, P1-T1-DNA1"):
        part.scenario("P1-T1-DNA1")


def test_preprocess_reads_a_chunk_of_the_candidates_and_the_alignment_properties():
    part = _part(scenario="somatic")
    wildcards = Wildcards(
        fromdict={"library_name": "P1-T1-DNA1", "member": "P1-N1-DNA1", "scatteritem": "2-of-2"}
    )

    inputs = part._get_input_files_preprocess(wildcards)

    assert inputs["candidates"] == "work/P1-T1-DNA1/out/P1-T1-DNA1/candidates/2-of-2.bcf"
    assert (
        inputs["alignment_properties"] == "work/P1-N1-DNA1/out/P1-N1-DNA1.alignment_properties.json"
    )
    assert inputs["bam"] == "P1-N1-DNA1.bam"
    assert part._get_params_preprocess(wildcards) == {"max_depth": 200, "extra_arguments": []}


# --- the model and the contract ------------------------------------------------------------


def test_related_grouping_of_the_somatic_scenario_needs_a_template():
    with pytest.raises(ValidationError, match="group_by: related with the somatic scenario"):
        VariantCalling(
            depends_on={"alignments": "m", "reference": "g", "candidates": "c"},
            tool="varlociraptor",
            varlociraptor={"scenario": "somatic"},
            group_by="related",
        )
    VariantCalling(
        depends_on={"alignments": "m", "reference": "g", "candidates": "c"},
        tool="varlociraptor",
        varlociraptor={"scenario": "somatic", "template": "relapse.yaml"},
        group_by="related",
    )


def test_changed_settings_of_the_other_scenario_are_rejected_and_defaults_pass():
    with pytest.raises(ValidationError, match="tumor_normal_mode: not used by the germline"):
        Varlociraptor(scenario="germline", tumor_normal_mode="paired")
    with pytest.raises(ValidationError, match="species: not used by the somatic"):
        Varlociraptor(scenario="somatic", species={"heterozygosity": 0.01})
    # A written-out default config names every field
    Varlociraptor(scenario="somatic", species={}, tumor_normal_mode="automatic")
    Varlociraptor(scenario="germline", species={"ploidy": {"male": {"all": 2}}})


def _config(**settings):
    return VariantCalling.model_construct(
        tool="varlociraptor", varlociraptor=Varlociraptor(**settings)
    )


def test_varlociraptor_produces_the_calls_of_its_scenario_tagged_with_the_tool():
    candidates = {
        "candidates": DataSignature(DataType.VARIANTS, frozenset({"candidates", "snv", "indel"}))
    }
    germline = VariantCallingWorkflow.task_produces(_config(scenario="germline"), candidates)
    somatic = VariantCallingWorkflow.task_produces(_config(scenario="somatic"), candidates)

    assert germline == (
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel", "varlociraptor"})),
    )
    assert somatic == (
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "varlociraptor"})),
    )
    assert VariantCallingWorkflow.filter_ids(_config(scenario="somatic")) == frozenset()


def test_the_scenario_sets_the_grouping_and_the_libraries_taken():
    assert VariantCallingWorkflow.default_group_by(_config(scenario="germline")) == "related"
    assert VariantCallingWorkflow.default_group_by(_config(scenario="somatic")) is None
    assert VariantCallingWorkflow.task_library_requirement(_config(scenario="germline")) is None
    assert (
        VariantCallingWorkflow.task_library_requirement(_config(scenario="somatic"))
        == "role == 'tumor' and extraction_type == 'dna'"
    )
    assert (
        VariantCallingWorkflow.task_library_requirement(
            _config(scenario="somatic", tumor_normal_mode="paired")
        )
        == MATCHED_NORMAL_REQUIREMENT
    )


# --- the candidate chunks ------------------------------------------------------------------


def _chunker():
    import importlib.util
    from pathlib import Path

    path = (
        Path(__file__).resolve().parents[3]
        / "snappy_wrappers/wrappers/varlociraptor/split_candidates/wrapper.py"
    )
    spec = importlib.util.spec_from_file_location("split_candidates_wrapper", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_chunks_hold_equal_depth_and_keep_breakend_events_together():
    chunker = _chunker()
    weights = [10, 10, 10, 10, 10, 10]
    groups = [None, ("event", "e1"), None, None, ("event", "e1"), None]

    assert chunker.plan(weights, groups, 3) == [0, 0, 1, 1, 0, 2]


def test_a_candidate_without_depth_weighs_one_read():
    chunker = _chunker()

    assert chunker.plan([1, 1, 1, 1], [None] * 4, 2) == [0, 0, 1, 1]
    assert chunker.plan([], [], 2) == []


def test_the_chunks_are_written_in_input_order(tmp_path):
    pysam = pytest.importorskip("pysam")
    chunker = _chunker()
    vcf = tmp_path / "candidates.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        "##contig=<ID=1,length=1000>\n##contig=<ID=2,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t10\tbnd_a\tN\tN[2:500[\t.\t.\tDP=30;SVTYPE=BND;MATEID=bnd_b\n"
        "1\t20\t.\tA\tC\t.\t.\tDP=30\n"
        "1\t30\t.\tA\tG\t.\t.\tDP=30\n"
        "2\t500\tbnd_b\tN\t]1:10]N\t.\t.\tDP=30;SVTYPE=BND;MATEID=bnd_a\n"
    )
    outputs = [str(tmp_path / f"chunk{i}.bcf") for i in range(2)]

    chunker.main(str(vcf), outputs)

    chunks = [[(r.chrom, r.pos) for r in pysam.VariantFile(path)] for path in outputs]
    assert chunks == [[("1", 10), ("1", 20), ("2", 500)], [("1", 30)]]


def test_the_scatter_splits_the_candidates_of_the_entity_and_preprocess_reads_a_chunk():
    part = _part(scenario="somatic")
    part._scattergather = lambda kind: (
        lambda template: [template.format(scatteritem=f"{i}-of-2") for i in (1, 2)]
    )

    assert part._get_input_files_scatter(Wildcards(fromdict={"library_name": "P1-T1-DNA1"})) == {
        "candidates": "candidates/P1-T1-DNA1.vcf.gz"
    }
    assert part._get_output_files_scatter()["candidates"] == [
        "work/{library_name}/out/{library_name}/candidates/1-of-2.bcf",
        "work/{library_name}/out/{library_name}/candidates/2-of-2.bcf",
    ]
    wildcards = Wildcards(
        fromdict={"library_name": "P1-T1-DNA1", "member": "P1-N1-DNA1", "scatteritem": "2-of-2"}
    )
    inputs = part._get_input_files_preprocess(wildcards)
    assert inputs["candidates"] == "work/P1-T1-DNA1/out/P1-T1-DNA1/candidates/2-of-2.bcf"
    assert "region" not in inputs


def test_sv_candidates_give_sv_calls():
    candidates = DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"}))

    produced = VariantCallingWorkflow.task_produces(
        _config(scenario="germline"), {"candidates": candidates}
    )

    assert produced == (
        DataSignature(DataType.VARIANTS, frozenset({"germline", "sv", "varlociraptor"})),
    )


# --- unsequenced relatives -----------------------------------------------------------------


def _donors(*rows):
    """A donor table: (donor, father, mother, sex, library)"""
    return pd.DataFrame(
        [
            {
                "donor_name": donor,
                "father_name": father,
                "mother_name": mother,
                "sex": sex,
                "library_name": library,
                "pedigree": "index",
            }
            for donor, father, mother, sex, library in rows
        ]
    )


def _pedigree_part(donors, members):
    part = _part("related", scenario="germline")
    part.parent.pedigree_donors = lambda index: donors
    rows = {
        row["library_name"]: _row(row["library_name"], row["sex"])
        for row in donors.to_dict("records")
        if row["library_name"]
    }
    part.parent.library_row = lambda name: rows[name]
    part.members = lambda library_name: members
    return part


def test_a_sequenced_trio_adds_no_unobserved_sample():
    donors = _donors(
        ("child", "father", "mother", "male", "child-DNA1"),
        ("father", "0", "0", "male", "father-DNA1"),
        ("mother", "0", "0", "female", "mother-DNA1"),
    )
    part = _pedigree_part(donors, ["child-DNA1", "father-DNA1", "mother-DNA1"])

    samples = part.pedigree_samples("child-DNA1")

    assert [s["library_name"] for s in samples] == ["child-DNA1", "father-DNA1", "mother-DNA1"]
    assert (samples[0]["father_name"], samples[0]["mother_name"]) == ("father-DNA1", "mother-DNA1")
    assert all(s.get("observed", True) for s in samples)


def test_the_unsequenced_parent_of_a_duo_joins_as_an_unobserved_sample():
    donors = _donors(
        ("child", "father", "mother", "male", "child-DNA1"),
        ("father", "0", "0", "male", ""),
        ("mother", "0", "0", "female", "mother-DNA1"),
    )
    part = _pedigree_part(donors, ["child-DNA1", "mother-DNA1"])

    samples = part.pedigree_samples("child-DNA1")
    scenario = germline_scenario(samples, SPECIES)

    assert samples[-1] == {
        "library_name": "father",
        "sex": "male",
        "father_name": "0",
        "mother_name": "0",
        "observed": False,
    }
    assert scenario["samples"]["child-DNA1"]["inheritance"] == {
        "mendelian": {"from": ["mother-DNA1", "father"]}
    }
    assert scenario["samples"]["father"] == {"sex": "male"}
    assert scenario["events"] == {
        "denovo": "child-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father:0.0",
        "present": "(child-DNA1:]0.0,1.0] | mother-DNA1:]0.0,1.0]) & !$denovo",
    }


def test_cousins_are_linked_through_their_unsequenced_parents_and_grandparents():
    donors = _donors(
        ("gf", "0", "0", "male", ""),
        ("gm", "0", "0", "female", ""),
        ("a", "gf", "gm", "male", ""),
        ("a_wife", "0", "0", "female", ""),
        ("b", "gf", "gm", "female", ""),
        ("b_husband", "0", "0", "male", ""),
        ("c1", "a", "a_wife", "male", "c1-DNA1"),
        ("c2", "b_husband", "b", "female", "c2-DNA1"),
        ("stranger", "0", "0", "male", ""),
    )
    part = _pedigree_part(donors, ["c1-DNA1", "c2-DNA1"])

    samples = part.pedigree_samples("c1-DNA1")

    unobserved = sorted(s["library_name"] for s in samples if not s.get("observed", True))
    assert unobserved == ["a", "a_wife", "b", "b_husband", "gf", "gm"]
    by_name = {s["library_name"]: s for s in samples}
    assert (by_name["c1-DNA1"]["father_name"], by_name["c1-DNA1"]["mother_name"]) == ("a", "a_wife")
    assert (by_name["a"]["father_name"], by_name["a"]["mother_name"]) == ("gf", "gm")


def test_unsequenced_parents_of_a_singleton_are_left_out():
    donors = _donors(
        ("child", "father", "mother", "male", "child-DNA1"),
        ("father", "0", "0", "male", ""),
        ("mother", "0", "0", "female", ""),
    )
    part = _pedigree_part(donors, ["child-DNA1"])

    samples = part.pedigree_samples("child-DNA1")

    assert [s["library_name"] for s in samples] == ["child-DNA1"]
    assert (samples[0]["father_name"], samples[0]["mother_name"]) == ("0", "0")


def test_siblings_share_their_unsequenced_parents():
    donors = _donors(
        ("s1", "father", "mother", "male", "s1-DNA1"),
        ("s2", "father", "mother", "female", "s2-DNA1"),
        ("father", "0", "0", "male", ""),
        ("mother", "0", "0", "female", ""),
    )
    part = _pedigree_part(donors, ["s1-DNA1", "s2-DNA1"])

    scenario = germline_scenario(part.pedigree_samples("s1-DNA1"), SPECIES)

    assert set(scenario["samples"]) == {"s1-DNA1", "s2-DNA1", "father", "mother"}
    assert scenario["events"]["denovo"].count(" | ") == 1


def test_a_sequenced_parent_outside_the_task_is_still_rejected():
    donors = _donors(
        ("child", "father", "mother", "male", "child-DNA1"),
        ("father", "0", "0", "male", "father-DNA1"),
        ("mother", "0", "0", "female", "mother-DNA1"),
    )
    part = _pedigree_part(donors, ["child-DNA1", "mother-DNA1"])

    with pytest.raises(ValueError, match="parents father-DNA1 and mother-DNA1"):
        germline_scenario(part.pedigree_samples("child-DNA1"), SPECIES)
