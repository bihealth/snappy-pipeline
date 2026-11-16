"""Tests for varlociraptor calling in ``variant_calling``: the scenario presets, the libraries
observed, the settings passed to the wrappers and the contract of the tool."""

from types import SimpleNamespace

import pandas as pd
import pytest
import yaml
from pydantic import ValidationError
from snakemake.iocontainers import Wildcards

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
        **extra,
    }


#: A trio, the child's RNA library, a tumor with its normal (purity 0.8) and one without
ROWS = {
    "child-DNA1": _row("child-DNA1", "male", "father-DNA1", "mother-DNA1"),
    "father-DNA1": _row("father-DNA1", "male"),
    "mother-DNA1": _row("mother-DNA1", "female"),
    "child-RNA1": _row("child-RNA1", "male", "father-DNA1", "mother-DNA1"),
    "P1-T1-DNA1": _row("P1-T1-DNA1", matched_normal_lib="P1-N1-DNA1", purity=0.8),
    "P1-N1-DNA1": _row("P1-N1-DNA1"),
    "P2-T1-DNA1": _row("P2-T1-DNA1", purity=""),
}
TRIO = [ROWS["child-DNA1"], ROWS["father-DNA1"], ROWS["mother-DNA1"]]


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
        "denovo": "(child-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father-DNA1:0.0)",
        "present": "!($denovo | $absent)",
    }


def test_the_germline_preset_of_a_quad_ors_the_denovo_events_of_both_children():
    sibling = _row("sibling-DNA1", "female", "father-DNA1", "mother-DNA1")

    scenario = germline_scenario([*TRIO, sibling], SPECIES)

    assert scenario["events"]["denovo"] == (
        "(child-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father-DNA1:0.0)"
        " | (sibling-DNA1:]0.0,1.0] & mother-DNA1:0.0 & father-DNA1:0.0)"
    )


def test_the_germline_preset_without_a_child_has_only_present_and_unknown_sex_is_diploid():
    duo = [_row("a-DNA1", "unknown"), _row("b-DNA1", "female", "x-DNA1", "y-DNA1")]

    scenario = germline_scenario(duo, SPECIES)

    assert scenario["samples"] == {"a-DNA1": {"ploidy": 2}, "b-DNA1": {"sex": "female"}}
    assert scenario["events"] == {"present": "a-DNA1:]0.0,1.0] | b-DNA1:]0.0,1.0]"}


def test_the_somatic_preset_models_the_purity_as_contamination_by_the_normal():
    scenario = somatic_scenario("P1-T1-DNA1", "P1-N1-DNA1", 0.8)

    assert scenario["samples"]["P1-T1-DNA1"]["contamination"] == {
        "by": "P1-N1-DNA1",
        "fraction": 0.2,
    }
    assert scenario["events"]["somatic_tumor_high"] == "P1-T1-DNA1:[0.3,1.0] & P1-N1-DNA1:0.0"
    assert scenario["events"]["germline_het"] == "P1-N1-DNA1:0.5"


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


def _part(group_by="library", **settings):
    part = object.__new__(VarlociraptorStepPart)
    part.config = SimpleNamespace(varlociraptor=Varlociraptor(**settings))
    part.parent = SimpleNamespace(
        task_name="calling",
        effective_group_by=group_by,
        primary_dna_libraries=lambda: pd.DataFrame(
            {"library_name": ["child-DNA1", "father-DNA1", "mother-DNA1"]}
        ),
        get_related_libraries=lambda name: ["father-DNA1", "child-DNA1", "mother-DNA1"],
        library_row=lambda name: ROWS[name],
        get_upstream_paths=_upstream,
        reference_inputs=lambda: {"reference": "ref.fa"},
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


def test_the_purity_comes_from_the_sheet_column_and_is_checked():
    part = _part(scenario="somatic")

    assert part.purity("P1-T1-DNA1") == 0.8
    assert part.purity("P2-T1-DNA1") is None
    assert part.purity("P1-N1-DNA1") is None
    ROWS["bad"] = _row("bad", purity="high")
    with pytest.raises(ValueError, match="purity of bad is 'high'"):
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


def test_preprocess_reads_the_candidates_of_the_window_and_the_alignment_properties():
    part = _part(scenario="somatic")
    wildcards = Wildcards(
        fromdict={"library_name": "P1-T1-DNA1", "member": "P1-N1-DNA1", "scatteritem": "2-of-2"}
    )

    inputs = part._get_input_files_preprocess(wildcards)

    assert inputs["candidates"] == "candidates/P1-T1-DNA1.vcf.gz"
    assert inputs["region"] == "work/scatter/2-of-2.region.bed"
    assert (
        inputs["alignment_properties"] == "work/P1-N1-DNA1/out/P1-N1-DNA1.alignment_properties.json"
    )
    assert inputs["bam"] == "P1-N1-DNA1.bam"
    assert part._get_params_preprocess(wildcards) == {"max_depth": 200, "extra_arguments": []}


# --- the model and the contract ------------------------------------------------------------


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
    germline = VariantCallingWorkflow.task_produces(_config(scenario="germline"), {})
    somatic = VariantCallingWorkflow.task_produces(_config(scenario="somatic"), {})

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
