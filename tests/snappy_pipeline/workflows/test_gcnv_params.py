"""Tests for the tool settings that the gCNV step parts pass to their wrappers."""

from types import SimpleNamespace

import pytest

from snappy_pipeline.workflows.common.gcnv.gcnv_build_model import BuildGcnvModelStepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.gcnv_model.model import GcnvModel
from snappy_pipeline.workflows.sv_calling_wgs.model import Gcnv

WILDCARDS = SimpleNamespace(library_kit="kit", library_name="P001")


def _model_config(**kwargs):
    # The params only read the model settings, so skip the validation of ``depends_on``.
    return GcnvModel.model_construct(**kwargs)


def _build_part(**kwargs):
    part = object.__new__(BuildGcnvModelStepPart)
    part.config = _model_config(intervals="targeted", **kwargs)
    return part


def _run_part(**kwargs):
    part = object.__new__(RunGcnvStepPart)
    part.config = SimpleNamespace(gcnv=Gcnv(**kwargs))
    part.parent = SimpleNamespace(get_task_config=lambda name: _model_config(intervals="wgs"))
    return part


def test_model_defaults_keep_the_former_values():
    config = _model_config()
    assert str(config.interval_merging_rule) == "OVERLAPPING_ONLY"
    assert (config.bin_length, config.padding, config.intervals_per_shard) == (0, 0, 5000)
    assert (config.minimum_gc_content, config.maximum_gc_content) == (0.1, 0.9)
    assert (config.minimum_mappability, config.maximum_mappability) == (0.9, 1.0)
    priors = [config.autosome_ploidy_prior, config.x_ploidy_prior, config.y_ploidy_prior]
    assert ["\t".join(f"{p:g}" for p in prior) for prior in priors] == [
        "0\t0.01\t0.98\t0.01",
        "0.01\t0.49\t0.49\t0.01",
        "0.495\t0.495\t0.01\t0",
    ]


def test_preprocess_intervals_takes_bin_length_for_targeted_intervals():
    assert _build_part()._get_params_preprocess_intervals(WILDCARDS) == {
        "interval_merging_rule": "OVERLAPPING_ONLY",
        "bin_length": 0,
    }


def test_preprocess_intervals_takes_padding_for_wgs_intervals():
    part = _build_part()
    part.config.intervals = "wgs"
    assert part._get_params_preprocess_intervals(WILDCARDS) == {
        "interval_merging_rule": "OVERLAPPING_ONLY",
        "padding": 0,
    }


def test_filter_intervals_gets_the_bounds():
    assert _build_part(minimum_mappability=0.8)._get_params_filter_intervals(WILDCARDS) == {
        "interval_merging_rule": "OVERLAPPING_ONLY",
        "minimum_gc_content": 0.1,
        "maximum_gc_content": 0.9,
        "minimum_mappability": 0.8,
        "maximum_mappability": 1.0,
    }


def test_scatter_intervals_gets_the_shard_size():
    assert _build_part(intervals_per_shard=100)._get_params_scatter_intervals(WILDCARDS) == {
        "intervals_per_shard": 100
    }


def test_contig_ploidy_gets_the_priors():
    params = _build_part(y_ploidy_prior=[0.5, 0.5, 0, 0])._get_params_contig_ploidy(WILDCARDS)
    assert params["y_ploidy_prior"] == [0.5, 0.5, 0, 0]
    assert params["autosome_ploidy_prior"] == [0, 0.01, 0.98, 0.01]
    assert params["interval_merging_rule"] == "OVERLAPPING_ONLY"


@pytest.mark.parametrize("action", ["annotate_gc", "call_cnvs", "coverage"])
def test_model_actions_get_the_merging_rule(action):
    part = _build_part(interval_merging_rule="ALL")
    assert getattr(part, f"_get_params_{action}")(WILDCARDS) == {"interval_merging_rule": "ALL"}


def test_post_germline_calls_of_the_model_gets_the_contigs():
    assert _build_part(allosomal_contigs=["X"])._get_params_post_germline_calls(WILDCARDS) == {
        "autosomal_ref_copy_number": 2,
        "allosomal_contigs": ["X"],
    }


def test_callers_take_the_merging_rule_of_the_model_task():
    part = _run_part()
    assert part._get_params_coverage(WILDCARDS) == {"interval_merging_rule": "OVERLAPPING_ONLY"}


def test_post_germline_calls_of_the_callers_gets_the_contigs():
    assert _run_part(autosomal_ref_copy_number=2)._get_params_post_germline_calls(WILDCARDS) == {
        "autosomal_ref_copy_number": 2,
        "allosomal_contigs": ["X", "Y"],
    }
