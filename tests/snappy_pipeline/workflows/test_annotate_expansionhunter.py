# -*- coding: utf-8 -*-
"""Tests for the Annotate ExpansionHunter module code"""

import json
import os
from shutil import copyfile
import tempfile

import pytest

from snappy_pipeline.workflows.repeat_expansion.annotate_expansionhunter import (
    AnnotateExpansionHunter,
)


@pytest.fixture
def expansionhunter_expected():
    """:return: Returns dictionary with expected results for ExpansionHunter - NORMAL."""
    return {
        "AR": "Normal - repeat AR observed: 24/25. Normal range: 9 - 35. Expansion range: 38 - 62.",
        "FMR1": (
            "Normal - repeat FMR1 observed: 17/17. "
            "Normal range: 6 - 54. Expansion range: 200 - 1000."
        ),
        "RFC1": "Undefined - genotype not defined in ExpansionHunter for repeat RFC1.",
        "TBP": "Undefined - no annotation found for repeat TBP. Observed genotype: 35/36.",
    }


@pytest.fixture
def expansionhunter_json():
    """:return: Returns path to real simplified Repeat Expansion JSON output file - UNCLEAR. """
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.path.join("data", "expansionhunter_normal_result.json"),
    )


@pytest.fixture
def expansionhunter_expanded_expected():
    """:return: Returns dictionary with expected results for ExpansionHunter - EXPANDED."""
    return {
        "AR": (
            "Expanded - repeat AR observed: 38/41. "
            "Normal range: 9 - 35. Expansion range: 38 - 62."
        ),
        "FMR1": (
            "Expanded - repeat FMR1 observed: 201/210. "
            "Normal range: 6 - 54. Expansion range: 200 - 1000."
        ),
        "RFC1": "Undefined - genotype not defined in ExpansionHunter for repeat RFC1.",
        "TBP": "Undefined - no annotation found for repeat TBP. Observed genotype: 35/36.",
    }


@pytest.fixture
def expansionhunter_expanded_json():
    """:return: Returns path to real simplified Repeat Expansion JSON output file - EXPANDED. """
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.path.join("data", "expansionhunter_expanded_result.json"),
    )


@pytest.fixture
def expansionhunter_unclear_expected():
    """:return: Returns dictionary with expected results for ExpansionHunter - UNCLEAR."""
    return {
        "AR": (
            "Unclear - repeat AR observed outside defined ranges (34/37). "
            "Normal range: 9 - 35. Expansion range: 38 - 62."
        ),
        "FMR1": (
            "Unclear - repeat FMR1 observed outside defined ranges (17/60). "
            "Normal range: 6 - 54. Expansion range: 200 - 1000."
        ),
        "RFC1": "Undefined - genotype not defined in ExpansionHunter for repeat RFC1.",
        "TBP": "Undefined - no annotation found for repeat TBP. Observed genotype: 35/36.",
    }


@pytest.fixture
def expansionhunter_unclear_json():
    """:return: Returns path to real simplified Repeat Expansion JSON output file - UNCLEAR. """
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.path.join("data", "expansionhunter_unclear_result.json"),
    )


@pytest.fixture
def annotation_json():
    """:return: Returns path to real Repeat Expansion annotation JSON file."""
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.path.join("data", "repeat_expansion_annotation.json"),
    )


@pytest.fixture
def annotate_obj(expansionhunter_json, annotation_json):
    """:return: Returns AnnotateExpansionHunter object."""
    # Define temporary output file
    temporary_output = os.path.join(tempfile.gettempdir(), "output.json")
    # Return AnnotateExpansionHunter object
    return AnnotateExpansionHunter(
        eh_json=expansionhunter_json, annotation_json=annotation_json, output_path=temporary_output
    )


def test_class_constructor(expansionhunter_json):
    """Tests AnnotateExpansionHunter::__init__()"""
    # Define input
    fake_eh_json = "_not_real_eh_file.json"
    fake_annotation_json = "_not_real_annotation_file.json"
    temporary_output = os.path.join(tempfile.gettempdir(), "output.json")

    # All json files must exist
    with pytest.raises(ValueError):
        AnnotateExpansionHunter(
            eh_json=fake_eh_json,
            annotation_json=fake_annotation_json,
            output_path=temporary_output,
        )
    with pytest.raises(ValueError):
        AnnotateExpansionHunter(
            eh_json=expansionhunter_json,
            annotation_json=fake_annotation_json,
            output_path=temporary_output,
        )


def test_load_json(annotate_obj, expansionhunter_json, annotation_json):
    """Tests AnnotateExpansionHunter::load_json() """
    # Define expected
    expected_annotation_keys = ["AR", "FMR1", "RFC1"]
    expected_output_keys = ["LocusResults", "SampleParameters"]

    # Test output json
    actual = annotate_obj.load_json(path=expansionhunter_json)
    assert isinstance(actual, dict)
    assert len(actual) == len(expected_output_keys)
    assert all([key in expected_output_keys for key in actual])

    # Test annotation json
    actual = annotate_obj.load_json(path=annotation_json)
    assert isinstance(actual, dict)
    assert len(actual) == len(expected_annotation_keys)
    assert all([key in expected_annotation_keys for key in actual])


def test_explain_result(annotate_obj):
    """Tests AnnotateExpansionHunter::explain_result() """
    # Get expected
    path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.path.join("data", "explain_results_tests.json"),
    )
    with open(path) as json_file:
        expected = json.load(json_file)

    # Assert all
    for entry in expected.values():
        expected = entry[6]
        actual = annotate_obj.explain_result(
            repeat=entry[0],
            genotype=entry[1],
            normal_start=entry[2],
            normal_end=entry[3],
            expansion_start=entry[4],
            expansion_end=entry[5],
        )
        assert actual == expected


def test_annotate_normal(
    annotate_obj, expansionhunter_json, annotation_json, expansionhunter_expected
):
    """Tests AnnotateExpansionHunter::annotate() - NORMAL"""
    # Define input
    with open(expansionhunter_json) as json_file:
        output_dict = json.load(json_file)
        simplified_out_dict = output_dict.get("LocusResults")
    with open(annotation_json) as json_file:
        annotation_dict = json.load(json_file)
    # Run and assert
    run_dict = annotate_obj.annotate(simplified_out_dict, annotation_dict)
    for repeat in run_dict:
        assert run_dict.get(repeat) == expansionhunter_expected.get(repeat)


def test_annotate_unclear(
    annotate_obj, expansionhunter_unclear_json, annotation_json, expansionhunter_unclear_expected
):
    """Tests AnnotateExpansionHunter::annotate() - UNCLEAR"""
    # Define input
    with open(expansionhunter_unclear_json) as json_file:
        output_dict = json.load(json_file)
        simplified_out_dict = output_dict.get("LocusResults")
    with open(annotation_json) as json_file:
        annotation_dict = json.load(json_file)
    # Run and assert
    run_dict = annotate_obj.annotate(simplified_out_dict, annotation_dict)
    for repeat in run_dict:
        assert run_dict.get(repeat) == expansionhunter_unclear_expected.get(repeat)


def test_annotate_expanded(
    annotate_obj, expansionhunter_expanded_json, annotation_json, expansionhunter_expanded_expected
):
    """Tests AnnotateExpansionHunter::annotate() - EXPANDED"""
    # Define input
    with open(expansionhunter_expanded_json) as json_file:
        output_dict = json.load(json_file)
        simplified_out_dict = output_dict.get("LocusResults")
    with open(annotation_json) as json_file:
        annotation_dict = json.load(json_file)
    # Run and assert
    run_dict = annotate_obj.annotate(simplified_out_dict, annotation_dict)
    for repeat in run_dict:
        assert run_dict.get(repeat) == expansionhunter_expanded_expected.get(repeat)


def test_run_annotation(annotate_obj):
    """Tests AnnotateExpansionHunter::run() - full run"""
    # Define expected/
    expected_keys = ["annotation", "results", "raw", "sample_parameters"]
    # Temporary file path
    temporary_output = os.path.join(tempfile.gettempdir(), "output.json")
    # Run and assert
    annotate_obj.run()
    with open(temporary_output) as json_file:
        output_dict = json.load(json_file)
    assert len(output_dict) == len(expected_keys)
    assert all([key in expected_keys for key in output_dict])


def test_write_md5(annotate_obj):
    """Tests AnnotateExpansionHunter::write_md5() """

    # Define input
    destiny = os.path.join(tempfile.gettempdir(), "empty.json")
    source = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.path.join("data", "empty.json"),
    )
    copyfile(source, destiny)

    # Define expected
    expected_hash = "d41d8cd98f00b204e9800998ecf8427e"  # empty file
    expected_file = destiny + ".md5"

    # Call and assert
    annotate_obj.write_md5(path=destiny)
    msg_exists = "Method should create file: {0}".format(expected_file)
    assert os.path.exists(expected_file), msg_exists
    msg_hash = "Should return hash for empty file: {0}".format(expected_hash)
    assert open(expected_file).read() == expected_hash, msg_hash
