"""Tests for the step overview (the docs table and the Graphviz graph)."""

import shutil
import subprocess

import pytest

from snappy_pipeline import step_overview
from snappy_pipeline.orchestration import ON_HOLD_STEPS
from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY


@pytest.fixture(scope="module")
def contracts():
    return step_overview.step_contracts()


def test_every_registered_step_has_a_contract(contracts):
    assert {contract.name for contract in contracts} == set(WORKFLOW_REGISTRY) - ON_HOLD_STEPS


def test_flow_requirements_have_producers_and_resources_have_none(contracts):
    calling = next(contract for contract in contracts if contract.name == "variant_calling")
    requirements = {requirement.key: requirement for requirement in calling.requirements}
    assert "ngs_mapping" in step_overview.producers(
        requirements["alignments"], "variant_calling", contracts
    )
    # The genome comes from a resource step only
    assert step_overview.producers(requirements["reference"], "variant_calling", contracts) == []
    # The tool placeholder of an index matches any tool, but indexes are resources
    mapping = next(contract for contract in contracts if contract.name == "ngs_mapping")
    index = next(requirement for requirement in mapping.requirements if requirement.key == "index")
    assert step_overview.producers(index, "ngs_mapping", contracts) == []


def test_table_and_graph_name_every_step(contracts):
    table = step_overview.rst_table(contracts)
    graph = step_overview.dot(contracts)
    for contract in contracts:
        assert f":ref:`{contract.name} <step_{contract.name}>`" in table
        assert f'"{contract.name}" [' in graph
    assert "- ``reference``: ``sequences [dna, reference]``" in table
    assert "- ``dbsnp``: ``variants [dbsnp]`` (optional)" in table


@pytest.mark.skipif(shutil.which("dot") is None, reason="Graphviz is not installed")
def test_graph_renders(contracts):
    svg = subprocess.run(
        ["dot", "-Tsvg"],
        input=step_overview.dot(contracts),
        text=True,
        check=True,
        capture_output=True,
    ).stdout
    assert "<svg" in svg
