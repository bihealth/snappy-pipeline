"""Rule names: the task, part and action of ``<task>_<part>_<action>``."""

from types import SimpleNamespace

from snappy_pipeline.rule_names import RuleName, parse_rule_name, task_parts

PARTS = {
    "mapping": {"bwa": ("run",), "qc_summary": ("run",)},
    "mapping_qc": {"alfred": ("run",), "qc_summary": ("run",)},
}


def test_parse_rule_name_resolves_a_task_that_prefixes_another():
    assert parse_rule_name("mapping_bwa_run", PARTS) == RuleName("mapping", "bwa", "run")
    assert parse_rule_name("mapping_qc_summary_run", PARTS) == RuleName(
        "mapping", "qc_summary", "run"
    )
    assert parse_rule_name("mapping_qc_qc_summary_run", PARTS) == RuleName(
        "mapping_qc", "qc_summary", "run"
    )
    assert parse_rule_name("mapping_qc_alfred_run", PARTS) == RuleName(
        "mapping_qc", "alfred", "run"
    )


def test_parse_rule_name_keeps_the_task_of_other_rules():
    # the results job and a handwritten rule of no part: the longest task name
    assert parse_rule_name("mapping_results", PARTS) == RuleName("mapping", "", "")
    assert parse_rule_name("mapping_qc_results", PARTS) == RuleName("mapping_qc", "", "")
    assert parse_rule_name("mapping_bwa_index", PARTS) == RuleName("mapping", "", "")
    assert parse_rule_name("all", PARTS) is None


def test_task_parts_reads_the_step_instances():
    part = lambda name, *actions: SimpleNamespace(name=name, actions=actions)  # noqa: E731
    instances = {"mapping": SimpleNamespace(sub_steps={"bwa": part("bwa", "run", "index")})}
    assert task_parts(instances) == {"mapping": {"bwa": ("run", "index")}}
