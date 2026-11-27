"""Rule names: the task, part and action of ``<task>__<part>_<action>``."""

from types import SimpleNamespace

from snappy_pipeline.rule_names import RuleName, parse_rule_name, split_task, task_parts

PARTS = {"mapping": {"bwa": ("run",), "qc_summary": ("run",)}, "mapping_qc": {"alfred": ("run",)}}


def test_split_task_cuts_at_the_separator():
    assert split_task("mapping_qc__qc_summary_run") == ("mapping_qc", "qc_summary_run")
    assert split_task("mapping__results") == ("mapping", "results")
    assert split_task("all") is None


def test_parse_rule_name_finds_the_part_and_action():
    assert parse_rule_name("mapping__bwa_run", PARTS) == RuleName("mapping", "bwa", "run")
    assert parse_rule_name("mapping__qc_summary_run", PARTS) == RuleName(
        "mapping", "qc_summary", "run"
    )
    assert parse_rule_name("mapping_qc__alfred_run", PARTS) == RuleName(
        "mapping_qc", "alfred", "run"
    )


def test_parse_rule_name_keeps_the_task_of_other_rules():
    assert parse_rule_name("mapping__results", PARTS) == RuleName("mapping", "", "")
    assert parse_rule_name("mapping__bwa_index", PARTS) == RuleName("mapping", "", "")
    assert parse_rule_name("all", PARTS) is None


def test_task_parts_reads_the_step_instances():
    part = lambda name, *actions: SimpleNamespace(name=name, actions=actions)  # noqa: E731
    instances = {"mapping": SimpleNamespace(sub_steps={"bwa": part("bwa", "run", "index")})}
    assert task_parts(instances) == {"mapping": {"bwa": ("run", "index")}}
