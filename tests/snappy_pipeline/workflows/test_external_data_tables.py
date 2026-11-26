# -*- coding: utf-8 -*-
"""Tests for the files tables of ``external_data`` tasks (``library_files``)."""

import textwrap

import pytest

from snappy_pipeline.base import ConfigError
from snappy_pipeline.orchestration import create_task_instances, load_project, missing_files
from snappy_pipeline.sheet_summary import missing_library_files
from snappy_pipeline.workflows.external_data import ExternalDataWorkflow

SHEET = """\
[Data]
patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tfolderName\thpoTerms
P001\t.\t.\tF\tY\tWGS\tP001\t.
P002\t.\t.\tM\tN\tWGS\tP002\t.
"""

CONFIG = """\
project:
  name: tables
data_sets:
  first:
    file: sheet.tsv
    type: germline_variants
    naming_scheme: only_secondary_id
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: refs/genome.fa}
  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: refs/genome}
  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]
  - name: mapping
    step: ngs_mapping
    config:
      depends_on: {reads: %(reads)s, reference: genome, index: bwa_index}
      tool: bwa
      bwa: {}
"""


def _project(tmp_path, reads, table=None):
    (tmp_path / "config.yaml").write_text(CONFIG % {"reads": reads})
    (tmp_path / "sheet.tsv").write_text(SHEET)
    if table is not None:
        (tmp_path / "reads.tsv").write_text(textwrap.dedent(table))
    return load_project(
        __import__("yaml").safe_load((tmp_path / "config.yaml").read_text()), str(tmp_path)
    )


def _touch(root, *paths):
    for path in paths:
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        (root / path).touch()


def test_the_read_groups_of_a_table_are_ordered_by_the_left_path(tmp_path):
    table = """\
        library\treadgroup\tleft\tright
        P001-N1-DNA1-WGS1\tL2\traw/P001/L2_R1.fastq.gz\traw/P001/L2_R2.fastq.gz
        P002-N1-DNA1-WGS1\tx\traw/P002/x_R1.fastq.gz\traw/P002/x_R2.fastq.gz
        P001-N1-DNA1-WGS1\tL1\traw/P001/L1_R1.fastq.gz\traw/P001/L1_R2.fastq.gz
        """
    mapping = create_task_instances(None, _project(tmp_path, "reads", table=table))["mapping"]

    groups = mapping.read_groups("P001-N1-DNA1-WGS1")
    assert [group.name for group in groups] == ["L1", "L2"]
    assert groups[0].left == f"{tmp_path}/raw/P001/L1_R1.fastq.gz"
    assert mapping.reads_input("P002-N1-DNA1-WGS1") == [
        f"{tmp_path}/raw/P002/x_R1.fastq.gz",
        f"{tmp_path}/raw/P002/x_R2.fastq.gz",
    ]
    assert missing_library_files({"mapping": mapping}) == []


def test_a_library_without_rows_is_reported_with_the_way_to_add_it(tmp_path):
    table = (
        "library\treadgroup\tleft\tright\n"
        "P001-N1-DNA1-WGS1\tL1\traw/P001/L1_R1.fastq.gz\traw/P001/L1_R2.fastq.gz\n"
    )
    mapping = create_task_instances(None, _project(tmp_path, "reads", table=table))["mapping"]
    with pytest.raises(
        ConfigError, match="list no files of library 'P002-N1-DNA1-WGS1'; run `snappy scan-files`"
    ):
        mapping.read_groups("P002-N1-DNA1-WGS1")
    assert missing_library_files({"mapping": mapping}) == [
        (
            "mapping",
            "P002-N1-DNA1-WGS1",
            f"The files tables {tmp_path}/reads.tsv list no files of library 'P002-N1-DNA1-WGS1'; run "
            "`snappy scan-files` to add them",
        )
    ]


def test_a_missing_table_loads_and_is_reported_at_setup(tmp_path):
    project = _project(tmp_path, "reads")
    (tmp_path / "raw").mkdir()
    mapping = create_task_instances(None, project)["mapping"]
    with pytest.raises(ConfigError, match="Files table .*reads.tsv: No such file"):
        mapping.read_groups("P001-N1-DNA1-WGS1")
    assert missing_files(project) == {"reads": [f"{tmp_path}/reads.tsv"]} | {
        "genome": [
            f"{tmp_path}/refs/genome.fa",
            f"{tmp_path}/refs/genome.fa.fai",
            f"{tmp_path}/refs/genome.dict",
        ],
        "bwa_index": [f"{tmp_path}/refs/genome"],
    }


def test_rows_of_unknown_libraries_and_missing_files_are_reported(tmp_path):
    table = """\
        library\treadgroup\tleft\tright
        P001-N1-DNA1-WGS1\tL1\traw/P001/L1_R1.fastq.gz\traw/P001/L1_R2.fastq.gz
        P002-N1-DNA1-WGS1\tL1\traw/P002/L1_R1.fastq.gz\traw/P002/L1_R2.fastq.gz
        P009\tL1\traw/P009/L1_R1.fastq.gz\traw/P009/L1_R2.fastq.gz
        """
    project = _project(tmp_path, "reads", table=table)
    _touch(tmp_path, "raw/P001/L1_R1.fastq.gz", "raw/P001/L1_R2.fastq.gz")
    mapping = create_task_instances(None, project)["mapping"]
    assert missing_library_files({"mapping": mapping}) == [
        ("reads", "P009", "is in a files table but in no sample sheet")
    ]
    assert missing_files(project)["reads"] == [
        f"{tmp_path}/raw/P002/L1_R1.fastq.gz",
        f"{tmp_path}/raw/P002/L1_R2.fastq.gz",
        f"{tmp_path}/raw/P009/L1_R1.fastq.gz",
        f"{tmp_path}/raw/P009/L1_R2.fastq.gz",
    ]


def test_files_per_library_from_a_table(tmp_path):
    (tmp_path / "bams.tsv").write_text("library\tbam\tbai\nP001\tP001.bam\tP001.bam.bai\n")
    config = ExternalDataWorkflow.config_model_class.model_validate(
        {
            "produces": {"type": "alignments", "tags": ["dna"]},
            "library_files": [str(tmp_path / "bams.tsv")],
        }
    )
    project = load_project({"project": {"name": "x"}, "data_sets": {}, "tasks": []}, str(tmp_path))
    assert ExternalDataWorkflow.output_keys(config) == {"bam", "bai"}
    assert ExternalDataWorkflow.get_output_paths(config, library_name="P001", project=project) == {
        "bam": f"{tmp_path}/P001.bam",
        "bai": f"{tmp_path}/P001.bam.bai",
    }
    with pytest.raises(ConfigError, match="list no files of library 'P002'; add them to the table"):
        ExternalDataWorkflow.get_output_paths(config, library_name="P002", project=project)
    config = ExternalDataWorkflow.config_model_class.model_validate(
        {
            "produces": {"type": "alignments", "tags": ["dna"]},
            "library_files": [str(tmp_path / "none.tsv")],
        }
    )
    assert ExternalDataWorkflow.output_keys(config) is None


@pytest.mark.parametrize(
    "config, error",
    [
        ({"files": {"a": "x"}, "library_files": ["t.tsv"]}, "Set one of files"),
        (
            {"library_files": ["t.tsv"], "produces": {"type": "variants"}, "mixed_se_pe": True},
            "only for reads",
        ),
    ],
)
def test_model_errors(config, error):
    with pytest.raises(ValueError, match=error):
        ExternalDataWorkflow.config_model_class.model_validate(
            {"produces": {"type": "raw"}} | config
        )
