# -*- coding: utf-8 -*-
"""Code for testing the code in the "abstract" workflow
"""

import textwrap

from biomedsheets.shortcuts import GenericSampleSheet
import pytest
import ruamel.yaml as yaml
from snakemake.io import Wildcards

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    DataSetInfo,
    LinkInPathGenerator,
    LinkInStep,
    LinkOutStepPart,
)

from .conftest import patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# Tests for DataSetInfo ---------------------------------------------------------------------------


def test_data_set_info_load_germline_tsv(germline_sheet_fake_fs, config_lookup_paths, mocker):
    # Patch out the file system related things in the abstract workflow module
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Exercise the code under test
    info = DataSetInfo(
        "first_batch",
        "sheet.tsv",
        config_lookup_paths,
        ["/path"],
        [{"left": "*_R1.fastq.gz", "right": "*_R2.fastq.gz"}],
        "germline_variants",
        False,
        "secondary_id_pk",
        False,
        None,
        None,
        None,
    )
    # Check results
    assert info.name == "first_batch"
    assert info.sheet_path == "sheet.tsv"
    assert info.base_paths == config_lookup_paths
    assert info.search_paths == ["/path"]
    assert info.search_patterns == [{"left": "*_R1.fastq.gz", "right": "*_R2.fastq.gz"}]
    assert info.sheet_type == "germline_variants"
    assert not info.is_background
    assert info.naming_scheme == "secondary_id_pk"
    # -- in particular, check sheet
    assert info.sheet
    actual = sorted(list(info.sheet.bio_entities.keys()))
    expected = ["P001", "P002", "P003", "P004", "P005", "P006"]
    assert actual == expected


def test_data_set_info_load_cancer_tsv(cancer_sheet_fake_fs, config_lookup_paths, mocker):
    # Patch out the file system related things in the abstract workflow module
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    # Exercise the code under test
    info = DataSetInfo(
        "first_batch",
        "sheet.tsv",
        config_lookup_paths,
        ["/path"],
        [{"left": "*_R1.fastq.gz", "right": "*_R2.fastq.gz"}],
        "matched_cancer",
        False,
        "secondary_id_pk",
        False,
        None,
        None,
    )
    # Check results
    assert info.name == "first_batch"
    assert info.sheet_path == "sheet.tsv"
    assert info.base_paths == config_lookup_paths
    assert info.search_paths == ["/path"]
    assert info.search_patterns == [{"left": "*_R1.fastq.gz", "right": "*_R2.fastq.gz"}]
    assert info.sheet_type == "matched_cancer"
    assert not info.is_background
    assert info.naming_scheme == "secondary_id_pk"
    # -- in particular, check sheet
    assert info.sheet
    actual = sorted(list(info.sheet.bio_entities.keys()))
    expected = ["P001", "P002"]
    assert actual == expected


# Test LinkInPathGenerator ------------------------------------------------------------------------


def test_link_in_path_generator(germline_sheet_fake_fs, config_lookup_paths, work_dir, mocker):
    # Patch out the file system related things in the abstract workflow module
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Exercise the code under test
    info = DataSetInfo(
        "first_batch",
        "sheet.tsv",
        config_lookup_paths,
        ["/path"],
        [{"left": "*/*/*_R1.fastq.gz", "right": "*/*/*_R2.fastq.gz"}],
        "germline_variants",
        False,
        "secondary_id_pk",
        False,
        None,
        None,
    )
    generator = LinkInPathGenerator(work_dir, [info], [], cache_file_name="_cache_file")
    # Check results
    expected = [
        ("/path/P001/FCXXXXXX/L001", "FCXXXXXX/L001", "P001_R1.fastq.gz"),
        ("/path/P001/FCXXXXXX/L001", "FCXXXXXX/L001", "P001_R2.fastq.gz"),
    ]
    assert list(generator.run("P001")) == expected


# Test LinkInStep ------------------------------------------------------------------------


@pytest.fixture
def dummy_config():
    """Return dummy configuration OrderedDicts"""
    return yaml.round_trip_load(
        textwrap.dedent(
            r"""
        step_config: {}
        static_data_config: {}
        data_sets:
          first_batch:  # example for a matched cancer data set
            file: sheet.tsv
            search_patterns:
            # Note that currently only "left" and "right" key known
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: germline_variants
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def dummy_generic_step(
    dummy_workflow,
    dummy_config,
    dummy_cluster_config,
    config_lookup_paths,
    config_paths,
    work_dir,
    germline_sheet_fake_fs,
    mocker,
):
    """Return BaseStep sub class instance using generic sample sheets; for use in tests of the
    abstract workflow module
    """

    class DummyBaseStep(BaseStep):
        """Dummy BaseStep sub class; for use in tests"""

        name = "dummy"
        sheet_shortcut_class = GenericSampleSheet

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.register_sub_step_classes((LinkInStep, LinkOutStepPart))

        @classmethod
        def default_config_yaml(cls):
            """Return default config YAML"""
            return textwrap.dedent(
                r"""
                step_config:
                  dummy:
                    key: value
                """
            ).lstrip()

    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    return DummyBaseStep(
        dummy_workflow,
        dummy_config,
        dummy_cluster_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_link_in_step_part_get_input_files(dummy_generic_step):
    assert dummy_generic_step.get_input_files("link_in", "run") == []


def test_link_in_step_part_get_output_files(dummy_generic_step):
    expected = "work/input_links/{library_name}/.done"
    actual = dummy_generic_step.get_output_files("link_in", "run")
    assert actual == expected


def test_link_in_step_part_get_shell_cmd(
    germline_sheet_fake_fs, config_lookup_paths, work_dir, mocker, dummy_generic_step
):
    # Patch out file system--related stuff in abstract workflow
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Exercise code under test
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    actual = dummy_generic_step.get_shell_cmd("link_in", "run", wildcards)
    # Check results
    expected = textwrap.dedent(
        r"""
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz || ln -sr /path/P001/FCXXXXXX/L001/P001_R1.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz || ln -sr /path/P001/FCXXXXXX/L001/P001_R2.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001; }}
        """
    ).strip()
    assert actual == expected


def test_link_in_step_part_get_shell_cmd_double_link_in_regression(
    germline_sheet_fake_fs2, config_lookup_paths, work_dir, mocker, capsys, dummy_generic_step
):
    # Patch out file system--related stuff in abstract workflow
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs2, mocker)
    # Exercise code under test
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    actual = dummy_generic_step.get_shell_cmd("link_in", "run", wildcards)
    # Check stdout/stderr
    out, err = capsys.readouterr()
    assert out == ""
    assert err == ""
    # Check results
    expected = textwrap.dedent(
        r"""
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz || ln -sr /path/P001/FCXXXXXX/L001/P001_R1.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz || ln -sr /path/P001/FCXXXXXX/L001/P001_R2.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCYYYYYY/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCYYYYYY/L001/P001_R1.fastq.gz || ln -sr /path/P001/FCYYYYYY/L001/P001_R1.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCYYYYYY/L001; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCYYYYYY/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCYYYYYY/L001/P001_R2.fastq.gz || ln -sr /path/P001/FCYYYYYY/L001/P001_R2.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCYYYYYY/L001; }}
        """
    ).strip()
    assert actual == expected


# Tests for LinkOutStepPart -----------------------------------------------------------------------


def test_link_out_step_part_get_input_files(dummy_generic_step):
    func = dummy_generic_step.get_input_files("link_out", "run")
    assert callable(func)
    expected = "work/path/file.txt"
    wildcards = Wildcards(fromdict={"path": "path", "file": "file", "ext": "txt"})
    actual = func(wildcards)
    assert actual == expected


def test_link_out_step_part_get_output_files(dummy_generic_step):
    expected = "output/{path}/{file}.{ext}"
    actual = dummy_generic_step.get_output_files("link_out", "run")
    assert actual == expected


def test_link_out_step_part_get_shell_cmd(dummy_generic_step):
    # Define expected
    expected = (
        "test -h output/{wildcards.path}/{wildcards.file}.{wildcards.ext} || "
        "ln -sr work/{wildcards.path}/{wildcards.file}.{wildcards.ext} "
        "output/{wildcards.path}/{wildcards.file}.{wildcards.ext}"
    )
    # Get actual
    wildcards = Wildcards(fromdict={"path": "path", "file": "file", "ext": "txt"})
    actual = dummy_generic_step.get_shell_cmd("link_out", "run", wildcards)

    assert actual == expected


# Tests for BaseStep ------------------------------------------------------------------------------


def test_base_step_ensure_w_config(dummy_generic_step):
    dummy_generic_step.ensure_w_config(("step_config", "dummy", "key"), "should be OK")
    with pytest.raises(MissingConfiguration):
        dummy_generic_step.ensure_w_config(("step_config", "dummy", "foo"), "should fail")
