# -*- coding: utf-8 -*-
"""Code for testing the code in the "abstract" workflow
"""
from copy import deepcopy
import filecmp
from pathlib import Path
from tempfile import NamedTemporaryFile
import textwrap
from typing import TypedDict
from unittest.mock import MagicMock

from biomedsheets.shortcuts import GenericSampleSheet, GermlineCaseSheet
import pytest
import ruamel.yaml
import ruamel.yaml as ruamel_yaml
from snakemake.io import OutputFiles, Wildcards
import yaml

from snappy_pipeline.base import MissingConfiguration, merge_dictlikes
import snappy_pipeline.workflow_model
from snappy_pipeline.workflow_model import ConfigModel
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    DataSearchInfo,
    DataSetInfo,
    LinkInPathGenerator,
    LinkInStep,
    LinkInVcfExternalStepPart,
    LinkOutStepPart,
    WritePedigreeSampleNameStepPart,
    WritePedigreeStepPart,
)

from .conftest import DummyModel, patch_module_fs

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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


# Test LinkInPathGenerator -------------------------------------------------------------------------


def test_link_in_path_generator_data_set_info(
    germline_sheet_fake_fs, config_lookup_paths, work_dir, mocker
):
    """Tests LinkInPathGenerator.run() using ``DataSetInfo`` as input"""
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
        ("/path/P001/FCXXXXXX/L001", "FCXXXXXX/L001", "P001_R1.fastq.gz.md5"),
        ("/path/P001/FCXXXXXX/L001", "FCXXXXXX/L001", "P001_R2.fastq.gz.md5"),
    ]
    assert list(generator.run("P001")) == expected


def test_link_in_path_generator_data_search_info(
    germline_sheet_with_ext_vcf_fake_fs, config_lookup_paths, work_dir, mocker
):
    """Tests LinkInPathGenerator.run() using ``DataSearchInfo`` as input"""
    # Patch out the file system related things in the abstract workflow module
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_with_ext_vcf_fake_fs, mocker
    )
    # Define input
    info = DataSearchInfo(
        sheet_path="sheet.tsv",
        base_paths=config_lookup_paths,
        search_paths=["/vcf_path"],
        search_patterns=[{"vcf": "*_dragen.vcf.gz"}],
        mixed_se_pe=True,
    )
    # Define expected
    root_path = "/vcf_path/220911_A00000_0000_BH7MHCDMXY/P001-N1-DNA1-WGS1"
    expected = [
        (root_path, ".", "P001_dragen.vcf.gz"),
        (root_path, ".", "P001_dragen.vcf.gz.md5"),
    ]
    # Check results
    actual = list(
        LinkInPathGenerator(work_dir, [info], [], cache_file_name="_cache_file").run(
            "P001-N1-DNA1-WGS1", ("vcf", "vcf_md5")
        )
    )
    assert actual == expected


def test_link_in_path_generator_path_link_in(
    germline_sheet_fake_fs_path_link_in, config_lookup_paths, work_dir, mocker
):
    # Patch out the file system related things in the abstract workflow module
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs_path_link_in, mocker
    )
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
    generator = LinkInPathGenerator(
        work_dir, [info], [], cache_file_name="_cache_file", preprocessed_path="/preprocess"
    )
    # Check results
    expected = [
        (
            "/preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out",
            "FCXXXXXX/L001/out",
            "P001_R1.fastq.gz",
        ),
        (
            "/preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out",
            "FCXXXXXX/L001/out",
            "P001_R2.fastq.gz",
        ),
        (
            "/preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out",
            "FCXXXXXX/L001/out",
            "P001_R1.fastq.gz.md5",
        ),
        (
            "/preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out",
            "FCXXXXXX/L001/out",
            "P001_R2.fastq.gz.md5",
        ),
    ]
    assert list(generator.run("P001-N1-DNA1-WGS1")) == expected


# Test LinkInStep ------------------------------------------------------------------------


@pytest.fixture
def dummy_config():
    """Return dummy configuration OrderedDicts"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        step_config: {}
        static_data_config:
          reference:
            path: /path/to/reference.fasta
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
                    path_link_in: ""
                    key: value
                """
            ).lstrip()

    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)

    class DummyStepConfig(TypedDict, total=False):
        dummy: DummyModel

    mocker.patch("snappy_pipeline.workflow_model.StepConfig", DummyStepConfig)

    config = deepcopy(dummy_config)
    yaml = ruamel_yaml.YAML()
    local_config = yaml.load(DummyBaseStep.default_config_yaml())
    dummy_model = DummyModel(**local_config["step_config"]["dummy"])
    dummy_config = merge_dictlikes(config, {"step_config": {"dummy": dummy_model}})

    return DummyBaseStep(
        dummy_workflow,
        dummy_config,
        config_lookup_paths,
        config_paths,
        work_dir,
        config_model_class=DummyModel,
    )


@pytest.fixture
def dummy_generic_step_path_link_in(
    dummy_workflow,
    dummy_config,
    config_lookup_paths,
    config_paths,
    work_dir,
    germline_sheet_fake_fs_path_link_in,
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
                    path_link_in: "/preprocess"
                    key: value
                """
            ).lstrip()

    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs_path_link_in, mocker
    )

    class DummyStepConfig(TypedDict, total=False):
        dummy: DummyModel

    mocker.patch("snappy_pipeline.workflow_model.StepConfig", DummyStepConfig)

    config = deepcopy(dummy_config)
    yaml = ruamel_yaml.YAML()
    local_config = yaml.load(DummyBaseStep.default_config_yaml())
    dummy_model = DummyModel(**local_config["step_config"]["dummy"])
    dummy_config = merge_dictlikes(config, {"step_config": {"dummy": dummy_model}})

    return DummyBaseStep(
        dummy_workflow,
        dummy_config,
        config_lookup_paths,
        config_paths,
        work_dir,
        config_model_class=DummyModel,
    )


def test_link_in_step_part_get_input_files(dummy_generic_step):
    assert dummy_generic_step.get_input_files("link_in", "run") == []


def test_link_in_step_part_get_output_files(dummy_generic_step):
    expected = "work/input_links/{library_name}/.done"
    actual = dummy_generic_step.get_output_files("link_in", "run")
    assert actual == expected


def test_link_in_step_part_get_shell_cmd(germline_sheet_fake_fs, mocker, dummy_generic_step):
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
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R1.fastq.gz.md5 || ln -sr /path/P001/FCXXXXXX/L001/P001_R1.fastq.gz.md5 work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001 && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/P001_R2.fastq.gz.md5 || ln -sr /path/P001/FCXXXXXX/L001/P001_R2.fastq.gz.md5 work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001; }}
        """
    ).strip()
    assert actual == expected


def test_link_in_step_part_get_shell_cmd_path_link_in(
    germline_sheet_fake_fs_path_link_in,
    config_lookup_paths,
    work_dir,
    mocker,
    dummy_generic_step_path_link_in,
):
    # Patch out file system--related stuff in abstract workflow
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_fake_fs_path_link_in, mocker
    )
    # Exercise code under test
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    actual = dummy_generic_step_path_link_in.get_shell_cmd("link_in", "run", wildcards)
    # Check results
    expected = textwrap.dedent(
        r"""
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R1.fastq.gz || ln -sr /preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R1.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R2.fastq.gz || ln -sr /preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R2.fastq.gz work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R1.fastq.gz.md5 || ln -sr /preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R1.fastq.gz.md5 work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out && {{ test -h work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R2.fastq.gz.md5 || ln -sr /preprocess/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out/P001_R2.fastq.gz.md5 work/input_links/P001-N1-DNA1-WGS1/FCXXXXXX/L001/out; }}
        """
    ).strip()
    assert actual == expected


def test_link_in_step_part_get_shell_cmd_double_link_in_regression(
    germline_sheet_fake_fs2, mocker, capsys, dummy_generic_step
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


# Test LinkInExternalStepPart ----------------------------------------------------------------------


@pytest.fixture
def vcf_dummy_config():
    """Return dummy configuration OrderedDicts"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        step_config: {}
        static_data_config:
          reference:
            path: /path/to/reference.fasta
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
def vcf_dummy_generic_step(
    dummy_workflow,
    vcf_dummy_config,
    config_lookup_paths,
    config_paths,
    work_dir,
    germline_sheet_with_ext_vcf_fake_fs,
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
            self.register_sub_step_classes((LinkInVcfExternalStepPart, LinkOutStepPart))
            self.data_search_infos = list(self._load_data_search_infos())

        @classmethod
        def default_config_yaml(cls):
            """Return default config YAML"""
            return textwrap.dedent(
                r"""
                step_config:
                  dummy:
                    key: value
                    search_paths: ["/vcf_path"]  # Path to all VCF files.
                    search_patterns: [{"vcf": "*_dragen.vcf.gz"}] # List of search pattern.
                """
            ).lstrip()

    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_with_ext_vcf_fake_fs, mocker
    )

    class DummyStepConfig(TypedDict, total=False):
        dummy: DummyModel

    mocker.patch("snappy_pipeline.workflow_model.StepConfig", DummyStepConfig)

    config = deepcopy(vcf_dummy_config)
    yaml = ruamel_yaml.YAML()
    local_config = yaml.load(DummyBaseStep.default_config_yaml())
    dummy_model = DummyModel(**local_config["step_config"]["dummy"])
    dummy_config = merge_dictlikes(config, {"step_config": {"dummy": dummy_model}})

    return DummyBaseStep(
        dummy_workflow,
        dummy_config,
        config_lookup_paths,
        config_paths,
        work_dir,
        config_model_class=DummyModel,
    )


def test_link_in_external_step_part_get_input_files(vcf_dummy_generic_step):
    """Tests LinkInExternalStepPart.get_input_files()"""
    actual = vcf_dummy_generic_step.get_input_files("link_in_vcf_external", "run")
    assert len(actual) == 0


def test_link_in_external_step_part_get_output_files(vcf_dummy_generic_step):
    expected = "work/input_links/{library_name}/.done"
    actual = vcf_dummy_generic_step.get_output_files("link_in_vcf_external", "run")
    assert actual == expected


def test_link_in_external_step_part_get_shell_cmd(
    germline_sheet_with_ext_vcf_fake_fs,
    mocker,
    vcf_dummy_generic_step,
):
    # Patch out file system--related stuff in abstract workflow
    patch_module_fs(
        "snappy_pipeline.workflows.abstract", germline_sheet_with_ext_vcf_fake_fs, mocker
    )
    # Define input
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected = textwrap.dedent(
        r"""
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/. && {{ test -h work/input_links/P001-N1-DNA1-WGS1/./P001_dragen.vcf.gz || ln -sr /vcf_path/220911_A00000_0000_BH7MHCDMXY/P001-N1-DNA1-WGS1/P001_dragen.vcf.gz work/input_links/P001-N1-DNA1-WGS1/.; }}
        mkdir -p work/input_links/P001-N1-DNA1-WGS1/. && {{ test -h work/input_links/P001-N1-DNA1-WGS1/./P001_dragen.vcf.gz.md5 || ln -sr /vcf_path/220911_A00000_0000_BH7MHCDMXY/P001-N1-DNA1-WGS1/P001_dragen.vcf.gz.md5 work/input_links/P001-N1-DNA1-WGS1/.; }}
        """
    ).strip()
    # Check results
    actual = vcf_dummy_generic_step.get_shell_cmd("link_in_vcf_external", "run", wildcards)
    assert actual == expected


# Tests for LinkOutStepPart ------------------------------------------------------------------------


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


# Tests for WritePedigreeStepPart* -----------------------------------------------------------------


@pytest.fixture
def dummy_generic_step_w_write_pedigree(
    dummy_workflow,
    dummy_config,
    config_lookup_paths,
    config_paths,
    work_dir,
    germline_sheet_fake_fs,
    mocker,
):
    """Return BaseStep sub class instance using generic sample sheets; for use in tests of
    classes ``WritePedigreeStepPart`` and ``WritePedigreeSampleNameStepPart``.
    """

    class DummyBaseStep(BaseStep):
        """Dummy BaseStep sub class; for use in tests"""

        name = "dummy"
        sheet_shortcut_class = GermlineCaseSheet

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.register_sub_step_classes((WritePedigreeStepPart, WritePedigreeSampleNameStepPart))

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

    class DummyStepConfig(TypedDict, total=False):
        dummy: DummyModel

    mocker.patch("snappy_pipeline.workflow_model.StepConfig", DummyStepConfig)

    config = deepcopy(dummy_config)
    yaml = ruamel_yaml.YAML()
    local_config = yaml.load(DummyBaseStep.default_config_yaml())
    dummy_model = DummyModel(**local_config["step_config"]["dummy"])
    dummy_config = merge_dictlikes(config, {"step_config": {"dummy": dummy_model}})

    return DummyBaseStep(
        dummy_workflow,
        dummy_config,
        config_lookup_paths,
        config_paths,
        work_dir,
        config_model_class=DummyModel,
    )


def test_write_pedigree_step_part_get_input_files(dummy_generic_step_w_write_pedigree):
    """Tests WritePedigreeStepPart.get_input_files()"""
    wildcards = Wildcards(fromdict={"library_name": "P001-N1-DNA1-WGS1"})
    expected = []  # as 'ngs_mapping' is not in self.parent.sub_workflows
    actual = dummy_generic_step_w_write_pedigree.get_input_files("write_pedigree", "run")(wildcards)
    assert actual == expected


def test_write_pedigree_step_part_get_output_files(dummy_generic_step_w_write_pedigree):
    """Tests WritePedigreeStepPart.get_output_files()"""
    expected = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
    actual = dummy_generic_step_w_write_pedigree.get_output_files("write_pedigree", "run")
    assert actual == expected


def test_write_pedigree_step_part_run_for_library(dummy_generic_step_w_write_pedigree):
    """Tests WritePedigreeStepPart.run() - for library 'P001-N1-DNA1-WGS1'"""
    # Define input
    observed_file = NamedTemporaryFile()
    output_ = OutputFiles([observed_file.name])
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected_file = (Path(__file__).parent / "data/write_pedigree.P001-N1-DNA1-WGS1.ped").resolve()
    # Call and compare files
    dummy_generic_step_w_write_pedigree.sub_steps["write_pedigree"].run(wildcards, output_)
    assert filecmp.cmp(observed_file.name, expected_file)


def test_write_pedigree_step_part_run_for_whole_cohort(dummy_generic_step_w_write_pedigree):
    """Tests WritePedigreeStepPart.run() - for whole cohort"""
    # Define input
    observed_file = NamedTemporaryFile()
    output_ = OutputFiles([observed_file.name])
    wildcards = Wildcards(fromdict={"index_ngs_library": "whole_cohort"})
    # Define expected
    expected_file = (Path(__file__).parent / "data/write_pedigree.whole_cohort.ped").resolve()
    # Call and compare files
    dummy_generic_step_w_write_pedigree.sub_steps["write_pedigree"].run(wildcards, output_)
    assert filecmp.cmp(observed_file.name, expected_file)


def test_write_pedigree_sample_name_step_part_run_for_library(dummy_generic_step_w_write_pedigree):
    """Tests WritePedigreeStepPart.run() - for library 'P001-N1-DNA1-WGS1'"""
    # Define input
    observed_file = NamedTemporaryFile(delete=False)
    output_ = OutputFiles([observed_file.name])
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    # Define expected
    expected_file = (Path(__file__).parent / "data/write_pedigree.P001.ped").resolve()
    # Call and compare files
    dummy_generic_step_w_write_pedigree.sub_steps["write_pedigree_with_sample_name"].run(
        wildcards, output_
    )
    assert filecmp.cmp(observed_file.name, expected_file)


# Tests for BaseStep ------------------------------------------------------------------------------


def test_base_step_ensure_w_config(dummy_generic_step):
    dummy_generic_step.ensure_w_config(("step_config", "dummy", "key"), "should be OK")
    with pytest.raises(MissingConfiguration):
        dummy_generic_step.ensure_w_config(("step_config", "dummy", "foo"), "should fail")
