# -*- coding: utf-8 -*-
"""Shared fixtures for the workflows unit tests"""

from collections import namedtuple
import textwrap
from unittest.mock import MagicMock

from biomedsheets.shortcuts import GenericSampleSheet
from pyfakefs import fake_filesystem
import pytest
from ruamel.yaml.comments import CommentedMap

from snappy_pipeline.workflows.abstract import BaseStep


@pytest.fixture
def dummy_config():
    """Return dummy configuration OrderedDicts"""
    return CommentedMap([("data_sets", CommentedMap())])


@pytest.fixture
def dummy_cluster_config():
    """Return dummy cluster configuration OrderedDicts"""
    return CommentedMap()


@pytest.fixture
def dummy_workflow():
    """Return dummy Snakemake workflow object"""
    mock_workflow = MagicMock()
    return mock_workflow


@pytest.fixture
def work_dir():
    """Return work directory string for overall consistency in tests"""
    return "/work"


@pytest.fixture
def config_lookup_paths(fake_fs):
    """
    :return: Returns configuration lookup paths list for overall consistency in tests.
    Method also create paths in the fake file system.
    """
    lookup_paths = ["/decoy/config", "/work/config"]
    for l_path in lookup_paths:
        fake_fs.fs.makedirs(l_path, exist_ok=True)
    return lookup_paths


@pytest.fixture
def config_paths():
    """Return configuration paths list for overall consistency in tests"""
    return []


@pytest.fixture
def dummy_generic_step(
    dummy_workflow, dummy_config, dummy_cluster_config, work_dir, config_lookup_paths
):
    """Return BaseStep sub class instance using generic sample sheets; for use in tests"""

    class DummyBaseStep(BaseStep):
        """Dummy BaseStep sub class; for use in tests"""

        name = "dummy"
        sheet_shortcut_class = GenericSampleSheet

        @classmethod
        def default_config_yaml(cls):
            """Return default config YAML"""
            return "step_config:\n  dummy:\n    key: value"

    return DummyBaseStep(
        dummy_workflow, dummy_config, dummy_cluster_config, config_lookup_paths, work_dir
    )


@pytest.fixture
def germline_sheet_tsv():
    """Return contents for germline TSV file"""
    return textwrap.dedent(
        """
        [Custom Fields]
        key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
        libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.

        [Data]
        patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tlibraryKit\tfolderName\thpoTerms
        P001\tP002\tP003\tF\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP001\t.
        P002\t.\t.\tM\tN\tWGS\tAgilent SureSelect Human All Exon V6\tP002\t.
        P003\t.\t.\tF\tN\tWGS\tAgilent SureSelect Human All Exon V6\tP003\t.
        P004\tP005\tP006\tM\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP004\t.
        P005\t.\t.\tM\tN\tWGS\tAgilent SureSelect Human All Exon V6\tP005\t.
        P006\t.\t.\tF\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP006\t.
        """
    ).lstrip()


@pytest.fixture
def germline_trio_plus_sheet_tsv():
    """Return contents for germline trio plus TSV file"""
    return textwrap.dedent(
        """
        [Custom Fields]
        key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
        familyId\tbioEntity\tFamily\tstring\t.\t.\t.\t.\t.
        libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.

        [Data]
        familyId\tpatientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tlibraryKit\tfolderName\thpoTerms
        family1\tP001\tP002\tP003\tF\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP011\t.
        family1\tP002\t.\t.\tM\tN\tWGS\tAgilent SureSelect Human All Exon V6\tP002\t.
        family1\tP003\t.\t.\tF\tN\tWGS\tAgilent SureSelect Human All Exon V6\tP003\t.
        family2\tP004\tP005\tP006\tM\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP004\t.
        family2\tP005\t.\t.\tM\tN\tWGS\tAgilent SureSelect Human All Exon V6\tP005\t.
        family2\tP006\t.\t.\tF\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP006\t.
        family2\tP007\t.\t.\tF\tY\tWGS\tAgilent SureSelect Human All Exon V6\tP007\t.
        """
    ).lstrip()


@pytest.fixture
def cancer_sheet_tsv():
    """Return contents for cancer TSV file"""
    return textwrap.dedent(
        """
        patientName\tsampleName\tisTumor\tlibraryType\tfolderName
        P001\tN1\tN\tWGS\tP001-N1-DNA1-WGS1
        P001\tT1\tY\tWGS\tP001-T1-DNA1-WGS1
        P001\tT1\tY\tmRNA_seq\tP001-T1-RNA1-mRNAseq1
        P002\tN1\tN\tWGS\tP002-N1-DNA1-WGS1
        P002\tT1\tY\tWGS\tP002-T1-DNA1-WGS1
        P002\tT1\tY\tWGS\tP002-T1-DNA1-WGS2
        P002\tT2\tY\tWGS\tP002-T2-DNA1-WGS1
        P002\tT2\tY\tmRNA_seq\tP002-T2-RNA1-mRNAseq1
        """
    ).lstrip()


@pytest.fixture
def generic_rna_sheet_tsv():
    """Return contents for generic RNA TSV file"""
    return textwrap.dedent(
        """
        [Metadata]
        schema\tgeneric
        schema_version\tv1
        title\tExample generic RNA study
        description\tSimple example of a generic study sample sheet.

        [Data]
        bioEntity\tbioSample\ttestSample\tngsLibrary\textractionType\tlibraryType\tfolderName
        E001\tBS1\tTS1\tLIB1\tRNA\ttotal_RNA_seq\tE001-BS1-TS1-LIB1
        E001\tBS2\tTS1\tLIB1\tRNA\ttotal_RNA_seq\tE001-BS2-TS1-LIB1
        E002\tBS1\tTS1\tLIB1\tRNA\ttotal_RNA_seq\tE002-BS1-TS1-LIB1
        E002\tBS1\tTS1\tLIB2\tRNA\ttotal_RNA_seq\tE002-BS1-TS1-LIB2
        """
    ).lstrip()


@pytest.fixture
def generic_mix_extraction_sheet_tsv():
    """Return contents for generic DNA|RNA TSV file"""
    return textwrap.dedent(
        """
        [Metadata]
        schema\tgeneric
        schema_version\tv1
        title\tExample generic mix data study
        description\tSimple example of a generic study sample sheet.

        [Data]
        bioEntity\tbioSample\ttestSample\tngsLibrary\textractionType\tlibraryType\tfolderName
        E001\tBS1\tTS1\tLIB1\tRNA\ttotal_RNA_seq\tE001-BS1-TS1-LIB1
        E001\tBS2\tTS1\tLIB1\tRNA\ttotal_RNA_seq\tE001-BS2-TS1-LIB1
        E002\tBS1\tTS1\tLIB1\tDNA\tWES\tE002-BS1-TS1-LIB1
        E002\tBS1\tTS1\tLIB2\tDNA\tWES\tE002-BS1-TS1-LIB2
        """
    ).lstrip()


@pytest.fixture
def fake_fs():
    """Return ``namedtuple`` with fake file system objects"""
    klass = namedtuple("FakeFsBundle", "fs os open inter_process_lock")
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    fake_lock = MagicMock()
    return klass(fs=fake_fs, os=fake_os, open=fake_open, inter_process_lock=fake_lock)


@pytest.fixture
def fake_fs2():
    """Return ``namedtuple`` with fake file system objects"""
    klass = namedtuple("FakeFsBundle", "fs os open inter_process_lock")
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    fake_lock = MagicMock()
    return klass(fs=fake_fs, os=fake_os, open=fake_open, inter_process_lock=fake_lock)


@pytest.fixture
def sample_cache_dict():
    return {
        "cache_version": 1,
        "root_dirs": {
            "/path": [
                ["P001", "flowcell", "lane", "P001_R1.fastq.gz"],
                ["P001", "flowcell", "lane", "P001_R2.fastq.gz"],
            ]
        },
    }


@pytest.fixture
def germline_sheet_fake_fs(fake_fs, germline_sheet_tsv):
    """Return fake file system setup with files for the germline_sheet_tsv"""
    # Create work directory
    fake_fs.fs.makedirs("/work", exist_ok=True)
    # Create FASTQ read files for the samples
    tpl = "/path/{donor}/FCXXXXXX/L001/{donor}_R{i}.fastq.gz"
    for line in germline_sheet_tsv.splitlines()[1:]:
        donor = line.split("\t")[0]
        fake_fs.fs.create_file(tpl.format(donor=donor, i=1))
        fake_fs.fs.create_file(tpl.format(donor=donor, i=2))
    # Create the sample TSV file
    fake_fs.fs.create_file(
        "/work/config/sheet.tsv", contents=germline_sheet_tsv, create_missing_dirs=True
    )
    return fake_fs


@pytest.fixture
def germline_sheet_fake_fs2(fake_fs2, germline_sheet_tsv):
    """Return fake file system setup with files for the germline_sheet_tsv"""
    # Create work directory
    fake_fs2.fs.makedirs("/work", exist_ok=True)
    # Create FASTQ read files for the samples
    tpl = "/path/{donor}/{flowcell}/L001/{donor}_R{i}.fastq.gz"
    for line in germline_sheet_tsv.splitlines()[1:]:
        for fc in ("FCXXXXXX", "FCYYYYYY"):
            donor = line.split("\t")[0]
            fake_fs2.fs.create_file(tpl.format(donor=donor, flowcell=fc, i=1))
            fake_fs2.fs.create_file(tpl.format(donor=donor, flowcell=fc, i=2))
    # Create the sample TSV file
    fake_fs2.fs.create_file(
        "/work/config/sheet.tsv", contents=germline_sheet_tsv, create_missing_dirs=True
    )
    return fake_fs2


@pytest.fixture
def germline_trio_plus_sheet_fake_fs(fake_fs, germline_trio_plus_sheet_tsv):
    """Return fake file system setup with files for the germline_trio_plus_sheet_tsv"""
    # Create work directory
    fake_fs.fs.makedirs("/work", exist_ok=True)
    # Create the sample TSV file
    fake_fs.fs.create_file(
        "/work/config/sheet_trio_plus.tsv",
        contents=germline_trio_plus_sheet_tsv,
        create_missing_dirs=True,
    )
    return fake_fs


@pytest.fixture
def cancer_sheet_fake_fs(fake_fs, cancer_sheet_tsv):
    """Return fake file system setup with files for the cancer_sheet_tsv"""
    # Create work directory
    fake_fs.fs.makedirs("/work", exist_ok=True)
    # Create FASTQ read files for the samples
    tpl = "/path/{folder}/FCXXXXXX/L001/{folder}_R{i}.fastq.gz"
    for line in cancer_sheet_tsv.splitlines()[1:]:
        folder = line.split("\t")[4]
        fake_fs.fs.create_file(tpl.format(folder=folder, i=1))
        fake_fs.fs.create_file(tpl.format(folder=folder, i=2))
    # Create the sample TSV file
    fake_fs.fs.create_file(
        "/work/config/sheet.tsv", contents=cancer_sheet_tsv, create_missing_dirs=True
    )
    return fake_fs


@pytest.fixture
def aligner_indices_fake_fs(fake_fs):
    """Return fake file system setup with files for aligner indices"""
    d = {
        "bwa": (".fasta.amb", ".fasta.ann", ".fasta.bwt", ".fasta.pac", ".fasta.sa"),
        "star": ("/Genome", "/SA", "/SAindex"),
    }
    for aligner, suffixes in d.items():
        for suffix in suffixes:
            fake_fs.fs.create_file("/path/to/{}/index{}".format(aligner, suffix))
    return fake_fs


def patch_module_fs(module_name, fake_fs, mocker):
    """Helper function to mock out the file-system related things in the module with the given
    name using the given fake_fs and pytest-mock mocker
    """
    mocker.patch("{}.os".format(module_name), fake_fs.os)
    mocker.patch("{}.open".format(module_name), fake_fs.open, create=True)
    mocker.patch("{}.os".format(module_name), fake_fs.os)
    mocker.patch("snappy_pipeline.find_file.InterProcessLock", fake_fs.inter_process_lock)
    mocker.patch("snappy_pipeline.find_file.open", fake_fs.open, create=True)
    mocker.patch("snappy_pipeline.find_file.os", fake_fs.os)
