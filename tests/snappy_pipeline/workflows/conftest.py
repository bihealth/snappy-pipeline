# -*- coding: utf-8 -*-
"""Shared fixtures for the workflows unit tests"""

from collections import namedtuple
import random
import textwrap
from unittest.mock import MagicMock

from biomedsheets.shortcuts import GenericSampleSheet
from pyfakefs import fake_filesystem
import pytest
from ruamel.yaml.comments import CommentedMap

from snappy_pipeline.workflows.abstract import BaseStep


def random_gender():
    """Random gender.
    :return: Returns either 'F' or 'M' by change.
    """
    genders = ["F", "M"]
    # pseudo-random generator is unsafe for crypto use, but OK here
    return random.choice(genders)  # nosec


def random_affected_status():
    """Random affected status.
    :return: Returns either 'Y' or 'N' by change.
    """
    affected = ["Y", "N"]
    # pseudo-random generator is unsafe for crypto use, but OK here
    return random.choice(affected)  # nosec


def set_trio_ids(index_i):
    """Define trio identifiers.
    :param index_i: Index identifier. Example if '1', final identifier will be 'P001'.
    :type index_i: int
    :return: Returns tuple with trio identifiers (index, mother, father).
    Example: ('P001', 'P002', 'P003').
    """
    # Initialise variable
    sampleid_pattern = "P{i}"
    # Set sample ids for trio
    index_sampleid = sampleid_pattern.format(i=str(index_i).zfill(3))
    mother_i = index_i + 1
    mother_sampleid = sampleid_pattern.format(i=str(mother_i).zfill(3))
    father_i = index_i + 2
    father_sampleid = sampleid_pattern.format(i=str(father_i).zfill(3))
    # Return
    return index_sampleid, mother_sampleid, father_sampleid


def get_entry_for_sample_sheet(index_i, entry_type="trio"):
    """Sample sheet entry text.
    :param index_i: Index id. Example if '1', final id will be 'P001'
    :type index_i: int
    :param entry_type: Type of entry: 'solo', 'duo', 'trio'. Default: 'trio'.
    :type entry_type: str, optional
    :return: Returns text used to define a trio case in sample sheet.
    """
    # Initialise variables
    entry_pattern = (
        "{sampleid}\t{father_sampleid}\t{mother_sampleid}\t{gender}\t{affected}"
        "\tWGS\tAgilent SureSelect Human All Exon V6\t{sampleid}\t.\n"
    )

    # Validate entry option
    valid_options = ["solo", "duo", "trio"]
    if entry_type not in valid_options:
        valid_options_str = ", ".join(valid_options)
        err_msg = "Options '{in_}' is not valid. Valid entry types options: {valid}".format(
            in_=entry_type, valid=valid_options_str
        )
        raise ValueError(err_msg)

    # Set sample ids for trio
    index_sampleid, mother_sampleid, father_sampleid = set_trio_ids(index_i=index_i)

    # Set entries for trio
    index_gender = random_gender()
    index_entry = entry_pattern.format(
        sampleid=index_sampleid,
        father_sampleid=father_sampleid,
        mother_sampleid=mother_sampleid,
        gender=index_gender,
        affected="Y",
    )
    mother_affected = random_affected_status()
    mother_entry = entry_pattern.format(
        sampleid=mother_sampleid,
        father_sampleid=".",
        mother_sampleid=".",
        gender="F",
        affected=mother_affected,
    )
    father_affected = random_affected_status()
    father_entry = entry_pattern.format(
        sampleid=father_sampleid,
        father_sampleid=".",
        mother_sampleid=".",
        gender="M",
        affected=father_affected,
    )

    # Return accordingly
    if entry_type == "solo":
        index_entry = index_entry.replace(mother_sampleid, ".").replace(father_sampleid, ".")
        return index_entry
    elif entry_type == "duo":
        index_entry = index_entry.replace(father_sampleid, ".")
        return index_entry + mother_entry
    elif entry_type == "trio":
        return index_entry + mother_entry + father_entry
    return None


@pytest.fixture
def header_germline_sheet():
    """Returns germline TSV file header with wildcard '{entries}'"""
    return textwrap.dedent(
        """
        [Custom Fields]
        key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
        libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.
        [Data]
        patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tlibraryKit\tfolderName\thpoTerms
        {entries}
        """
    ).lstrip()


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
def text_large_cohort_trios_only():
    """Defines arbitrary large cohort entries for sample sheet - trio cases only"""
    # Initialise variables
    n_patients_cohorts = 501  # 167 indexes
    full_cohort_text = ""
    # Create cohort
    index_i = 1
    while index_i < n_patients_cohorts:
        new_trio = get_entry_for_sample_sheet(index_i=index_i)
        # Append
        full_cohort_text = full_cohort_text + new_trio
        # Update index counter
        index_i += 3
    # Return
    # Remove last newline so biomedsheet doesn't interpret it as an empty entry
    return full_cohort_text.rstrip("\n")


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
def tsv_large_cohort_trios_only_germline_sheet(header_germline_sheet, text_large_cohort_trios_only):
    """Returns contents for large cohort germline TSV file - trio cases only"""
    return header_germline_sheet.format(entries=text_large_cohort_trios_only)


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
def germline_sheet_fake_fs2(
    fake_fs2, germline_sheet_tsv, tsv_large_cohort_trios_only_germline_sheet
):
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
    # Create the sample TSV file - Small cohort
    fake_fs2.fs.create_file(
        "/work/config/sheet.tsv", contents=germline_sheet_tsv, create_missing_dirs=True
    )
    # Create the sample TSV file - Large cohort trio cases only
    fake_fs2.fs.create_file(
        "/work/config/sheet_large_cohort_trio.tsv",
        contents=tsv_large_cohort_trios_only_germline_sheet,
        create_missing_dirs=True,
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
