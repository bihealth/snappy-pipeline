# -*- coding: utf-8 -*-
"""Code for testing the code in ``Chunk`` class.
"""

from collections import defaultdict
import io
import random
import textwrap

from biomedsheets.io_tsv import read_germline_tsv_sheet
from biomedsheets.shortcuts import GermlineCaseSheet, is_background, is_not_background
import pytest

from snappy_pipeline.chunk import BatchInsufficientSpaceException, Chunk


def build_pedigree_size_dictionary(all_pedigrees):
    """Pedigree to size dictionary.

    :param all_pedigrees: List with all pedigrees in cohort.
    :type all_pedigrees: list

    :return: Returns dictionary with size of each pedigree. Key: index secondary identifier
    (e.g., 'P001' [string]); Value: size of pedigree/family (int).
    """
    # Initialise variables
    pedigree_count_per_index_dict = defaultdict(int)
    # Iterate over donors to populate dict
    for pedigree in all_pedigrees:
        size = len(pedigree.donors)
        index_secondary_id = pedigree.index.wrapped.secondary_id
        pedigree_count_per_index_dict[index_secondary_id] = size
    # Return
    return pedigree_count_per_index_dict


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


def get_entry_for_sample_sheet_with_custom_fields(index_i, batch_number, entry_type="trio"):
    """Sample sheet entry text with custom fields

    :param index_i: Index id. Example if '1', final id will be 'P001'
    :type index_i: int

    :param batch_number: Batch number.
    :type batch_number: int

    :param entry_type: Type of entry: 'solo', 'duo', 'trio'. Default: 'trio'.
    :type entry_type: str, optional

    :return: Returns text used to define a trio case in sample sheet with custom fields 'familyId'
    and 'batchNo'.
    """
    # Initialise variables
    entry_pattern = (
        "{familyid}\t{sampleid}\t{father_sampleid}\t{mother_sampleid}\t{gender}\t{affected}"
        "\t{batch}\tWGS\tAgilent SureSelect Human All Exon V6\t{sampleid}\t.\n"
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

    # Set family id
    familyid = "FAM_" + index_sampleid

    # Set entries for trio
    index_gender = random_gender()
    index_entry = entry_pattern.format(
        familyid=familyid,
        sampleid=index_sampleid,
        father_sampleid=father_sampleid,
        mother_sampleid=mother_sampleid,
        gender=index_gender,
        affected="Y",
        batch=batch_number,
    )
    mother_affected = random_affected_status()
    mother_entry = entry_pattern.format(
        familyid=familyid,
        sampleid=mother_sampleid,
        father_sampleid=".",
        mother_sampleid=".",
        gender="F",
        affected=mother_affected,
        batch=batch_number,
    )
    father_affected = random_affected_status()
    father_entry = entry_pattern.format(
        familyid=familyid,
        sampleid=father_sampleid,
        father_sampleid=".",
        mother_sampleid=".",
        gender="M",
        affected=father_affected,
        batch=batch_number,
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
def large_cohort_trios_only_text():
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
def large_cohort_diverse_text():
    """Defines arbitrary large cohort entries for sample sheet - diverse cases"""
    # Initialise variables
    full_cohort_text = ""
    index_i = 1
    # Create cohort trios - 30 cases, 90 samples
    for _ in range(30):
        new_trio = get_entry_for_sample_sheet(index_i=index_i)
        # Append
        full_cohort_text = full_cohort_text + new_trio
        # Update index counter
        index_i += 3
    # Create cohort duo - 30 cases, 60 samples
    for _ in range(30):
        new_duo = get_entry_for_sample_sheet(index_i=index_i, entry_type="duo")
        # Append
        full_cohort_text = full_cohort_text + new_duo
        # Update index counter
        index_i += 2
    # Create cohort solo - 50 cases, 50 samples
    for _ in range(50):
        new_solo = get_entry_for_sample_sheet(index_i=index_i, entry_type="solo")
        # Append
        full_cohort_text = full_cohort_text + new_solo
        # Update index counter
        index_i += 1
    # Return
    # Remove last newline so biomedsheet doesn't interpret it as an empty entry
    return full_cohort_text.rstrip("\n")


@pytest.fixture
def medium_cohort_solo_only_text():
    """Defines arbitrary medium cohort entries for sample sheet - solo case only"""
    # Initialise variables
    full_cohort_text = ""
    index_i = 900
    # Create cohort solo - 99
    for _ in range(99):
        new_solo = get_entry_for_sample_sheet(index_i=index_i, entry_type="solo")
        # Append
        full_cohort_text = full_cohort_text + new_solo
        # Update index counter
        index_i += 1
    # Return
    # Remove last newline so biomedsheet doesn't interpret it as an empty entry
    return full_cohort_text.rstrip("\n")


@pytest.fixture
def medium_cohort_diverse_with_custom_fields_text():
    """Defines arbitrary medium cohort entries for sample sheet with custom fields -
    diverse cases"""
    # Initialise variables
    full_cohort_text = ""

    # Create index list
    indexes_list = list(range(1, 31, 3))  # 10 trio cases, 30 samples
    indexes_list += list(range(31, 51, 2))  # 10 duo cases, 20 samples
    indexes_list += list(range(51, 101, 1))  # 50 solo cases, 50 samples

    # Create case type list - reflects `indexes_list` structure
    case_list = 10 * ["trio"]
    case_list += 10 * ["duo"]
    case_list += 50 * ["solo"]

    # Create batch list
    batch_list = []
    for i in range(1, 8, 1):  # 7 batches, 10 cases in each
        batch_list.extend(10 * [i])
    batch_list = batch_list[::-1]  # reverse order -> [start:stop:step]

    # Create cohort
    for index, batch, case in zip(indexes_list, batch_list, case_list):
        case = get_entry_for_sample_sheet_with_custom_fields(
            index_i=index, batch_number=batch, entry_type=case
        )
        # Append
        full_cohort_text = full_cohort_text + case

    # Return
    # Remove last newline so biomedsheet doesn't interpret it as an empty entry
    return full_cohort_text.rstrip("\n")


@pytest.fixture
def germline_sheet_header():
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
def germline_sheet_header_with_custom_fields():
    """Returns germline TSV file header with custom fields and with wildcard '{entries}'"""
    return textwrap.dedent(
        """
        [Custom Fields]
        key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
        batchNo\tbioEntity\tBatch No.\tinteger\t.\t.\t.\t.\t.
        familyId\tbioEntity\tFamily\tstring\t.\t.\t.\t.\t.
        libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.

        [Data]
        familyId\tpatientName\tfatherName\tmotherName\tsex\tisAffected\tbatchNo\tlibraryType\tlibraryKit\tfolderName\thpoTerms
        {entries}
        """
    ).lstrip()


@pytest.fixture
def medium_cohort_solo_only_tsv(germline_sheet_header, medium_cohort_solo_only_text):
    """Returns contents for medium cohort germline TSV file - solo cases only"""
    return germline_sheet_header.format(entries=medium_cohort_solo_only_text)


@pytest.fixture
def medium_cohort_diverse_with_custom_features_tsv(
    germline_sheet_header_with_custom_fields, medium_cohort_diverse_with_custom_fields_text
):
    """Returns contents for medium cohort germline with custom fields TSV file - diverse cases"""
    return germline_sheet_header_with_custom_fields.format(
        entries=medium_cohort_diverse_with_custom_fields_text
    )


@pytest.fixture
def large_cohort_trios_only_germline_sheet_tsv(germline_sheet_header, large_cohort_trios_only_text):
    """Returns contents for large cohort germline TSV file - trio cases only"""
    return germline_sheet_header.format(entries=large_cohort_trios_only_text)


@pytest.fixture
def large_cohort_diverse_cases_germline_sheet_tsv(germline_sheet_header, large_cohort_diverse_text):
    """Returns contents for large cohort germline TSV file - diverse cases"""
    return germline_sheet_header.format(entries=large_cohort_diverse_text)


@pytest.fixture
def germline_sample_sheet_object(germline_sheet_tsv):
    """Returns GermlineCaseSheet object with small cohort"""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(germline_sheet_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def germline_sample_sheet_object_large_cohort_trios_only(
    large_cohort_trios_only_germline_sheet_tsv,
):
    """Returns GermlineCaseSheet object with large cohort - trio cases only"""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(large_cohort_trios_only_germline_sheet_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def germline_sample_sheet_object_large_cohort_diverse_cases(
    large_cohort_diverse_cases_germline_sheet_tsv,
):
    """Returns GermlineCaseSheet object with large cohort - trio cases onl."""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(large_cohort_diverse_cases_germline_sheet_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def germline_sample_sheet_object_medium_cohort_solo_cases_background(medium_cohort_solo_only_tsv):
    """Returns GermlineCaseSheet object with medium background cohort - solo cases only"""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(medium_cohort_solo_only_tsv)
    sheet = GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))
    sheet.sheet.json_data["extraInfoDefs"]["is_background"] = {"type": "boolean", "default": False}
    sheet.sheet.extra_infos["is_background"] = True
    return sheet


@pytest.fixture
def germline_sample_sheet_object_medium_cohort_diverse_with_custom_features(
    medium_cohort_diverse_with_custom_features_tsv,
):
    """Returns GermlineCaseSheet object with medium cohort with custom features - diverse cases"""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(medium_cohort_diverse_with_custom_features_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def chunk_object_single_run(germline_sample_sheet_object):
    """Returns Chunk object for single run"""
    # Get germline sheet
    sheet = germline_sample_sheet_object
    # Run single
    return Chunk(method="single", sheet_list=[sheet])


def test_background_sheet_fixture(germline_sample_sheet_object_medium_cohort_solo_cases_background):
    """Tests fixture germline_sample_sheet_object_medium_cohort_solo_cases_background()"""
    sheet = germline_sample_sheet_object_medium_cohort_solo_cases_background
    assert len(list(filter(is_background, [sheet]))) == 1
    assert len(list(filter(is_not_background, [sheet]))) == 0


def test_chunk_constructor(germline_sample_sheet_object):
    """Tests Chunk::__init__()"""
    # Initialise variables
    valid_methods = ["single", "evenly", "incremental"]
    arbitrary_max_size = 100
    silly_sheet_list = [1, 2, 3]
    silly_dict = {1: "one"}
    sheet = germline_sample_sheet_object

    # Test only valid methods
    for method in valid_methods:
        c = Chunk(method=method, sheet_list=silly_sheet_list, maximal_chunk_size=arbitrary_max_size)
        assert isinstance(c, Chunk), "Expected instance of class `Chunk`."

    with pytest.raises(ValueError):
        Chunk(method="_invalid_method_", sheet_list=silly_sheet_list)

    # Test sample sheet list
    with pytest.raises(ValueError):
        Chunk(method="single", sheet_list=silly_dict)  # Expects list

    with pytest.raises(ValueError):
        Chunk(method="single", sheet_list=[])  # Expects list with at least one entry

    # Tests minimum chunk size - depends on the tool being used
    # Expects value greater than zero, cannot be None
    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], minimum_chunk_size=-2)

    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], minimum_chunk_size=0)

    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], minimum_chunk_size="number")

    # Tests maximum chunk size - used in 'evenly' and 'incremental'
    # Expects value greater than zero, cannot be None
    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], maximal_chunk_size=-2)

    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], maximal_chunk_size=0)

    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], maximal_chunk_size="_zero_")

    with pytest.raises(Exception):
        # For 'evenly' max size must be set
        Chunk(method="evenly", sheet_list=[sheet])

    with pytest.raises(ValueError):
        Chunk(method="incremental", sheet_list=[sheet], maximal_chunk_size=-2)

    with pytest.raises(ValueError):
        Chunk(method="incremental", sheet_list=[sheet], maximal_chunk_size=0)

    with pytest.raises(Exception):
        # For 'incremental' max size must be set
        Chunk(method="incremental", sheet_list=[sheet])


def test_order_donors_by_field(
    chunk_object_single_run,
    germline_sample_sheet_object_medium_cohort_diverse_with_custom_features,
):
    """Tests Chunk::order_donors_by_field()"""
    # Define expected
    # Possible point of confusion: the batches were assigned in the reversed order of family ids
    # sampleid: P001; familyid: FAM_P001; batch: 7 - first index of batch number 7
    first_familyid = "FAM_P001"
    first_batch = 7
    # sampleid: P100; familyid: FAM_P100; batch: 1 - last index of batch number 1
    last_familyid = "FAM_P100"
    last_batch = 1

    # Get sheet
    sheet = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features
    all_pedigrees = sheet.cohort.pedigrees

    # Call method - batch number
    reordered_pedigrees = chunk_object_single_run.order_donors_by_field(
        all_pedigrees=all_pedigrees, fields=["batchNo"]
    )
    # Get actual
    actual_first_batch = reordered_pedigrees[0].index.extra_infos.get("batchNo")
    actual_last_batch = reordered_pedigrees[-1].index.extra_infos.get("batchNo")
    assert actual_first_batch == last_batch
    assert actual_last_batch == first_batch
    # Expects increase in batch number
    batch_no = None
    first_iteration = True
    for pedigree in reordered_pedigrees:
        # Initialise batch number on first iteration
        if first_iteration:
            batch_no = pedigree.index.extra_infos.get("batchNo")
            first_iteration = False
            continue
        # Validate
        assert_msg = "Previous batch number should be smaller or equal to current batch number."
        assert batch_no <= pedigree.index.extra_infos.get("batchNo"), assert_msg
        batch_no = pedigree.index.extra_infos.get("batchNo")

    # Call method - family id
    reordered_pedigrees = chunk_object_single_run.order_donors_by_field(
        all_pedigrees=all_pedigrees, fields=["familyId"]
    )
    # Get actual
    actual_first_familyid = reordered_pedigrees[0].index.extra_infos.get("familyId")
    actual_last_familyid = reordered_pedigrees[-1].index.extra_infos.get("familyId")
    assert actual_first_familyid == first_familyid
    assert actual_last_familyid == last_familyid


def test_get_random_cohort_sample(
    chunk_object_single_run,
    germline_sample_sheet_object,
    germline_sample_sheet_object_medium_cohort_solo_cases_background,
):
    """Tests Chunk::get_random_cohort_sample()"""
    # Initialise variables
    expected_donor_type = "<class 'biomedsheets.shortcuts.germline.GermlineDonor'>"
    assertion_type_msg = "Should return '{0}' object.".format(expected_donor_type)

    # ----------------- #
    # Test small cohort #
    # ----------------- #
    s_sheet = germline_sample_sheet_object
    s_all_donors = s_sheet.cohort.pedigrees
    # Expects Value error for request larger than amount of donors
    with pytest.raises(ValueError):
        chunk_object_single_run.get_random_cohort_sample(all_donors=s_all_donors, n_selections=7)
    # Expects output equal input
    actual = chunk_object_single_run.get_random_cohort_sample(
        all_donors=s_all_donors, n_selections=6
    )
    assert set(actual) == {d.index for d in s_all_donors}, "Should return same input."
    assert all([str(type(donor)) == expected_donor_type for donor in actual])
    # Expects output equal input
    actual = chunk_object_single_run.get_random_cohort_sample(
        all_donors=s_all_donors, n_selections=3
    )
    assert len(set(actual)) == 1, "Should return one index, i.e., three donors."
    assert all([str(type(donor)) == expected_donor_type for donor in actual]), assertion_type_msg
    # Expects output equal input
    actual = chunk_object_single_run.get_random_cohort_sample(
        all_donors=s_all_donors, n_selections=5
    )
    assert_msg = "Should return one index - i.e., three donors - as it won't break the pedigrees."
    assert len(set(actual)) == 1, assert_msg
    assert all([str(type(donor)) == expected_donor_type for donor in actual]), assertion_type_msg

    # ---------------------- #
    # Test background cohort #
    # ---------------------- #
    b_sheet = germline_sample_sheet_object_medium_cohort_solo_cases_background
    b_all_donors = b_sheet.cohort.pedigrees
    # Expects ValueError for request larger than amount of donors
    with pytest.raises(ValueError):
        chunk_object_single_run.get_random_cohort_sample(all_donors=b_all_donors, n_selections=200)
    # Expects output equal input
    actual = chunk_object_single_run.get_random_cohort_sample(
        all_donors=b_all_donors, n_selections=99
    )
    assert len(set(actual)) == 99, "Should return 99 indexes, i.e., all solo donors."
    assert all([str(type(donor)) == expected_donor_type for donor in actual]), assertion_type_msg
    # Expects output equal input
    actual = chunk_object_single_run.get_random_cohort_sample(
        all_donors=b_all_donors, n_selections=50
    )
    assert len(set(actual)) == 50, "Should return 50 unique donors."
    assert all([str(type(donor)) == expected_donor_type for donor in actual]), assertion_type_msg


def test_get_cohort_size_from_sheet(
    chunk_object_single_run,
    germline_sample_sheet_object,
    germline_sample_sheet_object_large_cohort_trios_only,
    germline_sample_sheet_object_large_cohort_diverse_cases,
):
    """Tests Chunk::get_cohort_size_from_sheet()"""
    # Define expected
    expected_small_cohort = 6
    expected_large_cohort_trio_only = 501
    expected_large_cohort_diverse = 200
    # Get actual - small cohort
    actual = chunk_object_single_run.get_cohort_size_from_sheet(sheet=germline_sample_sheet_object)
    assert actual == expected_small_cohort
    # Get actual - large cohort, trio cases only
    actual = chunk_object_single_run.get_cohort_size_from_sheet(
        sheet=germline_sample_sheet_object_large_cohort_trios_only
    )
    assert actual == expected_large_cohort_trio_only
    # Get actual - large cohort, diverse cases
    actual = chunk_object_single_run.get_cohort_size_from_sheet(
        sheet=germline_sample_sheet_object_large_cohort_diverse_cases
    )
    assert actual == expected_large_cohort_diverse


def test_get_cohort_size_from_donors_list(
    chunk_object_single_run,
    germline_sample_sheet_object,
    germline_sample_sheet_object_large_cohort_trios_only,
    germline_sample_sheet_object_large_cohort_diverse_cases,
):
    """Tests Chunk::get_cohort_size_from_donors_list()"""
    # Test empty list
    expected_empty_list = 0
    actual = chunk_object_single_run.get_cohort_size_from_donors_list(donors_list=[])
    assert actual == expected_empty_list
    # Test small cohort
    expected_small_cohort = 6
    s_sheet = germline_sample_sheet_object
    s_all_donors = s_sheet.cohort.pedigrees
    actual = chunk_object_single_run.get_cohort_size_from_donors_list(donors_list=s_all_donors)
    assert actual == expected_small_cohort
    # Test large cohort, trio cases only
    expected_large_cohort_trio_only = 501
    lt_sheet = germline_sample_sheet_object_large_cohort_trios_only
    lt_all_donors = lt_sheet.cohort.pedigrees
    actual = chunk_object_single_run.get_cohort_size_from_donors_list(donors_list=lt_all_donors)
    assert actual == expected_large_cohort_trio_only
    # Test large cohort, diverse cases
    expected_large_cohort_diverse = 200
    ld_sheet = germline_sample_sheet_object_large_cohort_diverse_cases
    ld_all_donors = ld_sheet.cohort.pedigrees
    actual = chunk_object_single_run.get_cohort_size_from_donors_list(donors_list=ld_all_donors)
    assert actual == expected_large_cohort_diverse


def test_build_pedigree_count_per_index_dict(
    chunk_object_single_run,
    germline_sample_sheet_object,
    germline_sample_sheet_object_large_cohort_trios_only,
):
    """Tests Chunk::build_pedigree_count_per_index_dict()"""
    # Initialise variables
    small_sheet = germline_sample_sheet_object
    trio_only_sheet = germline_sample_sheet_object_large_cohort_trios_only

    # Test small cohort
    size, out_dict = chunk_object_single_run.build_pedigree_count_per_index_dict(
        all_pedigrees=small_sheet.cohort.pedigrees
    )
    assert size == 6, "Should return 6, individuals P001 - P006."
    assert all([k == 3 for k in out_dict.keys()]), "Only trio cases in cohort."
    assert all([len(v) == 2 for v in out_dict.values()]), "Two indexes: P001, P004."

    # Test large cohort
    size, out_dict = chunk_object_single_run.build_pedigree_count_per_index_dict(
        all_pedigrees=trio_only_sheet.cohort.pedigrees
    )
    assert size == 501, "Should return 501; individuals: P001 - P501."
    assert all([k == 3 for k in out_dict.keys()]), "Only trio cases in cohort."
    assert all([len(v) == 167 for v in out_dict.values()]), "Contains 167 indexes."


def test_insufficient_space_exception():
    """Tests BatchInsufficientSpaceException"""
    # Initialise variables
    batch_counter_dict = {0: 29, 1: 29, 2: 28}
    # Define expected
    expected_msg = (
        "Iterated over all batches and none has enough space. Max batch size: 30; "
        "current count: 3.\nBatch structure:\n"
        "batch 0 count: 29\n"
        "batch 1 count: 29\n"
        "batch 2 count: 28\n"
    )
    # Get actual
    with pytest.raises(BatchInsufficientSpaceException) as e:
        raise BatchInsufficientSpaceException(
            max_batch_size=30, size_key=3, counter_dict=batch_counter_dict
        )
    actual_msg = e.value.args[0]
    assert actual_msg == expected_msg


def test_batchlist_to_dictionary(chunk_object_single_run):
    """Tests Chunk::batchlist_to_dictionary()"""
    # Define expected
    silly_list1 = [0, 1, 2]
    silly_list2 = [0, 3, 6, 9]
    silly_list3 = [11, 22, 33, 44]
    input_ = [silly_list1, silly_list2, silly_list3]
    expected = {1: silly_list1, 2: silly_list2, 3: silly_list3}
    # Get actual
    actual = chunk_object_single_run.batchlist_to_dictionary(input_)
    assert actual == expected


def test_split_cohort_into_balanced_chunks(
    chunk_object_single_run, germline_sample_sheet_object_large_cohort_trios_only
):
    """Tests Chunk::split_cohort_into_balanced_chunks()"""
    # -- Test balanced families, all trio cases -- #
    # Initialise variables
    max_batch_size = 200
    balanced_sheet = germline_sample_sheet_object_large_cohort_trios_only
    balanced_pedigree_list = balanced_sheet.cohort.pedigrees
    # Create dict with size of each pedigree
    pedigree_size_dict = build_pedigree_size_dictionary(
        all_pedigrees=balanced_sheet.cohort.pedigrees
    )
    # Get actual
    actual = chunk_object_single_run.split_cohort_into_balanced_chunks(
        all_pedigrees=balanced_pedigree_list, max_batch_size=max_batch_size
    )
    # Expects three lists: ceiling(501 / 200)
    assert len(actual) == 3
    # Expects no batch larger than max
    for batch in actual:
        size_pedigrees_list = []
        for donor in batch:
            id_ = donor.wrapped.secondary_id
            size_pedigrees_list.append(pedigree_size_dict.get(id_))
        assert sum(size_pedigrees_list) < max_batch_size


def test_is_appropriate_batch_size(chunk_object_single_run):
    """Tests Chunk::is_appropriate_batch_size()"""
    # Initialise variables
    fake_batch_0 = 0
    fake_batch_50 = 50
    valid_assert_msg = "Should return 'True'. A lot of space available in batch."
    invalid_assert_msg = "Should return 'False'. Although batch is empty pedigree is too large."

    # Test empty batch
    assert chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_0, max_batch_size=50, pedigree_count=5
    ), valid_assert_msg
    assert chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_0, max_batch_size=25, pedigree_count=5
    ), valid_assert_msg
    assert chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_0, max_batch_size=50, pedigree_count=50
    ), valid_assert_msg
    assert not chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_0, max_batch_size=50, pedigree_count=100
    ), invalid_assert_msg

    # Test batch with elements
    assert chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_50, max_batch_size=100, pedigree_count=5
    ), valid_assert_msg
    assert chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_50, max_batch_size=100, pedigree_count=5
    ), valid_assert_msg
    assert chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_50, max_batch_size=51, pedigree_count=1
    ), valid_assert_msg
    assert not chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_50, max_batch_size=50, pedigree_count=2
    ), invalid_assert_msg
    assert not chunk_object_single_run.is_appropriate_batch_size(
        batch_count=fake_batch_50, max_batch_size=51, pedigree_count=2
    ), invalid_assert_msg


def test_split_cohort_into_chunks(chunk_object_single_run):
    """Tests Chunk::split_cohort_into_chunks()"""
    # -- Test Even Sized input -- #
    # Define method input
    full_list = list(range(100))
    # Get actual
    actual = chunk_object_single_run.split_cohort_into_chunks(
        all_donors_list=full_list, max_batch_size=25
    )
    # Should return a list
    actual_type_str = str(type(actual))
    assert isinstance(actual, list), "Expected 'list', instead got {t}".format(t=actual_type_str)
    # Expects four lists: 100 / 25
    assert len(actual) == 4
    # All lists should have the same size
    assert all([len(chunk) == 25 for chunk in actual])

    # -- Test Uneven Sized input -- #
    # Define method input
    full_list = list(range(60))
    expected_chunk_size_dict = {1: 25, 2: 25, 3: 10}
    # Get actual
    actual = chunk_object_single_run.split_cohort_into_chunks(
        all_donors_list=full_list, max_batch_size=25
    )
    # Should return a list
    actual_type_str = str(type(actual))
    assert isinstance(actual, list), "Expected 'list', instead got {t}".format(t=actual_type_str)
    # Expects three lists: ceiling( 60 / 25 )
    assert len(actual) == 3
    # All lists should have the same size
    assert all([len(chunk) == expected_chunk_size_dict.get(i) for i, chunk in enumerate(actual, 1)])


def test_method_single_chunk(germline_sample_sheet_object):
    """Tests Chunk::single_chunk() for small cohort"""
    # Expected donors
    expected_donors = ["P001", "P004"]

    # Get germline sheet
    sheet = germline_sample_sheet_object

    # Run method
    single_chunk_run_out = Chunk(method="single", sheet_list=[sheet]).run()

    # Expects a single key
    assert len(single_chunk_run_out) == 1
    for donors_list in single_chunk_run_out.values():
        # Expects: P001, P004
        assert len(donors_list) == 2
        for donor in donors_list:
            assert donor.wrapped.secondary_id in expected_donors


def test_method_single_chunk_large_cohort(germline_sample_sheet_object_large_cohort_trios_only):
    """Tests Chunk::single_chunk() for large cohort"""
    # Expected donors: P001, P004, ... P496, P499
    expected_donors = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 501, 3)]

    # Get germline sheet
    sheet = germline_sample_sheet_object_large_cohort_trios_only

    # Run method
    single_chunk_run_out = Chunk(method="single", sheet_list=[sheet]).run()

    # Expects a single key
    assert len(single_chunk_run_out) == 1
    for donors_list in single_chunk_run_out.values():
        # Expects 167 trios with 501 patients
        assert len(donors_list) == 167
        for donor in donors_list:
            assert donor.wrapped.secondary_id in expected_donors


def test_method_evenly_chunk(germline_sample_sheet_object):
    """Tests Chunk::_evenly_chunk() for small cohort"""
    # Expected donors
    expected_donors = ["P001", "P004"]

    # Get germline sheet
    sheet = germline_sample_sheet_object
    with pytest.warns(UserWarning):
        evenly_chunk_run_out = Chunk(
            method="evenly", sheet_list=[sheet], maximal_chunk_size=5
        ).run()
    # Expects a dictionary with a single entry
    assert len(evenly_chunk_run_out) == 1
    for donors_list in evenly_chunk_run_out.values():
        # Expects: P001, P004
        assert len(donors_list) == 2
        for donor in donors_list:
            assert donor.wrapped.secondary_id in expected_donors


def test_method_evenly_chunk_large_cohort_trio_only(
    germline_sample_sheet_object_large_cohort_trios_only,
):
    """Tests Chunk::_evenly_chunk() for large cohort - trio cases only"""
    # Initialise variable
    max_chunk = 200
    # Expected donors: P001, P004, ... P496, P499
    expected_donors = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 501, 3)]
    # Donor seen count dictionary
    expected_donors_seen_count_dict = {}
    for donor in expected_donors:
        expected_donors_seen_count_dict[donor] = 0

    # Get germline sheet
    sheet = germline_sample_sheet_object_large_cohort_trios_only

    # Create dict with size of each pedigree
    pedigree_size_dict = build_pedigree_size_dictionary(all_pedigrees=sheet.cohort.pedigrees)

    # Run method
    evenly_chunk_run_out = Chunk(
        method="evenly", sheet_list=[sheet], maximal_chunk_size=max_chunk
    ).run()

    # Expects three keys/batches: ceiling( 501 / 200 )
    assert len(evenly_chunk_run_out) == 3
    for donors_list in evenly_chunk_run_out.values():
        size_pedigrees_list = []
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            expected_donors_seen_count_dict[id_] = expected_donors_seen_count_dict.get(id_) + 1
            size_pedigrees_list.append(pedigree_size_dict.get(id_))
            # Expects all indexes names to be valid
            assert id_ in expected_donors
        # Expects no batch larger than max
        assert sum(size_pedigrees_list) <= max_chunk
    # Expects that donors only seen once
    assert all([count == 1 for count in expected_donors_seen_count_dict.values()])


def test_method_evenly_chunk_large_cohort_diverse_cases(
    germline_sample_sheet_object_large_cohort_diverse_cases,
):
    """Tests Chunk::evenly_chunk() for large cohort - diverse cases"""
    # Initialise variable
    max_chunk = 50
    # Expected donors: P001, P004, ... P199, P200
    expected_donors = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 91, 3)]  # trios
    expected_donors += ["P{i}".format(i=str(i).zfill(3)) for i in range(91, 151, 2)]  # duos
    expected_donors += ["P{i}".format(i=str(i).zfill(3)) for i in range(151, 201, 1)]  # solos
    # Donor seen count dictionary
    expected_donors_seen_count_dict = {}
    for donor in expected_donors:
        expected_donors_seen_count_dict[donor] = 0

    # Get germline sheet
    sheet = germline_sample_sheet_object_large_cohort_diverse_cases

    # Create dict with size of each pedigree
    pedigree_size_dict = build_pedigree_size_dictionary(all_pedigrees=sheet.cohort.pedigrees)

    # Run method
    evenly_chunk_run_out = Chunk(
        method="evenly", sheet_list=[sheet], maximal_chunk_size=max_chunk
    ).run()

    # Expects four keys/batches: 200 / 50
    assert len(evenly_chunk_run_out) == 4
    for donors_list in evenly_chunk_run_out.values():
        size_pedigrees_list = []
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            expected_donors_seen_count_dict[id_] = expected_donors_seen_count_dict.get(id_) + 1
            size_pedigrees_list.append(pedigree_size_dict.get(id_))
            # Expects all indexes names to be valid
            assert id_ in expected_donors
        # Expects no batch larger than max
        assert sum(size_pedigrees_list) <= max_chunk
    # Expects that donors only seen once
    assert all([count == 1 for count in expected_donors_seen_count_dict.values()])


def test_method_incremental_chunk(germline_sample_sheet_object):
    """Tests Chunk::_incremental_chunk() for small cohort"""
    # Expected donors
    expected_donors = ["P001", "P004"]

    # Get germline sheet
    sheet = germline_sample_sheet_object
    with pytest.warns(UserWarning):
        evenly_chunk_run_out = Chunk(
            method="incremental", sheet_list=[sheet], maximal_chunk_size=5
        ).run()
    # Expects a dictionary with a single entry
    assert len(evenly_chunk_run_out) == 1
    for donors_list in evenly_chunk_run_out.values():
        # Expects: P001, P004
        assert len(donors_list) == 2
        for donor in donors_list:
            assert donor.wrapped.secondary_id in expected_donors


def test_method_incremental_chunk_large_cohort_trio_only(
    germline_sample_sheet_object_large_cohort_trios_only,
):
    """Tests Chunk::_incremental_chunk() for large cohort - trio cases only"""
    # Initialise variable
    batch_to_donors_id_dict = defaultdict(list)
    max_chunk = 200
    # Expected donors: P001, P004, ... P496, P499
    expected_donors_all = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 501, 3)]
    # Donor seen count dictionary
    expected_donors_seen_count_dict = {}
    for donor in expected_donors_all:
        expected_donors_seen_count_dict[donor] = 0

    # Get germline sheet
    sheet = germline_sample_sheet_object_large_cohort_trios_only

    # Create dict with size of each pedigree
    pedigree_size_dict = build_pedigree_size_dictionary(all_pedigrees=sheet.cohort.pedigrees)

    # Call method
    chunk_obj = Chunk(method="incremental", sheet_list=[sheet], maximal_chunk_size=max_chunk)
    evenly_chunk_run_out = chunk_obj.run()

    # Expects three keys/batches: ceiling( 501 / 200 )
    assert len(evenly_chunk_run_out) == 3
    for i_batch in evenly_chunk_run_out.keys():
        donors_list = evenly_chunk_run_out.get(i_batch)
        size_pedigrees_list = []
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            expected_donors_seen_count_dict[id_] = expected_donors_seen_count_dict.get(id_) + 1
            size_pedigrees_list.append(pedigree_size_dict.get(id_))
            batch_to_donors_id_dict[i_batch].append(int(id_.replace("P", "")))
            # Expects all indexes names to be valid
            assert id_ in expected_donors_all
        # Expects no batch larger than max
        assert sum(size_pedigrees_list) <= max_chunk
    # Expects that donors only seen once
    assert all([count == 1 for count in expected_donors_seen_count_dict.values()])

    # Expects that all secondary ids in first batches smaller than the subsequent ones
    last_batch = len(batch_to_donors_id_dict)  # 1-based index
    for i_batch in sorted(batch_to_donors_id_dict.keys()):
        # Nothing to compare last batch to
        if i_batch == last_batch:
            break
        curr_max_id = max(batch_to_donors_id_dict.get(i_batch))
        next_batch = i_batch + 1
        next_min_id = min(batch_to_donors_id_dict.get(next_batch))
        assert curr_max_id < next_min_id


def test_method_incremental_chunk_large_cohort_diverse_cases(
    germline_sample_sheet_object_large_cohort_diverse_cases,
):
    """Tests Chunk::_incremental_chunk() for large cohort - diverse cases"""
    # Initialise variable
    max_chunk = 90

    # Get germline sheet
    sheet = germline_sample_sheet_object_large_cohort_diverse_cases

    # Expected donors: P001, P004, ... P199, P200
    batch_1 = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 91, 3)]  # trios (90 samples)
    batch_2 = ["P{i}".format(i=str(i).zfill(3)) for i in range(91, 151, 2)]  # duos (60 samples)
    batch_2 += ["P{i}".format(i=str(i).zfill(3)) for i in range(151, 181, 1)]  # solos (30 samples)
    batch_3 = ["P{i}".format(i=str(i).zfill(3)) for i in range(181, 201, 1)]  # solos (20 samples)
    expected_donors_id_dict = {
        1: batch_1,
        2: batch_2,
        3: batch_3,
    }

    # Create dict with size of each pedigree
    pedigree_size_dict = build_pedigree_size_dictionary(all_pedigrees=sheet.cohort.pedigrees)

    # Call method
    chunk_obj = Chunk(method="incremental", sheet_list=[sheet], maximal_chunk_size=max_chunk)
    incremental_chunk_run_out = chunk_obj.run()

    # Expects three keys/batches: ceiling( 200 / 90 )
    assert len(incremental_chunk_run_out) == 3

    # Iterate over batches
    for i_batch in sorted(incremental_chunk_run_out.keys()):
        donors_list = incremental_chunk_run_out.get(i_batch)
        size_pedigrees_list = []
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            size = pedigree_size_dict.get(id_)
            # Expects that:
            #   batch 1 -> size 3 (trios only)
            #   batch 2 -> size 1 or 2 (duo and solo)
            #   batch 3 -> size 1 (solo only)
            if i_batch == 1:
                assert size == 3, "Expects that batch 1 only contains trio cases."
            elif i_batch == 2:
                assert size in (1, 2), "Expects that batch 2 contains duo and solo cases."
            elif i_batch == 3:
                assert size == 1, "Expects that batch 3 contains only solo cases."
            size_pedigrees_list.append(size)
            # Expects all indexes names to be valid
            i_batch_donors_list = expected_donors_id_dict.get(i_batch)
            assert id_ in i_batch_donors_list
        # Expects no batch larger than max
        assert sum(size_pedigrees_list) <= max_chunk


def test_method_incremental_chunk_background_and_foreground_cohorts(
    germline_sample_sheet_object_large_cohort_diverse_cases,
    germline_sample_sheet_object_medium_cohort_solo_cases_background,
):
    """Tests Chunk::_incremental_chunk() for foreground and background samples"""
    # Initialise variable
    max_chunk = 90

    # Get germline sheets
    f_sheet = germline_sample_sheet_object_large_cohort_diverse_cases
    b_sheet = germline_sample_sheet_object_medium_cohort_solo_cases_background

    # Expected donors: P001, P004, ... P199, P200
    f_batch_1 = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 91, 3)]  # trios (90 samples)
    f_batch_2 = ["P{i}".format(i=str(i).zfill(3)) for i in range(91, 151, 2)]  # duos (60 samples)
    f_batch_2 += ["P{i}".format(i=str(i).zfill(3)) for i in range(151, 181, 1)]  # solos(30 samples)
    f_batch_3 = ["P{i}".format(i=str(i).zfill(3)) for i in range(181, 201, 1)]  # solos (20 samples)
    b_batch = ["P{i}".format(i=str(i).zfill(3)) for i in range(900, 1000, 1)]  # solos (99 samples)

    # Create dict with size of each pedigree
    pedigree_size_dict = {
        **build_pedigree_size_dictionary(all_pedigrees=f_sheet.cohort.pedigrees),
        **build_pedigree_size_dictionary(all_pedigrees=b_sheet.cohort.pedigrees),
    }

    # Call method
    chunk_obj = Chunk(
        method="incremental", sheet_list=[f_sheet, b_sheet], maximal_chunk_size=max_chunk
    )
    incremental_chunk_run_out = chunk_obj.run()

    # Expects three keys/batches: ceiling( 200 / 90 )
    assert len(incremental_chunk_run_out) == 3

    # Iterate over batches
    for i_batch in sorted(incremental_chunk_run_out.keys()):
        donors_list = incremental_chunk_run_out.get(i_batch)
        batch_size = 0
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            batch_size += pedigree_size_dict.get(id_)
            # Expects that:
            #   batch 1 -> size 3 (trios only - foreground)
            #   batch 2 -> size 1 or 2 (duo and solo - foreground)
            #   batch 3 -> size 1 (solo only - foreground and background)
            if i_batch == 1:
                assert_msg = "Expects that batch 1 only contains trio cases, foreground only."
                assert id_ in f_batch_1, assert_msg
            elif i_batch == 2:
                assert_msg = "Expects that batch 2 contains duo and solo cases, foreground only."
                assert id_ in f_batch_2, assert_msg
            elif i_batch == 3:
                assert_msg = (
                    "Expects that batch 3 contains only solo cases, both foreground and background."
                )
                assert id_ in f_batch_3 or id_ in b_batch, assert_msg
        # Expects each batch exactly same value as max_chunk - used background samples
        assert batch_size == max_chunk


def test_method_incremental_chunk_sort_by_custom_field(
    germline_sample_sheet_object_medium_cohort_diverse_with_custom_features,
):
    """Tests Chunk::_incremental_chunk() for sample sheet with custom field"""
    # Initialise variable
    max_chunk = 50
    trio_cases_list = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 31, 3)]  # 10 trio cases
    duo_cases_list = ["P{i}".format(i=str(i).zfill(3)) for i in range(31, 51, 2)]  # 10 duo cases
    solo_cases_list = ["P{i}".format(i=str(i).zfill(3)) for i in range(51, 101, 1)]  # 50 solo cases

    # Get germline sheet
    sheet = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features

    # Call method - order by family id first
    chunk_obj = Chunk(
        method="incremental",
        sheet_list=[sheet],
        maximal_chunk_size=max_chunk,
        order_by_custom_field=["familyId", "batchNo"],
    )
    incremental_chunk_run_out = chunk_obj.run()
    assert len(incremental_chunk_run_out) == 2, "Expects two batches, 50 samples in each."
    # Iterate over batches
    for i_batch in sorted(incremental_chunk_run_out.keys()):
        donors_list = incremental_chunk_run_out.get(i_batch)
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            # Expects when ordering by family id:
            #   batch 1 -> duo and trios cases
            #   batch 2 -> solo cases only
            if i_batch == 1:
                assert_msg = "Expects that batch 1 contains trio and duo cases."
                assert id_ in trio_cases_list or id_ in duo_cases_list, assert_msg
            elif i_batch == 2:
                assert_msg = "Expects that batch 2 contains only solo cases."
                assert id_ in solo_cases_list, assert_msg

    # Call method - order by batch number first
    chunk_obj = Chunk(
        method="incremental",
        sheet_list=[sheet],
        maximal_chunk_size=max_chunk,
        order_by_custom_field=["batchNo", "familyId"],
    )
    incremental_chunk_run_out = chunk_obj.run()
    assert len(incremental_chunk_run_out) == 2, "Expects two batches, 50 samples in each."
    # Iterate over batches
    for i_batch in sorted(incremental_chunk_run_out.keys()):
        donors_list = incremental_chunk_run_out.get(i_batch)
        for donor in donors_list:
            id_ = donor.wrapped.secondary_id
            # Expects when ordering by batch number:
            #   batch 1 -> solo cases only
            #   batch 2 -> duo and trios cases
            if i_batch == 1:
                assert_msg = "Expects that batch 1 contains only solo cases."
                assert id_ in solo_cases_list, assert_msg
            elif i_batch == 2:
                assert_msg = "Expects that batch 2 contains duo and trio cases."
                assert id_ in trio_cases_list or id_ in duo_cases_list, assert_msg
