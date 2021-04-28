# -*- coding: utf-8 -*-
"""Code for testing the code in ``Chunk`` class.
"""

from collections import defaultdict
import io
import random
import textwrap

from biomedsheets.io_tsv import read_germline_tsv_sheet
from biomedsheets.shortcuts import GermlineCaseSheet
import pytest

from snappy_pipeline.chunk import BatchInsufficientSpaceException, Chunk


def build_pedigree_size_dictionary(all_pedigrees):
    """Pedigree to size dictionary.

    :param all_pedigrees: List with all donors in cohort.
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


def get_entry_for_sample_sheet(index_i, entry_type="trio"):
    """Sample sheet entry text.

    :param index_i: Index id. Example if '1', final id will be 'P001'
    :type index_i: int

    :param entry_type: Type of entry: 'solo', 'duo', 'trio'. Default: 'trio'.
    :type entry_type: str, optional

    :return: Returns text used to define a trio case in sample sheet.
    """
    # Initialise variables
    sampleid_pattern = "P{i}"
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
    index_sampleid = sampleid_pattern.format(i=str(index_i).zfill(3))
    mother_i = index_i + 1
    mother_sampleid = sampleid_pattern.format(i=str(mother_i).zfill(3))
    father_i = index_i + 2
    father_sampleid = sampleid_pattern.format(i=str(father_i).zfill(3))

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
def large_cohort_trios_only_text():
    """Defines arbitrary large cohort entries for sample sheet - trio cases only"""
    # Initialise variables
    n_patients_cohorts = 501
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
def germline_sheet_header():
    """Returns germline TSV file header with wildcard {entries}"""
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
def large_cohort_trios_only_germline_sheet_tsv(germline_sheet_header, large_cohort_trios_only_text):
    """Returns contents for large cohort germline TSV file - trio cases only"""
    return germline_sheet_header.format(entries=large_cohort_trios_only_text)


@pytest.fixture
def large_cohort_diverse_cases_germline_sheet_tsv(germline_sheet_header, large_cohort_diverse_text):
    """Returns contents for large cohort germline TSV file - diverse cases"""
    return germline_sheet_header.format(entries=large_cohort_diverse_text)


@pytest.fixture
def germline_sample_sheet_object(germline_sheet_tsv):
    """Returns GermlineCaseSheet object with small cohort."""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(germline_sheet_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def germline_sample_sheet_object_large_cohort_trios_only(
    large_cohort_trios_only_germline_sheet_tsv,
):
    """Returns GermlineCaseSheet object with large cohort - trio cases only."""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(large_cohort_trios_only_germline_sheet_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def germline_sample_sheet_object_large_cohort_diverse_cases(
    large_cohort_diverse_cases_germline_sheet_tsv,
):
    """Returns GermlineCaseSheet object with large cohort - trio cases only."""
    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(large_cohort_diverse_cases_germline_sheet_tsv)
    return GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))


@pytest.fixture
def chunk_object_single_run(germline_sample_sheet_object):
    """Returns Chunk object for single run."""
    # Get germline sheet
    sheet = germline_sample_sheet_object
    # Run single
    return Chunk(method="single", sheet_list=[sheet])


def test_chunk_constructor(germline_sample_sheet_object):
    """Tests Chunk::__init__()"""
    # Initialise variables
    valid_methods = ["single", "evenly", "incremental"]
    silly_sheet_list = [1, 2, 3]
    silly_dict = {1: "one"}
    sheet = germline_sample_sheet_object

    # Test only valid methods
    for method in valid_methods:
        c = Chunk(method=method, sheet_list=silly_sheet_list)
        assert isinstance(c, Chunk), "Expected instance of class `Chunk`."

    with pytest.raises(ValueError):
        Chunk(method="_invalid_method_", sheet_list=silly_sheet_list)

    # Test sample sheet list
    with pytest.raises(ValueError):
        Chunk(method="single", sheet_list=silly_dict)  # Expects list

    with pytest.raises(ValueError):
        Chunk(method="single", sheet_list=[])  # Expects list with at least one entry

    # Tests maximum chunk size - used in 'evenly' and 'incremental'
    # Expects value greater or equal to zero
    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], maximal_chunk_size=-2)

    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], maximal_chunk_size=0)

    with pytest.raises(ValueError):
        Chunk(method="evenly", sheet_list=[sheet], maximal_chunk_size="_zero_")

    with pytest.raises(ValueError):
        Chunk(method="incremental", sheet_list=[sheet], maximal_chunk_size=-2)

    with pytest.raises(ValueError):
        Chunk(method="incremental", sheet_list=[sheet], maximal_chunk_size=0)


def test_call_insufficient_space_exception(chunk_object_single_run):
    """Tests Chunk::call_insufficient_space_exception"""
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
        chunk_object_single_run.call_insufficient_space_exception(
            max_batch_size=30, i_size=3, counter_dict=batch_counter_dict
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
    pedigree_size_dict = build_pedigree_size_dictionary(all_pedigrees=balanced_sheet.cohort.pedigrees)
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
    actual_type = str(type(actual))
    assert isinstance(actual, list), "Expected 'list', instead got {t}".format(t=actual_type)
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
    actual_type = str(type(actual))
    assert isinstance(actual, list), "Expected 'list', instead got {t}".format(t=actual_type)
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

    # Run single
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

    # Run single
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
    """Tests Chunk::_evenly_chunk() for large cohort"""
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

    # Run single
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

    # Run single
    evenly_chunk_run_out = Chunk(
        method="evenly", sheet_list=[sheet], maximal_chunk_size=max_chunk
    ).run()

    # Expects three keys/batches: 200 / 50
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
