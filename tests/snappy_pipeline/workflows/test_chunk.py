# -*- coding: utf-8 -*-
"""Code for testing the code in ``Chunk`` class.
"""

from collections import defaultdict
from operator import attrgetter

from biomedsheets.shortcuts import GermlineCaseSheet, is_background, is_not_background
import pytest

from snappy_pipeline.chunk import (
    BatchInsufficientSpaceException,
    Chunk,
    ChunkAuxiliary,
    ChunkHistory,
)


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


# Test Chunk =======================================================================================


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


def test_split_pedigress_by_sequencing_kit(
    chunk_object_single_run, germline_sample_sheet_object_medium_cohort_solo_only_diverse_kits
):
    """TestsChunk::split_pedigrees_by_sequencing_kit()"""
    # Define expected
    expected_kit_list = ["Agilent SureSelect Human All Exon V6", "Illumina TruSeq PCR-free"]
    expected_n_pedigrees_split_list = [49, 50]

    # Get sheet
    sheet = germline_sample_sheet_object_medium_cohort_solo_only_diverse_kits
    all_pedigrees = sheet.cohort.pedigrees

    # Call method
    split_pedigrees_list = chunk_object_single_run.split_pedigrees_by_sequencing_kit(
        all_pedigrees=all_pedigrees
    )

    # Expects two lists with number of elements in [49, 50]
    assert len(split_pedigrees_list) == 2
    assert len(split_pedigrees_list[0]) != len(split_pedigrees_list[1])
    assert len(split_pedigrees_list[0]) in expected_n_pedigrees_split_list
    assert len(split_pedigrees_list[1]) in expected_n_pedigrees_split_list

    # Iterate over subsets of pedigrees
    for pedigrees in split_pedigrees_list:
        first_entry = True
        subset_library_kit_name = None
        for pedigree in pedigrees:
            if first_entry:
                # Get first kit name
                subset_library_kit_name = (
                    pedigree.index.dna_ngs_library.ngs_library.extra_infos.get("libraryKit")
                )
                continue
            # Expects that library kits in lists are the same across the subset
            i_lib_name = pedigree.index.dna_ngs_library.ngs_library.extra_infos.get("libraryKit")
            assert i_lib_name == subset_library_kit_name
            # Expects that library name in expected list (XOR)
            assert (i_lib_name == expected_kit_list[0]) != (i_lib_name == expected_kit_list[1])


def test_order_pedigrees_by_sampleid(chunk_object_single_run, germline_sample_sheet_object):
    """Tests Chunk::order_pedigrees_by_sampleid()"""
    # Get sheet
    sheet = germline_sample_sheet_object
    all_pedigrees = sheet.cohort.pedigrees

    # Reverse order and sanity check
    all_pedigrees_reversed = all_pedigrees[::-1]
    assert not all(map(lambda x, y: x == y, all_pedigrees, all_pedigrees_reversed))

    # Call method
    restored_order = chunk_object_single_run.order_pedigrees_by_sampleid(all_pedigrees_reversed)

    # Expects that output in the same order as originally set in sample sheet
    assert all(map(lambda x, y: x == y, all_pedigrees, restored_order))
    assert not all(map(lambda x, y: x == y, all_pedigrees, all_pedigrees_reversed))


def test_get_max_field_value_from_donors(
    chunk_object_single_run, germline_sample_sheet_object_small_cohort_with_custom_features
):
    """Tests Chunk::get_max_field_value_from_donors()"""
    # Initialise variables
    field = "batchNo"
    # Define expected
    expected_dict = {"P001": 3, "P004": 2}
    # Get sheet and iterate over pedigrees
    sheet = germline_sample_sheet_object_small_cohort_with_custom_features
    all_pedigrees = sheet.cohort.pedigrees
    for pedigree in all_pedigrees:
        # Get actual
        sampleid = pedigree.index.wrapped.secondary_id
        max_value = chunk_object_single_run.get_max_field_value_from_donors(
            pedigree=pedigree, field=field
        )
        assert expected_dict.get(sampleid) == max_value


def test_order_pedigrees_by_field(
    chunk_object_single_run,
    germline_sample_sheet_object_medium_cohort_diverse_with_custom_features,
):
    """Tests Chunk::order_pedigrees_by_field()"""
    # Define expected
    # Possible point of confusion: the batches were assigned in the reversed order
    # of family ids for the purpose of testing.
    # sampleid: P001; familyid: FAM_P001; batch: 7 - first index of batch number 7
    first_familyid = "FAM_P001"
    first_batch = 7
    # sampleid: P028; familyid: FAM_P028; batch: 7 - last index of batch number 7
    last_fam_of_first_batch = "FAM_P028"
    # sampleid: P091; familyid: FAM_P091; batch: 1 - first index of batch number 1
    first_fam_of_last_batch = "FAM_P091"
    # sampleid: P100; familyid: FAM_P100; batch: 1 - last index of batch number 1
    last_familyid = "FAM_P100"
    last_batch = 1

    # Get sheet
    sheet = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features
    all_pedigrees = sheet.cohort.pedigrees

    # Call method - batch number
    reordered_pedigrees = chunk_object_single_run.order_pedigrees_by_field(
        all_pedigrees=all_pedigrees, fields=["batchNo"]
    )
    # Get actual
    actual_first_batch = reordered_pedigrees[0].index.extra_infos.get("batchNo")
    actual_first_familyid = reordered_pedigrees[0].index.extra_infos.get("familyId")
    actual_last_batch = reordered_pedigrees[-1].index.extra_infos.get("batchNo")
    actual_last_familyid = reordered_pedigrees[-1].index.extra_infos.get("familyId")
    assert actual_first_batch == last_batch
    assert actual_first_familyid == first_fam_of_last_batch
    assert actual_last_batch == first_batch
    assert actual_last_familyid == last_fam_of_first_batch
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
    reordered_pedigrees = chunk_object_single_run.order_pedigrees_by_field(
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
    actual = chunk_object_single_run.get_cohort_size_from_pedigrees_list(all_pedigrees=[])
    assert actual == expected_empty_list
    # Test small cohort
    expected_small_cohort = 6
    s_sheet = germline_sample_sheet_object
    s_all_donors = s_sheet.cohort.pedigrees
    actual = chunk_object_single_run.get_cohort_size_from_pedigrees_list(all_pedigrees=s_all_donors)
    assert actual == expected_small_cohort
    # Test large cohort, trio cases only
    expected_large_cohort_trio_only = 501
    lt_sheet = germline_sample_sheet_object_large_cohort_trios_only
    lt_all_donors = lt_sheet.cohort.pedigrees
    actual = chunk_object_single_run.get_cohort_size_from_pedigrees_list(
        all_pedigrees=lt_all_donors
    )
    assert actual == expected_large_cohort_trio_only
    # Test large cohort, diverse cases
    expected_large_cohort_diverse = 200
    ld_sheet = germline_sample_sheet_object_large_cohort_diverse_cases
    ld_all_donors = ld_sheet.cohort.pedigrees
    actual = chunk_object_single_run.get_cohort_size_from_pedigrees_list(
        all_pedigrees=ld_all_donors
    )
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


def test_pedigree_to_index_list(chunk_object_single_run, germline_sample_sheet_object):
    """Tests Chunk::pedigree_to_index_list()"""
    # Expected donors
    expected_donors_id_list = ["P001", "P004"]
    # Get germline sheet
    all_pedigrees = germline_sample_sheet_object.cohort.pedigrees
    # Run method
    index_list = chunk_object_single_run.pedigree_to_index_list(all_pedigrees=all_pedigrees)
    # Expects list with two indexes
    assert isinstance(index_list, list)
    assert len(index_list) == 2
    for donor in index_list:
        assert donor.wrapped.secondary_id in expected_donors_id_list


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


def test_method_incremental_chunk_mixed_sequencing_kits(
    germline_sample_sheet_object_medium_cohort_solo_only_diverse_kits,
    germline_sample_sheet_object_large_cohort_trios_only,
):
    """Tests Chunk::_incremental_chunk() for mix of sequencing kits in sample sheets"""
    # Initialise variables
    expected_kit_list = ["Agilent SureSelect Human All Exon V6", "Illumina TruSeq PCR-free"]
    max_chunk = 50

    # Get sheets
    l_sheet = germline_sample_sheet_object_large_cohort_trios_only
    m_sheet = germline_sample_sheet_object_medium_cohort_solo_only_diverse_kits

    # Run method
    incremental_run_out = Chunk(
        method="incremental", sheet_list=[l_sheet, m_sheet], maximal_chunk_size=max_chunk
    ).run()

    # Expects 12 chunks with ~50 Agilent samples each, and 1 chunk with 50 TruSeq
    # Agilent batches: (16 trios), (16 trios), ..., (16 trios), (10 trios, 20 solos), (30 solos)
    # TruSeq batch: (50 solos)
    assert len(incremental_run_out) == 13
    # Iterate over batches
    for i_batch in sorted(incremental_run_out.keys()):
        donors = incremental_run_out.get(i_batch)
        first_entry = True
        subset_library_kit_name = None
        for donor in donors:
            if first_entry:
                # Get first kit name
                subset_library_kit_name = donor.dna_ngs_library.ngs_library.extra_infos.get(
                    "libraryKit"
                )
                first_entry = False
                continue
            # Expects that library kits in lists are the same across for all subset
            i_lib = donor.dna_ngs_library.ngs_library.extra_infos.get("libraryKit")
            assert i_lib == subset_library_kit_name
            # Expects that library name in expected list
            assert i_lib in expected_kit_list


def test_method_incremental_chunk_mixed_sequencing_kits_and_background(
    germline_sample_sheet_object_medium_cohort_solo_only_diverse_kits,
    germline_sample_sheet_object_large_cohort_trios_only,
    germline_sample_sheet_object_medium_cohort_solo_cases_background,
):
    """Tests Chunk::_incremental_chunk() for mix of sequencing kits in sample sheets"""
    # Initialise variables
    max_chunk = 50

    # Get sheets
    l_sheet = germline_sample_sheet_object_large_cohort_trios_only
    m_sheet = germline_sample_sheet_object_medium_cohort_solo_only_diverse_kits
    b_sheet = germline_sample_sheet_object_medium_cohort_solo_cases_background

    # Create dict with size of each pedigree
    pedigree_size_dict = {
        **build_pedigree_size_dictionary(all_pedigrees=l_sheet.cohort.pedigrees),
        **build_pedigree_size_dictionary(all_pedigrees=m_sheet.cohort.pedigrees),
        **build_pedigree_size_dictionary(all_pedigrees=b_sheet.cohort.pedigrees),
    }

    # Run method
    incremental_run_out = Chunk(
        method="incremental", sheet_list=[l_sheet, m_sheet, b_sheet], maximal_chunk_size=max_chunk
    ).run()

    # Iterate over batches
    for i_batch in sorted(incremental_run_out.keys()):
        donors = incremental_run_out.get(i_batch)
        batch_size = 0
        for donor in donors:
            id_ = donor.wrapped.secondary_id
            batch_size += pedigree_size_dict.get(id_)
        # Expects all chunks with max given that background samples were provided
        assert batch_size == max_chunk


# Test ChunkHistory ================================================================================


def test_summarize_sample_sheet(
    germline_sample_sheet_object,
    germline_sample_sheet_object_small_cohort_with_custom_features,
    germline_sample_sheet_object_medium_cohort_diverse_with_custom_features,
):
    """Tests ChunkHistory::summarize_sample_sheet()"""
    # Initialise variables
    c_hist = ChunkHistory(sheet_list=[])

    # Rename fixtures
    sheet_exception = germline_sample_sheet_object
    sheet_small = germline_sample_sheet_object_small_cohort_with_custom_features
    sheet_medium = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features

    # Define expected
    expected_small_cohort = {1: 2, 3: 1, 2: 3}
    expected_medium_cohort = {1: 10, 2: 10, 3: 10, 4: 10, 5: 10, 6: 20, 7: 30}

    # Test exception - sheet must contain field 'batchNo'
    with pytest.raises(ValueError):
        c_hist.summarize_sample_sheet(sheet=sheet_exception)

    # Test small cohort
    actual = c_hist.summarize_sample_sheet(sheet=sheet_small)
    assert actual == expected_small_cohort

    # Test medium cohort
    actual = c_hist.summarize_sample_sheet(sheet=sheet_medium)
    assert actual == expected_medium_cohort


def test_cumulative_batch_count():
    """Tests ChunkHistory::cumulative_batch_count()"""
    # Initialise variables
    c_hist = ChunkHistory(sheet_list=[])

    # Define input
    in_dict_ones = defaultdict(lambda: 0, {1: 1, 2: 1, 3: 1, 4: 1})
    in_dict_twos = defaultdict(lambda: 0, {1: 2, 2: 2, 3: 2, 4: 2})

    # Define expected
    expected_one = {1: 1, 2: 2, 3: 3, 4: 4}
    expected_two = {1: 2, 2: 4, 3: 6, 4: 8}

    # All 1s
    actual = c_hist.cumulative_batch_count(batch_counter_dict=in_dict_ones)
    assert actual == expected_one

    # All 2s
    actual = c_hist.cumulative_batch_count(batch_counter_dict=in_dict_twos)
    assert actual == expected_two


def test_filter_sample_sheet(
    germline_sample_sheet_object_small_cohort_with_custom_features,
    germline_sample_sheet_object_medium_cohort_diverse_with_custom_features,
):
    """Tests ChunkHistory::filter_sample_sheet()"""
    # Initialise variables
    c_hist = ChunkHistory(sheet_list=[])

    # Rename fixture
    sheet_small = germline_sample_sheet_object_small_cohort_with_custom_features
    sheet_medium = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features

    # Batch smaller than two - small cohort
    max_batch_no = 2
    out_sheet = c_hist.filter_sample_sheet(sheet=sheet_small, max_batch_no=max_batch_no)
    for pedigree in out_sheet.cohort.pedigrees:
        for donor in pedigree.donors:
            batch_value = attrgetter("extra_infos")(donor).get("batchNo")
            assert batch_value <= max_batch_no

    # Batch smaller than four - medium cohort
    max_batch_no = 4
    out_sheet = c_hist.filter_sample_sheet(sheet=sheet_medium, max_batch_no=max_batch_no)
    for pedigree in out_sheet.cohort.pedigrees:
        for donor in pedigree.donors:
            batch_value = attrgetter("extra_infos")(donor).get("batchNo")
            assert batch_value <= max_batch_no

    # Batch smaller than two and not in list - medium cohort
    filter_donors = ["P{i}".format(i=str(i).zfill(3)) for i in range(91, 101)]  # remove 10 solos
    max_batch_no = 2
    out_sheet = c_hist.filter_sample_sheet(
        sheet=sheet_medium, max_batch_no=max_batch_no, secondary_ids=filter_donors
    )
    for pedigree in out_sheet.cohort.pedigrees:
        for donor in pedigree.donors:
            batch_value = attrgetter("extra_infos")(donor).get("batchNo")
            assert batch_value <= max_batch_no
            assert donor.wrapped.secondary_id not in filter_donors


def test_chunk_count(germline_sample_sheet_object):
    """Tests ChunkHistory::chunk_count()"""
    # Rename fixture
    sheet = germline_sample_sheet_object
    # Initialise variable
    c_hist = ChunkHistory(sheet_list=[sheet])
    # Run method - simple way to get results as they should look
    single_chunk = Chunk(method="single", sheet_list=[sheet]).run()
    # Define expected
    expected = 6  # number of samples in small germline sheet; indexes: P001, P004
    # Get actual
    actual = c_hist.chunk_count(chunk=list(single_chunk.values())[0])
    assert actual == expected


def test_neg_max_batch(germline_sample_sheet_object_small_cohort_with_custom_features):
    """Tests ChunkHistory::neg_max_batch()"""
    # Rename fixture
    sheet = germline_sample_sheet_object_small_cohort_with_custom_features
    # Initialise variable
    c_hist = ChunkHistory(sheet_list=[sheet])
    # Run method - simple way to get results as they should look
    single_chunk = Chunk(method="single", sheet_list=[sheet]).run()
    # Define expected
    expected = -2
    # Get actual
    actual = c_hist.neg_max_batch(chunk=list(single_chunk.values())[0])
    msg = "Max batch in cohort is '2', hence method should return '-2'."
    assert actual == expected, msg


def test_identify_closed_chunks(
    germline_sample_sheet_object_small_cohort_with_custom_features,
    germline_sample_sheet_object_medium_cohort_diverse_with_custom_features,
):
    """Tests ChunkHistory::identify_closed_chunks()"""
    # Rename fixture
    sheet_small = germline_sample_sheet_object_small_cohort_with_custom_features
    sheet_medium = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features

    # Initialise history object
    c_hist = ChunkHistory(sheet_list=[sheet_small, sheet_medium])

    # Run incremental method - simple way to get results as they should look
    incremental_small = Chunk(
        method="incremental", sheet_list=[sheet_small], maximal_chunk_size=40
    ).run()
    actual = c_hist.identify_closed_chunks(chunks=list(incremental_small.values()))
    # Given that there will be only one chunk for the small sheet,
    # it will return that none is closed.
    assert len(actual) == 0

    # Run incremental method - simple way to get results as they should look
    incremental_medium = Chunk(
        method="incremental", sheet_list=[sheet_medium], maximal_chunk_size=40
    ).run()
    actual = c_hist.identify_closed_chunks(chunks=list(incremental_medium.values()))
    # Given the medium sheet composition, it should return 2 closed chunks with max chunk size = 40.
    # Composition:
    #   - 10 trio cases, 30 samples
    #   - 10 duo cases, 20 samples
    #   - 50 solo cases, 50 samples
    assert len(actual) == 2


def test_history(germline_sample_sheet_object_medium_cohort_diverse_with_custom_features):
    """Tests ChunkHistory::history()"""
    # Initialise variables
    trio_cases_list = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 31, 3)]  # 10 trio cases
    duo_cases_list = ["P{i}".format(i=str(i).zfill(3)) for i in range(31, 51, 2)]  # 10 duo cases
    solo_cases_list = ["P{i}".format(i=str(i).zfill(3)) for i in range(51, 101, 1)]  # 50 solo cases
    index_count_dict = {}
    for index in trio_cases_list + duo_cases_list + solo_cases_list:
        index_count_dict[index] = 0

    # Rename fixture
    sheet = germline_sample_sheet_object_medium_cohort_diverse_with_custom_features

    # Initialise history object
    c_hist = ChunkHistory(
        sheet_list=[sheet],
        order_by_custom_field=["batchNo", "familyId"],
        maximal_chunk_size=30,
    )

    # Get history
    actual = c_hist.history()
    counter = 0
    for chunk in actual:
        counter += 1
        for index in chunk:
            index_count_dict[index.wrapped.secondary_id] += 1

    # Expected that all samples are present exactly once
    assert all([count == 1 for count in index_count_dict.values()])
