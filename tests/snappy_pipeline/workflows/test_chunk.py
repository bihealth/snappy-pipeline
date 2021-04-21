# -*- coding: utf-8 -*-
"""Code for testing the code in ``Chunk`` class.
"""

import io

from biomedsheets.io_tsv import read_germline_tsv_sheet
from biomedsheets.shortcuts import GermlineCaseSheet
import pytest

from snappy_pipeline.workflows.targeted_seq_cnv_calling.chunk import Chunk


def test_chunk_constructor():
    """Tests Chunk::__init__()"""
    # Initialise variables
    valid_methods = ["single", "evenly", "incremental"]
    silly_sheet_list = [1, 2, 3]
    silly_dict = {1: "one"}

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


def test_method_single_chunk(germline_sheet_tsv):
    """Tests Chunk::single_chunk() for small cohort"""
    # Expected donors
    expected_donors = ["P001", "P004"]

    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(germline_sheet_tsv)
    sheet = GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))
    # Run single
    single_chunk_run_out = Chunk(method="single", sheet_list=[sheet]).run()

    # Expects a single key
    assert len(single_chunk_run_out) == 1
    for donors_list in single_chunk_run_out.values():
        # Expects: P001, P004
        assert len(donors_list) == 2
        for donor in donors_list:
            assert donor.wrapped.secondary_id in expected_donors


def test_method_single_chunk_large_cohort(large_cohort_germline_sheet_tsv):
    """Tests Chunk::single_chunk() for large cohort"""
    # Expected donors: P001, P004, ... P496, P499
    expected_donors = ["P{i}".format(i=str(i).zfill(3)) for i in range(1, 501, 3)]

    # Create dna sample sheet based on germline sheet
    germline_sheet_io = io.StringIO(large_cohort_germline_sheet_tsv)
    sheet = GermlineCaseSheet(sheet=read_germline_tsv_sheet(germline_sheet_io))
    # Run single
    single_chunk_run_out = Chunk(method="single", sheet_list=[sheet]).run()

    # Expects a single key
    assert len(single_chunk_run_out) == 1
    for donors_list in single_chunk_run_out.values():
        # Expects 167 trios with 501 patients
        assert len(donors_list) == 167
        for donor in donors_list:
            assert donor.wrapped.secondary_id in expected_donors
