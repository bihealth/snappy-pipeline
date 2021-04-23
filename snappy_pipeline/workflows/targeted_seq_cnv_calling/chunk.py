# -*- coding: utf-8 -*-
""" Implementation of the ``chunk`` procedure

The cohort being analyzed will be splitted into chunks depending on the selected method:
- ``single``: create one large chunk for all.
- ``evenly``: create evenly sized chunks of maximal size.
- ``incremental``: iterate samples and create a new chunk after the maximal
                   size has been reached; should be combined with background data
                   to guarantee a minimal number of samples.
"""
from collections import defaultdict
from math import ceil
import warnings

from biomedsheets.shortcuts import is_not_background

from snappy_pipeline.utils import dictify, listify

# List of implemented split methods
ACTIVE_METHODS = ["single", "evenly", "incremental"]

# TODO: Figure out a good value for MIN_CHUNK_SIZE. Should it be a class variable and changeable?
# Minimum allowed chunk
MIN_CHUNK_SIZE = 50


class BatchInsufficientSpaceException(Exception):
    """Raised when none of the batches has enough space to accommodate pedigree"""


class Chunk:
    """Class contains methods to split cohort into chunks."""

    #: Selected chunk method
    selected_method = None

    #: List of Sample sheets
    sheet_list = None

    #: Maximum number of chunks - used in 'evenly' method
    maximal_number_of_chunks = None

    #: Maximum chunk size - used in 'evenly' method
    maximal_chunk_size = None

    def __init__(self, sheet_list, method, maximal_number_of_chunks=None, maximal_chunk_size=None):
        """Constructor

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :param method: Chunk method. Available implementations: 'single', 'evenly', 'incremental'.
        :type method: str

        :param maximal_number_of_chunks: Maximum number of chunks. Used in 'evenly' method to
        define how to split the cohort.
        :type maximal_number_of_chunks: int, optional

        :param maximal_chunk_size: Maximum chunk size - used in 'evenly' method.
        :type maximal_chunk_size: int, optional
        """
        # Validate selected method
        if method not in ACTIVE_METHODS:
            valid_methods = ", ".join(ACTIVE_METHODS)
            raise ValueError(
                "Method '{m}' is not implemented. Available methods: {am}".format(
                    m=method, am=valid_methods
                )
            )
        self.selected_method = method

        # Validate sample sheet list
        if len(sheet_list) == 0:
            raise ValueError("Sample Sheets list cannot be empty.")
        if not isinstance(sheet_list, list):
            raise ValueError(
                "Method was expecting 'list', instead got '{type_}'".format(type_=type(sheet_list))
            )
        self.sheet_list = sheet_list

        # Set number of chunks
        self.maximal_number_of_chunks = maximal_number_of_chunks

        # Set chunk size
        self.maximal_chunk_size = maximal_chunk_size

    def run(self):
        """Runs split in chunks selected method."""
        # Single Chunk
        if self.selected_method == "single":
            return self._single_chunk(sheet_list=self.sheet_list)
        # Evenly Chunk
        if self.selected_method == "evenly":
            return self._evenly_chunk(sheet_list=self.sheet_list, max_chunk=self.maximal_chunk_size)

        return None

    def _evenly_chunk(self, sheet_list, max_chunk):
        """Define a evenly divided chunks based on sample sheets list.

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :param max_chunk: Maximum number of samples per chunk.
        :type max_chunk: int

        :return: Returns dictionary with several chunks indices (keys) associated with the list of
        index donors from the cohorts defined in the sample sheets (values).
        """
        # Initialise variables
        donors_list = []

        # Validate max_chunk
        if (not isinstance(max_chunk, int)) or max_chunk <= 0:
            err_msg = (
                "Max number of samples per chunk must be an integer greater than zero. "
                "Invalid input: '{0}'.".format(max_chunk)
            )
            raise ValueError(err_msg)

        # Iterate over sheets
        sheets = list(filter(is_not_background, sheet_list))
        for sheet in sheets:
            if len(sheet.cohort.pedigrees) < MIN_CHUNK_SIZE:
                warn_msg = (
                    "Cohort in sample sheet contains only {count_} (<{min_}) "
                    "and won't be split into batches.".format(
                        count_=str(len(sheet.cohort.pedigrees)), min_=MIN_CHUNK_SIZE
                    )
                )
                warnings.warn(warn_msg)
                single_run_dict = self._single_chunk(sheet_list=[sheet])
                balanced_chunks_list = list(single_run_dict.values())
            else:
                balanced_chunks_list = self.split_cohort_into_balanced_chunks(
                    all_pedigrees=sheet.cohort.pedigrees, max_batch_size=max_chunk
                )
            for batch in balanced_chunks_list:
                donors_list.append(batch)
        # Return dictionary
        return self.batchlist_to_dictionary(batch_lists=donors_list)

    def _single_chunk(self, sheet_list):
        """Define a single chunk based on sample sheets list.

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :return: Returns dictionary with a single key ('1') associated with the list of index
        donors from the cohorts defined in the sample sheets (value).
        """
        # Initialise variables
        donors_list = []
        # Iterate over sheets
        sheets = list(filter(is_not_background, sheet_list))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                # Append index only
                index_donor = pedigree.index
                donors_list.append(index_donor)
        # Return simple dictionary
        return self.batchlist_to_dictionary(batch_lists=[donors_list])

    @staticmethod
    @dictify
    def batchlist_to_dictionary(batch_lists):
        """Turns a batch list of donors into a dictionary.

        :param batch_lists: List of lists with donors associated with each batch.
        :type batch_lists: list

        :return: Returns a dictionary with batch number/index (key) and list of donors (value).
        """
        for i, donors in enumerate(batch_lists, 1):  # 1-based index
            yield i, donors

    @staticmethod
    @listify
    def split_cohort_into_chunks(all_donors_list, max_batch_size):
        """Splits cohorts into chunks regardless of pedigree structure.

        :param all_donors_list: List with all donors in cohort.
        :type all_donors_list: list

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :return: Returns subset of donors (i.e., chunks) with maximum number of individuals.
        """
        # looping till length l
        for i in range(0, len(all_donors_list), max_batch_size):
            yield all_donors_list[i : i + max_batch_size]  # noqa: E203

    def split_cohort_into_balanced_chunks(self, all_pedigrees, max_batch_size):
        """Splits cohorts into chunks taking pedigree structure into account.

        :param all_pedigrees: List with all donors in cohort.
        :type all_pedigrees: list

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :return: Returns subset of donors (i.e., chunks) with maximum number of individuals.

        :raises: BatchInsufficientSpaceException: when none of the batches has enough space for
        a pedigree.
        """
        # Initialise variables
        cohort_size = 0
        pedigree_count_per_index_dict = defaultdict(list)  # key: size of family; value: donors list

        # Iterate over donors to populate size-by-donor dict
        for pedigree in all_pedigrees:
            size = len(pedigree.donors)
            index = pedigree.index
            pedigree_count_per_index_dict[size].append(index)
            # Update counter
            cohort_size += size

        # Max number of batches
        number_of_batches = ceil(cohort_size / max_batch_size)
        # Define batch dictionaries
        batch_donors_dict = {}  # key: tmp batch index; value: list of donors
        batch_counter_dict = {}  # key: tmp batch index; value: number of samples
        for i in range(number_of_batches):
            batch_donors_dict[i] = []
            batch_counter_dict[i] = 0

        # Iterate over pedigrees sizes and populate batch-donor dict
        # -> from larger sized families to smaller
        i_batch = 0
        for size_key in sorted(pedigree_count_per_index_dict.keys(), reverse=True):
            # Iterate over donors
            for donor in pedigree_count_per_index_dict.get(size_key):
                appropriate_batch_found = False
                attempts = 1
                while not appropriate_batch_found:
                    # Raise exception if maxed out
                    if attempts > number_of_batches:
                        self.call_insufficient_space_exception(
                            max_batch_size=max_batch_size,
                            i_size=size_key,
                            counter_dict=batch_counter_dict,
                        )

                    # Check if iterate all batches once
                    if i_batch not in batch_donors_dict:
                        i_batch = 0

                    # Test if current batch can receive pedigree
                    appropriate_batch_found = self.is_appropriate_batch_size(
                        batch_count=batch_counter_dict.get(i_batch),
                        max_batch_size=max_batch_size,
                        pedigree_count=size_key,
                    )
                    # Push donor to batch and update counter
                    if appropriate_batch_found:
                        curr_list = batch_donors_dict.get(i_batch)
                        curr_list.append(donor)
                        batch_donors_dict[i_batch] = curr_list
                        batch_counter_dict[i_batch] = batch_counter_dict.get(i_batch) + size_key

                    # Update counters
                    i_batch += 1
                    attempts += 1

        # Return list of batches (list of lists)
        return list(batch_donors_dict.values())

    @staticmethod
    def is_appropriate_batch_size(batch_count, max_batch_size, pedigree_count):
        """Checks current batch size exceeds maximum number of elements, i.e.,
        if it is appropriate for the number of elements in pedigree.

        :param batch_count: Number of samples in batch.
        :type batch_count: int

        :param max_batch_size: Maximum allowed batch size.
        :type max_batch_size: int

        :param pedigree_count: Number of elements in pedigree.
        :type pedigree_count: int

        :return: Returns True if number of donors (index and parents) in batch do not exceed the
        maximum set batch size; otherwise, False.
        """
        count = batch_count + pedigree_count
        if count > max_batch_size:
            return False
        else:
            return True

    @staticmethod
    def call_insufficient_space_exception(max_batch_size, i_size, counter_dict):
        """Creates error and raises batch insufficient space exception.

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :param i_size: Current pedigree size, i.e., number of samples.
        :type i_size:  int

        :param counter_dict: Dictionary with number of counters. Key: batch index; Value: number of
        samples in batch (int).
        :type counter_dict: dict
        """
        # Initialise variable
        err_msg = (
            "Iterated over all batches and none has enough space. "
            "Max batch size: {m}; current count: {c}.\n"
            "Batch structure:\n{s}"
        )
        # Dictionary structure to string
        dict_structure_str = ""
        for k, v in counter_dict.items():
            tmp_str = "batch {k} count: {v}\n".format(k=k, v=v)
            dict_structure_str += tmp_str
        err_msg = err_msg.format(m=max_batch_size, c=i_size, s=dict_structure_str)
        raise BatchInsufficientSpaceException(err_msg)
