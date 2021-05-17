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
from operator import attrgetter
import random
import warnings

from biomedsheets.shortcuts import is_background, is_not_background

from snappy_pipeline.utils import dictify, listify


class BatchInsufficientSpaceException(Exception):
    """Raised when none of the batches has enough space to accommodate pedigree"""

    def __init__(self, max_batch_size, size_key, counter_dict):
        """Constructor.

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :param size_key: Current pedigree size, i.e., number of samples.
        :type size_key:  int

        :param counter_dict: Dictionary with number of counters. Key: batch index; Value: number of
        samples in batch (int).
        :type counter_dict: dict
        """
        # Define custom message
        message = self.custom_message(
            max_batch_size=max_batch_size,
            size_key=size_key,
            counter_dict=counter_dict,
        )
        super().__init__(message)

    @staticmethod
    def custom_message(max_batch_size, size_key, counter_dict):
        """Creates custom error message.

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :param size_key: Current pedigree size, i.e., number of samples.
        :type size_key:  int

        :param counter_dict: Dictionary with number of counters. Key: batch index; Value: number of
        samples in batch (int).
        :type counter_dict: dict

        :return: Returns custom error messaged based on provided input.
        """
        # Initialise message
        message = (
            "Iterated over all batches and none has enough space. "
            "Max batch size: {m}; current count: {c}.\n"
            "Batch structure:\n{s}"
        )
        # Dictionary structure to string
        dict_structure_str = ""
        for k, v in counter_dict.items():
            tmp_str = "batch {k} count: {v}\n".format(k=k, v=v)
            dict_structure_str += tmp_str
        message = message.format(m=max_batch_size, c=size_key, s=dict_structure_str)
        # Return
        return message


class Chunk:
    """Class contains methods to split cohort into chunks."""

    #: List of implemented methods.
    active_methods = ["single", "evenly", "incremental"]

    #: Selected chunk method.
    selected_method = None

    #: List of Sample sheets.
    sheet_list = None

    #: Maximum chunk size - used in 'evenly' and 'incremental' methods.
    maximal_chunk_size = None

    #: List of custom fields in sample sheet - used in 'incremental' to order samples.
    order_by_custom_field = None

    #: Minimum allowed chunk - different for each tool.
    minimum_chunk_size = None

    def __init__(
        self,
        sheet_list,
        method,
        maximal_chunk_size=None,
        order_by_custom_field=None,
        minimum_chunk_size=50,
    ):
        """Constructor

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :param method: Chunk method. Available implementations: 'single', 'evenly', 'incremental'.
        :type method: str

        :param maximal_chunk_size: Maximum chunk size - used in 'evenly' method.
        :type maximal_chunk_size: int, optional

        :param order_by_custom_field: List of sample sheets custom fields to be used to order the
        samples in the 'incremental' method.
        :type order_by_custom_field: list, optional

        :param minimum_chunk_size: Minimum allowed chunk - different for each tool. Default: 50.
        :type minimum_chunk_size: int, optional
        """
        # Validate selected method
        if method not in self.active_methods:
            valid_methods = ", ".join(self.active_methods)
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

        # Validate min chunk size
        min_err_msg = (
            "Min number of samples per chunk must be an integer greater than zero. "
            "Invalid input: '{0}'.".format(minimum_chunk_size)
        )
        try:
            int(minimum_chunk_size)
            if minimum_chunk_size <= 1:
                raise ValueError(min_err_msg)
        except ValueError as e:
            raise type(e)(str(e) + min_err_msg)
        self.minimum_chunk_size = minimum_chunk_size

        # Validate max chunk size
        if maximal_chunk_size is not None:
            max_err_msg = (
                "Max number of samples per chunk must be an integer greater than zero. "
                "Invalid input: '{0}'.".format(maximal_chunk_size)
            )
            try:
                int(maximal_chunk_size)
                if maximal_chunk_size <= 1:
                    raise ValueError(max_err_msg)
            except ValueError as e:
                raise type(e)(str(e) + max_err_msg)
        self.maximal_chunk_size = maximal_chunk_size

        # Max chunk must be defined for these methods
        if method in ["evenly", "incremental"] and maximal_chunk_size is None:
            err_msg = (
                "For methods 'evenly' and 'incremental' argument `maximal_chunk_size` "
                "must be defined."
            )
            raise Exception(err_msg)

        # Set order custom field
        self.order_by_custom_field = order_by_custom_field

    def run(self):
        """Runs split in chunks selected method."""
        # Single Chunk
        if self.selected_method == "single":
            return self._single_chunk(sheet_list=self.sheet_list)
        # Evenly Chunk
        if self.selected_method == "evenly":
            return self._evenly_chunk(sheet_list=self.sheet_list, max_chunk=self.maximal_chunk_size)
        # Incremental Chunk
        if self.selected_method == "incremental":
            return self._incremental_chunk(
                sheet_list=self.sheet_list,
                max_chunk=self.maximal_chunk_size,
                fields=self.order_by_custom_field,
            )
        return None

    def _incremental_chunk(self, sheet_list, max_chunk, fields=None):
        """Define an evenly divided incremental chunk based on sample sheet lists.

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :param max_chunk: Number of samples per chunk, if not enough samples it will complement
        with background samples (when available).
        :type max_chunk: int

        :param fields: List of custom fields from donor to be used to order the pedigrees.
        :type fields: list, optional

        :return: Returns dictionary with several chunks indices (keys) associated with the list of
        index donors from the cohorts (and background samples if necessary) defined in the sample
        sheets (values).
        """
        # Initialise variables
        donors_list = []
        all_foreground_pedigrees = []
        all_background_pedigrees = []

        # Define foreground samples
        sheets = list(filter(is_not_background, sheet_list))
        for f_sheet in sheets:
            all_foreground_pedigrees.extend(f_sheet.cohort.pedigrees)

        # Define background samples
        background_sheets = list(filter(is_background, sheet_list))
        for b_sheet in background_sheets:
            all_background_pedigrees.extend(b_sheet.cohort.pedigrees)

        # Subset pedigrees by sequencing kit
        pedigrees_by_kit_lists = self.split_pedigrees_by_sequencing_kit(
            all_pedigrees=all_foreground_pedigrees
        )

        # Iterate over sheets
        for i_pedigrees in pedigrees_by_kit_lists:
            # Check if makes sense to split cohort
            cohort_size = self.get_cohort_size_from_pedigrees_list(all_pedigrees=i_pedigrees)
            if cohort_size < self.minimum_chunk_size:
                self.warn_user_cohort_smaller_than_min(
                    cohort_size=cohort_size, min_chunk=self.minimum_chunk_size
                )
                balanced_chunks_list = [self.pedigree_to_index_list(all_pedigrees=i_pedigrees)]
            else:
                balanced_chunks_list = self.split_cohort_into_incremental_chunks(
                    all_pedigrees=i_pedigrees,
                    all_background_pedigrees=all_background_pedigrees,
                    max_batch_size=max_chunk,
                    fields=fields,
                )
            for batch in balanced_chunks_list:
                donors_list.append(batch)
        # Return dictionary
        return self.batchlist_to_dictionary(batch_lists=donors_list)

    def _evenly_chunk(self, sheet_list, max_chunk):
        """Define an evenly divided chunks based on sample sheets list.

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :param max_chunk: Maximum number of samples per chunk.
        :type max_chunk: int

        :return: Returns dictionary with several chunks indices (keys) associated with the list of
        index donors from the cohorts defined in the sample sheets (values).
        """
        # Initialise variables
        donors_list = []
        all_foreground_pedigrees = []

        # Define foreground samples
        sheets = list(filter(is_not_background, sheet_list))
        for f_sheet in sheets:
            all_foreground_pedigrees.extend(f_sheet.cohort.pedigrees)

        # Subset pedigrees by sequencing kit
        pedigrees_by_kit_lists = self.split_pedigrees_by_sequencing_kit(
            all_pedigrees=all_foreground_pedigrees
        )

        # Iterate over sheets
        for i_pedigrees in pedigrees_by_kit_lists:
            # Check if makes sense to split cohort
            cohort_size = self.get_cohort_size_from_pedigrees_list(all_pedigrees=i_pedigrees)
            print(cohort_size)
            if cohort_size < self.minimum_chunk_size:
                self.warn_user_cohort_smaller_than_min(
                    cohort_size=cohort_size, min_chunk=self.minimum_chunk_size
                )
                balanced_chunks_list = [self.pedigree_to_index_list(all_pedigrees=i_pedigrees)]
                print(balanced_chunks_list)
            else:
                balanced_chunks_list = self.split_cohort_into_balanced_chunks(
                    all_pedigrees=i_pedigrees, max_batch_size=max_chunk
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
    @listify
    def pedigree_to_index_list(all_pedigrees):
        """Create list of indexes from pedigrees.

        :param all_pedigrees: List of all Pedigrees in cohort.
        :type all_pedigrees: list

        :return: Returns list of indexes from the inputted pedigrees.
        """
        for pedigree in all_pedigrees:
            yield pedigree.index

    @staticmethod
    @dictify
    def batchlist_to_dictionary(batch_lists):
        """Turns a batch list of donors into a dictionary.

        :param batch_lists: List of lists with donors associated with each batch.
        :type batch_lists: list

        :return: Returns a dictionary with batch number/index (key) and list of donors (value).
        Used to uninform the output of the different methods.
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
        for i in range(0, len(all_donors_list), max_batch_size):
            yield all_donors_list[i : i + max_batch_size]  # noqa: E203

    def split_cohort_into_incremental_chunks(
        self, all_pedigrees, all_background_pedigrees, max_batch_size, fields=None
    ):
        """Splits cohorts into incremental chunks taking pedigree order and structure into account.

        :param all_pedigrees: List with all donors in cohort.
        :type all_pedigrees: list

        :param all_background_pedigrees: List with all donors in background cohort.
        :type all_background_pedigrees: list

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :param fields: List of custom fields from donor to be used to order the pedigrees.
        :type fields: list, optional

        :return: Returns subset of donors (i.e., chunks) with maximum number of individuals and
        background samples if applicable.

        :raises: Exception: when final observed number of batches is smaller than expected.
        Expected: ceiling(cohort_size / max_batch_size). It can be larger as it expected that it
        won't always be possible to split the cohort evenly.
        """
        # Initialise variables
        cohort_size = 0
        batch_donors_dict = defaultdict(list)  # key: batch index; value: list of donors
        batch_counter_dict = defaultdict(lambda: 0)  # key: batch index; value: number of samples

        # Check if background samples available
        use_background = len(all_background_pedigrees) > 0
        n_background_samples = self.get_cohort_size_from_pedigrees_list(
            all_pedigrees=all_background_pedigrees
        )

        # Iterate over pedigrees sizes and populate batch-donor dict
        # -> use order of pedigrees list
        if fields is not None:
            all_pedigrees = self.order_pedigrees_by_field(
                all_pedigrees=all_pedigrees, fields=fields
            )
        i_batch = 0
        for pedigree in all_pedigrees:
            # Get donor info
            size = len(pedigree.donors)
            index = pedigree.index

            # Check if current batch has enough space
            current_batch_has_enough_space = self.is_appropriate_batch_size(
                batch_count=batch_counter_dict[i_batch],
                max_batch_size=max_batch_size,
                pedigree_count=size,
            )
            if current_batch_has_enough_space:
                # Append to current batch
                batch_donors_dict[i_batch].append(index)
            else:
                # Append to next batch
                i_batch += 1
                batch_donors_dict[i_batch] = [index]
            # Update size of batch dict and cohort size
            cohort_size += size
            batch_counter_dict[i_batch] += size

        # If applicable, use background samples to complete cohort
        if use_background:
            for i_batch in batch_counter_dict.keys():
                required_samples = max_batch_size - batch_counter_dict[i_batch]
                if required_samples >= n_background_samples:
                    batch_donors_dict[i_batch].extend(
                        [b_pedigree.index for b_pedigree in all_background_pedigrees]
                    )
                else:
                    random_background_sample = self.get_random_cohort_sample(
                        all_donors=all_background_pedigrees, n_selections=required_samples
                    )
                    batch_donors_dict[i_batch].extend(random_background_sample)

        # Sanity Check
        # It won't necessarily be possible to split the cohort exactly, hence
        # it should return at least the number of expected batches.
        expected_number_of_batches = ceil(cohort_size / max_batch_size)
        observed_number_of_batches = len(batch_donors_dict.keys())
        if expected_number_of_batches > observed_number_of_batches:
            err_msg = (
                "Observed number of batches ({onb}) is smaller than expected ({enb}). Expected "
                "logic: ceiling( {size} / {max_} )".format(
                    onb=observed_number_of_batches,
                    enb=expected_number_of_batches,
                    size=cohort_size,
                    max_=max_batch_size,
                )
            )
            raise Exception(err_msg)

        # Return list of batches (list of lists)
        return list(batch_donors_dict.values())

    def split_cohort_into_balanced_chunks(self, all_pedigrees, max_batch_size):
        """Splits cohorts into evenly divided chunks taking pedigree structure into account.

        :param all_pedigrees: List with all donors in cohort.
        :type all_pedigrees: list

        :param max_batch_size: Maximum batch size
        :type max_batch_size: int

        :return: Returns subset of donors (i.e., chunks) with maximum number of individuals.

        :raises: BatchInsufficientSpaceException: when none of the batches has enough space for
        a pedigree.
        """
        # Get cohort size and structure
        cohort_size, pedigree_count_per_index_dict = self.build_pedigree_count_per_index_dict(
            all_pedigrees
        )

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
                        raise BatchInsufficientSpaceException(
                            max_batch_size=max_batch_size,
                            size_key=size_key,
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

    def build_pedigree_count_per_index_dict(self, all_pedigrees):
        """Build dictionary with pedigree count per index donor.

        :param all_pedigrees: List with all donors in cohort.
        :type all_pedigrees: list

        :return: Returns tuple: number of samples in cohort; dictionary with pedigree count per
        index donor (key: size of family; value: donors list ordered by secondary id).
        """
        # Initialise variables
        cohort_size = 0
        pedigree_count_per_index_dict = defaultdict(list)  # key: size of family; value: donors list
        # Order pedigrees by secondary id
        ordered_pedigrees = self.order_pedigrees_by_sampleid(all_pedigrees=all_pedigrees)
        # Iterate over donors to populate size-by-donor dict
        for pedigree in ordered_pedigrees:
            size = len(pedigree.donors)
            index = pedigree.index
            pedigree_count_per_index_dict[size].append(index)
            # Update counter
            cohort_size += size
        # Return
        return cohort_size, pedigree_count_per_index_dict

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
    def warn_user_cohort_smaller_than_min(cohort_size, min_chunk):
        """Warns user that cohort is smaller than minimum and it won't be split.

        :param cohort_size: Number of samples in cohort.
        :type cohort_size: int

        :param min_chunk: Minimum allowed chunk.
        :type min_chunk: int
        """
        warn_msg = (
            "Cohort in sample sheet contains only {count_} (<{min_}) "
            "and won't be split into batches.".format(count_=str(cohort_size), min_=min_chunk)
        )
        warnings.warn(warn_msg)

    @staticmethod
    def get_cohort_size_from_sheet(sheet):
        """Get cohort size from sample sheet.

        :param sheet: Sample sheet.
        :type sheet: biomedsheets.shortcuts.GermlineCaseSheet

        :return: Returns number of samples in cohort represented by sample sheet.
        """
        # Initialise variable
        cohort_size = 0
        # Iterate over pedigrees
        for pedigree in sheet.cohort.pedigrees:
            cohort_size += len(pedigree.donors)
        # Return
        return cohort_size

    @staticmethod
    def get_cohort_size_from_pedigrees_list(all_pedigrees):
        """Get cohort size from pedigrees list.

        :param all_pedigrees: List of all pedigrees.
        :type all_pedigrees: list

        :return: Returns number of samples in pedigree list.
        """
        # Initialise variable
        cohort_size = 0
        # Iterate over pedigrees
        for d in all_pedigrees:
            cohort_size += len(d.donors)
        # Return
        return cohort_size

    def get_random_cohort_sample(self, all_donors, n_selections):
        """Random sample of donors.

        :param all_donors: List with all donors in cohort.
        :type all_donors: list

        :param n_selections: Number of selections from sample.
        :type n_selections: int

        :return: Returns list with random subset of donors without repetition. The number of
        selected samples is defined by the number of samples associated with a pedigree. Hence, the
        result might be non-intuitive. For example: `n_selection` equal 3 might return three
        different sets: one pedigree with three members (e.g., index, mother, father - output
        length is one); two pedigrees, one with just index and other with index and one parent
        (output length is two); or three pedigrees, one index for each (output length is three).

        TODO: Is it acceptable to return less donors than requested (current implementation)?
        """
        # Initialise variable
        random_donors = []
        # Validate input
        length_all_donors_list = self.get_cohort_size_from_pedigrees_list(all_pedigrees=all_donors)
        if n_selections > length_all_donors_list:
            err_msg = (
                "Cannot provide sample larger than there are donors. "
                "Amount of donors: {l}; number of requested donors: {n}.".format(
                    l=length_all_donors_list, n=n_selections
                )
            )
            raise ValueError(err_msg)
        # Create random list
        if n_selections == length_all_donors_list:
            random_donors = [donor.index for donor in all_donors]  # trivial case
        else:
            max_n_attempts = length_all_donors_list
            attempts = 0
            selection_random_donors = 0
            while selection_random_donors < n_selections:
                # Update counter and check if maxed out
                attempts += 1
                if attempts >= max_n_attempts:
                    break
                # Randomly select a donor from input list
                i_donor = random.choice(all_donors)  # nosec
                i_pedigree_size = len(i_donor.donors)
                i_index = i_donor.index
                tmp_selection_size = selection_random_donors + i_pedigree_size
                if i_index not in random_donors and tmp_selection_size <= n_selections:
                    random_donors.append(i_index)
                    selection_random_donors = tmp_selection_size
                    attempts = 0
        # Return
        return random_donors

    def order_pedigrees_by_field(self, all_pedigrees, fields):
        """Order pedigrees by fields in index.

        :param all_pedigrees: List with all pedigrees in cohort.
        :type all_pedigrees: list

        :param fields: List of custom fields from donor to be used to order the pedigrees.
        :type fields: list

        :return: Returns the pedigree list ordered by the fields. Ties are solved using the
        secondary_id (sample id).

        :raises: AttributeError: when Pedigree in provided list doesn't have nested attributes
        'index.extra_infos'.
        :raises: ValueError: when custom field selected for ordering pedigrees is None.
        """
        # Validate input
        if len(fields) == 0:
            raise ValueError("Argument 'fields' cannot be empty. Please review input.")
        # Return ordered pedigree list
        return sorted(
            all_pedigrees,
            key=lambda pedigree: tuple(
                self.get_max_field_value_from_donors(pedigree, field) for field in fields
            )
            + (pedigree.index.wrapped.secondary_id,),
        )

    @staticmethod
    def get_max_field_value_from_donors(pedigree, field):
        """
        :param pedigree: Pedigree object.
        :type pedigree: biomedsheets.shortcuts.germline.Pedigree

        :param field: Custom fields from donor.
        :type field: str

        :return: Max value of a field for a set of donors in a pedigree.
        """
        # Initialise variables
        values_list = []
        # Initialise exception messages
        attribute_error_msg = "\nExpects that Donor object have attribute 'extra_infos'."
        value_error_msg = "Values used for ordering cannot be None. Field '{f}'; donor: {p}"
        # Collect all values
        for donor in pedigree.donors:
            try:
                field_value = attrgetter("extra_infos")(donor).get(field)
                values_list.append(field_value)
                if field_value is None:
                    value_error_msg = value_error_msg.format(f=field, p=donor)
                    raise ValueError(value_error_msg)
            except AttributeError as ae:
                raise type(ae)(ae.message + attribute_error_msg)
        # Return max value
        return max(values_list)

    @staticmethod
    def order_pedigrees_by_sampleid(all_pedigrees):
        """Order pedigrees by sample id of index (secondary_id).

        :param all_pedigrees: List with all pedigrees in cohort.
        :type all_pedigrees: list

        :return: Returns the pedigree list ordered by sample identifier, i.e., `secondary_id` field
        from the `BioEntity` object.
        """
        return sorted(all_pedigrees, key=lambda pedigree: pedigree.index.wrapped.secondary_id)

    @staticmethod
    def split_pedigrees_by_sequencing_kit(all_pedigrees):
        """Split pedigrees based on index sequencing kit (dna_ngs_library).

        :param all_pedigrees: List with all pedigrees in cohort.
        :type all_pedigrees: list

        :return: Returns list of list of pedigrees segregated based on sequencing kit name. It
        assumes that all donors in a pedigree were sequenced with the same kit as the index.
        """
        # Initialise variable
        kit_pedigree_dict = defaultdict(list)  # key: sequencing kit name; value: list of pedigrees
        # Iterate over pedigrees
        for pedigree in all_pedigrees:
            library_kit = pedigree.index.dna_ngs_library.ngs_library.extra_infos.get("libraryKit")
            kit_pedigree_dict[library_kit].append(pedigree)
        # Return list of lists
        return list(kit_pedigree_dict.values())
