# -*- coding: utf-8 -*-
""" Implementation of the ``chunk`` procedure

The cohort being analyzed will be splitted into chunks depending on the selected method:
- ``single``: create one large chunk for all.
- ``evenly``: create evenly sized chunks of maximal size.
- ``incremental``: iterate samples and create a new chunk after the maximal
                   size has been reached; should be combined with background data
                   to guarantee a minimal number of samples.
"""

from biomedsheets.shortcuts import is_not_background

# List of implemented split methods
ACTIVE_METHODS = ["single", "evenly", "incremental"]


class Chunk:
    """Class contains methods to split cohort into chunks."""

    # Selected chunk method
    selected_method = None

    # List of Sample sheets
    sheet_list = None

    def __init__(self, method, sheet_list):
        """Constructor

        :param method: Chunk method. Available implementations: 'single', 'evenly', 'incremental'.
        :type method: str

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list
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

    def run(self):
        """Runs split in chunks.

        :return:
        """
        # Single Chunk
        if self.selected_method == "single":
            return self.single_chunk(sheet_list=self.sheet_list)

        return None

    @staticmethod
    def single_chunk(sheet_list):
        """Define a single chunk based on sample sheets list.

        :param sheet_list: List of Sample Sheets.
        :type sheet_list: list

        :return: Returns dictionary with a single key ('1') associated with the list of index
        donors from the cohorts defined in the sample sheets.
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
        return {1: donors_list}
