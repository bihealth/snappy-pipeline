# -*- coding: utf-8 -*-
"""Shared fixtures for the workflows unit tests"""

import os
import textwrap
from unittest.mock import patch

import pytest


@pytest.fixture(autouse=True)
def mock_settings_env_vars():
    """For tests, we want to have fake use medium partition."""
    with patch.dict(os.environ, {"SNAPPY_PIPELINE_PARTITION": "medium"}):
        yield


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
        [Custom Fields]
        key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
        isTumor\tbioSample\tnormal/tumor\tstring\t.\t.\t.\t.\t.
        libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.
        extractionType\ttestSample\textraction type\tstring\t.\t.\t.\t.\t.

        [Data]
        patientName\tsampleName\tisTumor\tlibraryType\tfolderName\tlibraryKit\textractionType
        P001\tN1\tN\tWGS\tP001_N1_DNA1_WGS1\tAgilent SureSelect Human All Exon V6\tDNA
        P001\tT1\tY\tWGS\tP001_T1_DNA1_WGS1\tAgilent SureSelect Human All Exon V6\tDNA
        P001\tT1\tY\tmRNA_seq\tP001_T1_RNA1_mRNA_seq1\tNone\tRNA
        P002\tN1\tN\tWGS\tP002_N1_DNA1_WGS1\tAgilent SureSelect Human All Exon V6\tDNA
        P002\tT1\tY\tWGS\tP002_T1_DNA1_WGS1\tAgilent SureSelect Human All Exon V6\tDNA
        P002\tT2\tY\tWGS\tP002_T2_DNA1_WGS1\tAgilent SureSelect Human All Exon V6\tDNA
        P002\tT2\tY\tmRNA_seq\tP002_T2-RNA1_mRNA_seq1\tNone\tRNA
        # P003\tT1\tY\tWGS\tP003_T1_DNA1_WGS1\tNone\tDNA
        """
    ).lstrip()
