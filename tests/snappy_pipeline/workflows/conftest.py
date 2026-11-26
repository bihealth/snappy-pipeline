# -*- coding: utf-8 -*-
"""Shared fixtures for the workflows unit tests"""

import os
from unittest.mock import patch

import pytest


@pytest.fixture(autouse=True)
def mock_settings_env_vars():
    """For tests, we want to have fake use medium partition."""
    with patch.dict(os.environ, {"SNAPPY_PIPELINE_PARTITION": "medium"}):
        yield


@pytest.fixture
def germline_sheet_tsv():
    """Two trios, the families named after their index donors"""
    return (
        "donor_name\tfamily_id\tfather\tmother\tsex\tis_affected\tis_index\tsample_name\tlibrary_name\tlibrary_type\tlibrary_kit\tfolder_name\n"
        "P001\tP001\tP002\tP003\tfemale\ttrue\ttrue\tP001-N1\tP001-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP001\n"
        "P002\tP001\t\t\tmale\tfalse\t\tP002-N1\tP002-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP002\n"
        "P003\tP001\t\t\tfemale\tfalse\t\tP003-N1\tP003-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP003\n"
        "P004\tP004\tP005\tP006\tmale\ttrue\ttrue\tP004-N1\tP004-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP004\n"
        "P005\tP004\t\t\tmale\tfalse\t\tP005-N1\tP005-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP005\n"
        "P006\tP004\t\t\tfemale\ttrue\t\tP006-N1\tP006-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP006\n"
    )


@pytest.fixture
def germline_trio_plus_sheet_tsv():
    """Two families with ids; P007 is unrelated by the parent columns but in family2, P001's folder is P011"""
    return (
        "donor_name\tfamily_id\tfather\tmother\tsex\tis_affected\tis_index\tsample_name\tlibrary_name\tlibrary_type\tlibrary_kit\tfolder_name\n"
        "P001\tfamily1\tP002\tP003\tfemale\ttrue\ttrue\tP001-N1\tP001-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP011\n"
        "P002\tfamily1\t\t\tmale\tfalse\t\tP002-N1\tP002-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP002\n"
        "P003\tfamily1\t\t\tfemale\tfalse\t\tP003-N1\tP003-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP003\n"
        "P004\tfamily2\tP005\tP006\tmale\ttrue\ttrue\tP004-N1\tP004-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP004\n"
        "P005\tfamily2\t\t\tmale\tfalse\t\tP005-N1\tP005-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP005\n"
        "P006\tfamily2\t\t\tfemale\ttrue\t\tP006-N1\tP006-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP006\n"
        "P007\tfamily2\t\t\tfemale\ttrue\t\tP007-N1\tP007-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP007\n"
    )


@pytest.fixture
def cancer_sheet_tsv():
    """Two cancer donors with tumors, a relapse and RNA libraries"""
    return (
        "donor_name\tsample_name\tis_tumor\tlibrary_name\tlibrary_type\tlibrary_kit\tfolder_name\tis_primary\n"
        "P001\tP001-N1\tfalse\tP001-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP001_N1_DNA1_WGS1\t\n"
        "P001\tP001-T1\ttrue\tP001-T1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP001_T1_DNA1_WGS1\t\n"
        "P001\tP001-T1\ttrue\tP001-T1-RNA1-mRNA_seq1\tmRNA_seq\t\tP001_T1_RNA1_mRNA_seq1\t\n"
        "P002\tP002-N1\tfalse\tP002-N1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP002_N1_DNA1_WGS1\t\n"
        "P002\tP002-T1\ttrue\tP002-T1-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP002_T1_DNA1_WGS1\ttrue\n"
        "P002\tP002-T2\ttrue\tP002-T2-DNA1-WGS1\tWGS\tAgilent SureSelect Human All Exon V6\tP002_T2_DNA1_WGS1\tfalse\n"
        "P002\tP002-T2\ttrue\tP002-T2-RNA1-mRNA_seq1\tmRNA_seq\t\tP002_T2-RNA1_mRNA_seq1\t\n"
    )
