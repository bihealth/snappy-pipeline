# -*- coding: utf-8 -*-
"""Unit tests for the library DataFrame, relationship resolution, and entity helpers."""

import io
import textwrap
from types import SimpleNamespace

import pandas as pd
import pytest
from biomedsheets.io_tsv import read_germline_tsv_sheet
from biomedsheets.shortcuts import GermlineCaseSheet
from biomedsheets.shortcuts.cancer import CancerCaseSheet, CancerCaseSheetOptions
from biomedsheets.io_tsv import read_cancer_tsv_sheet

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.workflows.abstract import (
    apply_library_selection,
    build_library_dataframe,
    cohort_members,
    output_entity_names,
    resolve_relationships,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_data_set_info(sheet_type="germline_variants", is_background=False):
    return SimpleNamespace(sheet_type=sheet_type, is_background=is_background)


def _parse_cancer_sheet(tsv_str):
    sheet = read_cancer_tsv_sheet(io.StringIO(tsv_str))
    return CancerCaseSheet(
        sheet,
        options=CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True),
    )


def _parse_germline_sheet(tsv_str):
    sheet = read_germline_tsv_sheet(io.StringIO(tsv_str))
    return GermlineCaseSheet(sheet=sheet)


def _get_lib_names(df):
    """Return set of library names from a DataFrame."""
    return set(df["library_name"].tolist())


def _get_donor_libs(df, donor_name):
    """Return library names for a given donor."""
    return set(df[df["donor_name"] == donor_name]["library_name"].tolist())


# ---------------------------------------------------------------------------
# A. build_library_dataframe
# ---------------------------------------------------------------------------


class TestBuildLibraryDataframe:
    def test_cancer_sheet_correct_columns(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])

        expected_cols = {
            "library_name",
            "extraction_type",
            "kind",
            "role",
            "is_primary",
            "sex",
            "donor_name",
            "sample_name",
            "cohort_name",
            "father_name",
            "mother_name",
            "disease_state",
            "tissue_type",
        }
        extra_sheet_cols = {
            "folder_name",
            "library_kit",
            "library_type",
            "ncbi_taxon",
            "seq_platform",
        }
        assert expected_cols | extra_sheet_cols == set(df.columns)

    def test_cancer_sheet_roles_and_tissue_type(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])

        tumor_rows = df[df["role"] == "tumor"]
        assert (tumor_rows["tissue_type"] == "tumor").all()

        normal_rows = df[df["role"] == "normal"]
        assert (normal_rows["tissue_type"] == "normal").all()

        assert (df["kind"] == "cancer").all()

    def test_cancer_sheet_primary_flags(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])

        normal_rows = df[df["role"] == "normal"]
        assert normal_rows["is_primary"].all()

        for donor in df["donor_name"].unique():
            donor_tumor = df[(df["donor_name"] == donor) & (df["role"] == "tumor")]
            assert donor_tumor["is_primary"].any(), f"Donor {donor} has no primary tumor"

    def test_cancer_donor_count(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])
        # Two donors: P001 and P002 (P003 is commented out)
        donors = df["donor_name"].unique()
        assert len(donors) == 2
        assert all("P001" in d or "P002" in d for d in donors)

    def test_germline_trio_sheet_roles(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])

        assert (df["kind"] == "germline").all()

        # Two trios → 2 index, 2 father, 2 mother
        roles = df["role"].value_counts().to_dict()
        assert roles["index"] == 2
        assert roles["father"] == 2
        assert roles["mother"] == 2

    def test_germline_trio_pedigree_fields(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])

        # Index P001 should have father and mother library names
        idx_p001 = df[(df["role"] == "index") & (df["donor_name"].str.startswith("P001"))]
        assert len(idx_p001) == 1
        assert idx_p001.iloc[0]["father_name"] != "0"
        assert idx_p001.iloc[0]["mother_name"] != "0"

    def test_germline_sex_and_disease(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])

        # P001 is female and affected
        p001 = df[df["donor_name"].str.startswith("P001")]
        assert (p001["sex"] == "female").all()
        assert (p001["disease_state"] == 2).all()

        # P002 is male and unaffected
        p002 = df[df["donor_name"].str.startswith("P002")]
        assert (p002["sex"] == "male").all()
        assert (p002["disease_state"] == 1).all()

    def test_germline_tissue_type_is_unknown(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])
        assert (df["tissue_type"] == "unknown").all()

    def test_germline_cohort_grouping(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])

        # Two cohorts (trios), each with 3 members
        cohorts = df.groupby("cohort_name")
        assert len(cohorts) == 2
        for name, group in cohorts:
            assert len(group) == 3

    def test_mixed_sheets(self, cancer_sheet_tsv, germline_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info_c = _make_data_set_info("matched_cancer")
        info_g = _make_data_set_info("germline_variants")
        df = build_library_dataframe(
            [info_c, info_g],
            [csheet.sheet, gsheet.sheet],
            [csheet, gsheet],
        )
        assert set(df["kind"].unique()) == {"cancer", "germline"}

    def test_empty_data_returns_empty_dataframe(self):
        df = build_library_dataframe([], [], [])
        assert df.empty
        assert "library_name" in df.columns
        assert "tissue_type" in df.columns

    def test_background_dataset_skipped(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer", is_background=True)
        df = build_library_dataframe([info], [csheet.sheet], [csheet])
        assert df.empty

    def test_cancer_rna_libraries_present(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])
        rna = df[df["extraction_type"] == "rna"]
        assert len(rna) > 0

    def test_extra_sheet_columns_values(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])

        tumor_dna = df[df["library_name"].str.startswith("P001-T1-DNA1")].iloc[0]
        assert tumor_dna["folder_name"] == "P001_T1_DNA1_WGS1"
        assert tumor_dna["library_kit"] == "Agilent SureSelect Human All Exon V6"
        assert tumor_dna["library_type"] == "WGS"

    def test_sheet_keys_of_standard_columns_are_not_duplicated(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])

        # sex, isAffected, fatherName, motherName and extractionType only feed standard columns.
        assert {"is_affected", "father_pk", "mother_pk", "is_tumor"}.isdisjoint(df.columns)
        assert (df["library_type"] == "WGS").all()

    def test_library_selection_on_extra_sheet_column(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        df = build_library_dataframe([info], [csheet.sheet], [csheet])

        selected = apply_library_selection(df, "library_type == 'mRNA_seq'", "cancer")
        assert len(selected) > 0
        assert (selected["extraction_type"] == "rna").all()

    @pytest.mark.parametrize(
        "custom_fields, values, error",
        [
            # Maps to the standard column ``donor_name``.
            ([("donorName", "ngsLibrary")], ["X"], "standard library column"),
            # ``batchNo`` and ``batch_no`` both map to ``batch_no``.
            ([("batchNo", "bioEntity"), ("batch_no", "ngsLibrary")], ["B1", "B2"], "both map"),
        ],
    )
    def test_extra_sheet_column_collisions_raise(self, custom_fields, values, error):
        with pytest.raises(ValueError, match=error):
            _build_germline_singleton_dataframe(custom_fields, values)

    def test_same_extra_value_at_two_levels_is_kept_once(self):
        df = _build_germline_singleton_dataframe(
            [("batchNo", "bioEntity"), ("batch_no", "ngsLibrary")], ["B1", "B1"]
        )
        assert df.iloc[0]["batch_no"] == "B1"


def _build_germline_singleton_dataframe(custom_fields, values):
    """Build the library dataframe of a one-donor germline sheet with extra custom fields."""
    field_lines = "".join(
        f"{key}\t{entity}\t.\tstring\t.\t.\t.\t.\t.\n" for key, entity in custom_fields
    )
    header = "\t".join(key for key, _ in custom_fields)
    row = "\t".join(values)
    tsv = (
        "[Custom Fields]\n"
        "key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern\n"
        f"{field_lines}"
        "\n"
        "[Data]\n"
        f"patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tfolderName\thpoTerms\t{header}\n"
        f"P001\t.\t.\tF\tY\tWGS\tP001\t.\t{row}\n"
    )
    gsheet = _parse_germline_sheet(tsv)
    return build_library_dataframe(
        [_make_data_set_info("germline_variants")], [gsheet.sheet], [gsheet]
    )


# ---------------------------------------------------------------------------
# B. apply_library_selection
# ---------------------------------------------------------------------------


class TestApplyLibrarySelection:
    def _cancer_df(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        return build_library_dataframe([info], [csheet.sheet], [csheet])

    def test_cancer_default(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        result = apply_library_selection(df, None, "cancer")
        assert (result["role"] == "tumor").all()
        assert (result["extraction_type"] == "dna").all()

    def test_germline_default(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        df = build_library_dataframe([info], [gsheet.sheet], [gsheet])
        result = apply_library_selection(df, None, "germline")
        assert (result["extraction_type"] == "dna").all()

    def test_custom_expression_normals_only(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        result = apply_library_selection(df, "role == 'normal'", "cancer")
        assert (result["role"] == "normal").all()

    def test_custom_expression_rna(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        result = apply_library_selection(df, "extraction_type == 'rna'", "cancer")
        assert (result["extraction_type"] == "rna").all()

    def test_invalid_expression_raises(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        with pytest.raises(ValueError, match="invalid"):
            apply_library_selection(df, "this is not valid sql !!!", "cancer")

    def test_empty_dataframe(self):
        empty = pd.DataFrame(
            columns=[
                "library_name",
                "extraction_type",
                "kind",
                "role",
                "is_primary",
                "sex",
                "donor_name",
                "sample_name",
                "cohort_name",
                "father_name",
                "mother_name",
                "disease_state",
                "tissue_type",
            ]
        )
        result = apply_library_selection(empty, "role == 'tumor'", "cancer")
        assert result.empty

    def test_tissue_type_in_query(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        result = apply_library_selection(df, "tissue_type == 'tumor'", "cancer")
        assert (result["tissue_type"] == "tumor").all()


# ---------------------------------------------------------------------------
# C. resolve_relationships
# ---------------------------------------------------------------------------


class TestResolveRelationships:
    def _cancer_df(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        return build_library_dataframe([info], [csheet.sheet], [csheet])

    def test_matched_normal_lib(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
            )
        }
        result = resolve_relationships(df, rels)
        assert "matched_normal_lib" in result.columns

        # Each tumor DNA row should have a matched normal
        tumors = result[(result["role"] == "tumor") & (result["extraction_type"] == "dna")]
        assert (tumors["matched_normal_lib"] != "").all()

    def test_tumor_only_no_match(self):
        data = [
            {
                "library_name": "T1",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "tumor",
                "is_primary": True,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "T1",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "tumor",
            }
        ]
        df = pd.DataFrame(data)
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
            )
        }
        result = resolve_relationships(df, rels)
        assert result.iloc[0]["matched_normal_lib"] == ""

    def test_many_relationship(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "all_tumor_libs": RelationshipDefinition(
                via="donor_name",
                target="role == 'tumor'",
                many=True,
            )
        }
        result = resolve_relationships(df, rels)
        assert "all_tumor_libs" in result.columns
        for _, row in result.iterrows():
            assert isinstance(row["all_tumor_libs"], list)

    def test_invalid_via_column_raises(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "bad": RelationshipDefinition(
                via="nonexistent_column",
                target="role == 'tumor'",
            )
        }
        with pytest.raises(ValueError, match="not found"):
            resolve_relationships(df, rels)

    def test_invalid_target_expression_raises(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "bad": RelationshipDefinition(
                via="donor_name",
                target="this is not valid sql !!!",
            )
        }
        with pytest.raises(ValueError, match="failed"):
            resolve_relationships(df, rels)

    def test_multiple_relationships_sequential(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
            ),
            "donor_tumor_libs": RelationshipDefinition(
                via="donor_name",
                target="role == 'tumor'",
                many=True,
            ),
        }
        result = resolve_relationships(df, rels)
        assert "matched_normal_lib" in result.columns
        assert "donor_tumor_libs" in result.columns

    def test_relationship_column_available_for_selection(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
            )
        }
        result = resolve_relationships(df, rels)
        has_normal = result[result["matched_normal_lib"] != ""]
        assert len(has_normal) > 0

    def test_custom_column_name(self, cancer_sheet_tsv):
        df = self._cancer_df(cancer_sheet_tsv)
        rels = {
            "my_rel": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
                column="custom_col",
            )
        }
        result = resolve_relationships(df, rels)
        assert "custom_col" in result.columns
        assert "my_rel" not in result.columns

    def test_singular_relationship_errors_on_multiple_matches(self):
        """When many=False and multiple related rows match, raise ValueError."""
        data = [
            {
                "library_name": "T1",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "tumor",
                "is_primary": True,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "T1",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "tumor",
            },
            {
                "library_name": "N1",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "normal",
                "is_primary": False,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "N1",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "normal",
            },
            {
                "library_name": "N2",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "normal",
                "is_primary": False,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "N2",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "normal",
            },
        ]
        df = pd.DataFrame(data)
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
                many=False,
            )
        }
        with pytest.raises(ValueError, match="expected single match.*got 2"):
            resolve_relationships(df, rels)

    def test_many_relationship_accepts_multiple_matches(self):
        """When many=True, multiple matches are returned as a list."""
        data = [
            {
                "library_name": "T1",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "tumor",
                "is_primary": True,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "T1",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "tumor",
            },
            {
                "library_name": "N1",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "normal",
                "is_primary": False,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "N1",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "normal",
            },
            {
                "library_name": "N2",
                "extraction_type": "dna",
                "kind": "cancer",
                "role": "normal",
                "is_primary": False,
                "sex": "unknown",
                "donor_name": "P1",
                "sample_name": "N2",
                "cohort_name": "P1",
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "normal",
            },
        ]
        df = pd.DataFrame(data)
        rels = {
            "matched_normal_libs": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
                many=True,
            )
        }
        result = resolve_relationships(df, rels)
        tumor_row = result[result["library_name"] == "T1"].iloc[0]
        assert isinstance(tumor_row["matched_normal_libs"], list)
        assert set(tumor_row["matched_normal_libs"]) == {"N1", "N2"}

    def test_empty_dataframe(self):
        df = pd.DataFrame(
            columns=[
                "library_name",
                "extraction_type",
                "kind",
                "role",
                "is_primary",
                "sex",
                "donor_name",
                "sample_name",
                "cohort_name",
                "father_name",
                "mother_name",
                "disease_state",
                "tissue_type",
            ]
        )
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal'",
            )
        }
        result = resolve_relationships(df, rels)
        assert result.empty


# ---------------------------------------------------------------------------
# D. output_entity_names / cohort_members
# ---------------------------------------------------------------------------


class TestOutputEntityNames:
    def test_germline_trio_cohort_grouping(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        entities = output_entity_names(
            [info], [gsheet.sheet], [gsheet], selection=None, group_by="cohort"
        )
        # Two trios → two entities
        assert len(entities) == 2

    def test_germline_no_grouping(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        entities = output_entity_names(
            [info], [gsheet.sheet], [gsheet], selection=None, group_by=None
        )
        # 6 donors × 1 library each = 6 entities
        assert len(entities) == 6

    def test_cancer_with_library_selection(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        entities = output_entity_names(
            [info],
            [csheet.sheet],
            [csheet],
            selection="role == 'tumor' and extraction_type == 'dna'",
        )
        assert len(entities) > 0

    def test_empty_data(self):
        entities = output_entity_names([], [], [], selection=None)
        assert entities == []

    def test_cohort_members_dict(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        members = cohort_members([info], [gsheet.sheet], [gsheet], selection=None)
        assert isinstance(members, dict)
        assert len(members) == 2  # Two trios
        for cohort, libs in members.items():
            assert len(libs) == 3  # Each trio has 3 members

    def test_cohort_members_with_selection(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")
        members = cohort_members(
            [info], [gsheet.sheet], [gsheet], selection="extraction_type == 'dna'"
        )
        for cohort, libs in members.items():
            assert len(libs) == 3

    def test_entity_names_with_relationship(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
            )
        }
        entities = output_entity_names(
            [info],
            [csheet.sheet],
            [csheet],
            selection="role == 'tumor' and extraction_type == 'dna'",
            relationships=rels,
        )
        assert len(entities) > 0


# ---------------------------------------------------------------------------
# E. Model validation
# ---------------------------------------------------------------------------


class TestModelValidation:
    def test_relationship_definition_valid(self):
        rd = RelationshipDefinition(
            via="donor_name",
            target="role == 'normal'",
            column="my_col",
            many=True,
        )
        assert rd.via == "donor_name"
        assert rd.target == "role == 'normal'"
        assert rd.column == "my_col"
        assert rd.many is True

    def test_relationship_definition_defaults(self):
        rd = RelationshipDefinition(via="x", target="y == 1")
        assert rd.column is None
        assert rd.many is False

    def test_relationship_definition_missing_required(self):
        with pytest.raises(Exception):
            RelationshipDefinition()

    def test_snappy_step_model_has_fields(self):
        from snappy_pipeline.workflows.somatic_msi_calling.model import SomaticMsiCalling

        fields = list(SomaticMsiCalling.model_fields.keys())
        assert "library_selection" in fields
        assert "group_by" in fields
        assert "relationships" in fields

    def test_config_with_relationships_roundtrips(self):
        from snappy_pipeline.workflows.variant_calling.model import VariantCalling

        config = VariantCalling(
            depends_on={"alignments": "mapping", "reference": "genome", "dbsnp": "dbsnp"},
            tool="gatk4_hc_gvcf",
            gatk4_hc_gvcf={},
            relationships={
                "matched_normal_lib": {
                    "via": "donor_name",
                    "target": "role == 'normal' and extraction_type == 'dna'",
                }
            },
        )
        dumped = config.model_dump()
        assert dumped["relationships"]["matched_normal_lib"]["via"] == "donor_name"
        assert "role == 'normal'" in dumped["relationships"]["matched_normal_lib"]["target"]

    def test_config_without_relationships(self):
        from snappy_pipeline.workflows.variant_calling.model import VariantCalling

        config = VariantCalling(
            depends_on={"alignments": "mapping", "reference": "genome", "dbsnp": "dbsnp"},
            tool="gatk4_hc_gvcf",
            gatk4_hc_gvcf={},
        )
        assert config.relationships is None
        assert config.library_selection is None
        assert config.group_by is None


# ---------------------------------------------------------------------------
# F. Integration: relationship + selection + entity roundtrip
# ---------------------------------------------------------------------------


class TestIntegrationRoundtrip:
    def test_cancer_tumor_with_matched_normal(self, cancer_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        info = _make_data_set_info("matched_cancer")
        rels = {
            "matched_normal_lib": RelationshipDefinition(
                via="donor_name",
                target="role == 'normal' and extraction_type == 'dna'",
            )
        }

        df = build_library_dataframe([info], [csheet.sheet], [csheet], relationships=rels)
        assert "matched_normal_lib" in df.columns

        selected = apply_library_selection(
            df, "role == 'tumor' and extraction_type == 'dna'", "cancer"
        )
        assert (selected["role"] == "tumor").all()
        assert (selected["extraction_type"] == "dna").all()
        assert (selected["matched_normal_lib"] != "").all()

    def test_germline_trio_entity_roundtrip(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")

        entities = output_entity_names(
            [info], [gsheet.sheet], [gsheet], selection=None, group_by="cohort"
        )
        members = cohort_members([info], [gsheet.sheet], [gsheet], selection=None)

        for entity in entities:
            found = False
            for cohort, libs in members.items():
                if entity in libs:
                    found = True
                    break
            assert found, f"Entity {entity} not found in any cohort members"


# ---------------------------------------------------------------------------
# G. Parametrized scenario tests
# ---------------------------------------------------------------------------


# -- Helper: build a minimal cancer DataFrame directly ---------------------


def _cancer_df_from_libs(libs):
    """Build a cancer DataFrame from a list of (name, role, ext, donor) tuples."""
    rows = []
    for name, role, ext, donor in libs:
        rows.append(
            {
                "library_name": name,
                "extraction_type": ext,
                "kind": "cancer",
                "role": role,
                "is_primary": True,
                "sex": "unknown",
                "donor_name": donor,
                "sample_name": name,
                "cohort_name": donor,
                "father_name": "0",
                "mother_name": "0",
                "disease_state": 0,
                "tissue_type": "tumor" if role == "tumor" else "normal",
            }
        )
    return pd.DataFrame(rows)


MATCHED_PAIRED_RELS = {
    "matched_normal_lib": RelationshipDefinition(
        via="donor_name",
        target="role == 'normal' and extraction_type == 'dna'",
    )
}


class TestScenarioSomaticTumorNormalPaired:
    """Scenario 1: Somatic tumor-normal paired (two donors, each with T+N)."""

    def test_relationship_resolves_normal(self):
        df = _cancer_df_from_libs(
            [
                ("T1", "tumor", "dna", "D1"),
                ("N1", "normal", "dna", "D1"),
                ("T2", "tumor", "dna", "D2"),
                ("N2", "normal", "dna", "D2"),
            ]
        )
        result = resolve_relationships(df, MATCHED_PAIRED_RELS)
        t1 = result[result["library_name"] == "T1"].iloc[0]
        assert t1["matched_normal_lib"] == "N1"
        t2 = result[result["library_name"] == "T2"].iloc[0]
        assert t2["matched_normal_lib"] == "N2"

    def test_selection_tumor_only(self):
        df = _cancer_df_from_libs(
            [
                ("T1", "tumor", "dna", "D1"),
                ("N1", "normal", "dna", "D1"),
            ]
        )
        result = apply_library_selection(
            df, "role == 'tumor' and extraction_type == 'dna'", "cancer"
        )
        assert len(result) == 1
        assert result.iloc[0]["library_name"] == "T1"

    def test_entity_count(self):
        df = _cancer_df_from_libs(
            [
                ("T1", "tumor", "dna", "D1"),
                ("N1", "normal", "dna", "D1"),
                ("T2", "tumor", "dna", "D2"),
                ("N2", "normal", "dna", "D2"),
            ]
        )
        result = apply_library_selection(
            df, "role == 'tumor' and extraction_type == 'dna'", "cancer"
        )
        assert len(result) == 2


class TestScenarioSomaticTumorOnly:
    """Scenario 2: Somatic tumor-only (no normals in data)."""

    def test_no_match_in_relationship(self):
        df = _cancer_df_from_libs(
            [
                ("T1", "tumor", "dna", "D1"),
                ("T2", "tumor", "dna", "D2"),
            ]
        )
        result = resolve_relationships(df, MATCHED_PAIRED_RELS)
        assert (result["matched_normal_lib"] == "").all()

    def test_selection(self):
        df = _cancer_df_from_libs(
            [
                ("T1", "tumor", "dna", "D1"),
            ]
        )
        result = apply_library_selection(
            df, "role == 'tumor' and extraction_type == 'dna'", "cancer"
        )
        assert len(result) == 1


class TestScenarioSomaticNormalOnly:
    """Scenario 3: Somatic normal-only (panel of normals)."""

    def test_selection_normals(self):
        df = _cancer_df_from_libs(
            [
                ("N1", "normal", "dna", "D1"),
                ("N2", "normal", "dna", "D2"),
                ("T1", "tumor", "dna", "D1"),
            ]
        )
        result = apply_library_selection(
            df, "role == 'normal' and extraction_type == 'dna'", "cancer"
        )
        assert len(result) == 2
        assert (result["role"] == "normal").all()


class TestScenarioGermlineTrio:
    """Scenario 4: Germline trio (three members in one cohort)."""

    def test_cohort_grouping(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")

        entities = output_entity_names(
            [info], [gsheet.sheet], [gsheet], selection=None, group_by="cohort"
        )
        members = cohort_members([info], [gsheet.sheet], [gsheet], selection=None)

        assert len(entities) == 2
        for cohort, libs in members.items():
            assert len(libs) == 3

    def test_entity_is_in_members(self, germline_sheet_tsv):
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info = _make_data_set_info("germline_variants")

        entities = output_entity_names(
            [info], [gsheet.sheet], [gsheet], selection=None, group_by="cohort"
        )
        members = cohort_members([info], [gsheet.sheet], [gsheet], selection=None)

        all_member_libs = set()
        for libs in members.values():
            all_member_libs.update(libs)
        for entity in entities:
            assert entity in all_member_libs


class TestScenarioGermlineDuo:
    """Scenario 5: Germline duo (two members)."""

    def test_two_members(self):
        from biomedsheets.io_tsv import read_germline_tsv_sheet as read_gtsv

        tsv = textwrap.dedent("""\
            [Custom Fields]
            key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
            libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.

            [Data]
            patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tlibraryKit\tfolderName\thpoTerms
            P001\tP002\t.\tM\tY\tWGS\tkit\tP001\t.
            P002\t.\t.\tM\tN\tWGS\tkit\tP002\t.
        """)
        import io as _io

        sheet = read_gtsv(_io.StringIO(tsv))
        gsheet = GermlineCaseSheet(sheet=sheet)
        info = _make_data_set_info("germline_variants")

        df = build_library_dataframe([info], [sheet], [gsheet])
        assert len(df) == 2

        members = cohort_members([info], [sheet], [gsheet], selection=None)
        assert len(members) == 1
        for cohort, libs in members.items():
            assert len(libs) == 2


class TestScenarioGermlineSingleton:
    """Scenario 6: Germline singleton (one sample, no parents)."""

    def test_single_member(self):
        from biomedsheets.io_tsv import read_germline_tsv_sheet as read_gtsv

        tsv = textwrap.dedent("""\
            [Custom Fields]
            key\tannotatedEntity\tdocs\ttype\tminimum\tmaximum\tunit\tchoices\tpattern
            libraryKit\tngsLibrary\tEnrichment kit\tstring\t.\t.\t.\t.\t.

            [Data]
            patientName\tfatherName\tmotherName\tsex\tisAffected\tlibraryType\tlibraryKit\tfolderName\thpoTerms
            P001\t.\t.\tM\tY\tWGS\tkit\tP001\t.
        """)
        import io as _io

        sheet = read_gtsv(_io.StringIO(tsv))
        gsheet = GermlineCaseSheet(sheet=sheet)
        info = _make_data_set_info("germline_variants")

        df = build_library_dataframe([info], [sheet], [gsheet])
        assert len(df) == 1
        assert df.iloc[0]["role"] == "index"

        entities = output_entity_names([info], [sheet], [gsheet], selection=None, group_by=None)
        assert len(entities) == 1


class TestScenarioMixedCancerAndGermline:
    """Scenario 7: Mixed cancer + germline in one project."""

    def test_both_kinds_present(self, cancer_sheet_tsv, germline_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info_c = _make_data_set_info("matched_cancer")
        info_g = _make_data_set_info("germline_variants")

        df = build_library_dataframe(
            [info_c, info_g], [csheet.sheet, gsheet.sheet], [csheet, gsheet]
        )
        assert set(df["kind"].unique()) == {"cancer", "germline"}

    def test_per_kind_defaults(self, cancer_sheet_tsv, germline_sheet_tsv):
        csheet = _parse_cancer_sheet(cancer_sheet_tsv)
        gsheet = _parse_germline_sheet(germline_sheet_tsv)
        info_c = _make_data_set_info("matched_cancer")
        info_g = _make_data_set_info("germline_variants")

        entities = output_entity_names(
            [info_c, info_g],
            [csheet.sheet, gsheet.sheet],
            [csheet, gsheet],
            selection=None,
        )
        # Should have entities from both cancer and germline
        assert len(entities) > 0
