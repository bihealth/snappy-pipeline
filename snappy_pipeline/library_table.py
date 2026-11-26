"""The library table of a project: one row per NGS library of the sample sheets.

``Project.data_set_infos`` loads each data set's sheet once (:class:`DataSetInfo`);
:func:`build_library_dataframe` turns the sheets into the table that every step reads through
``BaseStep`` (``library_table``, ``build_library_dataframe``, ``libraries``), with the pedigree
columns of the biomedsheets germline view, :func:`build_donor_dataframe` gives the donors of the
germline sheets, sequenced or not, :func:`resolve_relationships` adds the relationship columns
(``matched_normal_lib``), and :func:`select_libraries` applies a task's ``library_selection``.
"""

import os
import re
from collections import OrderedDict
from typing import Any

from biomedsheets import io_tsv
from biomedsheets.io import SheetBuilder, json_loads_ordered
from biomedsheets.io_tsv.base import LIBRARY_TYPES, SheetIOException
from biomedsheets.naming import NAMING_SCHEMES, name_generator_for_scheme
from biomedsheets.ref_resolver import RefResolver

from snappy_pipeline.base import ConfigError
from snappy_pipeline.workflow_model import DataSetType


#: The biomedsheets readers of the TSV sheets by sheet type
_TSV_READERS = {
    DataSetType.MATCHED_CANCER: io_tsv.read_cancer_tsv_sheet,
    DataSetType.GERMLINE_VARIANTS: io_tsv.read_germline_tsv_sheet,
    DataSetType.GENERIC: io_tsv.read_generic_tsv_sheet,
}


class DataSetInfo:
    """A data set of the project with its loaded sample sheet; ``Project.data_set_infos`` loads
    each once per run."""

    def __init__(
        self,
        name,
        sheet_path,
        sheet_type,
        is_background,
        naming_scheme,
        pedigree_field=None,
    ):
        """Constructor.

        :param name: Name of the data set.
        :type name: str

        :param sheet_path: Path to sheet file that should be loaded.
        :type sheet_path: str

        :param sheet_type: Explicit sheet type (e.g. "matched_cancer"), if any.  Otherwise, will
        attempt to load from sheet.
        :type sheet_type: str

        :param is_background: Whether the data set info is to be used only for background.
        :type is_background: bool

        :param naming_scheme: Selected naming schema: either 'secondary_id_pk' or
        'only_secondary_id'.
        :type naming_scheme: str

        :param pedigree_field: Custom field from sample sheet used to define pedigree, e.g.,
        'familyId'. If none defined, it will set pedigree based on sample sheet 'row'.
        Default: None.
        :type pedigree_field: str
        """
        #: Name of the data set
        self.name = name
        #: Path to the sheet file, for loading
        self.sheet_path = sheet_path
        #: Explicit sheet type (e.g. "matched_cancer"), if any.  Otherwise, will attempt to load
        # from sheet.
        self.sheet_type = sheet_type
        #: Whether the data set info is to be used only for background
        self.is_background = is_background
        #: Selected naming schema
        if naming_scheme not in NAMING_SCHEMES:
            raise ValueError("Invalid naming scheme: {}".format(naming_scheme))  # pragma: no cover
        self.naming_scheme = naming_scheme
        #: The BioMed SampleSheet
        self.sheet = self._load_sheet()
        #: The (optional) custom field used to define pedigree
        self.pedigree_field_kwargs = None
        if pedigree_field:
            self.pedigree_field_kwargs = {"join_by_field": pedigree_field}

    def _load_sheet(self):
        # ``load_project()`` resolved the path against the project directory
        fname = self.sheet_path
        if not os.path.exists(fname):
            raise ConfigError(f"Could not find sample sheet file {fname}")
        if self.sheet_path.endswith(".tsv"):
            if self.sheet_type not in _TSV_READERS:
                raise ConfigError("Invalid sheet type {}".format(self.sheet_type))
            try:
                with open(fname, "rt") as f:
                    sheet = _TSV_READERS[self.sheet_type](f, self.sheet_path, self.naming_scheme)
            except SheetIOException as e:
                raise ConfigError(f"Sample sheet {fname}: {e}") from None
            # biomedsheets skips rows without libraryType silently (germline and generic sheets)
            if not any(
                test_sample.ngs_libraries
                for donor in sheet.bio_entities.values()
                for sample in donor.bio_samples.values()
                for test_sample in sample.test_samples.values()
            ):
                raise ConfigError(
                    f"Sample sheet {fname}: no row has a libraryType ({', '.join(LIBRARY_TYPES)}), "
                    "so there are no libraries. A row without libraryType is a person without "
                    "sequencing data, such as a parent who was not sequenced."
                )
            return sheet
        elif self.sheet_path.endswith(".json"):
            with open(fname, "rt") as f:
                sheet_json = json_loads_ordered(f.read())
            resolver = RefResolver(
                lookup_paths=os.path.dirname(os.path.abspath(self.sheet_path)),
                dict_class=OrderedDict,
            )
            resolved_json = resolver.resolve("file://" + self.sheet_path, sheet_json)
            return SheetBuilder(resolved_json).run(
                name_generator=name_generator_for_scheme(self.naming_scheme)
            )
        else:
            raise ValueError(  # pragma: no cover
                "Invalid sheet file type of {}".format(self.sheet_path)
            )


# ---------------------------------------------------------------------------
# Library DataFrame helpers (used by variant_annotation, variant_filtration, …)
# ---------------------------------------------------------------------------

#: Default pandas query expressions applied when library_selection is None.
_LIBRARY_SELECTION_DEFAULTS: dict[str, str] = {
    "cancer": "role == 'tumor' and extraction_type == 'dna'",
    "germline": "extraction_type == 'dna'",
}

#: Library types whose ``libraryKit`` is an enrichment kit
TARGETED_LIBRARY_TYPES = ("WES", "Panel_seq")

#: Standard columns of the library dataframe built by :func:`build_library_dataframe`.
_LIBRARY_COLUMNS = [
    "library_name",
    "extraction_type",
    "kind",
    "role",
    "is_primary",
    "sex",
    "donor_name",
    "sample_name",
    "related_group",
    "father_name",
    "mother_name",
    "disease_state",
    "tissue_type",
    "is_background",
    "data_set",
    "secondary_id",
]

#: Sample sheet keys that already feed standard columns and are not added again.
_SHEET_KEYS_IN_LIBRARY_COLUMNS = frozenset(
    (
        "extractionType",
        "isTumor",
        "sex",
        "isAffected",
        "fatherName",
        "motherName",
        "fatherPk",
        "motherPk",
    )
)


def _snake_case(name: str) -> str:
    """Return ``name`` in snake_case, e.g. ``libraryKit`` -> ``library_kit``."""
    name = re.sub(r"\W+", "_", name)
    name = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)
    return re.sub(r"_+", "_", name).strip("_").lower()


def _library_extra_columns(raw_sheet) -> dict[str, dict[str, Any]]:
    """Return library name -> extra columns from the sample sheet's ``extra_infos``.

    Each library gets the ``extra_infos`` of its bio entity, bio sample, test sample and
    itself, with keys converted to snake_case. Keys that already feed standard columns
    (:data:`_SHEET_KEYS_IN_LIBRARY_COLUMNS`) are skipped.

    Raises ``ConfigError`` when a key maps to a standard column, or when two keys map to
    the same column for one library with different values.
    """
    result: dict[str, dict[str, Any]] = {}
    for bio_entity in raw_sheet.bio_entities.values():
        for bio_sample in bio_entity.bio_samples.values():
            for test_sample in bio_sample.test_samples.values():
                for library in test_sample.ngs_libraries.values():
                    columns: dict[str, Any] = {}
                    sources: dict[str, str] = {}
                    levels = (
                        ("bio entity", bio_entity),
                        ("bio sample", bio_sample),
                        ("test sample", test_sample),
                        ("NGS library", library),
                    )
                    for level, entity in levels:
                        for key, value in (entity.extra_infos or {}).items():
                            if key in _SHEET_KEYS_IN_LIBRARY_COLUMNS:
                                continue
                            column = _snake_case(key)
                            source = f"{key!r} of {level} {entity.name!r}"
                            if column in _LIBRARY_COLUMNS:
                                raise ConfigError(
                                    f"Sample sheet column {source} maps to {column!r}, "
                                    "which is a standard library column."
                                )
                            if column in columns and columns[column] != value:
                                raise ConfigError(
                                    f"Sample sheet columns {sources[column]} and {source} both "
                                    f"map to {column!r} for library {library.name!r}, with "
                                    "different values."
                                )
                            columns.setdefault(column, value)
                            sources.setdefault(column, source)
                    result[library.name] = columns
    return result


def build_library_dataframe(data_set_infos, relationships=None, include_background=False):
    """Build a tidy pandas DataFrame with one row per NGS library.

    The biomedsheets views (CancerCaseSheet for cancer sheets; GenericSampleSheet for the
    libraries and GermlineCaseSheet for the pedigrees of the other sheets) are used directly,
    so that information already computed by biomedsheets (primary pairs, pedigree roles,
    affected status, sex, …) is re-used rather than re-derived.

    Parameters
    ----------
    data_set_infos:
        Iterable of :class:`DataSetInfo` objects (``self.data_set_infos``), each with its
        loaded biomedsheets ``sheet``.

    Returns
    -------
    pandas.DataFrame
        One row per NGS library.  Columns:

        * ``library_name`` — NGS library identifier
        * ``extraction_type`` — ``"dna"``, ``"rna"``, … (always lower-cased)
        * ``kind`` — ``"cancer"`` or ``"germline"``
        * ``role`` — ``"tumor"`` / ``"normal"`` / ``"index"`` / ``"father"``
          / ``"mother"`` / ``"affected"`` / ``"unaffected"``
        * ``is_primary`` — ``True`` for the primary tumor (cancer) or the DNA library of the
          pedigree index (germline); ``False`` otherwise
        * ``sex`` — ``"male"``, ``"female"``, or ``"unknown"``
        * ``donor_name`` — patient / family identifier
        * ``sample_name`` — bio-sample name
        * ``tissue_type`` — ``"tumor"``, ``"normal"``, or ``"unknown"``
        * ``is_background`` — whether the library is from a background data set; these are
          only included with ``include_background``
        * ``data_set`` — the name of the library's data set in the config
        * ``secondary_id`` — the secondary ID of the donor in the sample sheet

        Further sample sheet columns (``extra_infos`` of the library and its test sample,
        bio sample and bio entity) are added in snake_case, e.g. ``libraryKit`` becomes
        ``library_kit``; see :func:`_library_extra_columns`.
    """
    import pandas as pd
    from biomedsheets.shortcuts import GenericSampleSheet
    from biomedsheets.shortcuts.cancer import CancerCaseSheet, CancerCaseSheetOptions

    rows: list[dict] = []

    for info in data_set_infos:
        if info.is_background and not include_background:
            continue
        _raw_sheet = info.sheet

        extra_columns = _library_extra_columns(_raw_sheet) if _raw_sheet is not None else {}

        sheet_type = _sheet_type(info)

        # ── matched_cancer ───────────────────────────────────────────────────
        if sheet_type == DataSetType.MATCHED_CANCER:
            csheet = CancerCaseSheet(
                _raw_sheet,
                options=CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True),
            )
            # Collect library names that belong to the primary tumor per donor
            # (first tumor in primary_sample_pairs).
            primary_tumor_lib_names: set[str] = set()
            for pair in csheet.primary_sample_pairs:
                if pair and pair.tumor_sample and pair.tumor_sample.dna_ngs_library:
                    primary_tumor_lib_names.add(pair.tumor_sample.dna_ngs_library.name)
                if pair and pair.tumor_sample and pair.tumor_sample.rna_ngs_library:
                    primary_tumor_lib_names.add(pair.tumor_sample.rna_ngs_library.name)

            for donor in csheet.donors:
                donor_name = donor.name
                for bio_sample in donor.bio_samples.values():
                    is_tumor = bio_sample.extra_infos.get("isTumor", False)
                    role = "tumor" if is_tumor else "normal"
                    sample_name = bio_sample.name
                    for ts in bio_sample.test_samples.values():
                        ext = ts.extra_infos.get("extractionType", "unknown").lower()
                        for lib in ts.ngs_libraries.values():
                            is_primary = lib.name in primary_tumor_lib_names if is_tumor else True
                            rows.append(
                                {
                                    "library_name": lib.name,
                                    "extraction_type": ext,
                                    "kind": "cancer",
                                    "role": role,
                                    "is_primary": is_primary,
                                    "sex": "unknown",
                                    "donor_name": donor_name,
                                    "sample_name": sample_name,
                                    "related_group": _group_prefix(info) + donor_name,
                                    "father_name": "0",
                                    "mother_name": "0",
                                    "disease_state": 0,
                                    "tissue_type": "tumor" if is_tumor else "normal",
                                    "is_background": info.is_background,
                                    "data_set": info.name,
                                    "secondary_id": donor.bio_entity.secondary_id,
                                    **extra_columns.get(lib.name, {}),
                                }
                            )

        # ── germline_variants (and generic) ──────────────────────────────────
        else:
            # All libraries, with the pedigrees of the germline view (joined by the data set's
            # ``pedigree_field``, if set)
            index_donor_names: set[str] = set()
            index_lib_names: set[str] = set()
            father_lib_names: set[str] = set()
            mother_lib_names: set[str] = set()
            lib_sex: dict[str, str] = {}
            donor_family: dict[str, str] = {}
            lib_father: dict[str, str] = {}
            lib_mother: dict[str, str] = {}
            lib_disease: dict[str, int] = {}

            all_ngs_libs = GenericSampleSheet(_raw_sheet).all_ngs_libraries
            germline_sheet = _germline_case_sheet(info, _raw_sheet)
            pedigrees = germline_sheet.cohort.pedigrees

            if pedigrees:
                # First pass: map donor PKs to their DNA library names
                pk_to_lib = {}
                for ped in pedigrees:
                    for d in ped.donors:
                        lib = getattr(d, "dna_ngs_library", None)
                        if lib:
                            pk_to_lib[getattr(d, "pk", None)] = lib.name

                families: dict[str, str] = {}
                for ped in pedigrees:
                    family = _family_name(info, ped)
                    if family in families:
                        raise ConfigError(
                            f"the pedigrees of {families[family]} and {_donor_name(ped.index)} in "
                            f"data set {info.name!r} are both named {family}; group them with "
                            f"pedigree_field: {FAMILY_FIELD} on the data set, or fix the sheet"
                        )
                    families[family] = _donor_name(ped.index)
                    # Index / proband: all libraries of the donor, so the RNA library too
                    if ped.index and ped.index.dna_ngs_library:
                        index_donor_names.add(_donor_name(ped.index))
                        index_lib_names.add(ped.index.dna_ngs_library.name)

                    # Build sets of PKs referenced as father / mother by children
                    pk_is_father: set = {
                        int(d.father_pk)
                        for d in ped.donors
                        if getattr(d, "father_pk", None) is not None
                    }
                    pk_is_mother: set = {
                        int(d.mother_pk)
                        for d in ped.donors
                        if getattr(d, "mother_pk", None) is not None
                    }

                    for d in ped.donors:
                        donor_family[_donor_name(d)] = family
                        lib = getattr(d, "dna_ngs_library", None)
                        if lib is None:
                            continue
                        pk = getattr(d, "pk", None)
                        sex_raw = (
                            d.extra_infos.get("sex", "unknown")
                            if hasattr(d, "extra_infos")
                            else "unknown"
                        )
                        lib_sex[lib.name] = sex_raw if sex_raw in ("male", "female") else "unknown"
                        if pk in pk_is_father:
                            father_lib_names.add(lib.name)
                        elif pk in pk_is_mother:
                            mother_lib_names.add(lib.name)

                        father_pk = getattr(d, "father_pk", None)
                        mother_pk = getattr(d, "mother_pk", None)
                        lib_father[lib.name] = pk_to_lib.get(
                            int(father_pk) if father_pk is not None else None, "0"
                        )
                        lib_mother[lib.name] = pk_to_lib.get(
                            int(mother_pk) if mother_pk is not None else None, "0"
                        )

                        is_affected = (
                            d.extra_infos.get("isAffected") if hasattr(d, "extra_infos") else None
                        )
                        if is_affected is True or str(is_affected).lower() == "affected":
                            lib_disease[lib.name] = 2
                        elif is_affected is False or str(is_affected).lower() == "unaffected":
                            lib_disease[lib.name] = 1
                        else:
                            lib_disease[lib.name] = 0

            for lib in all_ngs_libs:
                ts = lib.test_sample
                bs = getattr(ts, "bio_sample", None)
                bio_entity = getattr(bs, "bio_entity", None)

                # Handle both GenericSampleSheet (TestSample has extra_infos)
                # and GermlineCaseSheet (TestSampleShortcut wraps TestSample).
                ts_extra = getattr(ts, "extra_infos", None)
                if ts_extra is None and hasattr(ts, "test_sample"):
                    ts_extra = getattr(ts.test_sample, "extra_infos", {})
                else:
                    ts_extra = ts_extra or {}

                be_extra = getattr(bio_entity, "extra_infos", {}) if bio_entity else {}

                ext = ts_extra.get("extractionType", "unknown").lower()
                lib_name = lib.name
                donor_name = getattr(bio_entity, "name", "")
                sample_name = getattr(bs, "name", "")
                is_affected = be_extra.get("isAffected") or ts_extra.get("isAffected", None)
                # TSV sheets hold the words, JSON sheets booleans
                status = str(is_affected).lower() if is_affected is not None else "unknown"
                sex = lib_sex.get(lib_name, "unknown")

                # Only the index's DNA library is primary: the library that stands for the pedigree
                is_primary = lib_name in index_lib_names
                if donor_name in index_donor_names:
                    role = "index"
                elif lib_name in father_lib_names:
                    role = "father"
                elif lib_name in mother_lib_names:
                    role = "mother"
                elif is_affected is True or status == "affected":
                    role = "affected"
                elif is_affected is False or status == "unaffected":
                    role = "unaffected"
                elif not pedigrees:
                    role = "index"  # no pedigree information: each library stands alone
                    is_primary = True
                else:
                    role = "unknown"

                rows.append(
                    {
                        "library_name": lib_name,
                        "extraction_type": ext,
                        "kind": "germline",
                        "role": role,
                        "is_primary": is_primary,
                        "sex": sex,
                        "donor_name": donor_name,
                        "sample_name": sample_name,
                        "related_group": donor_family.get(
                            donor_name, _group_prefix(info) + donor_name
                        ),
                        "father_name": lib_father.get(lib_name, "0"),
                        "mother_name": lib_mother.get(lib_name, "0"),
                        "disease_state": lib_disease.get(lib_name, 0),
                        "tissue_type": "unknown",
                        "is_background": info.is_background,
                        "data_set": info.name,
                        "secondary_id": getattr(bio_entity, "secondary_id", ""),
                        **extra_columns.get(lib_name, {}),
                    }
                )

    extra_column_names = sorted({column for row in rows for column in row} - set(_LIBRARY_COLUMNS))
    columns = _LIBRARY_COLUMNS + extra_column_names
    df = pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)
    if relationships:
        df = resolve_relationships(df, relationships)
    return df


#: The columns of the donor table of :func:`build_donor_dataframe`
_DONOR_COLUMNS = [
    "donor_name",
    "father_name",
    "mother_name",
    "sex",
    "library_name",
    "pedigree",
    "data_set",
    "is_background",
]


def _donor_name(donor) -> str:
    """The name of a donor as the library table names it: that of its bio entity."""
    return getattr(getattr(donor, "bio_entity", None), "name", None) or donor.name


# ---------------------------------------------------------------------------
# Related groups: their names, entities and libraries
# ---------------------------------------------------------------------------


#: The custom field that names a family on germline sheets, such as those of ``snappy pull-sheet``
FAMILY_FIELD = "familyId"

#: The prefix of a related group that is named after a donor, by sheet type
_GROUP_PREFIXES = {
    DataSetType.MATCHED_CANCER: "PAT_",
    DataSetType.GERMLINE_VARIANTS: "FAM_",
    DataSetType.GENERIC: "",
}


def _sheet_type(info) -> str:
    return getattr(info, "sheet_type", "") or str(getattr(info, "type", ""))


def _germline_case_sheet(info, raw_sheet):
    """The biomedsheets germline view of a sheet, with the pedigrees joined by the data set's ``pedigree_field``.
    Raises ``ConfigError`` when a donor lacks the field."""
    from biomedsheets.shortcuts import GermlineCaseSheet
    from biomedsheets.shortcuts.germline import UndefinedFieldException

    try:
        return GermlineCaseSheet(raw_sheet, **(getattr(info, "pedigree_field_kwargs", None) or {}))
    except UndefinedFieldException as exc:
        raise ConfigError(f"pedigree_field of data set {info.name!r}: {exc}") from exc


def _group_prefix(info) -> str:
    """The prefix of the related groups of a data set that are named after a donor."""
    return _GROUP_PREFIXES.get(_sheet_type(info), "")


def _family_name(info, pedigree) -> str:
    """Return the name of a pedigree, the ``related_group`` of its members.
    It is the value of the data set's ``pedigree_field``, else of the ``familyId`` custom field, over the donors that have one.
    Without any, it is the index donor's name with the prefix of the sheet type: ``FAM_`` on germline sheets, none on generic ones.
    Raises ``ConfigError`` when the donors disagree: the pedigree then comes from the parent columns, and the field names another family."""
    field = (getattr(info, "pedigree_field_kwargs", None) or {}).get("join_by_field")
    name_field = field or FAMILY_FIELD
    values = {
        str(value)
        for donor in pedigree.donors
        if (value := getattr(donor, "extra_infos", {}).get(name_field)) not in (None, "", ".")
    }
    index = _donor_name(pedigree.index)
    if len(values) > 1:
        raise ConfigError(
            f"the pedigree of {index} has several values of {name_field}: "
            f"{', '.join(sorted(values))}; group the pedigrees by the field with pedigree_field: "
            f"{name_field} on the data set, or fix the sheet"
        )
    if field and not values:
        raise ConfigError(
            f"{field} of data set {info.name!r} has no value for the pedigree of {index}"
        )
    return values.pop() if values else _group_prefix(info) + index


def entity_libraries(df, group_by: str | None) -> dict[str, str]:
    """Return the output entities of the libraries ``df`` with their primary library.
    With ``group_by == "related"`` each ``related_group`` is an entity, named like the group, and its primary library is the ``is_primary`` row, else the first one; otherwise every library is an entity of its own."""
    if df.empty:
        return {}
    if group_by == "related":
        entities = {}
        for name, group in df.groupby("related_group"):
            primary = group[group["is_primary"]]
            entities[name] = (primary if not primary.empty else group).iloc[0]["library_name"]
        return entities
    return {name: name for name in df["library_name"].unique()}


def _check_groups(df) -> None:
    """Raise ``ConfigError`` for a related group whose libraries ``df`` span several data sets,
    since the group's files are named after it."""
    if df.empty or not {"related_group", "data_set"} <= set(df.columns):
        return
    data_sets = df.groupby("related_group")["data_set"].agg(lambda names: sorted(set(names)))
    if colliding := {name: sets for name, sets in data_sets.items() if len(sets) > 1}:
        names = ", ".join(f"{name} ({', '.join(sets)})" for name, sets in colliding.items())
        first_data_set = next(iter(colliding.values()))[0]
        raise ConfigError(
            f"related groups of the same name in several data sets: {names}. A task names the "
            "files of a related group after it, so it can only take one group of a name. A family "
            "sequenced in batches belongs into one sheet: one data set, whose libraries the files "
            "tables list from all batches. For unrelated groups, select the libraries of one data "
            f'set with library_selection, for example library_selection: "data_set == '
            f"'{first_data_set}'\", or use one task per data set."
        )


def build_donor_dataframe(data_set_infos, include_background=False):
    """Build a DataFrame with one row per donor of the germline sample sheets, sequenced or
    not: the parents as donor names (``"0"`` for none), ``sex``, the donor's primary DNA
    library (``""`` without one) and ``pedigree``, the name of its pedigree (the ``related_group``
    of its members). Cancer sheets have no rows: their donors are all sequenced."""
    import pandas as pd

    rows: list[dict] = []
    for info in data_set_infos:
        if info.is_background and not include_background:
            continue
        sheet_type = _sheet_type(info)
        raw_sheet = info.sheet
        if sheet_type == DataSetType.MATCHED_CANCER or raw_sheet is None:
            continue
        germline_sheet = _germline_case_sheet(info, raw_sheet)
        for pedigree in germline_sheet.cohort.pedigrees:
            family = _family_name(info, pedigree)
            for donor in pedigree.donors:
                parents = []
                for pk in (donor.father_pk, donor.mother_pk):
                    parent = pedigree.pk_to_donor.get(str(pk)) if pk is not None else None
                    parents.append(_donor_name(parent) if parent else "0")
                sex = (
                    donor.extra_infos.get("sex", "unknown")
                    if hasattr(donor, "extra_infos")
                    else "unknown"
                )
                rows.append(
                    {
                        "donor_name": _donor_name(donor),
                        "father_name": parents[0],
                        "mother_name": parents[1],
                        "sex": sex if sex in ("male", "female") else "unknown",
                        "library_name": donor.dna_ngs_library.name if donor.dna_ngs_library else "",
                        "pedigree": family,
                        "data_set": info.name,
                        "is_background": info.is_background,
                    }
                )
    return pd.DataFrame(rows, columns=_DONOR_COLUMNS)


def resolve_relationships(df, relationships: dict):
    """Add relationship columns to the library DataFrame.

    For each relationship definition, a new column is added whose value is
    looked up from related rows sharing the same ``via`` column value and
    matching the ``target`` query expression.

    Parameters
    ----------
    df:
        Library DataFrame produced by :func:`build_library_dataframe`.
    relationships:
        Dict mapping relationship names to :class:`RelationshipDefinition`
        instances (or dicts with ``via``, ``target``, ``column``, ``many``
        keys).

    Returns
    -------
    pandas.DataFrame
        Copy of *df* with additional relationship columns.
    """
    if df.empty or not relationships:
        return df

    result = df.copy()

    for rel_name, rel_def in relationships.items():
        # Accept both model instances and plain dicts.
        via = rel_def.via if hasattr(rel_def, "via") else rel_def["via"]
        target = rel_def.target if hasattr(rel_def, "target") else rel_def["target"]
        column = (
            rel_def.column if hasattr(rel_def, "column") else rel_def.get("column")
        ) or rel_name
        many = rel_def.many if hasattr(rel_def, "many") else rel_def.get("many", False)

        via_columns = [via] if isinstance(via, str) else list(via)
        for via_column in via_columns:
            if via_column not in result.columns:
                raise ConfigError(
                    f"Relationship {rel_name!r}: join column {via_column!r} not found in "
                    f"DataFrame. Available columns: {list(result.columns)}"
                )
        # One join key per row; several ``via`` columns join on all of them
        via_keys = result[via_columns].astype(str).agg("/".join, axis=1)

        # Build a lookup: for each unique ``via`` value, find matching rows.
        # We query the full DataFrame for each group so that the target
        # expression can reference any column.
        matched: dict[str, list[str]] = {}
        for via_val, group in result.groupby(via_keys):
            try:
                hits = group.query(target)
            except Exception as exc:
                raise ConfigError(
                    f"Relationship {rel_name!r}: target query {target!r} failed: {exc}"
                ) from exc
            matched[str(via_val)] = hits["library_name"].tolist()

        if many:
            result[column] = via_keys.map(lambda v: matched.get(str(v), []))
        else:

            def _single_match(via_value):
                matches = matched.get(str(via_value), [])
                if len(matches) > 1:
                    raise ConfigError(
                        f"Relationship {rel_name!r}: expected single match for "
                        f"via value {via_value!r}, got {len(matches)}: {matches}. "
                        f"Use 'many: true' if multiple matches are expected."
                    )
                return matches[0] if matches else ""

            result[column] = via_keys.map(_single_match)

    return result


def select_libraries(
    df, selection: str | None, default: str | None = None, requirement: str | None = None
):
    """Return the libraries of ``df`` that a task takes.

    These are the rows of the step's ``requirement`` (what it can process at all) that are
    also rows of ``selection`` (the task's ``library_selection``), else those of the
    step's ``default``, else those of the default of each row's sheet kind
    (``_LIBRARY_SELECTION_DEFAULTS``); an empty expression takes all rows. Raises ``ConfigError``
    for an invalid expression and for libraries of the same name in several data sets.
    """
    import pandas as pd

    if df.empty:
        return df
    if requirement:
        df = _query_libraries(df, requirement)
    expr = selection if selection is not None else default
    if expr is None:
        df = pd.concat(
            [
                subset if subset.empty else _query_libraries(subset, kind_default)
                for kind, kind_default in _LIBRARY_SELECTION_DEFAULTS.items()
                for subset in [df[df["kind"] == kind]]
            ]
        )
    elif expr:
        df = _query_libraries(df, expr)
    if colliding := _colliding_libraries(df):
        raise ConfigError(colliding_libraries_message(colliding))
    _check_groups(df)
    return df


def _query_libraries(df, expr: str):
    try:
        return df.query(expr)
    except Exception as exc:
        raise ConfigError(f"library_selection expression {expr!r} is invalid: {exc}") from exc


def _colliding_libraries(df) -> dict[str, list[str]]:
    """Return library name -> data sets, for the names of several libraries in ``df``."""
    duplicated = df[df["library_name"].duplicated(keep=False)]
    return {name: list(group["data_set"]) for name, group in duplicated.groupby("library_name")}


def colliding_libraries_message(colliding: dict[str, list[str]]) -> str:
    """Return the error for libraries of the same name from several data sets, with the fix."""
    names = ", ".join(f"{name} ({', '.join(data_sets)})" for name, data_sets in colliding.items())
    first_data_set = next(iter(colliding.values()))[0]
    return (
        f"libraries of the same name in several data sets: {names}. A task names its files by "
        "library, so it can only take one library of a name. Select the libraries of one data "
        f'set with library_selection, for example library_selection: "data_set == '
        f"'{first_data_set}'\", or use one task per data set."
    )
