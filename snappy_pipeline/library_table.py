"""The library table of a project: one row per library of the sample sheets.

``Project.sheets`` reads the sample sheets once (:func:`snappy_pipeline.sample_sheet.read_sheets`);
:func:`build_library_dataframe` turns them into the table that every step reads through ``BaseStep``
(``library_table``, ``build_library_dataframe``, ``libraries``), :func:`build_donor_dataframe`
gives the donors, sequenced or not, :func:`resolve_relationships` adds the relationship columns
(``matched_normal_lib``), :func:`select_libraries` applies a task's ``library_selection`` and
:func:`entity_libraries` groups the libraries of a task by its ``group_by``.
"""

from __future__ import annotations

from snappy_pipeline.base import ConfigError
from snappy_pipeline.sample_sheet import (
    COLUMNS,
    TARGETED_LIBRARY_TYPES,
    Row,
    Sheets,
    primary_library,
)

__all__ = [
    "TARGETED_LIBRARY_TYPES",
    "DEFAULT_SELECTION",
    "DERIVED_COLUMNS",
    "GROUP_COLUMNS",
    "build_library_dataframe",
    "build_donor_dataframe",
    "entity_libraries",
    "group_names",
    "resolve_related",
    "resolve_relationships",
    "select_libraries",
]

#: The libraries a task takes without ``library_selection`` and without a step default: the
#: DNA libraries, within the step's requirement
DEFAULT_SELECTION = "extraction_type == 'dna'"

#: The columns the table derives from the sheet columns
DERIVED_COLUMNS = [
    "extraction_type",
    "role",
    "father_name",
    "mother_name",
    "disease_state",
]

#: ``group_by`` level -> the column whose value names the group; ``related`` resolves to
#: ``family`` or ``donor`` per task (:func:`resolve_related`)
GROUP_COLUMNS = {"library": "library_name", "donor": "donor_name", "family": "family_id"}

#: The columns of the table, in order: the names, the derived ones, then every sheet column
_COLUMN_ORDER = ["library_name", "donor_name", "sample_name", *DERIVED_COLUMNS] + [
    name for name in COLUMNS if name not in ("library_name", "donor_name", "sample_name")
]


def build_library_dataframe(sheets: Sheets, relationships=None, include_background=False):
    """Build a pandas DataFrame with one row per library of ``sheets``.

    Columns: every known sheet column (:data:`snappy_pipeline.sample_sheet.COLUMNS`, ``None``
    where the sheets have no value), the extra columns of the sheets as strings, and the derived
    columns:

    * ``extraction_type``: ``dna`` or ``rna``, from ``library_type``
    * ``role``: ``tumor`` for a library of a tumor sample (``is_tumor``), else ``normal``
    * ``father_name``, ``mother_name``: the primary DNA library of the parent, ``"0"`` without
      one
    * ``disease_state``: 2 (``is_affected`` true), 1 (false) or 0 (unknown)
    * ``sex``: ``unknown`` where the sheet has none
    * ``is_primary``: the primary library among a donor's libraries of one role and extraction
      type, marked or alone
    * ``is_index``: the index donor of its family, marked or alone
    * ``is_background``: false where the sheet has no value

    Background libraries (``is_background``) are left out unless ``include_background`` is set.
    ``relationships`` adds the columns of :func:`resolve_relationships`.
    """
    import pandas as pd

    rows = [_library_row(row, sheets) for row in sheets.libraries]
    extras = [name for name in sheets.columns if name not in COLUMNS]
    df = pd.DataFrame(rows, columns=_COLUMN_ORDER + extras)
    if not include_background:
        df = df[~df["is_background"]].reset_index(drop=True)
    if relationships:
        df = resolve_relationships(df, relationships)
    return df


def _library_row(row: Row, sheets: Sheets) -> dict:
    values = row.values
    record = {name: getattr(values, name) for name in COLUMNS}
    record.update(values.model_extra or {})
    record["extraction_type"] = values.extraction_type
    record["role"] = values.role
    record["sex"] = values.sex or "unknown"
    record["is_background"] = bool(values.is_background)
    record["is_primary"] = _is_primary(row, sheets)
    record["is_index"] = _is_index(values.donor_name, sheets)
    record["disease_state"] = {True: 2, False: 1}.get(values.is_affected, 0)
    for column in ("father", "mother"):
        parent = getattr(values, column)
        record[f"{column}_name"] = _primary_dna_library(parent, sheets) if parent else "0"
    return record


def _is_primary(row: Row, sheets: Sheets) -> bool:
    """Whether ``row`` is the primary library of its donor, role and extraction type: the one
    marked ``is_primary``, or the only one."""
    values = row.values
    group = [
        other
        for other in sheets.libraries
        if other.values.donor_name == values.donor_name
        and other.values.role == values.role
        and other.values.extraction_type == values.extraction_type
    ]
    return primary_library(group) is row


def _family_of(donor: Row) -> str:
    """The family of a donor: its ``family_id``, else the donor alone."""
    return donor.values.family_id or donor.values.donor_name


def _donors(sheets: Sheets) -> dict[str, Row]:
    donors: dict[str, Row] = {}
    for row in sheets.rows:
        donors.setdefault(row.values.donor_name, row)
    return donors


def _is_index(donor_name: str, sheets: Sheets) -> bool:
    """Whether the donor is the index of its family: marked ``is_index``, or the only donor of the
    family with a normal DNA library."""
    donors = _donors(sheets)
    donor = donors[donor_name]
    if donor.values.is_index is not None:
        return donor.values.is_index
    family = _family_of(donor)
    germline = {
        r.values.donor_name
        for r in sheets.libraries
        if r.values.role == "normal" and r.values.extraction_type == "dna"
    }
    members = [
        d for d in donors.values() if _family_of(d) == family and d.values.donor_name in germline
    ]
    return members == [donor]


def _primary_dna_library(donor_name: str, sheets: Sheets) -> str:
    """The primary normal DNA library of a donor, ``"0"`` without one."""
    group = [
        row
        for row in sheets.libraries
        if row.values.donor_name == donor_name
        and row.values.role == "normal"
        and row.values.extraction_type == "dna"
    ]
    return primary_library(group).values.library_name if group else "0"


#: The columns of the donor table of :func:`build_donor_dataframe`
_DONOR_COLUMNS = [
    "donor_name",
    "father_name",
    "mother_name",
    "sex",
    "library_name",
    "pedigree",
    "is_background",
]


def build_donor_dataframe(sheets: Sheets, include_background=False):
    """Build a DataFrame with one row per donor of the sheets, sequenced or not: the parents as
    donor names (``"0"`` for none), ``sex``, the donor's primary normal DNA library (``""``
    without one) and ``pedigree``, its ``family_id`` (``""`` without one). A donor is a
    background donor when all its libraries are."""
    import pandas as pd

    rows = []
    for donor in _donors(sheets).values():
        values = donor.values
        libraries = [r for r in sheets.libraries if r.values.donor_name == values.donor_name]
        background = bool(libraries) and all(r.values.is_background for r in libraries)
        if background and not include_background:
            continue
        library = _primary_dna_library(values.donor_name, sheets)
        rows.append(
            {
                "donor_name": values.donor_name,
                "father_name": values.father or "0",
                "mother_name": values.mother or "0",
                "sex": values.sex or "unknown",
                "library_name": "" if library == "0" else library,
                "pedigree": values.family_id or "",
                "is_background": background,
            }
        )
    return pd.DataFrame(rows, columns=_DONOR_COLUMNS)


def resolve_related(df, group_by: str) -> str:
    """Return ``group_by`` with ``related`` resolved for the libraries ``df``: ``family`` when
    every library has a ``family_id``, ``donor`` when none has; raises ``ConfigError`` when they
    mix, since the level is then a choice to make explicit."""
    if group_by != "related" or df.empty:
        return group_by
    with_family = df["family_id"].notna()
    if with_family.all():
        return "family"
    if not with_family.any():
        return "donor"
    raise ConfigError(
        "group_by: related is the family of libraries with a family_id and the donor of the "
        "others, but this task has both: "
        + ", ".join(df[with_family]["library_name"].head(3))
        + " have one, "
        + ", ".join(df[~with_family]["library_name"].head(3))
        + " have none; set group_by: family or group_by: donor"
    )


def group_names(df, group_by: str):
    """The name of the group of each library of ``df`` at the level ``group_by``: the library,
    its donor or its family (``related``: see :func:`resolve_related`). Raises ``ConfigError``
    for libraries without a ``family_id`` at the family level."""
    column = GROUP_COLUMNS[resolve_related(df, group_by)]
    names = df[column]
    if column == "family_id" and (missing := df[names.isna()]["library_name"].tolist()):
        raise ConfigError(
            "group_by: family needs a family_id in the sample sheet for every library; without: "
            + ", ".join(missing)
        )
    return names


def entity_libraries(df, group_by: str | None) -> dict[str, str]:
    """Return the output entities of the libraries ``df`` with their primary library.

    With ``group_by: library`` every library is an entity of its own. With ``donor``, ``family``
    or ``related`` (:func:`resolve_related`) each group is an entity, named like the group, and
    its primary library is the primary normal DNA library of the group's index donor, else the
    first primary library, else the first library.
    """
    if df.empty:
        return {}
    if not group_by or group_by == "library":
        return {name: name for name in df["library_name"].unique()}
    entities = {}
    for name, group in df.groupby(group_names(df, group_by), sort=False):
        index = group[group["is_index"] & group["is_primary"] & (group["role"] == "normal")]
        index = index[index["extraction_type"] == "dna"]
        primary = group[group["is_primary"]]
        chosen = index if not index.empty else primary if not primary.empty else group
        entities[name] = chosen.iloc[0]["library_name"]
    return entities


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
    df,
    selection: str | None,
    default: str | None = None,
    requirement: str | None = None,
    upstream_libraries=None,
):
    """Return the libraries of ``df`` that a task takes.

    These are the rows of the step's ``requirement`` (what it can process at all) that are
    also rows of ``selection`` (the task's ``library_selection``), else those of the step's
    ``default``, else the libraries of the task it reads (``upstream_libraries``, the names; see
    ``BaseStep.upstream_libraries``), else the DNA libraries (:data:`DEFAULT_SELECTION`); an
    empty expression takes all rows.
    Raises ``ConfigError`` for an invalid expression.
    """
    if df.empty:
        return df
    if requirement:
        df = _query_libraries(df, requirement)
    expr = selection if selection is not None else default
    if expr is None and upstream_libraries is not None:
        return df[df["library_name"].isin(set(upstream_libraries))]
    if expr is None:
        expr = DEFAULT_SELECTION
    if expr:
        df = _query_libraries(df, expr)
    return df


def _query_libraries(df, expr: str):
    try:
        return df.query(expr)
    except Exception as exc:
        raise ConfigError(f"library_selection expression {expr!r} is invalid: {exc}") from exc
