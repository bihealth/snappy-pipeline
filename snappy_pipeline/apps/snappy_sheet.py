"""``snappy sheet``: check plain TSV sample sheets, and convert biomedsheets sheets to them.

``convert`` loads a biomedsheets TSV the way snappy did (``biomedsheets`` must be installed) and
writes the plain layout with the same names: the library, donor and sample names of the library
table, the pedigree names as ``family_id``, the parents as donor names, the index that biomedsheets
chose as ``is_index`` and, where a donor has several libraries of one role and extraction type, the
one biomedsheets used as ``is_primary``. Custom fields become columns of the same name in snake
case; ``extractionType`` and ``ncbiTaxon`` are dropped, and ``seqPlatform`` when every library is
Illumina.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

from snappy_pipeline.base import ConfigError
from snappy_pipeline.sample_sheet import COLUMNS, Sheets, read_sheets

#: Sheet keys that the plain sheet derives or does not carry
DROPPED = frozenset({"extractionType", "ncbiTaxon", "isTumor", "fatherPk", "motherPk"})
#: biomedsheets' codes of the germline sheet -> the plain values
SEX = {"male": "male", "female": "female"}
AFFECTED = {"affected": True, "unaffected": False}


def convert(path: str, sheet_type: str, naming_scheme: str, pedigree_field: str | None) -> str:
    """Return the plain TSV of the biomedsheets sheet ``path``; raises ``ConfigError`` when
    biomedsheets cannot be imported or the sheet does not load."""
    try:
        from snappy_pipeline.library_table import (
            DataSetInfo,
            build_donor_dataframe,
            build_library_dataframe,
        )
    except ImportError as e:  # pragma: no cover - depends on the installed extras
        raise ConfigError(
            "snappy sheet convert needs the biomedsheets package (pip install biomedsheets)"
        ) from e
    info = DataSetInfo("converted", path, sheet_type, False, naming_scheme, pedigree_field)
    libraries = build_library_dataframe([info], include_background=True)
    donors = build_donor_dataframe([info], include_background=True)
    entities = _entities(info.sheet)
    rows = [_library_row(row, entities, libraries) for row in libraries.to_dict("records")]
    rows += [
        _donor_row(row, entities) for row in donors.to_dict("records") if not row["library_name"]
    ]
    return _write(rows)


def _entities(raw_sheet) -> dict[str, dict]:
    """Donor name (as the library table names it) -> the bio entity's facts and the names of its
    parents as donors."""
    by_secondary_id = {entity.secondary_id: entity for entity in raw_sheet.bio_entities.values()}
    entities = {}
    for entity in raw_sheet.bio_entities.values():
        infos = dict(entity.extra_infos or {})
        facts = {"secondary_id": entity.secondary_id}
        for column in ("father", "mother"):
            parent = infos.pop(f"{column}Name", None)
            parent = by_secondary_id.get(parent) if parent not in (None, "0") else None
            facts[column] = parent.name if parent is not None else None
        facts["sex"] = SEX.get(infos.pop("sex", None))
        facts["is_affected"] = AFFECTED.get(infos.pop("isAffected", None))
        hpo = infos.pop("hpoTerms", None)
        facts["hpo_terms"] = ",".join(hpo) if hpo else None
        facts["extras"] = infos
        entities[entity.name] = facts
    return entities


def _snake_case(name: str) -> str:
    name = re.sub(r"\W+", "_", name)
    name = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)
    return re.sub(r"_+", "_", name).strip("_").lower()


def _library_row(row: dict, entities: dict, libraries) -> dict:
    """The plain row of a library from its library table row."""
    donor = entities[row["donor_name"]]
    cancer = row["kind"] == "cancer"
    plain = {
        "donor_name": row["donor_name"],
        "family_id": None if cancer else row["related_group"],
        "father": donor["father"],
        "mother": donor["mother"],
        "sex": None if cancer else donor["sex"],
        "is_affected": None if cancer else donor["is_affected"],
        "is_index": True if row["role"] == "index" and not cancer else None,
        "hpo_terms": donor["hpo_terms"],
        "sample_name": row["sample_name"],
        "is_tumor": (row["role"] == "tumor") if cancer else None,
        "library_name": row["library_name"],
        "is_primary": _is_primary(row, libraries),
    }
    for column in COLUMNS:
        if column not in plain and column in row and not _missing(row[column]):
            plain[column] = row[column]
    if not plain.get("is_background"):
        plain.pop("is_background", None)  # the data set's flag, not a fact of the library
    for key, value in donor["extras"].items():
        if key not in DROPPED and not _missing(value):
            plain.setdefault(_snake_case(key), value)
    for column, value in row.items():
        if column not in plain and column not in _TABLE_ONLY and not _missing(value):
            plain[column] = value
    return plain


#: Library table columns that are no sheet columns
_TABLE_ONLY = frozenset(
    {
        "extraction_type",
        "kind",
        "role",
        "related_group",
        "father_name",
        "mother_name",
        "disease_state",
        "tissue_type",
        "is_background",
        "data_set",
        "secondary_id",
        "ncbi_taxon",
    }
)


def _is_primary(row: dict, libraries) -> bool | None:
    """``True`` for the library biomedsheets used among several of one donor, role and
    extraction type (the first it marked primary, else the first); ``None`` when alone."""
    same = libraries[
        (libraries["donor_name"] == row["donor_name"])
        & (libraries["role"] == row["role"])
        & (libraries["extraction_type"] == row["extraction_type"])
    ]
    if len(same) < 2:
        return None
    marked = same[same["is_primary"]]
    chosen = (marked if not marked.empty else same).iloc[0]["library_name"]
    return row["library_name"] == chosen


def _donor_row(row: dict, entities: dict) -> dict:
    donor = entities[row["donor_name"]]
    plain = {
        "donor_name": row["donor_name"],
        "family_id": row["pedigree"],
        "father": donor["father"],
        "mother": donor["mother"],
        "sex": donor["sex"],
        "is_affected": donor["is_affected"],
        "hpo_terms": donor["hpo_terms"],
    }
    for key, value in donor["extras"].items():
        if key not in DROPPED and not _missing(value):
            plain[_snake_case(key)] = value
    return plain


def _missing(value) -> bool:
    """Whether a table or sheet value is no value: empty, ``None`` or the string biomedsheets
    sheets wrote for it."""
    if value is None or value in ("", "None", []):
        return True
    try:
        return value != value  # NaN
    except Exception:  # pragma: no cover
        return False


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        return ",".join(map(str, value))
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _write(rows: list[dict]) -> str:
    """The TSV text: the known columns with a value in model order, then the extras."""
    present = {column for row in rows for column, value in row.items() if not _missing(value)}
    if "seq_platform" in present and {row.get("seq_platform") for row in rows} <= {
        "Illumina",
        None,
    }:
        present.discard("seq_platform")
    columns = [name for name in COLUMNS if name in present]
    columns += sorted(name for name in present if name not in COLUMNS)
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(_cell(row.get(column)) for column in columns))
    return "\n".join(lines) + "\n"


def summary(sheets: Sheets) -> str:
    """One line about the sheets: their libraries and donors."""
    donors = {row.values.donor_name for row in sheets.rows}
    files = ", ".join(sheet.path for sheet in sheets.sheets)
    return f"{len(sheets.libraries)} libraries of {len(donors)} donors in {files}"


def check(paths: Iterable[str]) -> Sheets:
    """Read the sheets; raises :class:`~snappy_pipeline.sample_sheet.SheetError`."""
    return read_sheets(paths)
