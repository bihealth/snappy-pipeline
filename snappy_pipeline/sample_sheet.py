"""Sample sheets: one plain TSV with one row per library, read and checked once per run.

A sheet has one header line and one row per library, or one row per donor without libraries (an
unsequenced relative). Columns are matched by name in any order; every column works as is in a
``library_selection`` query. Lines starting with ``#`` and blank lines are skipped. A missing value
is an empty cell: ``.`` is an error. Booleans are ``true`` or ``false``.

Sheets that collaborators typed are read as they are: the synonyms of :data:`SYNONYMS` (``F``,
``woman`` or ``XX`` for ``female``; ``yes``, ``Y`` or ``1`` for ``true``; ``exome`` for ``WES``;
``illumina`` for ``Illumina``) are read as the canonical value, and each such cell is a note, which
``snappy sheet check`` lists and ``snappy sheet check --fix`` rewrites. The PED codes of
``is_affected`` (``1``, ``2``) are rejected as ambiguous.

The core columns are ``donor_name``, ``library_name`` and ``library_type``. The other known columns
(:data:`COLUMNS`) are typed and checked at load; a step declares the columns it needs. Any other
snake_case column is kept as a string.

Each known column has a level, the key whose rows must agree on it (donor, sample, extract or
library), and the ISA-Tab header it maps to for SODAR. Choices are explicit: ``is_index`` marks the
index of a family and ``is_primary`` the primary library among a donor's libraries of one role and
extraction type; both are implied when there is one candidate and required when there are several.
"""

from __future__ import annotations

import csv
import enum
import keyword
import re
import types
import typing
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, ValidationError

from snappy_pipeline.base import ConfigError

NAME_PATTERN = r"^[A-Za-z0-9_.-]+$"
Name = typing.Annotated[str, StringConstraints(pattern=NAME_PATTERN)]
Sex = typing.Literal["male", "female", "unknown"]
SeqPlatform = typing.Literal["Illumina", "PacBio", "ONT", "other"]
LibraryType = typing.Literal["WES", "WGS", "Panel_seq", "mRNA_seq", "total_RNA_seq"]

#: Library type -> extraction type
EXTRACTION_TYPES = {
    "WES": "dna",
    "WGS": "dna",
    "Panel_seq": "dna",
    "mRNA_seq": "rna",
    "total_RNA_seq": "rna",
}
#: Library types whose ``library_kit`` is an enrichment kit
TARGETED_LIBRARY_TYPES = ("WES", "Panel_seq")
#: fgumi read structures of the FASTQ files of a read group, space separated
READ_STRUCTURE_PATTERN = r"^((\d+[TMBCS])*(\d+|\+)[TMBCS])( (\d+[TMBCS])*(\d+|\+)[TMBCS])*$"
HPO_TERMS_PATTERN = r"^HP:\d{7}(,HP:\d{7})*$"
#: Column names: usable as is in a pandas query
COLUMN_NAME = re.compile(r"^[a-z][a-z0-9_]*$")
#: Columns of the library table that are derived from the sheet, so a sheet cannot have them
DERIVED_COLUMNS = frozenset(
    {"extraction_type", "role", "father_name", "mother_name", "disease_state", "is_background_of"}
)
#: The columns every library row has
CORE_COLUMNS = ("donor_name", "library_name", "library_type")
#: Problems reported at most
MAX_PROBLEMS = 50


#: Column -> lower-case synonym -> canonical value (:func:`canonical`); the canonical values
#: themselves are accepted in any case
SYNONYMS: Mapping[str, Mapping[str, str]] = {
    "sex": {
        "female": "female", "f": "female", "woman": "female", "xx": "female", "2": "female",
        "male": "male", "m": "male", "man": "male", "xy": "male", "1": "male",
        "unknown": "unknown", "u": "unknown", "?": "unknown", "na": "unknown", "n/a": "unknown",
        "none": "unknown", "0": "unknown", "-9": "unknown",
    },
    "library_type": {
        "wes": "WES", "exome": "WES", "wxs": "WES",
        "wgs": "WGS", "genome": "WGS",
        "panel_seq": "Panel_seq", "panel": "Panel_seq", "panel-seq": "Panel_seq",
        "targeted": "Panel_seq",
        "mrna_seq": "mRNA_seq", "mrna-seq": "mRNA_seq", "mrnaseq": "mRNA_seq", "mrna": "mRNA_seq",
        "total_rna_seq": "total_RNA_seq", "total-rna-seq": "total_RNA_seq",
        "totalrnaseq": "total_RNA_seq", "total_rna": "total_RNA_seq",
    },
    "seq_platform": {
        "illumina": "Illumina",
        "pacbio": "PacBio", "pacbio_smrt": "PacBio",
        "ont": "ONT", "oxford_nanopore": "ONT", "nanopore": "ONT",
        "other": "other",
    },
}  # fmt: skip
#: Lower-case synonyms of the booleans
BOOLEANS = {
    "true": "true", "t": "true", "yes": "true", "y": "true", "1": "true",
    "false": "false", "f": "false", "no": "false", "n": "false", "0": "false",
}  # fmt: skip
#: Boolean columns whose digits are PED codes with another meaning, so no digits are read
NO_DIGIT_BOOLEANS = frozenset({"is_affected"})


def canonical(column: str, cell: str) -> str | None:
    """The canonical value of ``cell`` in ``column``, ``None`` when it is no known value."""
    key = cell.strip().lower()
    if column in BOOL_COLUMNS:
        if column in NO_DIGIT_BOOLEANS and key.isdigit():
            return None
        return BOOLEANS.get(key)
    return SYNONYMS.get(column, {}).get(key)


class Level(enum.StrEnum):
    """The key whose rows must agree on a column: all rows of a donor, a sample, an extract or a
    library carry the same value."""

    DONOR = "donor"
    SAMPLE = "sample"
    EXTRACT = "extract"
    LIBRARY = "library"


@dataclass(frozen=True)
class Column:
    """What a known sheet column is: its level, and the ISA-Tab header it maps to (the node of its
    level: Source, Sample, the Extract of the extraction or the Extract of the library)."""

    level: Level
    isa: str = ""


class SheetRow(BaseModel):
    """One row of a sample sheet: a library with its sample and donor, or a donor without
    libraries (every sample, extract and library column empty)."""

    model_config = ConfigDict(extra="allow", frozen=True)

    donor_name: typing.Annotated[Name, Column(Level.DONOR, "Source Name")]
    family_id: typing.Annotated[Name | None, Column(Level.DONOR, "Characteristics[Family]")] = None
    father: typing.Annotated[Name | None, Column(Level.DONOR, "Characteristics[Father]")] = None
    mother: typing.Annotated[Name | None, Column(Level.DONOR, "Characteristics[Mother]")] = None
    sex: typing.Annotated[Sex | None, Column(Level.DONOR, "Characteristics[Sex]")] = None
    is_affected: typing.Annotated[
        bool | None, Column(Level.DONOR, "Characteristics[Disease status]")
    ] = None
    is_index: typing.Annotated[bool | None, Column(Level.DONOR, "Comment[Index]")] = None
    hpo_terms: typing.Annotated[
        str | None,
        StringConstraints(pattern=HPO_TERMS_PATTERN),
        Column(Level.DONOR, "Characteristics[HPO terms]"),
    ] = None
    sample_name: typing.Annotated[Name | None, Column(Level.SAMPLE, "Sample Name")] = None
    is_tumor: typing.Annotated[bool | None, Column(Level.SAMPLE, "Characteristics[Is tumor]")] = (
        None
    )
    purity: typing.Annotated[
        float | None, Field(ge=0, le=1), Column(Level.SAMPLE, "Characteristics[Tumor purity]")
    ] = None
    is_ffpe: typing.Annotated[bool | None, Column(Level.SAMPLE, "Comment[FFPE]")] = None
    extract_name: typing.Annotated[Name | None, Column(Level.EXTRACT, "Extract Name")] = None
    library_name: typing.Annotated[Name | None, Column(Level.LIBRARY, "Extract Name")] = None
    library_type: typing.Annotated[
        LibraryType | None, Column(Level.LIBRARY, "Parameter Value[Library strategy]")
    ] = None
    library_kit: typing.Annotated[
        str | None, Column(Level.LIBRARY, "Parameter Value[Library kit]")
    ] = None
    folder_name: typing.Annotated[
        str | None,
        StringConstraints(pattern=r"^[^/]+$"),
        Column(Level.LIBRARY, "Characteristics[Folder name]"),
    ] = None
    seq_platform: typing.Annotated[
        SeqPlatform | None, Column(Level.LIBRARY, "Parameter Value[Platform]")
    ] = None
    read_structure: typing.Annotated[
        str | None,
        StringConstraints(pattern=READ_STRUCTURE_PATTERN),
        Column(Level.LIBRARY, "Comment[Read structure]"),
    ] = None
    batch_no: typing.Annotated[str | None, Column(Level.LIBRARY, "Characteristics[Batch]")] = None
    is_primary: typing.Annotated[bool | None, Column(Level.LIBRARY, "Comment[Primary]")] = None
    is_background: typing.Annotated[bool | None, Column(Level.LIBRARY, "Comment[Background]")] = (
        None
    )

    @property
    def is_library(self) -> bool:
        """Whether the row is a library, not a donor without libraries."""
        return self.library_name is not None

    @property
    def extraction_type(self) -> str | None:
        return EXTRACTION_TYPES.get(self.library_type or "")

    @property
    def role(self) -> str | None:
        """``tumor`` or ``normal`` for a library row."""
        if not self.is_library:
            return None
        return "tumor" if self.is_tumor else "normal"

    def extra(self, name: str):
        """The value of an extra (unknown) column, ``None`` when empty."""
        return (self.model_extra or {}).get(name)


def _column_of(name: str) -> Column:
    return next(m for m in SheetRow.model_fields[name].metadata if isinstance(m, Column))


#: The known columns, in the order of the model
COLUMNS: Mapping[str, Column] = {name: _column_of(name) for name in SheetRow.model_fields}


def _is_bool(annotation) -> bool:
    if annotation is bool:
        return True
    if typing.get_origin(annotation) in (typing.Union, types.UnionType):
        return bool in typing.get_args(annotation)
    return False


BOOL_COLUMNS = frozenset(
    name for name, f in SheetRow.model_fields.items() if _is_bool(f.annotation)
)


def header(columns: Iterable[str] = ()) -> str:
    """The header line of a sheet: the core columns, then ``columns`` in model order."""
    wanted = set(columns) | set(CORE_COLUMNS)
    return "\t".join(name for name in COLUMNS if name in wanted) + "\n"


@dataclass(frozen=True)
class Problem:
    """One problem of a sheet: ``path:line: column: message``."""

    path: str
    line: int
    column: str
    message: str

    def __str__(self) -> str:
        where = f"{self.path}:{self.line}" if self.line else self.path
        return (
            f"{where}: {self.column}: {self.message}" if self.column else f"{where}: {self.message}"
        )


class SheetError(ConfigError):
    """The problems of the sample sheets, all at once."""

    def __init__(self, problems: list[Problem]):
        self.problems = sorted(problems, key=lambda p: (p.path, p.line, p.column))
        shown = self.problems[:MAX_PROBLEMS]
        lines = [str(problem) for problem in shown]
        if len(self.problems) > MAX_PROBLEMS:
            lines.append(f"... and {len(self.problems) - MAX_PROBLEMS} more")
        count = f"{len(self.problems)} problem" + ("s" if len(self.problems) != 1 else "")
        super().__init__(f"The sample sheets have {count}:\n  " + "\n  ".join(lines))


@dataclass(frozen=True)
class Row:
    """A validated row with where it came from."""

    path: str
    line: int
    values: SheetRow


@dataclass(frozen=True)
class Sheet:
    """A sample sheet: its header, in file order, its rows, and the notes about cells read as a
    synonym of their canonical value."""

    path: str
    columns: tuple[str, ...]
    rows: tuple[Row, ...] = ()
    notes: tuple[Problem, ...] = ()

    @property
    def libraries(self) -> list[Row]:
        return [row for row in self.rows if row.values.is_library]


@dataclass
class Sheets:
    """The sample sheets of a project, checked together."""

    sheets: tuple[Sheet, ...] = ()
    problems: list[Problem] = field(default_factory=list)

    @property
    def rows(self) -> list[Row]:
        return [row for sheet in self.sheets for row in sheet.rows]

    @property
    def libraries(self) -> list[Row]:
        return [row for row in self.rows if row.values.is_library]

    @property
    def columns(self) -> list[str]:
        """The columns of all sheets, known ones first in model order, then the extras."""
        present = {name for sheet in self.sheets for name in sheet.columns}
        known = [name for name in COLUMNS if name in present]
        return known + sorted(present - set(COLUMNS))

    @property
    def notes(self) -> list[Problem]:
        """The cells of all sheets read as a synonym (``snappy sheet check --fix`` rewrites them)."""
        return [note for sheet in self.sheets for note in sheet.notes]


def read_sheets(paths: Iterable[str]) -> Sheets:
    """Read and check the sample sheets ``paths`` together; raise :class:`SheetError` with every
    problem of every sheet."""
    sheets, problems = [], []
    for path in paths:
        sheet, sheet_problems = read_sheet(path)
        problems += sheet_problems
        if sheet is not None:
            sheets.append(sheet)
    result = Sheets(tuple(sheets), problems)
    if not problems:
        result.problems += check_rows(result.rows)
    if result.problems:
        raise SheetError(result.problems)
    return result


def read_sheet(path: str) -> tuple[Sheet | None, list[Problem]]:
    """Read one sheet: its rows, validated one by one; the problems of its header and rows.
    Returns no sheet for a header problem."""
    try:
        with open(path, encoding="utf-8-sig", newline="") as f:
            lines = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    except OSError as e:
        return None, [Problem(path, 0, "", f"cannot read the sample sheet: {e.strerror}")]
    return parse_sheet(lines, path)


def parse_sheet(lines: list[list[str]], path: str) -> tuple[Sheet | None, list[Problem]]:
    """Parse the cells of the lines of a sheet (see :func:`read_sheet`)."""
    numbered = [
        (number, cells)
        for number, cells in enumerate(lines, start=1)
        if cells and any(cells) and not cells[0].lstrip().startswith("#")
    ]
    if not numbered:
        return None, [Problem(path, 0, "", "the sample sheet has no header line")]
    header_line, header = numbered[0]
    problems = _check_header(path, header_line, header)
    if problems:
        return None, problems
    rows, notes = [], []
    for number, cells in numbered[1:]:
        if len(cells) != len(header):
            problems.append(
                Problem(
                    path,
                    number,
                    "",
                    f"{len(cells)} cells, the header has {len(header)}; a missing value is an "
                    "empty cell (some editors drop trailing tabs)",
                )
            )
            continue
        values, cell_problems, cell_notes = _convert_cells(path, number, header, cells)
        problems += cell_problems
        notes += cell_notes
        if cell_problems:
            continue
        try:
            row = SheetRow.model_validate(values)
        except ValidationError as e:
            problems += [
                Problem(path, number, str(error["loc"][0]) if error["loc"] else "", error["msg"])
                for error in e.errors()
            ]
            continue
        problems += _check_row(path, number, row)
        rows.append(Row(path, number, row))
    return Sheet(path, tuple(header), tuple(rows), tuple(notes)), problems


def _check_header(path: str, line: int, header: list[str]) -> list[Problem]:
    problems = []
    for name in header:
        if not COLUMN_NAME.match(name) or keyword.iskeyword(name):
            problems.append(
                Problem(path, line, name, "column names are snake_case: ^[a-z][a-z0-9_]*$")
            )
        elif name in DERIVED_COLUMNS:
            problems.append(
                Problem(path, line, name, "the library table derives this column; rename it")
            )
    for name in sorted({name for name in header if header.count(name) > 1}):
        problems.append(Problem(path, line, name, "the column is repeated"))
    for name in CORE_COLUMNS:
        if name not in header:
            problems.append(Problem(path, line, name, "the column is required"))
    return problems


def _convert_cells(
    path: str, line: int, header: list[str], cells: list[str]
) -> tuple[dict[str, object], list[Problem], list[Problem]]:
    """The values of a row; its problems; its notes (cells read as a synonym)."""
    values: dict[str, object] = {}
    problems, notes = [], []
    for name, cell in zip(header, cells):
        if cell != cell.strip():
            problems.append(
                Problem(path, line, name, f"{cell!r} has leading or trailing whitespace")
            )
            continue
        if cell == ".":
            problems.append(Problem(path, line, name, "'.' is no value; leave the cell empty"))
            continue
        if cell == "":
            values[name] = None
            continue
        if name in BOOL_COLUMNS or name in SYNONYMS:
            value = canonical(name, cell)
            if value is None:
                hint = (
                    "; the PED codes 1 and 2 are ambiguous, write true or false"
                    if name in NO_DIGIT_BOOLEANS and cell.isdigit()
                    else ""
                )
                problems.append(Problem(path, line, name, f"{cell!r} is no known value{hint}"))
                continue
            if value != cell:
                notes.append(Problem(path, line, name, f"{cell!r} is read as {value!r}"))
            values[name] = value == "true" if name in BOOL_COLUMNS else value
        else:
            values[name] = cell
    return values, problems, notes


def fix_sheet(path: str) -> int:
    """Rewrite the cells of ``path`` that are synonyms with their canonical values; return how
    many. Everything else stays as it is, comments and blank lines included."""
    sheet, problems = read_sheet(path)
    if sheet is None or not sheet.notes:
        return 0
    with open(path, encoding="utf-8-sig", newline="") as f:
        text = f.read()
    newline = "\r\n" if "\r\n" in text else "\n"
    lines = text.split(newline)
    by_line: dict[int, list[Problem]] = {}
    for note in sheet.notes:
        by_line.setdefault(note.line, []).append(note)
    for number, notes in by_line.items():
        cells = lines[number - 1].split("\t")
        for note in notes:
            index = sheet.columns.index(note.column)
            cells[index] = canonical(note.column, cells[index])
        lines[number - 1] = "\t".join(cells)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(newline.join(lines))
    return len(sheet.notes)


#: The columns a donor-only row leaves empty: everything below the donor
_BELOW_DONOR = [name for name, column in COLUMNS.items() if column.level != Level.DONOR]


def _check_row(path: str, line: int, row: SheetRow) -> list[Problem]:
    problems = []
    if row.is_library:
        if row.library_type is None:
            problems.append(Problem(path, line, "library_type", "a library needs its type"))
        if row.purity is not None and not row.is_tumor:
            problems.append(Problem(path, line, "purity", "only a tumor sample has a purity"))
    else:
        filled = [name for name in _BELOW_DONOR if getattr(row, name) is not None]
        filled += sorted(name for name, value in (row.model_extra or {}).items() if value)
        if filled:
            problems.append(
                Problem(
                    path,
                    line,
                    "library_name",
                    "a row without a library describes a donor; it leaves "
                    f"{', '.join(filled)} empty",
                )
            )
    return problems


def check_rows(rows: list[Row]) -> list[Problem]:
    """The checks across rows (of all sheets together): unique names, consistent levels, one
    donor row or library rows per donor, parents, families and their index, the primary
    libraries, and at least one library."""
    problems: list[Problem] = []
    problems += _check_unique(rows, "library_name")
    problems += _check_unique(rows, "folder_name")
    problems += _check_unique(rows, "sample_name", key="donor_name", what="donor")
    problems += _check_unique(rows, "extract_name", key="sample_name", what="sample")
    problems += _check_levels(rows)
    problems += _check_donor_rows(rows)
    problems += _check_parents(rows)
    problems += _check_index(rows)
    problems += _check_primary(rows)
    if rows and not any(row.values.is_library for row in rows):
        problems.append(
            Problem(rows[0].path, 0, "library_name", "the sample sheets have no library")
        )
    return problems


def _check_unique(rows, column, key=None, what="") -> list[Problem]:
    """Each value of ``column`` names one thing: one row (no ``key``), or one value of ``key``
    (a sample belongs to one donor, an extract to one sample)."""
    problems = []
    seen: dict[str, Row] = {}
    for row in rows:
        value = getattr(row.values, column)
        if value is None:
            continue
        if value not in seen:
            seen[value] = row
            continue
        first = seen[value]
        if key is None:
            problems.append(
                Problem(row.path, row.line, column, f"{value!r} is also on line {first.line}")
            )
        elif getattr(row.values, key) != getattr(first.values, key):
            problems.append(
                Problem(
                    row.path,
                    row.line,
                    column,
                    f"{value!r} belongs to {what} {getattr(first.values, key)!r} on line "
                    f"{first.line}, here to {getattr(row.values, key)!r}",
                )
            )
    return problems


def _check_levels(rows) -> list[Problem]:
    """The rows of one donor, sample or extract agree on the columns of that level."""
    problems = []
    keys = {Level.DONOR: "donor_name", Level.SAMPLE: "sample_name", Level.EXTRACT: "extract_name"}
    for level, key in keys.items():
        columns = [name for name, column in COLUMNS.items() if column.level == level]
        first: dict[str, Row] = {}
        for row in rows:
            value = getattr(row.values, key)
            if value is None:
                continue
            if value not in first:
                first[value] = row
                continue
            for name in columns:
                if getattr(row.values, name) != getattr(first[value].values, name):
                    problems.append(
                        Problem(
                            row.path,
                            row.line,
                            name,
                            f"{getattr(row.values, name)!r} differs from "
                            f"{getattr(first[value].values, name)!r} on line {first[value].line} "
                            f"for {level} {value!r}",
                        )
                    )
    return problems


def _check_donor_rows(rows) -> list[Problem]:
    """A donor has library rows, or exactly one row without a library."""
    problems = []
    by_donor: dict[str, list[Row]] = {}
    for row in rows:
        by_donor.setdefault(row.values.donor_name, []).append(row)
    for donor, donor_rows in by_donor.items():
        donor_only = [row for row in donor_rows if not row.values.is_library]
        if donor_only and len(donor_rows) > 1:
            row = donor_only[-1]
            problems.append(
                Problem(
                    row.path,
                    row.line,
                    "library_name",
                    f"donor {donor!r} has a row without a library and other rows; a donor has "
                    "library rows or one row without",
                )
            )
    return problems


def _donors(rows) -> dict[str, Row]:
    """Donor -> its first row."""
    donors: dict[str, Row] = {}
    for row in rows:
        donors.setdefault(row.values.donor_name, row)
    return donors


def _check_parents(rows) -> list[Problem]:
    problems = []
    donors = _donors(rows)
    expected_sex = {"father": "female", "mother": "male"}  # the sex a parent must not have
    as_father = {d.values.father for d in donors.values() if d.values.father}
    as_mother = {d.values.mother for d in donors.values() if d.values.mother}
    for donor, row in donors.items():
        for column in ("father", "mother"):
            parent = getattr(row.values, column)
            if parent is None:
                continue
            if parent == donor:
                problems.append(
                    Problem(row.path, row.line, column, "a donor is not its own parent")
                )
                continue
            if parent not in donors:
                problems.append(
                    Problem(
                        row.path,
                        row.line,
                        column,
                        f"{parent!r} is no donor of the sheets; a parent needs a row, with or "
                        "without a library"
                        + (" (0 is no value; leave the cell empty)" if parent == "0" else ""),
                    )
                )
                continue
            parent_row = donors[parent]
            if parent_row.values.sex == expected_sex[column]:
                problems.append(
                    Problem(
                        row.path,
                        row.line,
                        column,
                        f"{parent!r} is {parent_row.values.sex} (line {parent_row.line})",
                    )
                )
            if row.values.family_id != parent_row.values.family_id:
                problems.append(
                    Problem(
                        row.path,
                        row.line,
                        column,
                        f"{parent!r} is in family {parent_row.values.family_id!r} (line "
                        f"{parent_row.line}), not in {row.values.family_id!r}",
                    )
                )
    for both in sorted(as_father & as_mother):
        row = donors[both]
        problems.append(
            Problem(row.path, row.line, "donor_name", f"{both!r} is a father and a mother")
        )

    # Cycles: a donor among its own ancestors
    def parents_of(name: str) -> list[str]:
        values = donors[name].values
        return [p for p in (values.father, values.mother) if p in donors and p != name]

    for donor, row in donors.items():
        ancestors, queue = set(), parents_of(donor)
        while queue:
            ancestor = queue.pop()
            if ancestor in ancestors:
                continue
            ancestors.add(ancestor)
            queue += parents_of(ancestor)
        if donor in ancestors:
            problems.append(
                Problem(row.path, row.line, "father", f"the parents of {donor!r} form a cycle")
            )
            break
    return problems


def _check_index(rows) -> list[Problem]:
    """One index per family where several donors have a normal DNA library; the index has one."""
    problems = []
    germline = {
        row.values.donor_name
        for row in rows
        if row.values.is_library
        and row.values.role == "normal"
        and row.values.extraction_type == "dna"
    }
    members: dict[str, list[Row]] = {}
    indexes: dict[str, list[Row]] = {}
    for row in _donors(rows).values():
        family = row.values.family_id or row.values.donor_name
        if row.values.donor_name in germline:
            members.setdefault(family, []).append(row)
        if row.values.is_index:
            indexes.setdefault(family, []).append(row)
            if row.values.donor_name not in germline:
                problems.append(
                    Problem(
                        row.path,
                        row.line,
                        "is_index",
                        "the index of a family has a normal DNA library",
                    )
                )
    for family, marked in indexes.items():
        if len(marked) > 1:
            row = marked[-1]
            problems.append(
                Problem(
                    row.path,
                    row.line,
                    "is_index",
                    f"family {family!r} has {len(marked)} index donors: "
                    + ", ".join(r.values.donor_name for r in marked),
                )
            )
    for family, rows_ in members.items():
        if len(rows_) > 1 and family not in indexes:
            row = rows_[0]
            problems.append(
                Problem(
                    row.path,
                    row.line,
                    "is_index",
                    f"family {family!r} has {len(rows_)} donors with a normal DNA library and no "
                    "index; mark one is_index",
                )
            )
    return problems


def _check_primary(rows) -> list[Problem]:
    """Among a donor's libraries of one role and extraction type, exactly one is primary when
    there are several: marked ``is_primary``."""
    problems = []
    groups: dict[tuple, list[Row]] = {}
    for row in rows:
        if row.values.is_library:
            key = (row.values.donor_name, row.values.role, row.values.extraction_type)
            groups.setdefault(key, []).append(row)
    for (donor, role, extraction), group in groups.items():
        marked = [row for row in group if row.values.is_primary]
        if len(group) > 1 and len(marked) != 1:
            row = group[-1]
            names = ", ".join(r.values.library_name for r in group)
            problems.append(
                Problem(
                    row.path,
                    row.line,
                    "is_primary",
                    f"donor {donor!r} has {len(group)} {role} {extraction} libraries ({names}); "
                    f"mark exactly one is_primary, {len(marked)} are marked",
                )
            )
    return problems


def primary_library(group: list[Row]) -> Row:
    """The primary library among ``group``, a donor's libraries of one role and extraction type
    (:func:`_check_primary` ensured one): the marked one, else the only one."""
    marked = [row for row in group if row.values.is_primary]
    return marked[0] if marked else group[0]
