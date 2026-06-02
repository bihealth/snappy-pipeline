from typing import Any
from collections.abc import Iterable

import pandas as pd

from biomedsheets.models import NGSLibrary, Sheet


def sample_sheets(sheets: list[Sheet], filter_background: bool = True) -> pd.DataFrame:
    """Creates a pandas data frame from snappy's list of samples sheets

    :param sheets: list of sheets provided by the abstract BaseStep class.
       The list **MUST** be ``self.sheets`` and **NOT** ``self.shortcut_sheets``,
       (because the latter diverges for cancer, germline & generic sheets).
    :returns: a pandas data frame containing all extra info columns.
       The data frame is guaranteed to have at least 4 columns with their entry names
       (``bio_entity``, ``bio_sample``, ``test_sample`` & ``ngs_library``).
       The data frame is indexed by the ngs library name.
       No duplicated rows nor duplicated ngs library names are allowed.
       Duplicate labels (column names) for extra info are also forbidden.
    """
    table: pd.DataFrame = None

    for sheet in sheets:
        if filter_background and sheet.extra_infos.get("is_background", False):
            continue
        for bio_entity in sheet.bio_entities.values():
            if bio_entity.disabled:
                continue
            for bio_sample in bio_entity.bio_samples.values():
                if bio_sample.disabled:
                    continue
                for test_sample in bio_sample.test_samples.values():
                    if test_sample.disabled:
                        continue
                    for ngs_library in test_sample.ngs_libraries.values():
                        if ngs_library.disabled:
                            continue
                        d = _ngs_library_to_df(ngs_library)
                        table = pd.concat([table, d], axis=0, ignore_index=True)

    assert not any(table.duplicated()), "Duplicated entries in sample sheets"
    assert not any(table["ngs_library"].duplicated()), "Duplicated NGS libraries"

    # table.set_index("ngs_library", drop=False, inplace=True)
    return table


def _ngs_library_to_df(ngs_library: NGSLibrary) -> pd.DataFrame:
    test_sample = ngs_library.test_sample
    bio_sample = test_sample.bio_sample
    bio_entity = bio_sample.bio_entity

    d = {
        "bio_entity": bio_entity.name,
        "bio_sample": bio_sample.name,
        "test_sample": test_sample.name,
        "ngs_library": ngs_library.name,
    }

    for o in (bio_entity, bio_sample, test_sample, ngs_library):
        extra_infos = getattr(o, "extra_infos")
        for k, v in extra_infos.items():
            assert k not in d, f"Extra info '{k}' already present elsewhere in {ngs_library.name}"
            d[k] = v

    return pd.DataFrame.from_dict({k: [v] for k, v in d.items()})


def filter_table(table: pd.DataFrame, filters: dict[str, Any]) -> pd.DataFrame:
    for expected_column, filter_value in filters.items():
        assert expected_column in table.columns, f"Missing mandatory column '{expected_column}'"
        if isinstance(filter_value, Iterable):
            table = table.loc[table[expected_column] in filter_value]
        else:
            table = table.loc[table[expected_column] == filter_value]
    return table


def filter_table_by_modality(table: pd.DataFrame, modality: str = "dna") -> pd.DataFrame:
    assert "extractionType" in table.columns, "Missing mandatory column 'extractionType'"
    return table.loc[table["extractionType"].apply(lambda x: x.lower() == modality)]


def tumor_to_normal_mapping(table: pd.DataFrame) -> dict[str, str]:
    assert "isTumor" in table.columns, "Missing mandatory column 'isTumor'"
    normals = table[~table["isTumor"]]
    assert all([n == 1 for n in normals.groupby(by="bio_entity").size()]), (
        "Multiple normals for at least one donor"
    )
    tumors = table[table["isTumor"]]
    tumor_normal_map = tumors[["ngs_library", "bio_entity"]].merge(
        normals[["ngs_library", "bio_entity"]], on="bio_entity"
    )
    return pd.Series(
        tumor_normal_map.ngs_library_y.values, index=tumor_normal_map.ngs_library_x.values
    ).to_dict()
