import pandas as pd

from biomedsheets.models import NGSLibrary, Sheet


def sample_sheets(sheets: list[Sheet]) -> pd.DataFrame:
    """Creates a pandas data frame from snappy's list of samples sheets

    :param sheets: list of sheets provided by the abstract BaseStep class.
       The list **MUST** be ``self.sheets`` and **NOT** ``self.shortcut_sheets``,
       (because the latter diverges for cancer, germline & generic sheets).
    :returns: a pandas data frame containing all extra info columns.
       The data frame is guaranteed to have at least 4 columns with their entry names
       (``bio_entity``, ``bio_sample``, ``test_sample`` & ``ngs_library``).
       The data frame is indexed by the ngs library name.
       No duplicated rows not duplicated ngs librari names are allowed.
       Duplicate lables for extra info are also forbidden.
    """
    table: pd.DataFrame = None

    for sheet in sheets:
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

    table.set_index("ngs_library", drop=False, inplace=True)
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
