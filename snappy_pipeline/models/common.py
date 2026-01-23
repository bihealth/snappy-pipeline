import enum

from pydantic import model_validator

from biomedsheets.io_tsv.base import EXTRACTION_TYPE_DNA, EXTRACTION_TYPE_RNA, LIBRARY_TYPES


class ExtractionType(enum.StrEnum):
    DNA = EXTRACTION_TYPE_DNA
    RNA = EXTRACTION_TYPE_RNA


class LibraryType(enum.StrEnum):
    WES = "WES"
    WGS = "WGS"
    Panel = "Panel_seq"
    mRNA = "mRNA_seq"

    @model_validator(mode="after")
    def ensure_types_match_biomedsheets(self):
        if not all([x in LIBRARY_TYPES for x in list(self)]):
            raise ValueError("Library types incompatible with biomedsheets")
        return self
