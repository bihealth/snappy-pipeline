import enum

from snappy_pipeline.models import SnappyModel


class VepTxFlag(enum.StrEnum):
    gencode_basic = "gencode_basic"
    refseq = "refseq"
    merged = "merged"


class Vep(SnappyModel):
    cache_dir: str = ""
    """Defaults to $HOME/.vep Not a good idea on the cluster"""

    species: str = "homo_sapiens"

    assembly: str = "GRCh38"

    cache_version: str = "102"
    """WARNING- this must match the wrapper's vep version!"""

    tx_flag: VepTxFlag = VepTxFlag.gencode_basic
    """The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged"."""

    pick_order: list[str] = [
        "biotype",
        "mane",
        "appris",
        "tsl",
        "ccds",
        "canonical",
        "rank",
        "length",
    ]
    num_threads: int = 8
    buffer_size: int = 1000
    output_options: list[str] = ["everything"]
