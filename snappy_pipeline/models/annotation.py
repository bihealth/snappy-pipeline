import enum

from typing import Annotated

from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel


class VepTxFlag(enum.StrEnum):
    gencode_basic = "gencode_basic"
    refseq = "refseq"
    merged = "merged"


class VepPlugin(SnappyModel):
    name: str
    path: str | None = None
    url: Annotated[
        str | None, Field(examples=["https://github.com/Ensembl/VEP_plugins/<plugin_name>.pm"])
    ] = None

    @model_validator(mode="after")
    def ensure_name_and_path_or_url(self):
        if not self.name:
            raise ValueError("Missing plugin name")
        if not (self.path or self.url):
            raise ValueError(f"Either path or URL must be defined for plugin {self.name}")
        return self


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
        "mane_select",
        "mane_plus_clinical",
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
    plugins: list[VepPlugin] = []
    plugins_dir: str = ""
