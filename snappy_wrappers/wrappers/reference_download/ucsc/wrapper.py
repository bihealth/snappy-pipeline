"""Download a genome FASTA from UCSC"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper
from snappy_wrappers.tools import reference_download

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


def main() -> None:
    args = snakemake.params.args
    url = reference_download.ucsc_url(args["url"], args["db"], args["build"], args["file_name"])
    reference_download.fetch(
        url, str(snakemake.output.fasta), args["contigs"], args["contigs_regex"]
    )


if __name__ == "__main__":
    PythonWrapper(snakemake).run(main)
