"""Download a genome FASTA from Ensembl"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper
from snappy_wrappers.tools import reference_download

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


def main() -> None:
    args = snakemake.params.args
    reference_download.fetch_ensembl(
        str(snakemake.output.fasta),
        contigs=args["contigs"],
        contigs_regex=args["contigs_regex"],
        url=args["url"],
        species=args["species"],
        version=args["version"],
        build=args["build"],
        datatype=args["datatype"],
        branch=args["branch"],
        chromosome=args["chromosome"],
        subset=args["subset"],
        assembly=args["assembly"],
        file_name=args["file_name"],
    )


if __name__ == "__main__":
    PythonWrapper(snakemake).run(main)
