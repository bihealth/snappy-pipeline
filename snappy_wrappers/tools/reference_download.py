"""Download a genome FASTA from Ensembl, RefSeq or UCSC, or copy it from a ``file://`` URL.

The ``reference_download`` wrappers call ``fetch`` or ``fetch_ensembl``,
which keep the contigs that ``contigs`` and ``contigs_regex`` select.
"""

import gzip
import os
import re
import shutil
import subprocess
import tempfile
from urllib.parse import urlparse
from urllib.request import url2pathname


def curl(url: str, path_out: str) -> subprocess.CompletedProcess:
    """Download ``url`` to ``path_out``, with three retries."""
    cmd = [
        "curl",
        "--location",
        "--fail",
        "--silent",
        "--show-error",
        "--retry",
        "3",
        "--retry-all-errors",
        "--output",
        path_out,
        url,
    ]
    return subprocess.run(cmd, check=False, text=True, capture_output=True)


def download(url: str, path_out: str) -> None:
    """Copy a ``file://`` URL or download ``url`` to ``path_out``."""
    parsed = urlparse(url)
    if parsed.scheme == "file":
        local_path = url2pathname(parsed.path)
        if not os.path.exists(local_path):
            raise FileNotFoundError(f"Local file not found: {local_path}")
        shutil.copy(local_path, path_out)
        return

    proc = curl(url, path_out)
    if proc.returncode != 0:
        raise RuntimeError(f"curl download failed for {url}: {proc.stderr.strip()}")


def url_exists(url: str) -> bool:
    parsed = urlparse(url)
    if parsed.scheme == "file":
        local_path = url2pathname(parsed.path)
        return os.path.exists(local_path)
    cmd = ["curl", "--location", "--head", "--silent", "--show-error", "--fail", url]
    return subprocess.run(cmd, check=False).returncode == 0


def download_with_fallback(url: str, path_out: str) -> str:
    """Copy or download the first of ``url`` and its ``ftp://`` variant that exists.

    Return the URL used.
    """
    candidates = [url]
    if url.startswith("https://"):
        candidates.append(url.replace("https://", "ftp://", 1))

    for candidate in candidates:
        if not url_exists(candidate):
            continue

        parsed = urlparse(candidate)
        if parsed.scheme == "file":
            local_path = url2pathname(parsed.path)
            shutil.copy(local_path, path_out)
            return candidate

        if curl(candidate, path_out).returncode == 0:
            return candidate

    candidate_list = "\n".join(f"- {u}" for u in candidates)
    raise RuntimeError(
        "Unable to download/copy reference from Ensembl. Tried:\n"
        f"{candidate_list}\n"
        "Check species/build/release/datatype and server availability."
    )


def append_fasta(path: str, fout) -> None:
    """Append the FASTA at ``path``, gunzipped if it ends with ``.gz``, to the text file ``fout``."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fin:
        shutil.copyfileobj(fin, fout)


def filter_fasta(path_in: str, path_out: str, contigs: list, contigs_regex: str) -> None:
    """Write the sequences of ``contigs`` (all if empty) that match ``contigs_regex``."""
    contig_set = set(contigs)
    rx = re.compile(contigs_regex) if contigs_regex else None

    with open(path_in, "rt") as fin, open(path_out, "wt") as fout:
        keep = True
        for line in fin:
            if line.startswith(">"):
                contig = line[1:].split()[0]
                keep = True
                if contig_set:
                    keep = contig in contig_set
                if keep and rx is not None:
                    keep = bool(rx.search(contig))
            if keep:
                fout.write(line)


def refseq_url(url: str, assembly_accession: str, assembly_name: str, file_name: str) -> str:
    if url:
        return url

    if not assembly_accession or not assembly_name:
        raise ValueError(
            "refseq download requires either explicit url or both assembly_accession and assembly_name"
        )

    prefix, rest = assembly_accession.split("_", 1)
    numeric = rest.split(".", 1)[0]
    chunks = [numeric[0:3], numeric[3:6], numeric[6:9]]
    dir_name = f"{assembly_accession}_{assembly_name}"
    file_name = file_name or f"{dir_name}_genomic.fna.gz"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/"
        f"{chunks[0]}/{chunks[1]}/{chunks[2]}/{dir_name}/{file_name}"
    )


def ucsc_url(url: str, db: str, file_name: str) -> str:
    if url:
        return url

    if not db:
        raise ValueError("ucsc download requires either explicit url or db")

    file_name = file_name or f"{db}.fa.gz"
    return f"https://hgdownload.soe.ucsc.edu/goldenPath/{db}/bigZips/{file_name}"


def ensembl_prefix(species: str, version: str, build: str, datatype: str, branch: str) -> str:
    """The URL of the release's FASTA files without the suffix of ``ensembl_suffixes``."""
    species = species.lower()
    release = int(version)

    if release >= 81 and build == "GRCh37":
        branch = "grch37/"
    elif branch:
        branch = branch.strip("/") + "/"
    else:
        branch = ""

    spec = build if release > 75 else f"{build}.{release}"
    # Ensembl FASTA files use Genus capitalized and remaining parts lowercase.
    species_cap = species.capitalize()
    base_url = "https://ftp.ensembl.org/pub"
    return f"{base_url}/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}"


def ensembl_suffixes(
    datatype: str, chromosome: list, subset: str, assembly: str, file_name: str
) -> list:
    """The file name suffixes to try in turn, or to concatenate for ``chromosome``."""
    subset = subset or assembly or "primary_assembly"

    if file_name:
        return [file_name]

    if datatype == "dna":
        if chromosome:
            return [f"dna.chromosome.{chrom}.fa.gz" for chrom in chromosome]
        return [f"dna.{subset}.fa.gz", "dna.toplevel.fa.gz"]
    if datatype == "cdna":
        return ["cdna.all.fa.gz"]
    if datatype == "cds":
        return ["cds.all.fa.gz"]
    if datatype == "ncrna":
        return ["ncrna.fa.gz"]
    if datatype == "pep":
        return ["pep.all.fa.gz"]

    raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")


def fetch(url: str, fasta: str, contigs: list, contigs_regex: str) -> None:
    """Download ``url`` (RefSeq, UCSC) and write the selected contigs to ``fasta``."""
    contigs = list(contigs or [])
    contigs_regex = contigs_regex or ""

    os.makedirs(os.path.dirname(fasta), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        path_dl = os.path.join(tmpdir, "download.fa.gz" if url.endswith(".gz") else "download.fa")
        path_raw = os.path.join(tmpdir, "raw.fa")
        download(url, path_dl)
        with open(path_raw, "wt") as fout:
            append_fasta(path_dl, fout)

        filter_fasta(path_raw, fasta, contigs, contigs_regex)

    print(f"Downloaded: {url}")
    print(f"Contigs: {contigs}")
    print(f"Regex: {contigs_regex}")


def fetch_ensembl(
    fasta: str,
    contigs: list,
    contigs_regex: str,
    url: str,
    species: str,
    version: str,
    build: str,
    datatype: str,
    branch: str,
    chromosome: list,
    subset: str,
    assembly: str,
    file_name: str,
) -> None:
    """Download from Ensembl and write the selected contigs to ``fasta``.

    An explicit ``url`` replaces the Ensembl URL.
    Without ``chromosome``, the first suffix of ``ensembl_suffixes`` that exists wins;
    with it, the files of all chromosomes are concatenated.
    """
    contigs = list(contigs or [])
    contigs_regex = contigs_regex or ""
    chromosomes = list(chromosome or [])

    if chromosomes and datatype != "dna":
        raise ValueError("Ensembl chromosome selection requires datatype='dna'")

    os.makedirs(os.path.dirname(fasta), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        path_raw = os.path.join(tmpdir, "raw.fa")
        downloaded_urls = []
        success = False

        if url:
            # Single processing path for explicit URL overrides (e.g., local files)
            path_dl = os.path.join(
                tmpdir, "download.fa.gz" if url.endswith(".gz") else "download.fa"
            )
            try:
                used_url = download_with_fallback(url, path_dl)
                success = True
                downloaded_urls.append(used_url)
                with open(path_raw, "wt") as out_raw:
                    append_fasta(path_dl, out_raw)
            except Exception as e:
                raise RuntimeError(f"Unable to copy or download custom URL: {url}. Error: {e}")
        else:
            # Construction-based loop
            prefix = ensembl_prefix(species, version, build, datatype, branch)
            suffixes = ensembl_suffixes(datatype, chromosomes, subset, assembly, file_name)
            with open(path_raw, "wt") as out_raw:
                for suffix in suffixes:
                    url = f"{prefix}.{suffix}"
                    path_dl = os.path.join(
                        tmpdir,
                        f"download_{len(downloaded_urls)}.fa.gz"
                        if url.endswith(".gz")
                        else "download.fa",
                    )
                    try:
                        used_url = download_with_fallback(url, path_dl)
                    except RuntimeError:
                        if chromosomes:
                            raise
                        continue

                    success = True
                    downloaded_urls.append(used_url)
                    append_fasta(path_dl, out_raw)

                    if not chromosomes:
                        break

        if not success:
            raise RuntimeError(
                "Unable to download requested Ensembl reference. "
                "Please verify species/build/release/datatype combination."
            )

        filter_fasta(path_raw, fasta, contigs, contigs_regex)

    print("Downloaded URLs:")
    for item in downloaded_urls:
        print(f"- {item}")
    print(f"Contigs: {contigs}")
    print(f"Regex: {contigs_regex}")
