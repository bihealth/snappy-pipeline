import gzip
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from typing import TYPE_CHECKING
from urllib.parse import urlparse
from urllib.request import url2pathname

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


def _md5sum(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as f:  # noqa: S324
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_md5(path: str, md5_path: str) -> None:
    base = os.path.basename(path)
    with open(md5_path, "wt") as f:
        f.write(f"{_md5sum(path)}  {base}\n")


def _run_conda_cmd(cmd: list[str], path_out: str) -> None:
    with open(path_out, "wt") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False)


def _url_exists(url: str) -> bool:
    parsed = urlparse(url)
    if parsed.scheme == "file":
        local_path = url2pathname(parsed.path)
        return os.path.exists(local_path)
    cmd = ["curl", "--location", "--head", "--silent", "--show-error", "--fail", url]
    return subprocess.run(cmd, check=False).returncode == 0


def _download_with_curl(url: str, path_out: str) -> str:
    candidates = [url]
    if url.startswith("https://"):
        candidates.append(url.replace("https://", "ftp://", 1))

    for candidate in candidates:
        if not _url_exists(candidate):
            continue

        parsed = urlparse(candidate)
        if parsed.scheme == "file":
            local_path = url2pathname(parsed.path)
            shutil.copy(local_path, path_out)
            return candidate

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
            candidate,
        ]
        proc = subprocess.run(cmd, check=False, text=True, capture_output=True)
        if proc.returncode == 0:
            return candidate

    candidate_list = "\n".join(f"- {u}" for u in candidates)
    raise RuntimeError(
        "Unable to download/copy reference from Ensembl. Tried:\n"
        f"{candidate_list}\n"
        "Check species/build/release/datatype and server availability."
    )


def _filter_fasta(path_in: str, path_out: str, contigs: list[str], contigs_regex: str) -> None:
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


def _resolve_ensembl_prefix(params: dict) -> str:
    if params.get("url"):
        return params["url"]

    species = params["species"].lower()
    release = int(params["version"])
    build = params["build"]
    datatype = params.get("datatype", "dna")

    branch = ""
    if release >= 81 and build == "GRCh37":
        branch = "grch37/"
    elif params.get("branch"):
        branch = params["branch"].strip("/") + "/"

    spec = build if release > 75 else f"{build}.{release}"
    # Ensembl FASTA files use Genus capitalized and remaining parts lowercase.
    species_cap = species.capitalize()
    base_url = "https://ftp.ensembl.org/pub"
    return f"{base_url}/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}"


def _suffixes_for_datatype(params: dict) -> list[str]:
    datatype = params.get("datatype", "dna")
    chromosomes = list(params.get("chromosome") or [])
    subset = params.get("subset") or params.get("assembly") or "primary_assembly"

    if params.get("file_name"):
        return [params["file_name"]]

    if datatype == "dna":
        if chromosomes:
            return [f"dna.chromosome.{chrom}.fa.gz" for chrom in chromosomes]
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


def _symlink_output(work_path: str, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if os.path.lexists(out_path):
        os.remove(out_path)
    target = os.path.relpath(work_path, start=os.path.dirname(out_path))
    os.symlink(target, out_path)


def main() -> None:
    params = dict(snakemake.params.args)
    log = snakemake.log

    os.makedirs(os.path.dirname(log.log), exist_ok=True)
    _run_conda_cmd(["conda", "info"], str(log.conda_info))
    _run_conda_cmd(["conda", "list"], str(log.conda_list))

    contigs = list(params.get("contigs") or [])
    contigs_regex = params.get("contigs_regex") or ""
    chromosomes = list(params.get("chromosome") or [])

    if chromosomes and params.get("datatype", "dna") != "dna":
        raise ValueError("Ensembl chromosome selection requires datatype='dna'")

    os.makedirs(os.path.dirname(str(snakemake.output.fasta)), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        path_raw = os.path.join(tmpdir, "raw.fa")
        downloaded_urls = []
        success = False

        if params.get("url"):
            # Single processing path for explicit URL overrides (e.g., local files)
            url = params["url"]
            path_dl = os.path.join(
                tmpdir, "download.fa.gz" if url.endswith(".gz") else "download.fa"
            )
            try:
                used_url = _download_with_curl(url, path_dl)
                success = True
                downloaded_urls.append(used_url)
                if path_dl.endswith(".gz"):
                    with gzip.open(path_dl, "rt") as fin, open(path_raw, "wt") as fout:
                        shutil.copyfileobj(fin, fout)
                else:
                    shutil.copy(path_dl, path_raw)
            except Exception as e:
                raise RuntimeError(f"Unable to copy or download custom URL: {url}. Error: {e}")
        else:
            # Construction-based loop
            prefix_or_url = _resolve_ensembl_prefix(params)
            suffixes = _suffixes_for_datatype(params)
            with open(path_raw, "wt") as out_raw:
                for suffix in suffixes:
                    url = f"{prefix_or_url}.{suffix}"
                    path_dl = os.path.join(
                        tmpdir,
                        f"download_{len(downloaded_urls)}.fa.gz"
                        if url.endswith(".gz")
                        else "download.fa",
                    )
                    try:
                        used_url = _download_with_curl(url, path_dl)
                    except RuntimeError:
                        if chromosomes:
                            raise
                        continue

                    success = True
                    downloaded_urls.append(used_url)

                    if path_dl.endswith(".gz"):
                        with gzip.open(path_dl, "rt") as fin:
                            shutil.copyfileobj(fin, out_raw)
                    else:
                        with open(path_dl, "rt") as fin:
                            shutil.copyfileobj(fin, out_raw)

                    if not chromosomes:
                        break

        if not success:
            raise RuntimeError(
                "Unable to download requested Ensembl reference. "
                "Please verify species/build/release/datatype combination."
            )

        _filter_fasta(path_raw, str(snakemake.output.fasta), contigs, contigs_regex)

    _write_md5(str(snakemake.output.fasta), str(snakemake.output.fasta_md5))
    with open(log.log, "wt") as f:
        f.write("Downloaded URLs:\n")
        for item in downloaded_urls:
            f.write(f"- {item}\n")
        f.write(f"Contigs: {contigs}\n")
        f.write(f"Regex: {contigs_regex}\n")
    _write_md5(str(log.log), str(log.log_md5))
    _write_md5(str(log.conda_info), str(log.conda_info_md5))
    _write_md5(str(log.conda_list), str(log.conda_list_md5))

    for dst in snakemake.output.output_links:
        src = str(dst).replace("output/", "work/", 1)
        _symlink_output(src, str(dst))


if __name__ == "__main__":
    main()
