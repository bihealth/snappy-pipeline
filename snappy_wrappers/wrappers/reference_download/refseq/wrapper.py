import gzip
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from typing import TYPE_CHECKING

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


def _download(url: str, path_out: str) -> None:
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
    proc = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"curl download failed for {url}: {proc.stderr.strip()}")


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


def _resolve_refseq_url(params: dict) -> str:
    if params.get("url"):
        return params["url"]

    accession = params.get("assembly_accession") or ""
    assembly_name = params.get("assembly_name") or ""
    if not accession or not assembly_name:
        raise ValueError(
            "refseq download requires either explicit url or both assembly_accession and assembly_name"
        )

    prefix, rest = accession.split("_", 1)
    numeric = rest.split(".", 1)[0]
    chunks = [numeric[0:3], numeric[3:6], numeric[6:9]]
    dir_name = f"{accession}_{assembly_name}"
    file_name = params.get("file_name") or f"{dir_name}_genomic.fna.gz"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/"
        f"{chunks[0]}/{chunks[1]}/{chunks[2]}/{dir_name}/{file_name}"
    )


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

    url = _resolve_refseq_url(params)
    contigs = list(params.get("contigs") or [])
    contigs_regex = params.get("contigs_regex") or ""

    os.makedirs(os.path.dirname(str(snakemake.output.fasta)), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        path_dl = os.path.join(tmpdir, "download.fa.gz" if url.endswith(".gz") else "download.fa")
        path_raw = os.path.join(tmpdir, "raw.fa")
        _download(url, path_dl)

        if path_dl.endswith(".gz"):
            with gzip.open(path_dl, "rt") as fin, open(path_raw, "wt") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            shutil.copy(path_dl, path_raw)

        _filter_fasta(path_raw, str(snakemake.output.fasta), contigs, contigs_regex)

    _write_md5(str(snakemake.output.fasta), str(snakemake.output.fasta_md5))
    with open(log.log, "wt") as f:
        f.write(f"Downloaded: {url}\n")
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

