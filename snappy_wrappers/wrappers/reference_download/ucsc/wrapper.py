import gzip
import os
import re
import shutil
import subprocess
import tempfile
from typing import TYPE_CHECKING
from urllib.parse import urlparse
from urllib.request import url2pathname

from snappy_wrappers.snappy_wrapper import PythonWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


def _download(url: str, path_out: str) -> None:
    parsed = urlparse(url)
    if parsed.scheme == "file":
        local_path = url2pathname(parsed.path)
        if not os.path.exists(local_path):
            raise FileNotFoundError(f"Local file not found: {local_path}")
        shutil.copy(local_path, path_out)
        return

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


def _resolve_ucsc_url(params: dict) -> str:
    if params.get("url"):
        return params["url"]

    db = params.get("db") or params.get("build")
    if not db:
        raise ValueError("ucsc download requires either explicit url or db/build")

    file_name = params.get("file_name") or f"{db}.fa.gz"
    return f"https://hgdownload.soe.ucsc.edu/goldenPath/{db}/bigZips/{file_name}"


def _symlink_output(work_path: str, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if os.path.lexists(out_path):
        os.remove(out_path)
    target = os.path.relpath(work_path, start=os.path.dirname(out_path))
    os.symlink(target, out_path)


def main() -> None:
    params = dict(snakemake.params.args)

    url = _resolve_ucsc_url(params)
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

    print(f"Downloaded: {url}")
    print(f"Contigs: {contigs}")
    print(f"Regex: {contigs_regex}")

    for dst in snakemake.output.output_links:
        src = str(dst).replace("output/", "work/", 1)
        _symlink_output(src, str(dst))


if __name__ == "__main__":
    PythonWrapper(snakemake).run(main)
