import argparse
import csv
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import time

from multiprocessing import Manager, Process, ProcessError, Queue, current_process
from pathlib import Path
from typing import Any

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

NETMHCSTABPAN_CMD = r"""
allele=$1
length=$2
fasta=$3
out=$4

export NETMHCpan={netMHCpan_path}
export NETMHCstabpan={netMHCstabpan_path}
export TMPDIR={tmpdir}

grep -q $allele $NETMHCstabpan/data/MHC_pseudo.dat
if [[ $? -eq 0 ]]
then
    $NETMHCstabpan/bin/netMHCstabpan -dirty \
        -a $allele -l $length \
        -xls -xlsfile $out \
        $fasta
else
    touch $out
fi
"""

HLA_CLEAN_PATTERN = re.compile(r"[:\*\-]+")
ALLELE_PROXIMITY_PATTERN = re.compile(
    r"^(.*?) : Distance to trai?ning data\s+(\d.\d+).*? nearest neighbor (.*?)\)$", re.MULTILINE
)
NETMHCSTABPAN_OUTPUT_PATTERN = re.compile(r"")
FILENAME_PATTERN = re.compile(r"^(?P<hla>[^_]+)_(?P<len>[0-9]+)\.(fasta|xls|head)$")
COMMA = re.compile(r",")


def _clean_hla(hla_type: str) -> str:
    return HLA_CLEAN_PATTERN.sub("", hla_type)


def _create_sequence_files(
    seq_dir: str,
    epitopes: list[dict[str, str]],
    hla_types: set[str],
    lengths: set[int],
    seq_name: str = "MT Epitope Seq",
    hla_name: str = "HLA Allele",
    len_name: str = "Peptide Length",
) -> list[str]:
    sequences = {}
    for i, epitope in enumerate(epitopes):
        hla_type = epitope[hla_name]
        epitope_length = len(epitope[seq_name])
        if hla_type in hla_types and epitope_length in lengths:
            k = f"{_clean_hla(epitope[hla_name])}_{str(epitope_length)}.fasta"
            if k not in sequences:
                sequences[k] = []
            sequences[k].append((i, epitope[seq_name]))

    file_list = []
    for fn, seqs in sequences.items():
        file_name = os.path.join(seq_dir, fn)
        with open(file_name, "wt") as f:
            for seq in seqs:
                f.write(f">seq_{seq[0]}\n")
                f.write(seq[1] + "\n")
        file_list.append(file_name)
    return file_list


def _parse_netMHCstabpan_output(out: io.StringIO):
    """Parses long netMHCstabpan output"""
    distance = 0.0
    replacement = ""
    for line in out:
        m = ALLELE_PROXIMITY_PATTERN.match(line.strip())
        if m:
            distance = float(m.group(2))
            replacement = m.group(3)
    return (replacement, distance)


def _parse_output_file(fn: str) -> list[dict[str, Any]]:
    results = []
    with open(fn, "rt") as f:
        hla_type = f.readline().strip()  # noqa: F841
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = {k: v.strip() for k, v in row.items()}
            results.append(row)
    return results


def _run_netMHCstabpan_command(cmd: list[str], worker_tmp: str, timeout: int = 3600):
    """
    Runs netMHCstabpan to find peptide stability score

    The sequences are first saved in a temp directory, and netMHCstabpan is run.
    The output is parsed and stability scores are stored in the epitope object.
    The site position and its score are retained.
    """
    logging.debug(f"Process {current_process().name}, cmd = {' '.join(cmd)}, timeout = {timeout}")

    # Compute stability scores
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    try:
        out, err = p.communicate(timeout=timeout)
    except TimeoutError:
        p.kill()
        raise (f"The command {' '.join(cmd)} has timed out")
    if p.returncode != 0:
        raise ChildProcessError(f"Command {' '.join(cmd)} failed with return code {p.returncode}")

    # Parse cleavage sites & put them into epitope object
    return _parse_netMHCstabpan_output(out.decode("utf-8").split("\n"))


def _worker(
    task_queue: Queue,
    error_queue: Queue,
    return_dict: dict[str, Any],
    cmd: list[str],
    worker_tmp: str,
    timeout: int = 3600,
):
    """Multi-processing intermediate for netchop"""
    logging.debug(
        f"Starting worker with cmd = {' '.join(cmd)}, tmpdir = {worker_tmp}, timeout = {timeout}"
    )
    while not task_queue.empty():
        hla_type, length, sequences, out = task_queue.get()
        logging.debug(
            f"Processing queued job with arguments {hla_type}, {length}, {sequences}, {out}"
        )
        try:
            return_dict[f"{_clean_hla(hla_type)}_{str(length)}.head"] = _run_netMHCstabpan_command(
                cmd + [hla_type, str(length), sequences, out], worker_tmp, timeout
            )
        except Exception as e:
            error_queue.put((e, hla_type, length))


def run_netMHCstabpan(
    epitopes: list[dict[str, str]],
    hla_types: set[str],
    lengths: set[int],
    base_path: str,
    container: str,
    netMHCstabpan_path: str,
    netMHCpan_path: str,
    epitope_seq_column_name: str = "MT Epitope Seq",
    n_workers: int = 1,
    timeout: int = 3600,
):
    """Runs netMHCstabpan within docker container"""
    base_path = os.path.realpath(base_path)
    temp_path = os.path.join(base_path, "tmp", "netMHCstabpan")
    container_path = "/base_path"
    if os.path.exists(temp_path):
        shutil.rmtree(temp_path)

    tmp_dir = "tmp"
    seq_dir = "seqs"
    out_dir = "out"
    for d in (tmp_dir, seq_dir, out_dir):
        os.makedirs(os.path.join(temp_path, d), mode=0o750, exist_ok=False)

    script = os.path.join(base_path, "scripts", "run_netMHCstabpan.sh")
    os.makedirs(os.path.dirname(script), mode=0o750, exist_ok=True)
    with open(script, "wt") as f:
        f.write(
            NETMHCSTABPAN_CMD.format(
                netMHCstabpan_path=os.path.dirname(os.path.dirname(netMHCstabpan_path)),
                netMHCpan_path=os.path.dirname(os.path.dirname(netMHCpan_path)),
                tmpdir=os.path.join(container_path, tmp_dir),
            )
        )

    cleaned_hla_types = {k: _clean_hla(k) for k in hla_types}
    assert len(set(cleaned_hla_types.values())) == len(hla_types), (
        f"HLA types {cleaned_hla_types} not unique after cleaning"
    )
    reverse_cleaned = {v: k for k, v in cleaned_hla_types.items()}

    file_list = _create_sequence_files(
        os.path.join(temp_path, seq_dir),
        epitopes,
        hla_types,
        lengths,
        seq_name=epitope_seq_column_name,
    )

    container = os.path.realpath(container)

    # Prepare netchop command
    cmd = [
        "apptainer",
        "run",
        "--no-home",
        "--bind",
        f"{temp_path}:{container_path}:rw",
        "--bind",
        f"{os.path.dirname(script)}:/scripts:ro",
        container,
        "bash",
        os.path.join("/scripts", os.path.basename(script)),
    ]
    logging.info(
        f"netMHCstabpan command: {' '.join(cmd + ['<hla_type>', '<length>', '<fasta>', '<xls>'])}"
    )

    # Prepare multiprocessing
    for i_worker in range(n_workers):
        worker_tmp = f"worker_{i_worker}"
        os.makedirs(os.path.join(temp_path, tmp_dir, worker_tmp), mode=0o700, exist_ok=False)

    task_queue = Queue()
    for file_name in file_list:
        m = FILENAME_PATTERN.match(os.path.basename(file_name))
        assert m, f"File {file_name} doesn't match expected pattern (hlatype_length)"
        hla_type = reverse_cleaned[m.group("hla")]
        length = int(m.group("len"))
        task_queue.put(
            (
                hla_type.replace("*", ""),
                int(length),
                os.path.join(container_path, seq_dir, os.path.basename(file_name)),
                os.path.join(
                    container_path,
                    out_dir,
                    cleaned_hla_types[hla_type] + "_" + str(length) + ".xls",
                ),
            )
        )
    error_queue = Queue()
    return_dict = Manager().dict()
    processes: list[Process] = []

    time.sleep(2.0)

    # Start the workers
    for i_worker in range(n_workers):
        worker_tmp = os.path.join(container_path, "tmp", f"worker_{i_worker}")
        p = Process(
            target=_worker, args=(task_queue, error_queue, return_dict, cmd, worker_tmp, timeout)
        )
        processes.append(p)
        logging.info(f"Starting worker {i_worker}")
        p.start()

    # Wait for completion
    for p in processes:
        p.join()
        logging.info(f"Worker {p.name} completed")

    # Check for errors
    error = False
    while not error_queue.empty():
        error = True
        e, hla_type, length = error_queue.get()
        logging.error(
            f"An error occurred during netMHCstabpan for {hla_type} epitopes of length {length} - message {e}"
        )
    if error:
        raise ProcessError("Error running one of netMHCstabpan sub-processes")

    for i in range(len(epitopes)):
        for k in ("Predicted Stability", "Half Life", "Stability Rank", "NetMHCstab allele"):
            epitopes[i][k] = "NA"

    # Rapatriate netMHCstabpan results into epitope
    for file_name in file_list:
        m = FILENAME_PATTERN.match(os.path.basename(file_name))
        assert m, f"File {file_name} doesn't match expected pattern (hlatype_length)"
        hla_type = reverse_cleaned[m.group("hla")]
        length = int(m.group("len"))
        k = f"{m.group('hla')}_{m.group('len')}"

        replacement = return_dict[k + ".head"]
        if replacement[0]:
            replacement = f"{replacement[0]} (distance: {float(replacement[1])})"
        else:
            replacement = f"{hla_type} (distance: {float(replacement[1])})"

        fn = os.path.join(temp_path, out_dir, k + ".xls")
        assert os.path.exists(fn), f"Can't find file {fn}"
        results = _parse_output_file(fn)
        for row in results:
            i = int(row["ID"][4:])
            for src, dest in (
                ("Pred", "Predicted Stability"),
                ("Thalf(h)", "Half Life"),
                ("Rank", "Stability Rank"),
            ):
                epitopes[i][dest] = row[src]
            epitopes[i]["NetMHCstab allele"] = replacement

    return epitopes


def read_epitopes_table(fn: str | Path) -> list[dict[str, str]]:
    """
    Reads neo-epitope table

    Produced by pVACseq (should work for pVACSplice and pVACfuse, but untested)
    The files are generally <sample>.<MHC class>.(all_epitopes|filtered).tsv.
    """
    records = []
    with open(fn, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert "Index" in reader.fieldnames, f"Index column missing from {fn}"
        for row in reader:
            records.append(row)
    return records


def write_output_table(
    f: io.TextIOBase,
    epitopes: list[dict[str, Any]],
):
    """Add 3 columns to the epitope table (Best cleavage pos & score, and digest of all cleavage sites)"""
    titles = epitopes[0].keys()
    f.write("\t".join(titles) + "\n")
    for epitope in epitopes:
        f.write("\t".join([epitope[k] for k in titles]) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="Wrapper around netchop",
        description="Runs netMHCstabpan on results from pVACtools modules",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verbosity level")
    parser.add_argument(
        "-w", "--workers", type=int, default=1, help="Number of threads for netMHCstabpan"
    )
    parser.add_argument(
        "--base-path", help="Base directory to create scripts & tmp/netMHCstabpan sub-directories"
    )

    parser.add_argument("--container", help="Path to the pVACtools container")
    parser.add_argument(
        "--netMHCstabpan",
        default="/opt/iedb/mhc_i/method/netmhcstabpan-1.0-executable/netmhcstabpan_1_0_executable/Linux_x86_64/bin/netMHCstabpan",
        help="Path to the netMHCstabpan binary in the docker container",
    )
    parser.add_argument(
        "--netMHCpan",
        default="/opt/iedb/mhc_i/method/netmhcpan-2.8-executable/netmhcpan_2_8_executable/Linux_x86_64/bin/netMHCpan",
        help="Path to the netMHCpan binary in the docker container",
    )

    parser.add_argument(
        "-t",
        "--tool",
        choices=("pvacseq", "pvacfuse", "pvacsplice"),
        default="pvacseq",
        help="pVACtool module (pvacseq, pvacfuse, pvacsplice)",
    )

    parser.add_argument("--hla-types", help="HLA types separated by a comma")
    parser.add_argument("--lengths", help="Epitope lengths separated by a comma")

    parser.add_argument("--timeout", type=int, default=3600, help="Netchop command timeout")
    parser.add_argument("-o", "--output", help="Output table filename (stdout if missing)")

    parser.add_argument(
        "epitopes", help="Neoepitope prediction results (*.all_epitopes.tsv or *.filtered.tsv)"
    )

    args = parser.parse_args()
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    if not args.base_path:
        base_path = os.getcwd()
    else:
        base_path = args.base_path

    records = read_epitopes_table(args.epitopes)
    if len(records) == 0:
        logging.info("No predicted neo-epitopes")
        Path.touch(args.output, mode=0o750)
        return 0
    logging.info(f"{len(records)} neo-epitope predictions have been read from file {args.epitopes}")

    colname = "MT Epitope Seq" if args.tool == "pvacseq" else "Epitope Seq"

    if args.hla_types:
        hla_types = COMMA.split(args.hla_types)
    else:
        hla_types = set(map(lambda record: record["HLA Allele"], records))
    if args.lengths:
        lengths = map(int, COMMA.split(args.lengths))
    else:
        lengths = set(map(lambda record: len(record[colname]), records))

    logging.info(
        f"Starting netMHCstabpan runs ({args.netMHCstabpan}) with {args.workers} processes"
    )
    epitopes = run_netMHCstabpan(
        records,
        hla_types=hla_types,
        lengths=lengths,
        base_path=base_path,
        container=args.container,
        netMHCstabpan_path=args.netMHCstabpan,
        netMHCpan_path=args.netMHCpan,
        epitope_seq_column_name=colname,
        n_workers=args.workers,
        timeout=args.timeout,
    )
    logging.info(f"{len(epitopes)} netMHCstabpan run completed")

    logging.info("Writing results")
    if args.output:
        f = open(args.output, "wt")
    else:
        f = sys.stdout
    write_output_table(f, epitopes)
    if args.output:
        f.flush()
        f.close()

    logging.info("Success - all done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
