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

from multiprocessing import Process, ProcessError, Queue, current_process
from pathlib import Path
from typing import Any, NamedTuple, Tuple, TextIO

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

logger = logging.getLogger(__name__)

NETMHCSTABPAN_CMD = r"""
allele=$1
length=$2
fasta=$3
out=$4

export TMPDIR=/base_temp
export NETMHCpan={netMHCpan_path}
export NETMHCstabpan={netMHCstabpan_path}

grep -q $allele $NETMHCstabpan/data/MHC_pseudo.dat
if [[ $? -eq 0 ]]
then
    $NETMHCstabpan/bin/netMHCstabpan \
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
COMMA = re.compile(r",")


class EpitopeSeq(NamedTuple):
    name: str
    seq: str
    indices: list[int] = []


def _epitope_length_filename(epitope_name: str, length: int) -> str:
    """Creates a filename for one epitope of a given length"""
    return f"{epitope_name}_{str(length)}"


def _create_sequence_files(
    seq_dir: str,
    epitopes: list[dict[str, str]],
    hla_types: set[str],
    lengths: set[int],
    seq_name: str = "MT Epitope Seq",
    hla_name: str = "HLA Allele",
) -> dict[str, Any]:
    """
    Creates fasta files containing all epitopes from one allele & one length.

    In each file, the sequences are unique (no duplicates).
    The returned file_list is a dict which has the file name (sometempdir/hlaname_epitopelength.fasta)
    as keys, and another dict as value. The second dict contains:
    - the hla type & length, and
    - the epitope sequence identifier, its AA sequence and the indices of such epitope in the
      epitope list.
      This triplet of values is indexed by the epitope sequence identifier, so that mapping the
      stability values back to the main epitope table is easy.
    """
    file_list = {}
    for iEpitope, epitope in enumerate(epitopes):
        hla_type = epitope[hla_name]
        epitope_length = len(epitope[seq_name])
        if hla_type not in hla_types or epitope_length not in lengths:
            continue

        fn = _epitope_length_filename(epitope[hla_name], epitope_length)
        if fn not in file_list:
            file_list[fn] = {
                "hla_type": hla_type,
                "length": epitope_length,
                "sequences": dict[str, EpitopeSeq](),
            }

        epitope_sequence = epitope[seq_name]
        if epitope_sequence not in file_list[fn]["sequences"]:
            n = len(file_list[fn]["sequences"].keys())
            file_list[fn]["sequences"][epitope_sequence] = EpitopeSeq(
                name=f"seq_{n}", seq=epitope_sequence, indices=[]
            )
        file_list[fn]["sequences"][epitope_sequence].indices.append(iEpitope)

    for fn, content in file_list.items():
        seq_by_name = {}
        for epitope_seq in content["sequences"].values():
            seq_by_name[epitope_seq.name] = epitope_seq

        with open(os.path.join(seq_dir, fn + ".fasta"), "wt") as f:
            for epitope_seq in seq_by_name.values():
                f.write(f">{epitope_seq.name}\n")
                f.write(epitope_seq.seq + "\n")

        file_list[fn]["sequences"] = seq_by_name

    return file_list


def _verify(
    file_list: dict[str, Any], epitopes: list[dict[str, str]], column_name: str = "MT Epitope Seq"
) -> bool:
    oks: list[bool | None] = [None] * len(epitopes)
    for fn, content in file_list.items():
        logger.debug(f"Starting check of {fn}")

        hla_type = content["hla_type"]
        length = content["length"]
        expected = f"{hla_type}_{str(length)}"
        if expected != fn:
            logger.error(f"Filename error: filename = {fn}, expected = {expected}")

        sequences: dict[str, EpitopeSeq] = content["sequences"]
        n = 0
        for name, epitope_seq in sequences.items():
            if name != epitope_seq.name:
                logger.error(f"Sequence ID error: key = {name}, sequence id = {epitope_seq.name}")
            for iEpitope in epitope_seq.indices:
                n += 1
                epitope = epitopes[iEpitope]
                if epitope["HLA Allele"] != hla_type:
                    logger.error(
                        f"HLA type mismatch for sequence {epitope_seq.name}, index {iEpitope}: {epitope['HLA Allele']} != {hla_type}"
                    )
                    oks[iEpitope] = False
                if epitope[column_name] != epitope_seq.seq:
                    logger.error(
                        f"Sequence mismatch for sequence {epitope_seq.name}, index {iEpitope}: {epitope[column_name]} != {epitope_seq.seq}"
                    )
                    oks[iEpitope] = False
                if len(epitope[column_name]) != length:
                    logger.error(
                        f"Sequence length mismatch for sequence {epitope_seq.name}, index {iEpitope}"
                    )
                    oks[iEpitope] = False
                if oks[iEpitope] is None:
                    oks[iEpitope] = True

        logger.debug(f"Check of {fn} complete, {n} epitopes assessed")

    n = 0
    for ok in oks:
        if ok is None:
            n += 1
    if n > 0:
        logger.error(f"{n} epitopes not assessed")

    return n == 0 and all(oks)


def _parse_netMHCstabpan_output(out: list[str]) -> Tuple[str, float]:
    """Parses long netMHCstabpan output"""
    distance = 0.0
    replacement = ""
    for line in out:
        m = ALLELE_PROXIMITY_PATTERN.match(line.strip())
        if m:
            distance = float(m.group(2))
            replacement = m.group(3)
            break
    return (replacement, distance)


def _parse_output_file(fn: str) -> list[dict[str, Any]]:
    """
    Parses the netMHCstabpan output file

    The file starts with one line containing the HLA allele,
    followed by a table with columns:

    Note that the file has generally a *.xls extension, but it is a text file.
    """
    results = []
    with open(fn, "rt") as f:
        hla_type = f.readline().strip()  # noqa: F841
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = {k: v.strip() for k, v in row.items()}
            results.append(row)
    return results


def _run_netMHCstabpan_command(cmd: list[str], worker_tmp: str, timeout: int = 43200):
    """
    Runs netMHCstabpan on one fasta file, to find peptide stability score

    netMHCstabpan creates a file with stability scores for each sequence, but
    stdout must also be parsed (in _parse_netMHCstabpan_output) to extract
    the actual allele which has been used to compute stability, and its closeness to
    the requested allele.

    The output file is separately parsed by _parse_output_file
    """
    logger.debug(f"Process {current_process().name}, cmd = {' '.join(cmd)}, timeout = {timeout}")

    # Compute stability scores
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    try:
        out, err = p.communicate(timeout=timeout)
    except TimeoutError:
        p.kill()
        raise TimeoutError(f"The command {' '.join(cmd)} has timed out")
    if p.returncode != 0:
        raise ChildProcessError(f"Command {' '.join(cmd)} failed with return code {p.returncode}")

    # Return parsed stability output (allele used and distance to requested allele)
    return _parse_netMHCstabpan_output(out.split("\n"))


def _worker(
    task_queue: Queue,
    output_queue: Queue,
    error_queue: Queue,
    cmd: list[str],
    timeout: int = 43200,
    log_level: int = logging.WARNING,
):
    """Multi-processing intermediate for netchop"""
    logger = logging.getLogger(current_process().name)
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=log_level,
    )

    logger.debug(f"Process {current_process().name} starts")
    while True:
        args = task_queue.get()
        if args is None:
            logger.debug(f"Process {current_process().name} received stop signal")
            break
        hla_type, length = (args[0], args[1])
        sequences = os.path.join(
            "/main_path", "seqs", _epitope_length_filename(hla_type, length) + ".fasta"
        )
        out = os.path.join("/main_path", "out", _epitope_length_filename(hla_type, length) + ".xls")
        logger.debug(
            f"Processing queued job with arguments {hla_type}, {length}, {sequences}, {out}"
        )
        try:
            output_queue.put(
                (
                    hla_type,
                    length,
                    _run_netMHCstabpan_command(
                        cmd + [hla_type.replace("*", ""), str(length), sequences, out],
                        "/base_temp",
                        timeout,
                    ),
                )
            )
        except Exception as e:
            error_queue.put((hla_type, length, e))


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
    timeout: int = 43200,
):
    """
    Runs netMHCstabpan within docker container

    Because computing stability is quite time-consuming, the function delegates
    the bulk of the work to a number of workers, each running netMHCstabpan on
    one allele/length combination. This adds some complexity to the function.

    The function prepares the jobs, creating the input data and separate temp
    for each worker. Then it creates and starts the worker processes, creates a
    task queue, and two queues to receive the results and possible errors.
    Finally, is waits for the workers to finish, and collects the results into
    the main epitope table.

    There are several directories that must be created and bound to the container.
    All those directories are temporary to a different degree. Some contain the
    input data prepared by the main process and must remain present until
    the job is finished. Some are private temporary directories used by the workers.
    Finally, a separate directory is created to hold the script that is executed
    in the container (to allow to keep in read-only during binding to the container).
    """
    # Real paths and bound paths
    base_path = os.path.realpath(base_path)
    base_path = os.path.join(base_path, "tmp", "netMHCstabpan")
    if os.path.exists(base_path):
        shutil.rmtree(base_path)

    main_path = os.path.join(base_path, "main")
    seq_dir = "seqs"
    out_dir = "out"
    for d in (seq_dir, out_dir):
        os.makedirs(os.path.join(main_path, d), mode=0o750, exist_ok=False)

    script = os.path.join(base_path, "scripts", "run_netMHCstabpan.sh")
    os.makedirs(os.path.dirname(script), mode=0o750, exist_ok=True)
    with open(script, "wt") as f:
        f.write(
            NETMHCSTABPAN_CMD.format(
                netMHCstabpan_path=os.path.dirname(os.path.dirname(netMHCstabpan_path)),
                netMHCpan_path=os.path.dirname(os.path.dirname(netMHCpan_path)),
            )
        )

    for i_worker in range(n_workers):
        os.makedirs(os.path.join(base_path, "workers", str(i_worker)), mode=0o700, exist_ok=False)

    # Prepare fasta files with peptide sequences for netMHCstabpan input
    file_list = _create_sequence_files(
        os.path.join(main_path, seq_dir),
        epitopes,
        hla_types,
        lengths,
        seq_name=epitope_seq_column_name,
    )
    assert _verify(file_list, epitopes), "Internal error: verification of sequence files failed"

    # Prepare netMHCstabpan command (binding, container path, ...)
    cmd_part1 = [
        "apptainer",
        "run",
        "--no-home",
        "--bind",
        f"{main_path}:/main_path:rw",
        "--bind",
        f"{os.path.dirname(script)}:/scripts:ro",
    ]
    cmd_part2 = [
        os.path.realpath(container),
        "bash",
        os.path.join("/scripts", os.path.basename(script)),
    ]
    logger.info(
        f"netMHCstabpan command: {' '.join(cmd_part1 + ['--bind', '<base_path>/workers/<i_worker>:/base_temp:rw'] + cmd_part2 + ['<hla_type>', '<length>', '<fasta>', '<xls>'])}"
    )

    # Prepare multiprocessing
    task_queue = Queue()
    output_queue = Queue()
    error_queue = Queue()
    processes: list[Process] = []

    # Fill all epitopes with stability default values
    for i in range(len(epitopes)):
        for k in ("Predicted Stability", "Half Life", "Stability Rank", "NetMHCstab allele"):
            epitopes[i][k] = "NA"

    # Start the workers
    for i_worker in range(n_workers):
        worker_temp = os.path.join(base_path, "workers", str(i_worker))
        cmd = (
            cmd_part1
            + [
                "--bind",
                worker_temp + ":/base_temp:rw",
            ]
            + cmd_part2
        )
        p = Process(
            target=_worker,
            args=(
                task_queue,
                output_queue,
                error_queue,
                cmd,
                timeout,
                logger.root.level,
            ),
        )
        processes.append(p)
        logger.info(f"Starting worker {i_worker}")
        p.start()

    for file_name, descr in file_list.items():
        task_queue.put((descr["hla_type"], descr["length"]))

    time.sleep(2.0)

    for i_worker in range(n_workers):
        task_queue.put(None)

    # Wait for completion
    for p in processes:
        p.join()
        logger.info(f"Worker {p.name} completed")

    # Fill the epitope tables with the output queue & the contents of output files
    while not output_queue.empty():
        hla_type, length, replacement = output_queue.get()

        if replacement[0]:
            replacement = f"{replacement[0]} (distance: {float(replacement[1])})"
        else:
            replacement = f"{hla_type} (distance: {float(replacement[1])})"

        results = _parse_output_file(
            os.path.join(
                main_path,
                "out",
                _epitope_length_filename(hla_type, length) + ".xls",
            )
        )

        fn = _epitope_length_filename(hla_type, length)
        content = file_list[fn]
        logger.debug(f"Inserting {len(results)} results from file {fn}")

        for row in results:
            iSeq = row["ID"]
            assert iSeq in content["sequences"], f"Can't find sequence {iSeq} in {fn}"
            for iEpitope in content["sequences"][iSeq].indices:
                for src, dest in (
                    ("Pred", "Predicted Stability"),
                    ("Thalf(h)", "Half Life"),
                    ("Rank", "Stability Rank"),
                ):
                    epitopes[iEpitope][dest] = row[src]
                epitopes[iEpitope]["NetMHCstab allele"] = replacement

    # Check for errors
    error = False
    while not error_queue.empty():
        error = True
        hla_type, length, e = error_queue.get()
        logger.error(
            f"An error occurred during netMHCstabpan for {hla_type} epitopes of length {length} - message {e}"
        )
    if error:
        raise ProcessError("Error running one of netMHCstabpan sub-processes")

    shutil.rmtree(base_path)

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
        assert isinstance(reader.fieldnames, list), (
            "Internal error: protects against pandas-created tables(?)"
        )
        assert "Index" in reader.fieldnames, f"Index column missing from {fn}"
        for row in reader:
            records.append(row)
    return records


def write_output_table(
    f: io.TextIOWrapper | TextIO,
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

    parser.add_argument("--timeout", type=int, default=43200, help="Netchop command timeout")
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
        logger.info("No predicted neo-epitopes")
        Path.touch(args.output, mode=0o750)
        return 0
    logger.info(f"{len(records)} neo-epitope predictions have been read from file {args.epitopes}")

    hla_types = set(map(lambda record: record["HLA Allele"], records))
    if args.hla_types:
        hla_types = hla_types.intersection(set(COMMA.split(args.hla_types)))

    colname = "MT Epitope Seq" if args.tool == "pvacseq" else "Epitope Seq"
    lengths = set(map(lambda record: len(record[colname]), records))
    if args.lengths:
        lengths = lengths.intersection(set(map(int, COMMA.split(args.lengths))))

    logger.info(f"Starting netMHCstabpan runs ({args.netMHCstabpan}) with {args.workers} processes")
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
    logger.info(f"{len(epitopes)} netMHCstabpan run completed")

    logger.info("Writing results")
    if args.output:
        f = open(args.output, "wt")
    else:
        f = sys.stdout
    write_output_table(f, epitopes)
    if args.output:
        f.flush()
        f.close()

    logger.info("Success - all done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
