import argparse
import csv
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

from operator import itemgetter
from multiprocessing import Manager, Process, ProcessError, Queue, current_process
from pathlib import Path
from typing import Any

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class Epitope:
    ALLOWED_AA = (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )

    def __init__(
        self,
        iRecord: int,
        epitope: str,
        mt_peptide: str,
        wt_peptide: str | None = None,
        flanking_sequence_length: int = 9,
    ):
        self.iRecord = iRecord
        self.epitope = epitope

        if wt_peptide:
            self.peptide = self._get_mutated_peptide_with_flanking_sequence(
                wt_peptide, mt_peptide, flanking_sequence_length
            )
            self.start_diff = flanking_sequence_length
        else:
            self.peptide, self.start_diff = self._extract_flanked_epitope(
                mt_peptide, self.epitope, flanking_sequence_length
            )

        self.seq_hash = str(hash(self.epitope + "\t" + self.peptide + "\t" + str(self.start_diff)))
        self.sites = []

    def _extract_flanked_epitope(
        self, full_peptide: str, epitope: str, flanking_sequence_length: int
    ) -> tuple[str, int]:
        assert epitope in full_peptide, (
            f"Epitope {epitope} not in {full_peptide} for {self.iRecord}"
        )
        ep_start = full_peptide.index(epitope)
        start = max(0, ep_start - flanking_sequence_length)
        start_diff = ep_start - start
        end = ep_start + len(self.epitope) + flanking_sequence_length
        return full_peptide[start:end], start_diff

    def _get_mutated_peptide_with_flanking_sequence(
        self, wt_peptide, mt_peptide, flanking_length
    ) -> str:
        wt_l = len(wt_peptide)
        mt_l = len(mt_peptide)

        n = min(wt_l, mt_l) - flanking_length + 1
        for start in range(n):
            if (
                wt_peptide[start : (start + flanking_length)]
                != mt_peptide[start : (start + flanking_length)]
            ):
                break

        n = max(min(wt_l, mt_l) - start - flanking_length, 1)
        for i in range(n):
            wt_i = wt_l - i - flanking_length + 1
            mt_i = mt_l - i - flanking_length + 1
            if (
                wt_peptide[wt_i : (wt_i + flanking_length)]
                != mt_peptide[mt_i : (mt_i + flanking_length)]
            ):
                break
        stop = min(mt_i + flanking_length, mt_l)

        mutant_subsequence = mt_peptide[start:stop]

        if mutant_subsequence[0] not in self.ALLOWED_AA:
            mutant_subsequence = mutant_subsequence[1:]
        if mutant_subsequence[-1] not in self.ALLOWED_AA:
            mutant_subsequence = mutant_subsequence[:-1]
        if not all([c in self.ALLOWED_AA for c in mutant_subsequence]):
            logging.warning(
                f"Mutant sequence contains unsupported amino acid. Skipping entry {self.iRecord}"
            )
            return ""
        return mutant_subsequence

    def set_sites(self, sites: list[tuple[int, float]]):
        for site in sites:
            pos = site[0]
            if self.start_diff <= pos and pos <= self.start_diff + len(self.epitope):
                self.sites.append((pos - self.start_diff, site[1]))
        self.sites.sort(key=itemgetter(1), reverse=True)

    def format_sites(self) -> str:
        sites = {
            "Best Cleavage Position": "NA",
            "Best Cleavage Score": "NA",
            "Cleavage Sites": "NA",
        }
        if self.sites:
            sites["Best Cleavage Position"] = self.sites[0][0]
            sites["Best Cleavage Score"] = self.sites[0][1]
            sites["Cleavage Sites"] = []
            for site in self.sites:
                sites["Cleavage Sites"].append(site)
            sites["Cleavage Sites"] = ",".join([f"{k}:{v}" for k, v in sites["Cleavage Sites"]])
        return "\t".join(map(str, sites.values()))


def _parse_netchop_output(
    out: io.StringIO,
    pattern: re.Pattern = re.compile(
        r"^ *(?P<pos>[0-9]+) +(?P<AA>[ACDEFGHIKLMNPQRSTVWY]) +(?P<site>[^ ]+) +(?P<score>[01]\.[0-9]+) +(?P<Ident>.+)"
    ),
):
    """Parses long netchop output"""
    sites = []
    for line in out:
        m = pattern.match(line.strip())
        if m and m.group("site") == "S":
            sites.append((int(m.group("pos")), float(m.group("score"))))
    return sites


def _run_netchop_command(
    epitope: Epitope, cmd: list[str], worker_tmp: str, timeout: int = 3600, line_length: int = 80
):
    """
    Runs netchop to find cleavage sites in on protein sequence

    The sequence is first saved in a temp directory, and netchop is run.
    The output is parsed and cleavage sites are stored in the epitope object.
    The site position and its score are retained.
    """
    logging.debug(
        f"Process {current_process().name} epitope {epitope.iRecord}, cmd = {' '.join(cmd)}, timeout = {timeout}"
    )

    # Write the mutated sequence
    fn = os.path.join(worker_tmp, f"sequence_{epitope.iRecord}.fasta")
    with open(fn, "wt") as f:
        f.write(f">seq_{epitope.iRecord}\n")
        i = 0
        while i < len(epitope.peptide):
            f.write(epitope.peptide[i : min(i + line_length, len(epitope.peptide))] + "\n")
            i += line_length

    # Find all cleavage sites
    # fn = "../../../resources/netchop-3.1/test/test.fsa"
    p = subprocess.Popen(cmd + [fn], stdout=subprocess.PIPE)
    try:
        out, err = p.communicate(timeout=timeout)
    except TimeoutError:
        p.kill()
        raise (f"The command {' '.join(cmd)} has timed out")
    if p.returncode != 0:
        raise ChildProcessError(
            f"Command {' '.join(cmd + [fn])} failed with return code {p.returncode}"
        )

    # Parse cleavage sites & put them into epitope object
    return _parse_netchop_output(out.decode("utf-8").split("\n"))


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
        epitope: Epitope = task_queue.get()
        try:
            return_dict[str(epitope.iRecord)] = _run_netchop_command(
                epitope, cmd, worker_tmp, timeout
            )
        except Exception as e:
            error_queue.put((e, epitope))


def _netchop_workaround(
    netchop: str, workaround_dir: str = os.path.join(os.getcwd(), "tmp"), force: bool = False
):
    """
    Workaround problems running netchop

    netchop doesn't run when the NETCHOP environment variable is set to an absolute path.
    (Perhaps because the path is too long?).
    The workaround creates a local temp directory, a symlink to the base netchop installation
    (netchop is in netchop installation/Linux_x86_64/bin), and sets the NETCHOP environment variable
    to netchop-3.1/Linux_x86_64 (valid from tmpdir).
    It seems that setting the temp directory somewhere within $TMPDIR doesn't work. Don't ask me why.
    Finally, all efforts trying to make the -tdir argument to work failed miserably.
    """
    netchop = os.path.realpath(netchop)
    netchop_dir = os.path.dirname(os.path.dirname(netchop))
    nmhome_dir = os.path.dirname(netchop_dir)

    # Workaround problems with running netchop in any directory
    nmhome_rel = os.path.basename(nmhome_dir)
    netchop_rel = os.path.join(nmhome_rel, os.path.basename(netchop_dir))

    # Create or use existing workaround_dir
    try:
        os.makedirs(workaround_dir, mode=0o755, exist_ok=False)
    except FileExistsError as e:
        if force:
            try:
                os.remove(workaround_dir)
                logging.warning(
                    f"File {workaround_dir} has been removed to make space for temp directory"
                )
            except OSError:
                shutil.rmtree(workaround_dir)
                logging.warning(
                    f"Directory {workaround_dir} has been removed to make space for temp directory"
                )
            os.makedirs(workaround_dir, mode=0o755, exist_ok=False)
        else:
            raise e

    # Create symlink to NMHOME in workaround_dir
    try:
        os.symlink(nmhome_dir, os.path.join(workaround_dir, nmhome_rel))
        logging.info(f"Created symlink {nmhome_rel} -> {nmhome_dir}")
    except FileExistsError as e:
        if force:
            try:
                os.remove(os.path.join(workaround_dir, nmhome_rel))
                logging.warning(
                    f"File {os.path.join(workaround_dir, nmhome_rel)} has been removed to make space for symlink"
                )
            except OSError:
                shutil.rmtree(nmhome_rel)
                logging.warning(
                    f"Directory {os.path.join(workaround_dir, nmhome_rel)} has been removed to make space for symlink"
                )
            os.symlink(nmhome_dir, os.path.join(workaround_dir, nmhome_rel))
            logging.info(f"Created symlink {nmhome_rel} -> {nmhome_dir}")
        else:
            raise e

    os.environ["NETCHOP"] = netchop_rel


def run_netchop(
    epitopes: dict[str, Epitope],
    netchop: str,
    args: dict[str, Any] = {},
    workaround_dir: str = os.path.join(os.getcwd(), "tmp"),
    clean: bool = True,
    force: bool = False,
    n_workers: int = 1,
    timeout: int = 3600,
):
    """
    Runs netchop on all variants creating neo-epitopes

    An ugly workaround is apparently necessary...
    """
    # Set variables as recommended
    netchop = os.path.realpath(netchop)
    os.environ["NETCHOP"] = os.path.dirname(os.path.dirname(netchop))
    os.environ["NMHOME"] = os.path.dirname(os.environ["NETCHOP"])

    # Workaround problems getting netchop to work with long paths(?)
    # The workaround creates temp directory, symlinks & alters environment variables
    current_dir = os.getcwd()
    _netchop_workaround(netchop, workaround_dir, force)
    os.chdir(workaround_dir)
    logging.info(f"Working from newly created directory {workaround_dir}")
    logging.info(f"NMHOME environment variable set to {os.environ['NMHOME']}")
    logging.info(f"NETCHOP environment variable set to {os.environ['NETCHOP']}")

    # Prepare netchop command
    cmd = [
        netchop,
        "-v",
        "1" if args.get("method", "cterm") == "20s" else "0",
        "-t",
        str(args.get("threshold", 0.5)),
    ]
    logging.info(f"netchop command: {' '.join(cmd + ['<fn>'])}")

    # Prepare multiprocessing
    tmpdir = tempfile.mkdtemp()
    worker_tmp_template = os.path.join(tmpdir, "worker_{i_worker}")
    for i_worker in range(n_workers):
        worker_tmp = worker_tmp_template.format(i_worker=i_worker)
        os.makedirs(worker_tmp, mode=0o700, exist_ok=False)

    # Avoid duplication: epitopes with same epitope, peptide & starting pos are just computed once
    task_queue = Queue()
    indices = {}
    for epitope in epitopes:
        h = epitope.seq_hash
        if h not in indices:
            indices[h] = []
            task_queue.put(epitope)
        indices[h].append(epitope.iRecord)
    error_queue = Queue()
    return_dict = Manager().dict()
    processes: list[Process] = []

    time.sleep(2.0)

    # Start the workers
    for i_worker in range(n_workers):
        worker_tmp = worker_tmp_template.format(i_worker=i_worker)
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
        e, epitope = error_queue.get()
        logging.error(
            f"An error occurred during netchop for epitope {epitope.iRecord} - message {e}"
        )
    if error:
        raise ProcessError("Error running one of netchop sub-processes")

    # Rapatriate netchop results into epitope
    for iRecord, sites in return_dict.items():
        iRecord = int(iRecord)
        h = epitopes[iRecord].seq_hash
        for i in indices[h]:
            epitopes[i].set_sites(sites)

    # Clean-up workaround
    os.chdir(current_dir)
    if clean:
        shutil.rmtree(workaround_dir)
        shutil.rmtree(tmpdir)


def read_fasta(fn: str | Path) -> dict[str, str]:
    sequences = {}
    seqname = None
    seq = ""

    with open(fn, "rt") as f:
        iLine = 0
        for line in f:
            iLine += 1
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if seqname:
                    if seqname in sequences:
                        logging.warning(
                            f"Duplicate sequence {seqname} on line {iLine} of {fn}, ignored"
                        )
                    else:
                        assert seq, f"Missing sequence for {seqname} of {fn}"
                        sequences[seqname] = seq
                seqname = line.strip("> ")
                seq = ""
            else:
                seq += line

    if not seqname:
        assert seq == "", f"{fn} is not FASTA format"
        logging.warning(f"{fn} contains no sequence")
    elif seqname in sequences:
        logging.warning(f"Duplicate sequence {seqname} on last line of {fn}, ignored")
    else:
        assert seq, f"Missing sequence for {seqname} of {fn}"
        sequences[seqname] = seq

    return sequences


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


def create_epitope_objects_pvacseq(
    sequences: dict[str, str], records: list[dict[str, str]]
) -> list[Epitope]:
    epitopes = []
    for iRecord, record in enumerate(records):
        epitopes.append(
            Epitope(iRecord, record["MT Epitope Seq"], sequences["MT." + record["Index"]])
        )
    return epitopes


def create_epitope_objects_pvacfuse(
    sequences: dict[str, str], records: list[dict[str, str]]
) -> list[Epitope]:
    epitopes = []
    for iRecord, record in enumerate(records):
        epitopes.append(Epitope(iRecord, record["Epitope Seq"], sequences[record["Index"]]))
    return epitopes


def create_epitope_objects_pvacsplice(
    sequences: dict[str, str], records: list[dict[str, str]]
) -> list[Epitope]:
    epitopes = []
    for iRecord, record in enumerate(records):
        epitopes.append(
            Epitope(
                iRecord,
                record["Epitope Seq"],
                sequences["ALT." + record["Index"]],
                sequences["WT." + record["Index"]],
            )
        )
    return epitopes


def write_output_table(
    f: io.TextIOBase,
    epitopes: list[Epitope],
    records: list[dict[str, Any]],
):
    """Add 3 columns to the epitope table (Best cleavage pos & score, and digest of all cleavage sites)"""
    assert len(records) == len(epitopes), (
        f"Internal error: input ({len(records)}) and output ({len(epitopes)}) out of sync"
    )
    titles = records[0].keys()
    f.write(
        "\t".join(
            list(titles) + ["Best Cleavage Position", "Best Cleavage Score", "Cleavage Sites"]
        )
        + "\n"
    )
    for iRecord in range(len(records)):
        record = records[iRecord]
        epitope = epitopes[iRecord]
        line = "\t".join([record[title] for title in titles]) + "\t" + epitope.format_sites()
        f.write(line + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="Wrapper around netchop",
        description="Runs netchop on results from pVACtools modules",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verbosity level")
    parser.add_argument(
        "-f", "--force", action="store_true", help="Force creation of temp directory and symlink"
    )
    parser.add_argument(
        "-T", "--tmpdir", default="tmp", help="Temp directory outside of $TMPDIR (don't ask why...)"
    )
    parser.add_argument(
        "-w", "--workers", type=int, default=1, help="Number of threads for netchop"
    )

    parser.add_argument(
        "-t",
        "--tool",
        choices=("pvacseq", "pvacfuse", "pvacsplice"),
        default="pvacseq",
        help="pVACtool module (pvacseq, pvacfuse, pvacsplice)",
    )
    parser.add_argument("-n", "--netchop", nargs=1, help="Path to the netchop binary")
    parser.add_argument("--timeout", type=int, default=3600, help="Netchop command timeout")
    parser.add_argument(
        "-m", "--method", choices=("cterm", "20s"), default="cterm", help="Netchop method"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Score threshold to filter cleaving sites",
    )
    parser.add_argument("-o", "--output", help="Output table filename (stdout if missing)")

    parser.add_argument(
        "epitopes", help="Neoepitope prediction results (*.all_epitopes.tsv or *.filtered.tsv)"
    )
    parser.add_argument(
        "sequences",
        help="Neighboring sequences near somatic variants in fasta format",
    )

    args = parser.parse_args()
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    records = read_epitopes_table(args.epitopes)
    if len(records) == 0:
        logging.info("No predicted neo-epitopes")
        Path.touch(args.output, mode=0o750)
        return 0
    logging.info(f"{len(records)} neo-epitope predictions have been read from file {args.epitopes}")
    sequences = read_fasta(args.sequences)
    logging.info(f"{len(sequences)} sequences have been read from file {args.sequences}")

    create_epitope_objects = globals().get(f"create_epitope_objects_{args.tool}", None)
    assert create_epitope_objects, f"Tool {args.tool} not implemented"
    epitopes = create_epitope_objects(sequences, records)

    logging.info(f"Starting netchop runs ({args.netchop}) with {args.workers} processes")
    run_netchop(
        epitopes,
        args.netchop[0],
        args={"method": args.method, "threshold": args.threshold},
        workaround_dir=args.tmpdir,
        clean=args.verbose,
        force=args.force,
        n_workers=args.workers,
        timeout=args.timeout,
    )
    logging.info(f"{len(epitopes)} netchop run completed")

    logging.info("Writing results")
    if args.output:
        f = open(args.output, "wt")
    else:
        f = sys.stdout
    write_output_table(f, epitopes, records)
    if args.output:
        f.flush()
        f.close()

    logging.info("Success - all done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
