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
from typing import Any, Self

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class Variant:
    """
    Slightly enhanced data class to store somatic variant

    The variant's HGVSp serves as unique identifier of the variant (at the protein level).
    If two variants share the same HGVSp, the second is ignored.
    Only protein sequence altering variants are stored.
    """

    AMINO_ACID: re.Pattern = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]+")
    NUCLEOTIDE: re.Pattern = re.compile(r"^[ACGT-]+")
    PROTEIN_POSITION: re.Pattern = re.compile(r"^([0-9]+)(-([0-9]+))?$")
    MUTATION: re.Pattern = re.compile(
        r"^([ACDEFGHIKLMNPQRSTVWY]*[X\*]|[ACDEFGHIKLMNPQRSTVWY]+[X\*]?)$"
    )

    def __init__(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        identifier: str,
        feature: str,
        sequence: str,
        wt_seq: str,
        mt_seq: str,
        start: int,
        end: int | None = None,
    ):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

        self.feature = feature
        self.identifier = identifier

        self.sequence = sequence
        self.wt_seq = wt_seq
        self.mt_seq = mt_seq

        self.start = start
        self.end = end if end is not None else self.start

        self.sites = [(-1, -1.0)]

        self._check_input()

    def _check_input(self):
        assert self.pos > 0, f"Negative variant position {self.pos}"
        assert self.NUCLEOTIDE.match(self.ref), f"Illegal reference allele {self.ref}"
        assert self.NUCLEOTIDE.match(self.alt), f"Illegal alt allele {self.alt}"
        assert self.AMINO_ACID.match(self.sequence), f"Illegal protein sequence {self.sequence}"
        assert self.start > 0 and self.end >= self.start and self.end <= len(self.sequence) + 1, (
            f"Mutation position {self.start}-{self.end} illegal or outside protein bounds (length {len(self.sequence)})"
        )
        assert self.MUTATION.match(self.wt_seq), f"Illegal wild-type sequence {self.wt_seq}"
        assert self.MUTATION.match(self.mt_seq), f"Illegal mutation sequence {self.mt_seq}"

        if self.end == len(self.sequence) + 1:
            assert self.wt_seq.endswith("*"), (
                f"Mutation position {self.start}-{self.end} illegal or outside protein bounds (length {len(self.sequence)})"
            )

    @staticmethod
    def _parse_table(out: str) -> list[Self]:
        """
        Parse the bcftools +split-vep output

        The table must have 10 columns:
        1. CHROM,
        2. POS,
        3. REF,
        4. ALT,
        5. HGVSp,
        6. Feature (transcript ID)
        7. Protein_position (can be from-to),
        8. Amino_acids,
        9. FrameshiftSequence, and
        10. WildtypeProtein

        Only 1 to 4 are always present, 7 is used as flag for the presence of a protein,
        and when empty, the variant is ignored. Silent variants are checked with 8, and also ignored.
        The sequence stored in the variant is 9 when not empty, and 10 otherwise.
        """
        variants = {}
        for line in out.split("\n"):
            if line == "":
                continue
            tokens = line.strip().split("\t")
            assert len(tokens) >= 10, f"Not enough fields in {line.strip()}"

            contig = tokens[0]
            pos = int(tokens[1])
            ref = tokens[2]
            alt = tokens[3]
            identifier = tokens[4]
            feature = tokens[5]
            assert identifier not in variants, f"Duplicated variant {identifier}"

            # Check that variant is in coding sequence
            m = Variant.PROTEIN_POSITION.match(tokens[6])
            if not m:
                continue
            start = int(m.group(1))
            end = int(m.group(3)) if m.group(2) is not None else None

            # Check that variant is not silent
            mutation = tokens[7].split("/")
            if len(mutation) == 1:
                continue
            assert len(mutation) == 2, f"Illegal mutation {tokens[7]}"
            wt_seq = mutation[0]
            mt_seq = mutation[1]

            sequence = tokens[9] if tokens[8] == "." or tokens[8] == "" else tokens[8]

            variant = Variant(
                contig,
                pos,
                ref,
                alt,
                identifier,
                feature,
                sequence,
                wt_seq,
                mt_seq,
                start,
                end,
            )
            variants[identifier] = variant

        return variants

    @staticmethod
    def parse_annotated_vcf(fn: str | Path, timeout: int = 300) -> list[tuple]:
        """Runs bcftools +split-vep & parses the output."""
        bcftools = shutil.which("bcftools")
        fmt = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "HGVSp",
            "Feature",
            "Protein_position",
            "Amino_acids",
            "FrameshiftSequence",
            "WildtypeProtein",
        ]
        cmd = [bcftools, "+split-vep", fn, "-f", "%{}".format("\t%".join(fmt))]
        logging.debug(f"VCF parsing command: {' '.join(cmd)}")
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        try:
            out, err = p.communicate(timeout=timeout)
        except TimeoutError:
            p.kill()
            raise (f"The command {' '.join(cmd)} has timed out")
        if p.returncode != 0:
            raise ChildProcessError(
                f"Command {' '.join(cmd)} failed with return code {p.returncode}"
            )
        variants = Variant._parse_table(out.decode("utf-8"))
        return variants


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
    return sorted(sites, key=lambda x: x[1])


def _run_netchop_command(
    variant: Variant, cmd: list[str], worker_tmp: str, timeout: int = 3600, line_length: int = 80
):
    """
    Runs netchop to find cleavage sites in on protein sequence

    The sequence is first saved in a temp directory, and netchop is run.
    The output is parsed and cleavage sites are stored in the variant object.
    The site position and its score are retained.
    """
    logging.debug(
        f"Process {current_process().name} variant {variant.identifier}, cmd = {' '.join(cmd)}, timeout = {timeout}"
    )

    # Write the mutated sequence
    fn = os.path.join(worker_tmp, "sequence.fasta")
    with open(fn, "wt") as f:
        f.write(f">{variant.feature}\n")
        i = 0
        while i < len(variant.sequence):
            f.write(variant.sequence[i : max(i + line_length, len(variant.sequence))] + "\n")
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

    # Parse cleavage sites & put them into variant object
    return _parse_netchop_output(out.decode("utf-8").split("\n"))


def _worker(
    task_queue: Queue,
    error_queue: Queue,
    return_dict: dict,
    cmd: list[str],
    worker_tmp: str,
    timeout: int = 3600,
):
    """Multi-processing intermediate for netchop"""
    logging.debug(
        f"Starting worker with cmd = {' '.join(cmd)}, tmpdir = {worker_tmp}, timeout = {timeout}"
    )
    while not task_queue.empty():
        variant: Variant = task_queue.get()
        try:
            return_dict[variant.identifier] = _run_netchop_command(
                variant, cmd, worker_tmp, timeout
            )
        except Exception as e:
            error_queue.put((e, variant))


def _netchop_workaround(
    netchop: str, tmpdir: str = os.path.join(os.getcwd(), "tmp"), force: bool = False
):
    """
    Workaround problems running netchop

    netchop doesn't run when the NETCHOP environment variable is set to an absolute path.
    (Perhaps because the path is too long?).
    The workaround create a local temp directory, creates a symlink to the base netchop installation
    (netchop is in netchop installation/Linux_x86_64/bin), and sets the NETCHOP environment variable
    to tmpdir/netchop-3.1/Linux_x86_64.
    It seems that setting the temp directory somewhere within $TMPDIR doesn't work. Don't ask me why.
    """
    netchop = os.path.realpath(netchop)
    netchop_dir = os.path.dirname(os.path.dirname(netchop))
    nmhome_dir = os.path.dirname(netchop_dir)

    # Workaround problems with running netchop in any directory
    nmhome_rel = os.path.basename(nmhome_dir)
    netchop_rel = os.path.join(nmhome_rel, os.path.basename(netchop_dir))

    # Create or use existing temp directory
    try:
        os.makedirs(tmpdir, mode=0o755, exist_ok=False)
    except FileExistsError as e:
        if force:
            try:
                os.remove(tmpdir)
                logging.warning(f"File {tmpdir} has been removed to make space for temp directory")
            except OSError:
                shutil.rmtree(tmpdir)
                logging.warning(
                    f"Directory {tmpdir} has been removed to make space for temp directory"
                )
            os.makedirs(tmpdir, mode=0o755, exist_ok=False)
        else:
            raise e

    # Create symlink to NMHOME in temp dir
    try:
        os.symlink(nmhome_dir, os.path.join(tmpdir, nmhome_rel))
        logging.info(f"Created symlink {nmhome_rel} -> {nmhome_dir}")
    except FileExistsError as e:
        if force:
            try:
                os.remove(os.path.join(tmpdir, nmhome_rel))
                logging.warning(
                    f"File {os.path.join(tmpdir, nmhome_rel)} has been removed to make space for symlink"
                )
            except OSError:
                shutil.rmtree(nmhome_rel)
                logging.warning(
                    f"Directory {os.path.join(tmpdir, nmhome_rel)} has been removed to make space for symlink"
                )
            os.symlink(nmhome_dir, os.path.join(tmpdir, nmhome_rel))
            logging.info(f"Created symlink {nmhome_rel} -> {nmhome_dir}")
        else:
            raise e

    os.environ["NETCHOP"] = netchop_rel


def run_netchop(
    variants: dict[str, Variant],
    netchop: str,
    args: dict[str, Any] = {},
    tmpdir: str = os.path.join(os.getcwd(), "tmp"),
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
    _netchop_workaround(netchop, tmpdir, force)
    os.chdir(tmpdir)
    logging.info(f"Working from newly created directory {tmpdir}")
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
    logging.debug(f"netchop command: {' '.join(cmd + ['<fn>'])}")

    # Prepare multiprocessing
    for i_worker in range(n_workers):
        worker_tmp = f"worker_{i_worker}"
        os.makedirs(worker_tmp, mode=0o700, exist_ok=False)

    task_queue = Queue()
    for variant in variants.values():
        task_queue.put(variant)
    error_queue = Queue()
    return_dict = Manager().dict()
    processes: list[Process] = []

    time.sleep(2.0)

    # Start the workers
    for i_worker in range(n_workers):
        worker_tmp = f"worker_{i_worker}"
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
        e, variant = error_queue.get()
        logging.error(
            f"An error occurred during netchop for variant {variant.identifier} - message {e}"
        )
    if error:
        raise ProcessError("Error running one of netchop sub-processes")

    # Rapatriate netchop results into variants
    for identifier, sites in return_dict.items():
        variants[identifier].sites = sites

    # Clean-up workaround
    os.chdir(current_dir)
    if clean:
        shutil.rmtree(tmpdir)


def read_epitopes_table(fn: str | Path, columns: list[str] = []) -> dict[str, dict[str, Any]]:
    """
    Reads neo-epitope table

    Produced by pVACseq (should work for pVACSplice and pVACfuse, but untested)
    The files are generally <sample>.<MHC class>.(all_epitopes|filtered).tsv.
    """
    epitopes = {}
    with open(fn, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if columns:
            assert all(column in reader.fieldnames for column in columns), (
                f"Not all requested column ({columns}) are present in file {fn}"
            )
        for row in reader:
            if columns:
                epitope = {}
                for column in columns:
                    epitope[column] = row[column]
            else:
                epitope = row
            identifier = f"{row['HGVSp']}_{row['HLA Allele']}_{row['MT Epitope Seq']}"
            assert identifier not in epitopes, f"Duplicated epitope {identifier}"
            epitopes[identifier] = epitope
    return epitopes


def match_variant_to_epitope(
    variants: dict[str, Variant], epitopes: dict[str, dict[str, Any]]
) -> dict[str, str]:
    """Creates a mapping table from epitope to corresponding variant, based on HGVSp identifiers"""
    mapping_table = {}
    for identifier in epitopes.keys():
        variant_id = identifier.split("_")[0]
        assert variant_id in variants, f"Variant identifer {variant_id} not in variant table"
        mapping_table[identifier] = variant_id
    return mapping_table


def _mutation_locus(epitope: dict[str, Any]) -> tuple[int, int]:
    """Extract the epitope position in the protein (may not be correct for default wildtype sequence)"""
    m = Variant.PROTEIN_POSITION.match(epitope["Protein Position"])
    assert m, f"Illegal protein position in epitope {epitope}"
    start_pos = int(m.group(1))
    if m.group(3):
        end_pos = int(m.group(3))
    else:
        end_pos = start_pos
    start = start_pos - int(epitope["Peptide Length"]) + int(epitope["Sub-peptide Position"])
    end = start + int(epitope["Peptide Length"]) + (end_pos - start_pos)
    return (start, end)


def _output_sites(epitope: dict[str, Any], variant: Variant) -> str:
    """Fill cleavage data for an epitope, using the sites in the variant object"""
    sites = {"Best Cleavage Position": "NA", "Best Cleavage Score": "NA", "Cleavage Sites": []}
    start, end = _mutation_locus(epitope)
    for site in variant.sites:
        pos = site[0]
        score = site[1]
        if start <= pos and pos <= end:
            sites["Best Cleavage Position"] = pos
            sites["Best Cleavage Score"] = score
            sites["Cleavage Sites"].append(site)
    sites["Cleavage Sites"] = (
        ",".join([f"{k}:{v}" for k, v in sites["Cleavage Sites"]])
        if sites["Cleavage Sites"]
        else "NA"
    )
    return "\t".join(map(str, sites.values()))


def write_output_table(
    f: io.TextIOBase,
    variants: dict[str, Variant],
    epitopes: dict[str, dict[str, Any]],
    mapping_table: dict[str, str],
):
    """Add 3 columns to the epitope table (Best cleavage pos & score, and digest of all cleavage sites)"""
    titles = list(epitopes.values())[0].keys()
    f.write(
        "\t".join(
            list(titles) + ["Best Cleavage Position", "Best Cleavage Score", "Cleavage Sites"]
        )
        + "\n"
    )
    for e_identifier, v_identifier in mapping_table.items():
        variant = variants[v_identifier]
        epitope = epitopes[e_identifier]
        line = (
            "\t".join([epitope[title] for title in titles]) + "\t" + _output_sites(epitope, variant)
        )
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

    parser.add_argument("-n", "--netchop", nargs=1, help="Path to the netchop binary")
    parser.add_argument("--timeout", type=int, default=3600, help="Netchop command timeout")
    parser.add_argument(
        "-m", "--method", choices=("cterm", "20s"), default="cterm", help="Netchop method"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.5,
        help="Score threshold to filter cleaving sites",
    )
    parser.add_argument("-o", "--output", help="Output table filename (stdout if missing)")

    parser.add_argument(
        "variants",
        help="Somatic variants vcf file annotated with VEP unsig Frameshift & Wildtype plugins",
    )
    parser.add_argument(
        "epitopes", help="Neoepitope prediction results (*.all_epitopes.tsv or *.filtered.tsv)"
    )

    args = parser.parse_args()
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    epitopes = read_epitopes_table(args.epitopes)
    logging.info(
        f"{len(epitopes)} neo-epitope predictions have been read from file {args.epitopes}"
    )
    variants = Variant.parse_annotated_vcf(args.variants)
    logging.info(f"{len(variants)} somatic variants have been read from file {args.variants}")

    mapping_table = match_variant_to_epitope(variants, epitopes)
    variants_in_epitopes = set(mapping_table.values())
    variants = {k: v for k, v in variants.items() if k in variants_in_epitopes}
    logging.info(f"{len(variants)} somatic variants produce {len(mapping_table)} neo-epitopes")

    logging.info(f"Starting netchop runs ({args.netchop}) with {args.workers} processes")
    run_netchop(
        variants,
        args.netchop[0],
        args={"method": args.method, "threshold": args.threshold},
        tmpdir=args.tmpdir,
        clean=args.verbose,
        force=args.force,
        n_workers=args.workers,
        timeout=args.timeout,
    )
    logging.info(f"{len(variants)} netchop run completed")

    logging.info("Writing results")
    if args.output:
        f = open(args.output, "wt")
    else:
        f = sys.stdout
    write_output_table(f, variants, epitopes, mapping_table)
    if args.output:
        f.flush()
        f.close()

    logging.info("Success - all done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
