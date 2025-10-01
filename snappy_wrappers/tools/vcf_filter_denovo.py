#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Program for filtering variants from a VCF file to *de novo* variants.

This program implements the filters and heuristics similar to the one by Wong et al. and
Besenbacher et al.
"""

import argparse
import collections
import datetime
import itertools
import logging
import os
import sys

import pysam
import vcfpy
from vcfpy.exceptions import InvalidRecordException

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers import genome_regions, sweep  # noqa: E402

# PED Parsing ======================================================================================

Individual = collections.namedtuple(
    "Individual", ["family", "name", "father", "mother", "gender", "affected"]
)

# gender/sex
Individual.UNKNOWN = "0"
Individual.MALE = "1"
Individual.FEMALE = "2"

# affection state
Individual.UNAFFECTED = "1"
Individual.AFFECTED = "2"


class Pedigree:
    """Representation of a pedigree."""

    @classmethod
    def load(klass, file):
        entries = []
        for line in file:
            entries.append(Individual(*line.strip().split("\t")[:6]))
        return Pedigree(entries)

    def __init__(self, entries):
        """Initialize with an iterable of entries."""
        self.entries = list(entries)
        self._father_of = {e.name: e.father for e in self.entries}
        self._mother_of = {e.name: e.mother for e in self.entries}
        self._by_name = {e.name: e for e in self.entries}
        self.by_family = {}
        for entry in self.entries:
            self.by_family.setdefault(entry.family, []).append(entry)

    def get_individual(self, name):
        """Return ``Individual`` object with the given name or ``None``."""
        return self._by_name.get(name)

    def get_father(self, name):
        """Return id of father, if any, otherwise ``None``."""
        result = self._father_of.get(name)
        if result == "0":
            return None
        return result

    def get_mother(self, name):
        """Return id of mother, if any, otherwise ``None``."""
        result = self._mother_of.get(name)
        if result == "0":
            return None
        return result

    def print_ped(self, file):
        for e in self.entries:
            print("\t".join(map(str, e)), file=file)

    def to_string(self):
        return "\n".join("\t".join(map(str, e)) for e in self.entries)


# Variant Annotation Related =======================================================================


class ClippedRegion(genome_regions.GenomeRegion):
    NONE = 0
    LEFT = 1
    RIGHT = 2
    BOTH = 3

    def __init__(self, chrom, begin, end, side):
        super().__init__(chrom, begin, end)
        self.side = side

    def is_left(self):
        return bool(self.side & 1)

    def is_right(self):
        return bool(self.side & 2)

    def is_both(self):
        return bool(self.side == 3)


class VcfProcessor:
    """Base class for VCF processors that augment VCF headers and modify records."""

    def __init__(self, file_in, file_out, args=None):
        #: file-like object for reading VCF file from
        self.file_in = file_in
        #: file-like object for writing VCF file to
        self.file_out = file_out
        #: ``vcfpy.Reader`` for reading VCF records
        self.reader = None
        #: ``vcfpy.Writer`` for writing VCF records
        self.writer = None
        #: Configuration
        self.args = args
        #: Number of written records
        self.count_written = 0
        #: Time of last block completion
        self.time_prev = None
        #: Position of last block completion
        self.pos_prev = (None, 0)

    def process_record(self, record):
        """Process record and write out record using ``write_record``.

        Records have to be written explicitely so it is possible to hold back records, e.g., in the
        case of BND records.
        """
        self.write_record(record)

    def done_processing(self):
        """Called after processing each record.

        This can be used for flushing functionality in case records have been held back.
        """

    def augment_header(self, header):
        """Augment header

        The default implementation fixes the meta information that was transformed into an
        ``OrderedDict``.
        """
        return header

    def write_record(self, record):
        """Write out processed record"""
        if self.count_written == 0:
            self.time_prev = datetime.datetime.now()
            self.pos_prev = (record.CHROM, record.POS)
        self.count_written += 1
        if self.count_written % 100000 == 0:
            this_time = datetime.datetime.now()
            if record.CHROM != self.pos_prev[0]:
                mbp = "?"
            else:
                mbp = (record.POS - self.pos_prev[1]) / 1000 / 1000
                mbp = "%.1f" % mbp
            spent = (this_time - self.time_prev).total_seconds()
            logging.info(
                "written %s records / %s Mbp / %s:%s in %.1fs",
                "{:,}".format(self.count_written),
                mbp,
                record.CHROM,
                "{:,}".format(record.POS),
                spent,
            )
            self.time_prev = this_time
            self.pos_prev = (record.CHROM, record.POS)
        self.writer.write_record(record)

    def run(self):
        """Perform the processing

        After processing, the written VCF file will have an augmented header and the records were
        processed with ``process_record``.
        """
        logging.info("Opening input VCF file...")
        self.reader = vcfpy.Reader.from_path(self.file_in)
        self.on_after_vcf_reader_creation()

        logging.info("Augmenting VCF header...")
        header = self.augment_header(self.reader.header)

        logging.info("Opening output VCF file...")
        self.writer = vcfpy.Writer.from_path(self.file_out, header)
        self.on_after_vcf_writer_creation()

        logging.info("Processing records...")
        self._process_records()
        self.done_processing()

    def _process_records(self):
        """Process records, either all or only from regions"""
        no = 0
        if not self.args or not self.args.regions:
            logging.info("Process whole file/genome")
            for no, record in enumerate(self._yield_skip_invalid(iter(self.reader))):
                self.process_record(record)
            logging.info("Processed %d records", no)
        else:
            for region in self.args.regions:
                logging.info("Processing region %s", region)
                region = "".join([x for x in region if x != ","])
                chrom, begin_end = region.split(":", 1)
                begin, end = begin_end.split("-", 1)
                begin, end = int(begin) - 1, int(end)
                try:
                    it = iter(self.reader.fetch(chrom, begin, end))
                    no += 1
                except ValueError as e:
                    it = []
                    logging.warning("WARNING: caught exception %s", e)
                for record in self._yield_skip_invalid(it):
                    if (record.POS - 1 >= end) or (
                        record.INFO.get("END", record.affected_end) < begin
                    ):
                        continue
                    self.process_record(record)
                logging.info("Processed %d records", no)

    def _yield_skip_invalid(self, it):
        """Wrapper that yields from iterator ``it``, errors ignored when configured so.

        We need to skip invalid records as for some reason HTSJDK sometimes writes
        out empty records. :(
        """
        if self.args and self.args.skip_invalid:
            yield from it
        else:
            while True:
                try:
                    yield next(it)
                except InvalidRecordException as e:  # bad but not fatal
                    logging.exception("Found invalid record, skipping. Exception was: %s", e)
                except StopIteration:
                    break  # stop iterating

    def on_after_vcf_reader_creation(self):
        """Hook called after VCF reader has been set"""

    def on_after_vcf_writer_creation(self):
        """Hook called after VCF writer has been set"""


class IsMnvAnnotationStep:
    """Flags MNV variants."""

    def __init__(self, processor):
        self.processor = processor
        self.args = self.processor.args
        self.pedigree = self.processor.pedigree

    def augment_header(self, header):
        """Add INFO line to ``header``."""
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "IsMNV"),
                    ("Type", "Flag"),
                    ("Number", "0"),
                    ("Description", "Valid MNV by itself"),
                ]
            )
        )
        return header

    def trim_common(self, a, b):
        while a and b and a[0] == b[0]:
            a = a[1:]
            b = b[1:]
        while a and b and a[-1] == b[-1]:
            a = a[:-1]
            b = b[:-1]
        return a, b

    def annotate(self, record):
        """Apply de novo filter for variant."""
        if len(record.ALT) != 1:
            return  # We currently ignore non-biallelic sites
        if len(record.REF) == 1 and len(record.ALT[0].value) == 1:
            return  # Nothing to do
        ref, alt = self.trim_common(record.REF, record.ALT[0].value)
        if (len(ref) >= 1 and len(alt) >= 2) or (len(ref) >= 2 and len(alt) >= 1):
            record.INFO["IsMNV"] = True


class DeNovoAnnotationStep:
    """Apply de novo filtering.

    A variant is marked as "de novo" for an individual if both parents are present and it is "0/0"
    in both parents.
    """

    def __init__(self, processor):
        self.processor = processor
        self.args = self.processor.args
        self.pedigree = self.processor.pedigree

    def augment_header(self, header):
        """Augment ``header`` with ``INFO`` lines."""
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "DeNovo"),
                    ("Type", "String"),
                    ("Number", "."),
                    ("Description", "Samples for which the variant is de novo"),
                ]
            )
        )
        return header

    def annotate(self, record):
        """Apply denovo filter for variant"""
        if len(record.ALT) != 1:
            return  # We currently ignore non-biallelic sites
        denovo_samples = []
        for sample in record.call_for_sample.keys():
            father = self.pedigree.get_father(sample)
            mother = self.pedigree.get_mother(sample)
            if not father or not mother:
                continue  # skip
            gt_sample = record.call_for_sample[sample]
            gt_father = record.call_for_sample[father]
            gt_mother = record.call_for_sample[mother]
            # Perform filtration based on genotypes alone
            if gt_sample.data["GT"] not in ("0/1", "0|1", "1/0", "1|0"):
                continue  # skip, is not HET
            if gt_father.data["GT"] not in ("0/0", "0|0"):
                continue  # skip, is not wild type
            if gt_mother.data["GT"] not in ("0/0", "0|0"):
                continue  # skip, is not wild type
            # If we reach here then the call is de novo for the child
            denovo_samples.append(sample)
        # Annotate de novo samples:
        if denovo_samples:
            record.INFO["DeNovo"] = list(sorted(denovo_samples))
            logging.info(
                "Variant at %s:%d is de novo for %s",
                record.CHROM,
                record.POS,
                record.INFO["DeNovo"],
            )


class HaplotypeProcessor:
    """Haplotype and phasing based processing

    Can be disabled with self.args.use_phase_info == False.
    """

    def __init__(self, soft_filter_proc, args):
        self.soft_filter_proc = soft_filter_proc
        self.args = args
        self.pedigree = self.soft_filter_proc.pedigree
        self.sweeper = sweep.OverlappingWindowSweeper(
            4 * self.soft_filter_proc.args.haplotype_window,
            self.soft_filter_proc.args.haplotype_window,
            self.process_record,
            vcf_record_to_genome_position,
        )
        self.father, self.mother = self._get_parents()

    def augment_header(self, header):
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "DeNovoOrigin"),
                    ("Type", "String"),
                    ("Number", "."),
                    ("Description", "Marker for de novo haplotype origin"),
                ]
            )
        )
        return header

    def _get_parents(self):
        """Return name of parents"""
        father, mother = None, None
        for entry in self.pedigree.entries:
            if entry.father != "0":
                if father is None:
                    father = entry.father
                else:
                    if father != entry.father:
                        raise Exception("Inconsistent pedigree, seen different fathers!")
            if entry.mother != "0":
                if mother is None:
                    mother = entry.mother
                else:
                    if mother != entry.mother:
                        raise Exception("Inconsistent pedigree, seen different mothers!")
        if father is None or mother is None:
            raise Exception("Missing father or mother in pedigree!")
        return father, mother

    def process_record(self, record, left, right):  # noqa: C901
        """Process record with neighborhood"""
        # TODO: much better would be to process all inner of a window at once
        # Check whether record is de novo at all
        if "DeNovo" not in record.INFO:
            return  # skip
        logging.info("Trying to phase de novo record at %s:%d", record.CHROM, record.POS)

        # Get index haplotype id
        index_hps = {}
        for name, sample in record.call_for_sample.items():
            if name == self.father and name == self.mother:
                continue  # only indices
            if not sample.data.get("HP"):
                continue  # no haplotype information
            # Store name of ALT haplotype here without "/" or "|" separator;
            # will be one of "00", "01", "10" or "11".
            sample_gt = "".join(x for x in sample.data.get("GT", "0/0") if x in "01")[:2]
            if "1" in sample_gt:
                index_hps[name] = sample.data["HP"][sample_gt.index("1")]

        if not index_hps:
            logging.info("=> halotype information is not available")
            return  # skip, no haplotype information

        # Look whether we can find the haplotype for the index in a phased
        # variant that is phase informative.  If this is the case then we can
        # use this to trace the haplotype to its origin.  We only look into
        # the haplotypes in the indices.
        origin = {}  # for each sample
        for rec in itertools.chain([record], left, right):
            logging.info("Considering %s:%s", rec.CHROM, rec.POS)
            for name, sample in rec.call_for_sample.items():
                if name == self.father or name == self.mother:
                    continue  # only indices here
                if not sample.is_phased:
                    logging.info("not phased on %s:%s", rec.CHROM, rec.POS)
                    continue  # no phasing for sample
                if not sample.data.get("HP", None):
                    logging.info("no haplotype in %s:%s", rec.CHROM, rec.POS)
                    continue  # no haplotype information for sample
                if name not in index_hps:
                    logging.info("%s not in index_hps", name)
                    continue  # no de novo HP for this index
                logging.info("-> record is phased and has haplotype")
                # If we reach here, we have a phased sample with a haplotype
                # in the current record.
                index_gt = sample.data.get("GT")
                if not index_gt:
                    continue  # no genotype, skip
                assert "/" not in index_gt, "Must not be phased"
                if index_hps[name] not in sample.data["HP"]:
                    continue  # haplotype info does not match
                # Try to infer paternal/maternal information from phasing
                try:
                    idx_paternal = 0 if self.args.phase_paternal_first else 1
                    if sample.data["HP"].index(index_hps[name]) == idx_paternal:
                        evidence = "paternal"
                    else:
                        evidence = "maternal"
                except ValueError:
                    logging.info("-> record is not in HP")
                    continue  # skip, not in HP

                logging.info("-> found evidence for %s origin", evidence)

                # Convert flags to annotation in file
                if name not in origin:
                    origin[name] = evidence
                elif origin[name] != evidence:
                    origin[name] = "inconsistent"
                logging.info(
                    "%s denovo=%s:%s-> %s",
                    name,
                    record.CHROM,
                    "{:,}".format(record.POS),
                    origin[name],
                )

        # Check whether index hps of record consistent with child
        info_value = []
        for index, rating in sorted(origin.items()):
            other = "."
            if rating == "maternal":
                other = self.mother
            elif rating == "paternal":
                other = self.father
            info_value.append("{}|{}|{}".format(index, rating, other))
        # print('info', info_value, file=sys.stderr)
        if info_value:
            record.INFO["DeNovoOrigin"] = info_value

    def push(self, record):
        if not self.args.use_phase_info:
            # Handle case of not using phase information and short-circuit
            self.soft_filter_proc.write(record)
            return
        for done_record in self.sweeper.push(record):
            self.soft_filter_proc.write(done_record)

    def finish(self):
        if not self.args.use_phase_info:
            return  # short-circuit, nothing to do
        # print('Finishing haplo...', file=sys.stderr)
        for done_record in self.sweeper.finish():
            self.soft_filter_proc.write(done_record)


def vcf_record_to_genome_position(record):
    return sweep.GenomePosition(record.CHROM, record.POS - 1)


class NeighborhoodProcessor:
    def __init__(self, soft_filter_proc, haplo_proc):
        self.soft_filter_proc = soft_filter_proc
        self.args = self.soft_filter_proc.args
        self.haplo_proc = haplo_proc
        self.sweeper = sweep.OverlappingWindowSweeper(
            3 * self.soft_filter_proc.args.exclusive_neighborhood,
            self.soft_filter_proc.args.exclusive_neighborhood,
            self.process_record,
            vcf_record_to_genome_position,
        )

    def augment_header(self, header):
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "Neighbor"),
                    ("Type", "String"),
                    ("Number", "."),
                    ("Description", "Another de novo variant in neighborhood"),
                ]
            )
        )
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "NeighborMNV"),
                    ("Type", "String"),
                    ("Number", "."),
                    ("Description", "Another de novo variant in MNV neighborhood"),
                ]
            )
        )
        return header

    def process_record(self, record, left, right):
        """Process record with tentative neighbors"""
        # print('nbr processing', file=sys.stderr)
        for other in itertools.chain(left, right):
            if self._is_neighbor(record, other, self.args.exclusive_neighborhood):
                self._flag(record, other, "Neighbor")
            if self._is_neighbor(record, other, self.args.mnv_neighborhood):
                self._flag(record, other, "NeighborMNV")

    def _is_neighbor(self, lhs, rhs, dist):
        if lhs.CHROM != rhs.CHROM:
            return False
        return abs(lhs.POS - rhs.POS) < dist

    def _flag(self, record, other, flag):
        if "DeNovo" in record.INFO and "DeNovo" in other.INFO:
            in_record = set(record.INFO["DeNovo"])
            in_other = set(other.INFO["DeNovo"])
            in_both = in_record & in_other
            if not in_both:
                return  # skip, not de novo in both
            record_n = set(record.INFO.get(flag, []))
            other_n = set(other.INFO.get(flag, []))
            record.INFO[flag] = list(sorted(record_n | in_both))
            other.INFO[flag] = list(sorted(other_n | in_both))

    def finish(self):
        logging.info("Finishing nbr...")
        for done_record in self.sweeper.finish():
            self.haplo_proc.push(done_record)
        self.haplo_proc.finish()

    def push(self, record):
        for done_record in self.sweeper.push(record):
            self.haplo_proc.push(done_record)


class BesenbacherFilterAnnotationStep:
    """Apply "Besenbacher" filter to de novo calls

    Must be after the DeNovoAnnotationStep
    """

    def __init__(self, processor):
        self.processor = processor
        self.args = self.processor.args
        self.pedigree = self.processor.pedigree
        self.warned = False

    def annotate(self, record):
        if not record.INFO.get("DeNovo"):
            return  # skip, no sample marked as de novo
        for child in record.INFO.get("DeNovo"):
            # shortcuts to names of mother and father
            father = self.pedigree.get_father(child)
            mother = self.pedigree.get_mother(child)
            # shortcuts to genotype calls
            gt_child = record.call_for_sample[child]
            gt_father = record.call_for_sample[father]
            gt_mother = record.call_for_sample[mother]
            if "AD" in gt_child.data:  # looks like GATK
                self._fix_gatk_record(record)
                if not self.annotate_gatk(record, gt_child, gt_father, gt_mother):
                    break  # added Besenbacher filter
            else:
                if not self.warned:
                    self.warned = True
                    logging.warning("WARNING: missing AD is not GATK; warning only shown once")
                    logging.warning("%s", gt_child.data)

    def _fix_gatk_record(self, record):
        """Fix GATK record in case of ``.`` occurenced, entries AD and DP."""
        for call in record.calls:
            if "AD" in call.data and call.data.get("AD") is None:
                call.data["AD"] = [0 for i in range(len(call.ALT or [])) + 1]
            if "DP" in call.data and call.data.get("DP") is None:
                call.data["DP"] = 0

    def augment_header(self, header):
        """Augment and return augmented ``header``."""
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "Besenbacher"),
                    (
                        "Description",
                        "Does not satisfy Besenbacher filter (only considered for tentative de novo calls",
                    ),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [("ID", "InParent"), ("Description", "Variant also seen in a parent")]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [("ID", "InFamily"), ("Description", "Variant also seen in another family member")]
            )
        )
        return header

    def annotate_gatk(self, record, gt_child, gt_father, gt_mother):
        besenbacher = False
        # perform Besenbacher filtering
        if not gt_child.data["DP"] or gt_child.data["GQ"] is None:
            return False
        if (
            gt_child.data["DP"] < self.args.min_dp
            or gt_child.data["DP"] > self.args.max_dp
            or gt_child.data["GQ"] < self.args.min_gq
        ):
            record.add_filter("Besenbacher")
            besenbacher = True
        assert gt_child.data["DP"] > 0  # the following might otherwise fail
        ab_child = gt_child.data["AD"][1] / sum(gt_child.data["AD"])
        if ab_child < self.args.min_ab or ab_child > self.args.max_ab:
            if not besenbacher:
                record.add_filter("Besenbacher")
                besenbacher = True
        for gt in record.calls:
            if gt.sample != gt_child.sample:
                if gt.data["AD"] and gt.data["AD"][1] > 0:
                    if record.FILTER and "InFamily" not in record.FILTER:
                        record.add_filter("InFamily")
                    break
        for gt in [gt_father, gt_mother]:
            if gt.data["AD"] and gt.data["AD"][1] > 0:
                if record.FILTER and "InParent" not in record.FILTER:
                    record.add_filter("InParent")
            if any(gt.data[key] is None for key in ("DP", "GQ", "AD")):
                if not besenbacher:
                    record.add_filter("Besenbacher")
                return False
            if (
                gt.data["DP"] < self.args.min_dp
                or gt.data["DP"] > self.args.max_dp
                or gt.data["GQ"] < self.args.min_gq
                or gt.data["AD"][1] > self.args.max_ad2
            ):
                if not besenbacher:
                    record.add_filter("Besenbacher")
                    besenbacher = True
        return not besenbacher


class SoftFilterProcessor(VcfProcessor):
    """VcfProcessor for the soft-annotation of filters.

    The following annotation/soft-filtration will be performed:

    ``INFO/DeNovo``
        List of individuals for which the variant is de novo based on their parent calls (both
        parent calls must be available.

    ``INFO/Neighbor``
        List of individuals for which there is a de novo variant within 1kb (default) of the
        current variant.

    ``INFO/NeighborMNV``
        List of individuals for which there is a de novo variant within 20bp (default) of the
        current variant.

    ``INFO/ClippedStack``
        for de novo indels that presumably are caused by a BWA-MEM clipped read stack, require
        clipping on both sides.

    ``INFO/HalfClippedStack``
        for de novo indels that presumably are caused by a BWA-MEM clipped read stack, require
        clipping on one side only.

    ``INFO/DeNovoOrigin``
        If the variant haplotype of a de novo variant was found in mother/father/both, marks
        the variant as maternal/paternal/inconsistent in the form of ``"${offspring}|${class}"``,
        e.g. ``"child1|maternal,child2|inconsistent,child3|paternal"``.

    ``FILTER``
        ``Besenbacher`` if Besenbacher filter fails.
    """

    def __init__(self, args):
        super().__init__(args.input_vcf, args.output_vcf, args)
        # Setup the logging.
        self.offspring_bam_file = pysam.AlignmentFile(args.offspring_bam, "rb")
        self._setup_logging()
        # Load pedigree.
        self.pedigree = self.load_pedigree()
        self.counter = 0
        self.steps = (
            IsMnvAnnotationStep(self),
            DeNovoAnnotationStep(self),
            BesenbacherFilterAnnotationStep(self),
        )
        #: Helper for haplotype determination, will write back to this
        #: processor.
        self.haplo_proc = HaplotypeProcessor(self, self.args)
        #: Helper for neighborhood determination, will write into
        #: self.haplotype_proc.
        self.neighbor_proc = NeighborhoodProcessor(self, self.haplo_proc)
        #: (chrom, pos) of haplotype sweeping window start
        self.haplo_sweep_chrom = None
        self.haplo_sweep_pos = None

    def _setup_logging(self):
        logging.basicConfig(
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
        )
        logger = logging.getLogger("")  # root logger
        if self.args.verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)

    def load_pedigree(self):
        logging.info("Loading pedigree...")
        with open(self.args.input_ped, "rt") as inputf:
            pedigree = Pedigree.load(inputf)
        logging.info("Pedigree:\n%s", pedigree.to_string())
        return pedigree

    def augment_header(self, header):
        """Augment meta/header information in reader"""
        # TODO: can we use pedigree from VCF in the future?
        header = super().augment_header(header)
        header = self.haplo_proc.augment_header(header)
        header = self.neighbor_proc.augment_header(header)
        header = self._add_infos(header)
        return header

    def _add_infos(self, header):
        """Add INFO header records"""
        # Register the INFO headers that we will add.
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "ClippedStack"),
                    ("Type", "Flag"),
                    ("Number", 0),
                    (
                        "Description",
                        "Presumably caused by BWA-MEM clipped read stack for index sample (two-sided "
                        "clipping), ONLY ANNOTATED FOR INDEX",
                    ),
                ]
            )
        )
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "HalfClippedStack"),
                    ("Type", "Flag"),
                    ("Number", 0),
                    (
                        "Description",
                        "Presumably caused by BWA-MEM clipped read stack for index sample (at least "
                        "one-sided clipping), ONLY ANNOTATED FOR INDEX",
                    ),
                ]
            )
        )
        # Register headers from sub steps.
        for step in self.steps:
            header = step.augment_header(header)
        return header

    def _add_pedigree(self, reader):
        """Add PEDIGREE header records, if not already present"""
        if reader.metadata.get("PEDIGREE"):
            return  # information already there
        # Extend META fields
        # TODO(holtgrewe): Don't append already existing one
        metas = [
            collections.OrderedDict(
                [
                    ("ID", "Sex"),
                    ("Type", "String"),
                    ("Number", "1"),
                    # ('Values', 'UNKNOWN,MALE,FEMALE'),
                ]
            ),
            collections.OrderedDict(
                [
                    ("ID", "Disease"),
                    ("Type", "String"),
                    ("Number", "1"),
                    # ('Values', 'UNKNOWN,UNAFFECTED,AFFECTED'),
                ]
            ),
            collections.OrderedDict([("ID", "Father"), ("Type", "String"), ("Number", "1")]),
            collections.OrderedDict([("ID", "Mother"), ("Type", "String"), ("Number", "1")]),
        ]
        reader.metadata.setdefault("META", [])
        reader.metadata["META"] += metas
        # Write out PEDIGREE fields
        entries = []
        for entry in self.pedigree.entries:
            val = collections.OrderedDict()
            val["ID"] = entry.name
            if entry.father != Individual.UNKNOWN:
                val["Father"] = entry.father
            if entry.mother != Individual.UNKNOWN:
                val["Mother"] = entry.mother
            if entry.gender == Individual.UNKNOWN:
                val["Sex"] = "UNKNOWN"
            elif entry.gender == Individual.MALE:
                val["Sex"] = "MALE"
            else:
                val["Sex"] = "FEMALE"
            if entry.affected == Individual.AFFECTED:
                val["Affected"] = "AFFECTED"
            elif entry.affected == Individual.UNAFFECTED:
                val["Affected"] = "UNAFFECTED"
            else:
                val["Affected"] = "UNKNOWN"
            entries.append(val)
        reader.metadata["PEDIGREE"] = entries

    def process_record(self, record):
        for step in self.steps:
            step.annotate(record)
        # Add ClippedStack info for de novos
        self._flag_read_stack(record)
        # Put into neighborhood and haplotype processing pipeline
        self.neighbor_proc.push(record)

    def write(self, obj):
        """Write to output through super class process_record()"""
        super().process_record(obj)

    def done_processing(self):
        self.neighbor_proc.finish()
        self.writer.close()

    def _flag_read_stack(self, record):
        if self.args.index_name in record.INFO.get("DeNovo", []):
            bam = self.offspring_bam_file
            # Collect interesting BAM records that have clipping
            RADIUS = 100
            FUZZ = 10
            bam_records = bam.fetch(record.CHROM, max(0, record.POS - RADIUS), record.POS + RADIUS)
            regions = self._bam_records_to_clipped_regions(bam_records)
            stacks = self._stack_regions(regions)
            # print('RECORD', record, file=sys.stderr)
            # print('STACKS', stacks, file=sys.stderr)
            # import sys; print(stacks, file=sys.stderr)
            for stack in stacks:
                num_clipped, num_half_clipped = 0, 0
                for region in stack:
                    # print('BOTH', region.is_both(), file=sys.stderr)
                    # print('LEFT', region.is_left(), file=sys.stderr)
                    # print('RIGHT', region.is_right(), file=sys.stderr)
                    if region.is_both():
                        num_half_clipped += 1
                        num_clipped += 1
                    elif region.is_left() and region.begin + FUZZ >= record.POS:
                        num_half_clipped += 1
                    elif region.is_right() and region.end < record.POS + FUZZ:
                        num_half_clipped += 1
                # print('CLIPPING', num_half_clipped, num_clipped,
                #       file=sys.stderr)
                # Add flag, depending on clipped count
                MIN_CLIPPED = 3
                if num_clipped >= MIN_CLIPPED:
                    record.INFO["ClippedStack"] = True
                if num_half_clipped >= MIN_CLIPPED:
                    record.INFO["HalfClippedStack"] = True

    def _stack_regions(self, clipped_regions):
        """Stack regions using a greedy algorithm"""
        current = []
        result = [current]
        for region in clipped_regions:
            if not current:
                current.append(region)
            elif region.jaccard(current[0]) > 0.5:
                current.append(region)
            else:
                current = [region]
                result.append(current)
        if not result[-1]:
            result.pop()
        return result

    def _bam_records_to_clipped_regions(self, bam_records):
        """Convert list of BAM records to ClippedRegion objects"""
        MASK_IGNORE = (
            pysam.FUNMAP | pysam.FSECONDARY | pysam.FQCFAIL | pysam.FDUP | pysam.FSUPPLEMENTARY
        )
        CLIPPING_THRESH = 20  # ignore smaller clippings
        result = []
        for bam_record in bam_records:
            if bam_record.flag & MASK_IGNORE:
                continue  # skip unmapped
            clippings = [
                (operation, length)
                for (operation, length) in bam_record.cigartuples
                if operation in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP) and length > CLIPPING_THRESH
            ]
            side = ClippedRegion.NONE
            if len(clippings) == 2:
                side = ClippedRegion.BOTH
            elif len(clippings) == 1:
                if bam_record.cigartuples[0][0] in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP):
                    side = ClippedRegion.LEFT
                else:
                    side = ClippedRegion.RIGHT
            if side != ClippedRegion.NONE:
                result.append(
                    ClippedRegion(
                        bam_record.reference_name, bam_record.pos, bam_record.reference_end, side
                    )
                )
        return result


def run(args):
    """Program main entry point after parsing command line arguments."""
    time_start = datetime.datetime.now()
    processor = SoftFilterProcessor(args)
    processor.run()
    time_end = datetime.datetime.now()
    spent = time_end - time_start
    logging.info("Spent %.1f s in annotation", spent.total_seconds())


def main(argv=None):
    """Program's main entry point (before parsing command line arguments)."""

    # Setup command line parser
    parser = argparse.ArgumentParser(
        description=(
            "Annotate VCF file with various of soft filters for sequence variant de novo filtration"
        )
    )

    parser.add_argument(
        "--verbose",
        "-v",
        dest="verbose",
        default=False,
        action="store_true",
        help="Enable verbose logging",
    )

    group = parser.add_argument_group("Input / Output")
    group.add_argument("--index-name", required=True, help="Name of expected index library")
    group.add_argument("--input-vcf", required=True, help="Input VCF file")
    group.add_argument("--input-ped", required=True, help="Input PED file")
    group.add_argument("--output-vcf", required=True, help="Output VCF file (NOT bgziped)")
    group.add_argument(
        "--offspring-bam",
        type=str,
        required=True,
        help="Path to BAM file for offspring, for read stack artifact filter of indels",
    )
    group.add_argument(
        "--region",
        dest="regions",
        type=str,
        action="append",
        default=[],
        nargs="+",
        help="(List) of region(s) to process",
    )
    group.add_argument(
        "--skip-invalid",
        default=False,
        action="store_true",
        help="Ignore invalid VCF records, (default: false)",
    )

    group = parser.add_argument_group("INFO/FORMAT-based annotation (Besenbacher)")
    group.add_argument(
        "--min-gq", default=50, type=float, help="Minimal GQ quality for Besenbacher filter"
    )
    group.add_argument("--min-dp", default=10, type=int, help="Minimal DP for Besenbacher filter")
    group.add_argument("--max-dp", default=120, type=int, help="Minimal DP for Besenbacher filter")
    group.add_argument(
        "--min-allele-balance",
        dest="min_ab",
        default=0.3,
        type=float,
        help="Minimal alelle balance value",
    )
    group.add_argument(
        "--max-allele-balance",
        dest="max_ab",
        default=0.7,
        type=float,
        help="Maximal allele balance value",
    )
    group.add_argument(
        "--max-ad2", default=0, type=int, help="Maximal alternative observation for parent"
    )

    group = parser.add_argument_group("Neighborhood-based filter")
    group.add_argument(
        "--exclusive-neighborhood",
        type=int,
        default=1000,
        help="If there is more than one variant in such a neighborhood then flag it with Neighbor",
    )
    group.add_argument(
        "--mnv-neighborhood",
        type=int,
        default=20,
        help="Flag if there is more than one in this neighbourhood with NeighborMNV",
    )

    group = parser.add_argument_group("Phasing-based annotation")
    group.add_argument(
        "--use-phase-info",
        action="store_true",
        default=False,
        help="Use phasing information, default is to not use phasing or haplotype information",
    )
    group.add_argument(
        "--haplotype-window",
        type=int,
        default=100000,
        help="Haplotype window size in bp, defaults to 100kbp",
    )
    group.add_argument(
        "--phase-maternal-first",
        default=False,
        dest="phase_paternal_first",
        action="store_true",
        help="Phased data shows maternal allele first default is paternal allele first",
    )

    # Parse arguments, postprocess and kick-off program.
    args = parser.parse_args(argv)
    args.regions = [item for sublist in args.regions for item in sublist]

    run(args)


if __name__ == "__main__":
    sys.exit(main())
