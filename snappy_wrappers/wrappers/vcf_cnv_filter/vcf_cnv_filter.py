#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Apply soft filters to CNV VCF file.

The following fields will be added to the ``INFO`` column:

- ``AFFECTED_CARRIERS``: Number of affected samples (from pedigree) that
        carry the variant
- ``UNAFFECTED_CARRIERS``: Number of unaffected samples (from pedigree) that
        carry the variant
- ``BACKGROUND_CARRIERS``: Number of samples outside family that carry the
        variant.

- ``INHERITANCE``: Information on compatible inheritance modes, comma-separated
                   list of strings. As the genotyping of SVs is usually
                   problematic, the precise genotype is ignored, only carrier
                   and non-carriers are differentiated.
    - ``DE_NOVO``: Individual has parents, individual carries variants and
                   parents do not.
    - ``DOMINANT``: Individual has parents, individual carries variants and
                    exactly one parent carries the variant.
"""

# TODO: consolidate with vcf_sv_filter.
# TODO: annotate with database overlap?
# XXX: affecteds can have no entry!

import argparse
import functools
import itertools
import logging
import sys
import warnings

import pysam
import vcfpy

# White-listed chromosomes.
_CHROMS = tuple(itertools.chain(map(str, range(1, 23)), ("X", "Y")))
CHROMS = tuple(itertools.chain(_CHROMS, ["chr" + c for c in _CHROMS]))


class PedigreeMember:
    """Representation of one PED file line"""

    UNKNOWN = "0"
    MALE = "1"
    FEMALE = "2"
    UNAFFECTED = "1"
    AFFECTED = "2"

    @classmethod
    def parse_line(klass, line):
        arr = line.strip().split("\t")
        return PedigreeMember(*arr[0:6], data=arr[:6])

    def __init__(self, family, name, father, mother, gender, disease, data=[]):
        self.family = family
        self.name = name
        self.father = father
        self.mother = mother
        self.gender = gender
        self.disease = disease
        self.data = list(data)

    def __str__(self):
        return "Pedigree({})".format(
            ", ".join(
                map(
                    str,
                    [self.family, self.name, self.father, self.mother, self.gender, self.disease],
                )
            )
        )


class Pedigree:
    """Representation of a pedigree"""

    @classmethod
    def parse(klass, f):
        """Parse from file-like object"""
        members = []
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            members.append(PedigreeMember.parse_line(line))
        return Pedigree(members)

    def __init__(self, members=[]):
        self.members = list(members)
        self.by_name = {m.name: m for m in self.members}

    @property
    @functools.lru_cache(maxsize=1)
    def affecteds(self):
        """Return list of affected individuals"""
        return [m for m in self.members if m.disease == PedigreeMember.AFFECTED]

    @property
    @functools.lru_cache(maxsize=1)
    def affecteds_names(self):
        """Return list of names of affected individuals"""
        return [m.name for m in self.affecteds]


def full_chromosomes(reader):
    """Return list of regions of all chromosomes of VCF reader."""
    for line in reader.header.get_lines("contig"):
        if line.id in CHROMS:
            name = line.id
            length = line.length or 1_000_000_000
            yield "{}:{}-{}".format(name, 1, length)


class FilterStep:
    """Base class for filter step."""

    def __init__(self, owner, args, inner=None):
        #: Owner App object
        self.owner = owner
        #: Command line arguments to pass in.
        self.args = args
        #: Inner filter to pass through, ``None`` if no inner filter.
        self.inner = inner
        #: Whether or not active by configuration.
        self.active = False
        #: Hook for implementing setup
        self._setup()

    def apply(self, record):
        if self.inner:
            for proc_record in self.inner.apply(record):
                yield from self._apply(proc_record)
        else:
            yield from self._apply(record)

    def _setup(self):
        """Perform extended setup."""

    def _apply(self, record):
        """Actual implementation."""
        raise NotImplementedError("Implement me!")


class InheritanceAnnoFilterStep(FilterStep):
    def _setup(self):
        self.active = True

    def _apply(self, record):
        for member in self.owner.pedigree.members:
            if member in record.call_for_sample:
                compatible = []
                if self._compatible_de_novo(record, member):
                    compatible.append("DE_NOVO")
                if self._compatible_dominant(record, member):
                    compatible.append("DOMINANT")
                if compatible:
                    record.add_format("INHERITANCE", [])
                    record.call_for_sample[member.name].data["INHERITANCE"] = compatible
        yield record

    def _compatible_de_novo(self, record, member):
        if member.father == "0" or member.mother == "0":
            return False
        if (
            member.father not in record.call_for_sample
            or member.mother not in record.call_for_sample
        ):
            return False
        parent_variants = {
            not record.call_for_sample[member.father].gt_type,
            not record.call_for_sample[member.mother].gt_type,
        }
        return parent_variants == set([False])

    def _compatible_dominant(self, record, member):
        if member.father == "0" or member.mother == "0":
            return False
        if (
            member.father not in record.call_for_sample
            or member.mother not in record.call_for_sample
        ):
            return False
        parent_variants = {
            not record.call_for_sample[member.father].gt_type,
            not record.call_for_sample[member.mother].gt_type,
        }
        return len(parent_variants) == 2


class AnnotateCarrierCountsFilterStep(FilterStep):
    def _setup(self):
        self.active = True
        self.tabixes = {
            anno_args["info"]: pysam.TabixFile(anno_args["path"])
            for anno_args in self.args.annotation_beds
        }

    def _apply(self, record):
        for key, tabix in self.tabixes.items():
            cnv_start = record.affected_start
            cnv_end = record.INFO.get("END", cnv_start)
            genes = []
            try:
                tbx_iter = tabix.fetch(record.CHROM, cnv_start, cnv_end)
            except ValueError:
                tbx_iter = []
            for line in tbx_iter:
                arr = line.rstrip().split("\t")
                start = int(arr[1])
                end = int(arr[2])
                if start < cnv_end and end > cnv_start:
                    genes.append(arr[3])
            record.INFO[key] = genes
        yield record


class AnnotateOverlappingGenesFilter(FilterStep):
    def _setup(self):
        self.active = True
        # Create shortcuts of samples in pedigree and samples in affected.
        self.in_pedigree = set()
        self.affected = set()
        for member in self.owner.pedigree.members:
            self.in_pedigree.add(member.name)
            if member.disease == PedigreeMember.AFFECTED:
                self.affected.add(member.name)

    def _apply(self, record):
        # Count carriers.
        ped_carriers = 0
        affected_carriers = 0
        all_carriers = 0
        for call in record:
            if not call.is_variant:
                continue  # skip
            all_carriers += 1
            if call.sample in self.in_pedigree:
                ped_carriers += 1
            if call.sample in self.affected:
                affected_carriers += 1
        # Write out counts
        record.INFO["AFFECTED_CARRIERS"] = affected_carriers
        record.INFO["UNAFFECTED_CARRIERS"] = ped_carriers - affected_carriers
        record.INFO["BACKGROUND_CARRIERS"] = all_carriers - ped_carriers
        yield record


#: The filters to apply.
FILTERS = (
    # Add annotation for compatibility with inheritance.
    InheritanceAnnoFilterStep,
    # Annotate with carrier counts.
    AnnotateCarrierCountsFilterStep,
    # Annotate with overlapping genes.
    AnnotateOverlappingGenesFilter,
)


class VcfFilterApp:
    """Container class for storing the application's state.

    We are using a class for the implementation in contrast to global
    functions because this gives us a way to store the arguments etc.
    implicitely without passing it as parameters.
    """

    def __init__(self, args, filters=FILTERS):
        #: Command line arguments.
        self.args = args
        #: Pedigree to use.
        self.pedigree = self._load_pedigree()
        # Setup the logging.
        self._setup_logging()
        #: Filter steps to run.
        self.filters = filters
        #: Setup the filter chain.
        self.filter_chain = self._build_filter_chain()

    def _load_pedigree(self):
        logging.info("Loading pedigree file %s", self.args.ped_file)
        with open(self.args.ped_file, "rt") as pedf:
            return Pedigree.parse(pedf)

    def _setup_logging(self):
        logging.basicConfig(
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
        )
        logger = logging.getLogger("")
        if self.args.verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)

    def _build_filter_chain(self):
        """Build the filter chain."""
        result = None
        for klass in self.filters:
            tmp = klass(self, self.args, result)
            logging.info("%s %s", klass, tmp.active)
            if tmp.active:
                result = tmp
        return result or (lambda x: x)

    def _print_header(self):
        logging.info("SV VCF filter")
        logging.info("Arguments: %s", self.args)

    def _process_region(self, region, writer):
        """Process a single region and write its result to the writer."""

    def _augment_header(self, header):
        """Augment header information"""
        header = self._augment_filter(header)
        header = self._augment_info(header)
        header = self._augment_format(header)
        return header

    def _augment_filter(self, header):
        """Augment header for FILTER column"""
        return header

    def _augment_info(self, header):
        """Augment header for INFO column"""
        # Information on carriers
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "AFFECTED_CARRIERS"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Number of affected samples from pedigree that are carriers"),
                ]
            )
        )
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "UNAFFECTED_CARRIERS"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Number of unaffected samples from pedigree that are carriers"),
                ]
            )
        )
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "BACKGROUND_CARRIERS"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Number of background samples that are carriers"),
                ]
            )
        )
        for anno_args in self.args.annotation_beds:
            header.add_info_line(
                vcfpy.OrderedDict(
                    [
                        ("ID", anno_args["info"]),
                        ("Number", "."),
                        ("Type", "String"),
                        ("Description", anno_args["description"]),
                    ]
                )
            )
        return header

    def _augment_format(self, header):
        """Augment header for FORMAT column"""
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "INHERITANCE"),
                    ("Number", "."),
                    ("Type", "String"),
                    ("Description", "Compatible modes of inheritance"),
                ]
            )
        )
        return header

    def run(self):
        self._print_header()
        with vcfpy.Reader.from_path(self.args.input_vcf) as reader:
            # If no regions are given, fall back to all chromosomes.
            regions = list(self.args.regions or full_chromosomes(reader))
            # Extend header with new lines.
            header = self._augment_header(reader.header)
            # Open the VCF writer for writing and process each region.
            with vcfpy.Writer.from_path(self.args.output_vcf, header) as writer:
                for region in regions:
                    logging.info("Processing %s", region)
                    try:
                        records = reader.fetch(region)
                    except ValueError:
                        logging.warning("Could not fetch records for %s", region)
                    else:
                        for record in records:
                            for proc_record in self.filter_chain.apply(record):
                                if proc_record.REF == ".":  # fix non-validating VCF
                                    proc_record.REF = "N"
                                writer.write_record(proc_record)


def main(argv=None):
    # Suppress warnings on missing "length" attribute of contig lines.
    warnings.filterwarnings(action="once")

    parser = argparse.ArgumentParser(description="Targeted CNV call annotation tool")

    # -----------------------------------------------------------------------
    group = parser.add_argument_group("General Options")
    group.add_argument("-v", "--verbose", default=0, action="count")

    group = parser.add_argument_group("Input / Output Options")
    group.add_argument("--input-vcf", required=True, help="input VCF file")
    group.add_argument("--output-vcf", help="output VCF file", default="/dev/stdout")
    group.add_argument("--ped-file", required=True, help="Path to PED file to use.")
    group.add_argument(
        "--region",
        type=str,
        required=False,
        default=[],
        action="append",
        dest="regions",
        nargs="+",
        help=("region(s) to limit analysis to"),
    )
    group.add_argument(
        "--annotation-bed",
        type=str,
        required=False,
        default=[],
        action="append",
        dest="annotation_beds",
        help="Annotate overlaps with these regions",
    )

    args = parser.parse_args(argv)

    for i, xs in enumerate(args.annotation_beds):
        keys = ("info", "description", "path")
        arr = xs.split("|")
        assert len(arr) == 3
        args.annotation_beds[i] = dict(zip(keys, arr))

    args.regions = [r for lst in args.regions for r in lst]
    return VcfFilterApp(args).run()


if __name__ == "__main__":
    sys.exit(main())
