#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Apply soft filters to SV caller file.

Currently supported callers are Delly2 and PopDel.

The following values will be added to the ``FILTER`` columns:

- ``ALU_OVL``: variant overlaps with ALU element
- ``DB_OVL``: variant overlaps with a variant from database

The following fields will be added to the ``INFO`` column:

- ``SIZE_CLASS``: size class of the SV
- ``BEST_DB_JACCARD``: Jaccard value for best overlap with DB
- ``BEST_ALU_JACCARD``: Jaccad value for the best overlap ALU DB
- ``AFFECTED_CARRIERS``: Number of affected samples (from pedigree) that
        carry the variant
- ``UNAFFECTED_CARRIERS``: Number of unaffected samples (from pedigree) that
        carry the variant
- ``BACKGROUND_CARRIERS``: Number of samples outside family that carry the
        variant.

The following fields will be added to the ``FORMAT`` column/genotypes:

- ``FT``: list of strings with filter values for the genotype:
    - ``HET_SNVS``: too may heterozygous SNVs in deletion (only applicable
            to deletions)
    - ``MIN_PE_COUNT``: does not pass minimal paired-end support count filter
    - ``MIN_PE_AAF``: does not pass minimal paired-end support AAF filter
    - ``MIN_SR_COUNT``: does not pass minimal split read support count filter
    - ``MIN_SR_AAF``: does not pass minimal split read support AAF filter
- ``PE_AAF``: alternate allele support in fraction of paired-end read pairs
- ``PE_COUNT``: alternate allele support in number of paired-end read pairs
- ``SR_AAF``: alternate allele support in fraction of split reads
- ``SR_COUNT``: alternate allele support in number of split reads
- ``HET_SNVS``: number of heterozygous SNVs in deletion (only applicable
                to deletions)
- ``INHERITANCE``: Information on compatible inheritance modes, comma-separated
                   list of strings. As the genotyping of SVs is usually
                   problematic, the precise genotype is ignored, only carrier
                   and non-carriers are differentiated.
    - ``DE_NOVO``: Individual has parents, individual carries variants and
                   parents do not.
    - ``DOMINANT``: Individual has parents, individual carries variants and
                    exactly one parent carries the variant.
"""

# XXX: affecteds can have no entry!

import argparse
import functools
import itertools
import logging
import re
import sys

import tabix
import vcfpy

# The tools supported by the filter.
TOOLS = ("delly", "manta")
# The size classes.
SIZE_CLASSES = ("SMALL", "MEDIUM", "LARGE")
# The SV types that are considered.
SV_TYPES = ("DEL", "DUP", "INV", "INS", "BND")
# White-listed chromosomes.
_CHROMS = tuple(itertools.chain(map(str, range(1, 23)), ("X", "Y")))
CHROMS = tuple(itertools.chain(_CHROMS, ["chr" + c for c in _CHROMS]))


def append_unique(lst, elem):
    """Append element to FILTER list if not present.

    If present, ``"PASS"`` will be removed.
    """
    # Remove PASS if present.
    if "PASS" in lst:
        lst.remove("PASS")
    # Add new element.
    if elem not in lst:
        lst.append(elem)


class GenomeRegion:
    """Genome region with half-open intervals"""

    @staticmethod
    def from_bed_line(line):
        """Return GenomeRegion from BED line"""
        chrom, begin, end = line.strip().split("\t")
        return GenomeRegion(chrom, int(begin), int(end))

    @staticmethod
    def parse(user_readable):
        user_readable = user_readable.replace(",", "")
        m = re.match("^(?P<chrom>.*?):(?P<start>[0-9]+)-(?P<end>[0-9]+)$", user_readable)
        if not m:
            raise ValueError("Invalid region string: {}".format(user_readable))
        else:
            return GenomeRegion(m.group("chrom"), int(m.group("start")), int(m.group("end")))

    def __init__(self, chrom, begin, end):
        self.chrom = chrom
        self.begin = begin
        self.end = end

    def as_bed(self):
        """Return BED reprentation (half-open intervals)"""
        tpl = "{}\t{}\t{}"
        return tpl.format(self.chrom, self.begin, self.end)

    def human_readable(self):
        """Return human readable string"""
        tpl = "{}:{:,}-{:,}"
        return tpl.format(self.chrom, self.begin + 1, self.end)

    @property
    def length(self):
        """Return length"""
        return self.end - self.begin

    def __str__(self):
        tpl = "GenomeRegion({}, {}, {})"
        return tpl.format(*list(map(repr, [self.chrom, self.begin, self.end])))

    def __repr__(self):
        return str(self)

    def overlaps(self, other):
        """Return whether the region overlas with ``other``"""
        # TODO: write test for me
        if self.chrom != other.chrom:
            return False
        if self.begin < other.end and other.begin < self.end:
            return True
        else:
            return False

    def extend(self, by):
        return GenomeRegion(self.chrom, max(0, self.begin - by), self.end + by)

    def jaccard(self, other):
        """Return jaccard index between ``self`` and ``other``"""
        if self.chrom != other.chrom:
            return 0
        pos_begin = max(self.begin, other.begin)
        pos_end = min(self.end, other.end)
        if pos_end >= pos_begin:
            overlap = pos_end - pos_begin
        else:
            overlap = 0
        union = self.length + other.length - overlap
        return overlap / union

    def __eq__(self, other):
        if self is other:
            return True
        else:
            return (self.chrom, self.begin, self.end) == (other.chrom, other.begin, other.end)


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
            length = line.length
            yield "{}:{}-{}".format(name, 1, length)


class GenotypeMetrics:
    """Unified description of a variant's metrics

    Required as the variant callers use different FORMAT entries and different
    conditions for including them.
    """

    def __init__(self):
        #: split reads supporting wild-type
        self.wt_split_reads = None
        #: paired reads supporting wild-type
        self.wt_paired_reads = None
        #: split reads supporting alternative allele
        self.alt_split_reads = None
        #: paired reads supporting alternative allele
        self.alt_paired_reads = None

    def sr_aaf(self):
        """Return alternative allele fraction for SR

        Return None if there is not enough information.
        """
        if self.alt_split_reads is None or self.wt_split_reads is None:
            return None
        if self.alt_split_reads + self.wt_split_reads == 0:
            return 0
        else:
            return self.alt_split_reads / (self.alt_split_reads + self.wt_split_reads)

    def pe_aaf(self):
        """Return alternative allele fraction for PR

        Return None if there is not enough information.
        """
        if self.alt_paired_reads is None or self.wt_paired_reads is None:
            return None
        if self.alt_paired_reads + self.wt_paired_reads == 0:
            return 0
        else:
            return self.alt_paired_reads / (self.alt_paired_reads + self.wt_paired_reads)


class GenotypeMetricsBuilder:
    """Base class for helper classes that generate ``GenotypeMetrics``
    objects from calls in VCF files
    """

    def build_call_metrics(self, record):
        raise NotImplementedError("Write me!")

    def get_length(self, record):
        """Return length of the the variant

        Return `None` for non-linear changes.
        """
        raise NotImplementedError("Write me!")

    def get_inner_region(self, record):
        """Return GenomeRegion inside the CI

        In the case of an empty region, the approximate positions 10% are
        removed from each side.  None if empty.
        """
        raise NotImplementedError("Write me!")


class DellyGenotypeMetricsBuilder(GenotypeMetricsBuilder):
    """Generate ``GenotypeMetrics`` for Delly VCF records"""

    def build_call_metrics(self, call):
        result = GenotypeMetrics()
        result.wt_split_reads = call.data["RR"]
        result.alt_split_reads = call.data["RV"]
        result.wt_paired_reads = call.data["DR"]
        result.alt_paired_reads = call.data["DV"]
        return result

    def get_length(self, record):
        if record.INFO.get("SVTYPE") == "INS":
            return record.INFO.get("INSLEN", 0)
        if record.INFO.get("CHR2") and record.INFO["CHR2"] != record.CHROM:
            return None
        elif not record.INFO.get("END"):
            return None
        else:
            return record.INFO.get("END") - record.POS + 1

    def get_inner_region(self, record):
        length = self.get_length(record)
        if length is None:
            return None
        pos_begin = record.POS - 1
        pos_end = pos_begin + length
        ci_pos = record.INFO.get("CIPOS")
        ci_end = record.INFO.get("CIEND")
        if abs(ci_pos[1]) + abs(ci_end[0]) >= length:
            # confidence intervals overlap
            pos_begin += int(0.1 * length)
            pos_end -= int(0.1 * length)
            if pos_end >= pos_begin:
                return None
            else:
                return GenomeRegion(record.CHROM, pos_begin, pos_end)
        else:
            # confidence intervals don't overlap
            pos_begin += ci_pos[1]
            pos_end += ci_pos[0]  # negative
            return GenomeRegion(record.CHROM, pos_begin, pos_end)


class PopDelGenotypeMetricsBuilder(GenotypeMetricsBuilder):
    """Generate ``GenotypeMetrics`` for PopDel VCF records"""

    def build_call_metrics(self, call):
        result = GenotypeMetrics()
        lad = call.data["LAD"]
        result.wt_paired_reads = lad[0]
        result.alt_paired_reads = lad[2]
        return result

    def get_length(self, record):
        return -int(record.INFO["SVLEN"])

    def get_inner_region(self, record):
        length = self.get_length(record)
        pos_begin = record.POS - 1
        pos_end = pos_begin + length
        return GenomeRegion(record.CHROM, pos_begin, pos_end)


class MantaGenotypeMetricsBuilder(GenotypeMetricsBuilder):
    """Generate ``GenotypeMetrics`` for Manta VCF records"""

    def build_call_metrics(self, call):
        result = GenotypeMetrics()
        srs = call.data.get("SR", [0, 0])
        result.wt_split_reads = srs[0]
        result.alt_split_reads = srs[1]
        prs = call.data.get("PR", [0, 0])
        result.wt_paired_reads = prs[0]
        result.alt_paired_reads = prs[1]
        return result

    def get_length(self, record):
        if record.INFO.get("SVTYPE") == "INS":
            return record.INFO.get("SVLEN", [0])[0]
        elif not record.INFO.get("END"):
            return None
        else:
            return record.INFO.get("END") - record.POS + 1

    def get_inner_region(self, record):
        # TODO: dupe from Delly
        length = self.get_length(record)
        if length is None:
            return None
        pos_begin = record.POS - 1
        pos_end = pos_begin + length
        ci_pos = record.INFO.get("CIPOS") or [0, 0]
        ci_end = record.INFO.get("CIEND") or [0, 0]
        if abs(ci_pos[1]) + abs(ci_end[0]) >= length:
            # confidence intervals overlap
            pos_begin += int(0.1 * length)
            pos_end -= int(0.1 * length)
            if pos_end >= pos_begin:
                return None
            else:
                return GenomeRegion(record.CHROM, pos_begin, pos_end)
        else:
            # confidence intervals don't overlap
            pos_begin += ci_pos[1]
            pos_end += ci_end[0]  # negative
            return GenomeRegion(record.CHROM, pos_begin, pos_end)


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
        #: Metric builder to use, built on first usage.
        self._metric_builder = None

    @property
    def metric_builder(self):
        if not self._metric_builder:
            if "RR" in self.owner.vcf_header.format_ids():
                self._metric_builder = DellyGenotypeMetricsBuilder()
            elif "PR" in self.owner.vcf_header.format_ids():
                self._metric_builder = MantaGenotypeMetricsBuilder()
            else:
                self._metric_builder = PopDelGenotypeMetricsBuilder()
        return self._metric_builder

    def apply(self, record):
        if self.inner:
            record = self.inner.apply(record)
        return self._apply(record)

    def _setup(self):
        """Perform extended setup."""

    def _apply(self, record):
        """Actual implementation."""
        raise NotImplementedError("Implement me!")


class AnnotateSizeClassFilterStep(FilterStep):
    """Add ``INFO/SIZE_CLASS`` entry."""

    def _setup(self):
        self.active = True

    def _apply(self, record):
        # Compute length and exit if no length given.
        length = self.metric_builder.get_length(record)
        if length is None:
            return record
        # Add INFO/SVLEN if not existing.
        if "SVLEN" not in record.INFO:
            record.INFO["SVLEN"] = [length]
        # Compute length by thresholds.
        if length <= self.args.small_sv_max_size:
            record.INFO["SIZE_CLASS"] = "SMALL"
        elif length <= self.args.medium_sv_max_size:
            record.INFO["SIZE_CLASS"] = "MEDIUM"
        else:
            record.INFO["SIZE_CLASS"] = "LARGE"
        return record


class HetSnvFilterStep(FilterStep):
    """Filter by existence of het. SNVs in deletions."""

    def _setup(self):
        if self.args.small_var_vcf:
            self.active = True
            self.reader = vcfpy.Reader.from_path(self.args.small_var_vcf)

    def _apply(self, record):
        # Handle cases of non-deletion or large deletions.
        assert self.reader
        if not record.INFO.get("SIZE_CLASS") in ("SMALL", "MEDIUM"):
            return record  # ignore, LARGE or non-linear
        elif record.INFO.get("SVTYPE") != "DEL":
            return record  # ignore, not deletion
        # Get region to process.
        inner_region = self.metric_builder.get_inner_region(record)
        if not inner_region:
            return record  # empty inner region, skip
        # Check whether the deletion calls pass the het. SNV filter
        counters = dict((a.name, 0) for a in self.owner.pedigree.affecteds)
        for vcf_record in self.reader.fetch(
            inner_region.chrom, inner_region.begin, inner_region.end
        ):
            if not vcf_record.is_snv():
                continue  # only consider SNVs
            for sample in self.owner.pedigree.affecteds:
                if sample.name not in vcf_record.call_for_sample:
                    continue  # skip, sample not in file
                call = vcf_record.call_for_sample[sample.name]
                if call.gt_type != 1:
                    continue  # only consider HET calls
                if (call.data.get("GQ", 0) or 0) <= self.args.snv_min_gq:
                    continue  # ignore, GQ too bad
                if (call.data.get("DP", 0) or 0) <= self.args.snv_min_dp:
                    continue  # ignore, DP too low
                # TODO: limit on alternative allele fraction?
                counters[call.sample] += 1  # have het. SNV for sample
        # Augment FORMAT and sample data
        record.add_format("HET_SNVS", 0)
        record.add_format("FT", ["PASS"])
        for sample, num in counters.items():
            if sample not in record.call_for_sample:
                continue
            call = record.call_for_sample[sample]
            call.data["HET_SNVS"] = num
            if num >= self.args.snv_min_count:
                call = record.call_for_sample[sample]
                if "FT" not in call.data:
                    call.data["FT"] = ["PASS"]
                append_unique(call.data["FT"], "HET_SNVS")
        return record


class AlternateAlleleSupportFilterStep(FilterStep):
    def _setup(self):
        self.active = True
        self.count_mults = {
            "DEL": self.args.sv_del_count_mult,
            "DUP": self.args.sv_dup_count_mult,
            "INV": self.args.sv_inv_count_mult,
            "INS": self.args.sv_ins_count_mult,
            "BND": self.args.sv_bnd_count_mult,
        }

    def _apply(self, record):
        # Precompute multiplicator to use
        count_mult = self.count_mults.get(record.INFO["SVTYPE"], 1)
        # Add defaults to FORMAT field if necessary
        record.add_format("PE_COUNT", 0)
        record.add_format("PE_AAF", 0)
        record.add_format("SR_COUNT", 0)
        record.add_format("SR_AAF", 0)
        record.add_format("FT", ["PASS"])
        for call in record:
            if "FT" not in call.data:
                call.data["FT"] = ["PASS"]
            metrics = self.metric_builder.build_call_metrics(call)
            if metrics.pe_aaf() is not None:
                call.data["PE_AAF"] = metrics.pe_aaf()
            if metrics.alt_paired_reads is not None:
                call.data["PE_COUNT"] = metrics.alt_paired_reads
            if metrics.sr_aaf() is not None:
                call.data["SR_AAF"] = metrics.sr_aaf()
            if metrics.alt_split_reads is not None:
                call.data["SR_COUNT"] = metrics.alt_split_reads
            if call.sample in self.owner.pedigree.affecteds_names:
                if (
                    metrics.pe_aaf() is not None
                    and metrics.pe_aaf() < self.args.affected_min_pe_aaf
                ):
                    append_unique(call.data["FT"], "MIN_PE_AAF")
                if (
                    metrics.sr_aaf() is not None
                    and metrics.sr_aaf() < self.args.affected_min_sr_aaf
                ):
                    append_unique(call.data["FT"], "MIN_SR_AAF")
                if (
                    metrics.alt_paired_reads is not None
                    and metrics.alt_paired_reads < self.args.affected_min_pe_count * count_mult
                ):
                    append_unique(call.data["FT"], "MIN_PE_COUNT")
                if (
                    metrics.alt_split_reads is not None
                    and metrics.alt_split_reads < self.args.affected_min_sr_count * count_mult
                ):
                    append_unique(call.data["FT"], "MIN_SR_COUNT")
            else:
                if (
                    metrics.pe_aaf() is not None
                    and metrics.pe_aaf() > self.args.unaffected_max_pe_aaf
                ):
                    append_unique(call.data["FT"], "MAX_PE_AAF")
                if (
                    metrics.sr_aaf() is not None
                    and metrics.sr_aaf() > self.args.unaffected_max_sr_aaf
                ):
                    append_unique(call.data["FT"], "MAX_SR_AAF")
                if (
                    metrics.alt_paired_reads is not None
                    and metrics.alt_paired_reads > self.args.unaffected_max_pe_count * count_mult
                ):
                    append_unique(call.data["FT"], "MAX_PE_COUNT")
                if (
                    metrics.alt_split_reads is not None
                    and metrics.alt_split_reads > self.args.unaffected_max_sr_count * count_mult
                ):
                    append_unique(call.data["FT"], "MAX_SR_COUNT")
        return record


class InheritanceAnnoFilterStep(FilterStep):
    def _setup(self):
        self.active = True

    def _apply(self, record):
        for member in self.owner.pedigree.members:
            compatible = []
            if self._compatible_de_novo(record, member):
                compatible.append("DE_NOVO")
            if self._compatible_dominant(record, member):
                compatible.append("DOMINANT")
            if compatible:
                record.add_format("INHERITANCE", [])
                record.call_for_sample[member.name].data["INHERITANCE"] = compatible
        return record

    def _compatible_de_novo(self, record, member):
        if member.father == "0" or member.mother == "0":
            return False
        if (
            member.father not in record.call_for_sample
            or member.mother not in record.call_for_sample
        ):
            return False
        parent_variants = {
            record.call_for_sample[member.father].is_variant,
            record.call_for_sample[member.mother].is_variant,
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
            record.call_for_sample[member.father].is_variant,
            record.call_for_sample[member.mother].is_variant,
        }
        return len(parent_variants) == 2


class AnnotateAluOverlapFilterStep(FilterStep):
    """Add ``ALU_OVL`` filter to SVs overlapping with ALU over a Jaccard
    threshold.

    Only ``SMALL`` and ``MEDIUM`` variants are considered in this step
    for performance reasons.
    """

    def _setup(self):
        if self.args.alu_bed:
            self.active = True
            self.tabix = tabix.open(self.args.alu_bed)

    def _apply(self, record):
        if not self.tabix:
            return record  # no ALU database
        if not record.INFO.get("SIZE_CLASS") in ("SMALL", "MEDIUM"):
            return record  # ignore, is LARGE or non-linear
        if not record.INFO.get("SVTYPE") == "DEL":
            return record  # ignore, is not deletion
        chrom = record.CHROM
        pos_begin = record.POS - 1
        pos_end = pos_begin + self.metric_builder.get_length(record)
        sv_region = GenomeRegion(chrom, pos_begin, pos_end)
        try:
            best_jaccard = -1
            for bed_record in self.tabix.query(chrom, pos_begin, pos_end):
                record_region = GenomeRegion(bed_record[0], int(bed_record[1]), int(bed_record[2]))
                if sv_region.jaccard(record_region) > best_jaccard:
                    best_jaccard = sv_region.jaccard(record_region)
            if best_jaccard >= self.args.alu_jaccard_threshold:
                record.add_filter("ALU_OVL")
                if "ALU_OVL" not in record.FILTER:
                    record.add_filter("ALU_OVL")
            if best_jaccard >= 0:
                record.INFO["BEST_ALU_JACCARD"] = round(best_jaccard, 4)
        except tabix.TabixError:
            pass  # swallow, probably unknown contig
        return record


class AnnotateDatabaseOverlapFilterStep(FilterStep):
    """Add DB_OVL filter to SVs overlapping with DB entry over a Jaccard
    treshold

    Currently, this is only implemented for deletions.

    Also, this is only activated for small and medium deletions because of
    performance.
    """

    def _setup(self):
        if self.args.db_bed:
            self.active = True
            self.tabix = [tabix.open(p) for p in self.args.db_bed]

    def _apply(self, record):
        if not record.INFO.get("SIZE_CLASS") in ("SMALL", "MEDIUM"):
            return record  # ignore, is LARGE or non-linear
        if not record.INFO.get("SVTYPE") == "DEL":
            return record  # ignore, is not deletion
        chrom = record.CHROM
        pos_begin = record.POS - 1
        pos_end = pos_begin + self.metric_builder.get_length(record)
        sv_region = GenomeRegion(chrom, pos_begin, pos_end)
        best_jaccard = -1
        for tbx in self.tabix:
            try:
                for bed_record in tbx.query(chrom, pos_begin, pos_end):
                    record_region = GenomeRegion(
                        bed_record[0], int(bed_record[1]), int(bed_record[2])
                    )
                    if sv_region.jaccard(record_region) >= best_jaccard:
                        best_jaccard = sv_region.jaccard(record_region)
            except tabix.TabixError:
                pass  # swallow, probably unknown contig
        if best_jaccard >= self.args.db_jaccard_threshold:
            if "DB_OVL" not in record.FILTER:
                record.add_filter("DB_OVL")
        if best_jaccard >= 0:
            record.INFO["BEST_DB_JACCARD"] = round(best_jaccard, 4)
        return record


class AnnotateCarrierCountsFilterStep(FilterStep):
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
        return record


#: The filters to apply.
FILTERS = (
    # The first step is always the annotation with the size class.
    AnnotateSizeClassFilterStep,
    # Filter based on het. SNV calls within the deletion.
    HetSnvFilterStep,
    # Filter based on alternate allele support.
    AlternateAlleleSupportFilterStep,
    # Add annotation for compatibility with inheritance.
    InheritanceAnnoFilterStep,
    # Overlaps with ALUs
    AnnotateAluOverlapFilterStep,
    # Overlaps with DBV gold standard database
    AnnotateDatabaseOverlapFilterStep,
    # Annotate with carrier counts.
    AnnotateCarrierCountsFilterStep,
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
        #: The VCF header, set after construction.
        self.vcf_header = None

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
        # Record-wise FILTER entries
        #
        # VCF header FILTER entry for ALU overlap
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "ALU_OVL"),
                    (
                        "Description",
                        "Jaccard index with an ALU element >= {}".format(
                            self.args.alu_jaccard_threshold
                        ),
                    ),
                ]
            )
        )
        # VCF header FILTER entry for overlap with database
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "DB_OVL"),
                    (
                        "Description",
                        "Jaccard index with an database element >= {}".format(
                            self.args.db_jaccard_threshold
                        ),
                    ),
                ]
            )
        )
        # Genotype-wise FILTER entries
        #
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "HET_SNVS"),
                    (
                        "Description",
                        "Overlaps with at least {} heterozygous SNVs".format(
                            self.args.snv_min_count
                        ),
                    ),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MIN_PE_COUNT"),
                    ("Description", "Does not pass minimal paired read support in affected"),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MIN_PE_AAF"),
                    ("Description", "Does not pass minimal paired read support in affected"),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MIN_SR_COUNT"),
                    ("Description", "Does not pass minimal split read support in affected"),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MIN_SR_AAF"),
                    ("Description", "Does not pass minimal split read support in affected"),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MAX_PE_COUNT"),
                    (
                        "Description",
                        "Does not pass maximal paired read support in unaffected/background",
                    ),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MAX_PE_AAF"),
                    (
                        "Description",
                        "Does not pass maximal paired read support in unaffected/background",
                    ),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MAX_SR_COUNT"),
                    (
                        "Description",
                        "Does not pass maximal split read support in unaffected/background",
                    ),
                ]
            )
        )
        header.add_filter_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "MAX_SR_AAF"),
                    (
                        "Description",
                        "Does not pass maximal split read support in unaffected/background",
                    ),
                ]
            )
        )
        return header

    def _augment_info(self, header):
        """Augment header for INFO column"""
        # VCF header SVLEN giving SV size.
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "SVLEN"),
                    ("Number", "."),
                    ("Type", "Integer"),
                    ("Description", "Difference in length between REF and ALT alleles"),
                ]
            )
        )
        # VCF header SIZE_CLASS entry for describing SV size
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "SIZE_CLASS"),
                    ("Number", 1),
                    ("Type", "String"),
                    (
                        "Description",
                        ("SV size class, one of SMALL (<= {}), MEDIUM (<= {}), " "LARGE").format(
                            self.args.small_sv_max_size, self.args.medium_sv_max_size
                        ),
                    ),
                ]
            )
        )
        # Highest Jaccard value for DB overlap
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "BEST_DB_JACCARD"),
                    ("Number", 1),
                    ("Type", "Float"),
                    ("Description", "Best Jaccard value for DB overlap"),
                ]
            )
        )
        # Highest Jaccard value for ALU overlap
        header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "BEST_ALU_JACCARD"),
                    ("Number", 1),
                    ("Type", "Float"),
                    ("Description", "Best Jaccard value for ALU overlap"),
                ]
            )
        )
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
        return header

    def _augment_format(self, header):
        """Augment header for FORMAT column"""
        header = vcfpy.header_without_lines(header, (("FORMAT", "FT"),))
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "FT"),
                    ("Number", "1"),
                    ("Type", "String"),
                    ("Description", "Semicolon-separated list of filters"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "INHERITANCE"),
                    ("Number", "."),
                    ("Type", "String"),
                    ("Description", "List of compatible inheritance modes (DE_NOVO, DOMINANT)"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "HET_SNVS"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Number of overlapping heterozygous SNVs"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "PE_AAF"),
                    ("Number", 1),
                    ("Type", "Float"),
                    ("Description", "Paired-end support of variant as alternate allele fraction"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "PE_COUNT"),
                    ("Number", 1),
                    ("Type", "Float"),
                    ("Description", "Paired-end support of variant as read pair count"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "SR_AAF"),
                    ("Number", 1),
                    ("Type", "Float"),
                    ("Description", "Split-end support of variant as alternate allele fraction"),
                ]
            )
        )
        header.add_format_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "SR_COUNT"),
                    ("Number", 1),
                    ("Type", "Float"),
                    ("Description", "Split-end support of variant as read pair count"),
                ]
            )
        )
        return header

    def run(self):
        self._print_header()
        with vcfpy.Reader.from_path(self.args.input_vcf) as reader:
            # If no regions are given, fall back to all chromosomes.
            regions = self.args.regions or full_chromosomes(reader)
            # Extend header with new lines.
            header = self._augment_header(reader.header)
            # Store header in ``self.vcf_header``.
            self.vcf_header = header
            # Open the VCF writer for writing and process each region.
            with vcfpy.Writer.from_path(self.args.output_vcf, header) as writer:
                for region in regions:
                    logging.info("Processing %s", region)
                    try:
                        records = reader.fetch(region)
                    except ValueError:
                        records = []
                        logging.warning("Could not fetch records for %s", region)
                    for record in records:
                        record = self.filter_chain.apply(record)
                        writer.write_record(record)


def main(argv=None):
    parser = argparse.ArgumentParser(description="(Delly) SV VCF soft-filter application tool")

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

    # -----------------------------------------------------------------------
    group = parser.add_argument_group("SR/PR Support Configuration")
    group.add_argument(
        "--unaffected-max-pe-count",
        type=int,
        default=4,
        help=("Maximal support for unaffected samples; " "default is 4"),
    )
    group.add_argument(
        "--unaffected-max-pe-aaf",
        type=float,
        default=0.2,
        help=("Maximal support for unaffected samples; " "default is 0.2"),
    )
    group.add_argument(
        "--unaffected-max-sr-count",
        type=int,
        default=4,
        help=("Maximal support for unaffected samples; " "default is 4"),
    )
    group.add_argument(
        "--unaffected-max-sr-aaf",
        type=float,
        default=0.2,
        help=("Maximal support for unaffected samples; " "default is 0.2"),
    )
    group.add_argument(
        "--affected-min-pe-count",
        type=int,
        default=4,
        help=("Minimal support for affected samples; " "default is 4"),
    )
    group.add_argument(
        "--affected-min-pe-aaf",
        type=float,
        default=0.2,
        help=("Minimal support for affected samples; " "default is 0.2"),
    )
    group.add_argument(
        "--affected-min-sr-count",
        type=int,
        default=4,
        help=("Minimal support for affected samples; " "default is 4"),
    )
    group.add_argument(
        "--affected-min-sr-aaf",
        type=float,
        default=0.2,
        help=("Minimal support for affected samples; " "default is 0.2"),
    )

    group = parser.add_argument_group("Multipliers for PE/SR counts")
    group.add_argument(
        "--sv-del-count-mult", type=int, default=1, help="Modifying multiplicator for counting"
    )
    group.add_argument(
        "--sv-dup-count-mult", type=int, default=2, help="Modifying multiplicator for counting"
    )
    group.add_argument(
        "--sv-inv-count-mult", type=int, default=2, help="Modifying multiplicator for counting"
    )
    group.add_argument(
        "--sv-ins-count-mult", type=int, default=1, help="Modifying multiplicator for counting"
    )
    group.add_argument(
        "--sv-bnd-count-mult", type=int, default=2, help="Modifying multiplicator for counting"
    )

    # -----------------------------------------------------------------------
    group = parser.add_argument_group("SV Size Configuration")
    group.add_argument(
        "--small-sv-max-size",
        type=int,
        default=400,
        help='Maximal SV length to be classified as "SMALL"',
    )
    group.add_argument(
        "--medium-sv-max-size",
        type=int,
        default=50000,
        help='Maximal SV length to be classified as "MEDIUM"',
    )

    # -----------------------------------------------------------------------
    group = parser.add_argument_group("ALU / Variant Database Overlap")
    group.add_argument(
        "--alu-jaccard-threshold",
        type=float,
        default=0.7,
        help=("Jaccard index threshold for flagging as ALU overlap; " "default: 0.7"),
    )
    group.add_argument("--alu-bed", type=str, help="Path to ALU BED file")
    group.add_argument(
        "--db-jaccard-threshold",
        type=float,
        default=0.7,
        help=("Jaccard index threshold for flagging as DB overlap; " "default: 0.7"),
    )
    group.add_argument(
        "--db-bed",
        type=str,
        default=[],
        action="append",
        help="Path to BED database with known variants",
    )

    # -----------------------------------------------------------------------
    group = parser.add_argument_group("Het. SNV filter")
    group.add_argument("--small-var-vcf", type=str, help="Path to small variant VCF file")
    group.add_argument(
        "--snv-min-count",
        type=int,
        default=2,
        help="Smallest number of SNVs to use for subtraction",
    )
    group.add_argument(
        "--snv-min-dp", type=int, default=10, help="Small variant minimal DP for subtraction"
    )
    group.add_argument(
        "--snv-min-gq", type=int, default=60, help="Small variant minimal GQ for subtraction"
    )

    args = parser.parse_args(argv)
    args.regions = [r for lst in args.regions for r in lst]
    return VcfFilterApp(args).run()


if __name__ == "__main__":
    sys.exit(main())
