# -*- coding: utf-8 -*-
"""Code for testing vcf_sv_filter wrapper"""

import vcfpy

from snappy_wrappers.wrappers.vcf_sv_filter.vcf_sv_filter import (
    GenomeRegion,
    MantaGenotypeMetricsBuilder,
)


def test_manta_genotype_metrics_builder_get_length_dup_no_end():
    """Tests MantaGenotypeMetricsBuilder.get_length() - DUP no END"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=500,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("DUP")],
        QUAL=None,
        FILTER=[],
        INFO=vcfpy.OrderedDict({"SVTYPE": "DUP"}),
    )
    actual = MantaGenotypeMetricsBuilder().get_length(record)
    assert actual is None


def test_manta_genotype_metrics_builder_get_length_del():
    """Tests MantaGenotypeMetricsBuilder.get_length() - DEL"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=500,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("DEL")],
        QUAL=None,
        FILTER=[],
        INFO=vcfpy.OrderedDict({"END": 1535, "SVTYPE": "DEL"}),
    )
    expected = 1036
    actual = MantaGenotypeMetricsBuilder().get_length(record)
    assert actual == expected


def test_manta_genotype_metrics_builder_get_length_ins():
    """Tests MantaGenotypeMetricsBuilder.get_length() - INS without SVLEN"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=500,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("INS")],
        QUAL=None,
        FILTER=[],
        INFO=vcfpy.OrderedDict({"SVTYPE": "INS", "SVLEN": [1500, 9999]}),
    )
    expected = 1500  # first coordinate in list
    actual = MantaGenotypeMetricsBuilder().get_length(record)
    assert actual == expected


def test_manta_genotype_metrics_builder_get_length_ins_no_svlen():
    """Tests MantaGenotypeMetricsBuilder.get_length() - INS without SVLEN"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=500,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("INS")],
        QUAL=None,
        FILTER=[],
        INFO=vcfpy.OrderedDict({"SVTYPE": "INS"}),
    )
    expected = 0  # default value
    actual = MantaGenotypeMetricsBuilder().get_length(record)
    assert actual == expected


def test_manta_genotype_metrics_builder_get_inner_region():
    """Tests MantaGenotypeMetricsBuilder.get_inner_region()"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=500,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("DEL")],
        QUAL=None,
        FILTER=[],
        INFO=vcfpy.OrderedDict(
            {
                "END": 1535,
                "SVTYPE": "DEL",
                "SVLEN": [-1036],
                "IMPRECISE": True,
                "CIPOS": [-538, 538],
                "CIEND": [-294, 294],
                "SIZE_CLASS": "MEDIUM",
            }
        ),
    )
    expected = GenomeRegion("chr1", 1037, 1241)
    actual = MantaGenotypeMetricsBuilder().get_inner_region(record)
    assert actual == expected


def test_manta_genotype_metrics_builder_get_inner_region_no_length():
    """Tests MantaGenotypeMetricsBuilder.get_inner_region() - Length is None"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=500,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("DEL")],
        QUAL=None,
        FILTER=[],
        INFO=vcfpy.OrderedDict({"SVTYPE": "DEL"}),
    )
    actual = MantaGenotypeMetricsBuilder().get_inner_region(record)
    assert actual is None
