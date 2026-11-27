# -*- coding: utf-8 -*-
"""Tests of the calls report:
the VCF reader, the varlociraptor and the generic sample columns,
the compact cells and the datavzrd configuration."""

import json
import re
from pathlib import Path

import fastjsonschema
import polars as pl
import pytest

from snappy_pipeline.workflows.variant_report.model import Linkout
from snappy_wrappers.tools.calls_report import Options, summarize, variant_type
from snappy_wrappers.tools.calls_report.varlociraptor import (
    biases,
    compact_afd,
    observations,
    phred_to_probability,
)
from snappy_wrappers.tools.calls_report.vcf import header_attributes, read_vcf
from snappy_wrappers.tools.calls_report import annotation
from snappy_wrappers.tools.qc_summary import (
    AGAINST,
    CELL_CONFIG,
    MUTED,
    SEQUENTIAL,
    STATUS,
    VEGA_LITE,
    config,
)

FIXTURES = Path(__file__).parent / "fixtures" / "qc_summary"
CALLS = FIXTURES / "calls.vcf"

#: A VCF of another caller: a multi-allelic record and a record without data for one sample
GENERIC_VCF = """##fileformat=VCFv4.2
##contig=<ID=1>
##FILTER=<ID=low_qual,Description="low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths, comma-separated">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\trs1\tG\tT,C\t50.5\tPASS\tDP=40\tGT:AD:DP:GQ\t1/2:10,20,10:40:99\t0/1:30,10,0:40:60
1\t200\t.\tA\tAT\t.\tlow_qual\tDP=10\tGT:AD:DP:GQ\t0/0:10,0:10:30\t./.:.:.:.
"""

#: Structural variants: a deletion with END, an insertion with SVLEN only, a breakend with neither
SV_VCF = """##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant (inclusive, 1-based).">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1500;SVLEN=-500\tGT\t0/1
1\t3000\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=50\tGT\t0/1
1\t5000\tbnd_a\tN\tN[2:100[\t.\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT\t0/1
"""


def _validate(summary):
    """The datavzrd configuration of ``summary`` is valid against the schema."""
    cfg = config([summary.table, *summary.views], "calls", "/data", summary.default_view)
    schema = json.loads((FIXTURES / "datavzrd-2.78.1.schema.json").read_text(encoding="utf-8"))
    fastjsonschema.validate(schema, cfg)
    return cfg


def test_read_vcf_takes_the_header_fields_and_one_record_per_alt(tmp_path):
    vcf = read_vcf(str(CALLS))

    assert vcf.samples == [
        "case001subregion-N1-DNA1-WES1",
        "case001subregionFather-N1-DNA1-WES1",
        "case001subregionMother-N1-DNA1-WES1",
    ]
    assert [name for name in vcf.info if name.startswith("PROB_")] == [
        "PROB_DENOVO",
        "PROB_PRESENT",
        "PROB_ARTIFACT",
        "PROB_ABSENT",
    ]
    assert "Description" in vcf.info["PROB_DENOVO"] and "AF" in vcf.formats
    assert len(vcf.records) == 6
    record = vcf.records[0]
    assert (record.chromosome, record.position, record.alt_index, record.alts) == (
        "12",
        23687354,
        0,
        1,
    )
    assert record.info["PROB_PRESENT"] == "0"
    assert record.samples[vcf.samples[0]]["GT"] == "0/1"

    path = tmp_path / "generic.vcf"
    path.write_text(GENERIC_VCF, encoding="utf-8")
    generic = read_vcf(str(path))

    assert generic.contigs == ["1"] and list(generic.formats) == ["GT", "AD", "DP", "GQ"]
    assert [(r.alt, r.alt_index, r.alts) for r in generic.records] == [
        ("T", 0, 2),
        ("C", 1, 2),
        ("AT", 0, 1),
    ]
    assert header_attributes('ID=AD,Number=R,Description="Allelic depths, comma-separated"') == {
        "ID": "AD",
        "Number": "R",
        "Description": "Allelic depths, comma-separated",
    }


def test_observations_split_alt_and_ref_reads_by_strength():
    records = observations("26V2B", "27V2b")

    assert records == [
        {
            "allele": "alt",
            "strength": "very strong",
            "mapping": "certain",
            "count": 26,
            "reads": 26,
        },
        {"allele": "alt", "strength": "barely", "mapping": "certain", "count": 2, "reads": 2},
        {
            "allele": "ref",
            "strength": "very strong",
            "mapping": "certain",
            "count": 27,
            "reads": -27,
        },
        {"allele": "ref", "strength": "barely", "mapping": "uncertain", "count": 2, "reads": -2},
    ]
    assert observations("", ".") == []


def test_compact_afd_thins_to_the_points_and_keeps_the_mode():
    assert compact_afd("0.000=924.70,0.500=-0.00,1.000=915.34") == "0=924.7,0.5=0.0,1=915.3"
    assert compact_afd("0.000=inf,0.500=0.00") == "0.5=0.0"
    assert compact_afd(".") == ""

    dense = ",".join(f"{i / 100:.3f}={abs(37 - i) * 3:.2f}" for i in range(101))
    entries = compact_afd(dense, 10).split(",")

    assert len(entries) == 10
    assert entries[0] == "0=111.0" and entries[-1] == "1=189.0"
    assert "0.37=0.0" in entries  # the mode survives the thinning


def test_biases_name_the_biases_found():
    assert biases({"SB": "+", "ROB": ".", "RPB": ".", "SCB": "$", "HE": ".", "ALB": "."}) == (
        "strand +, softclip $"
    )
    assert biases({"SB": ".", "ROB": "."}) == ""


def test_phred_to_probability_and_variant_type():
    assert phred_to_probability("0") == 1.0
    assert phred_to_probability("10") == pytest.approx(0.1)
    assert phred_to_probability("inf") == 0.0
    assert phred_to_probability(".") is None
    assert variant_type("A", "AT", {}) == "INS"
    assert variant_type("AT", "A", {}) == "DEL"
    assert variant_type("AT", "GC", {}) == "MNV"
    assert variant_type("N", "<DEL>", {"SVTYPE": "DEL"}) == "DEL"


def test_varlociraptor_calls_get_the_probabilities_and_the_per_sample_columns():
    summary = summarize(str(CALLS), Options(page_size=10, max_in_memory_rows=1000))
    table = summary.table
    frame = table.frame

    assert table.name == "all" and table.pinned == ("chromosome", "position")
    assert (table.page_size, table.max_in_memory_rows) == (10, 1000)
    assert frame.columns[:10] == [
        "chromosome",
        "position",
        "type",
        "ref",
        "alt",
        "event",
        "P(denovo)",
        "P(present)",
        "P(artifact)",
        "P(absent)",
    ]
    assert frame.height == 6
    row = frame.row(0, named=True)
    assert (row["chromosome"], row["position"], row["type"]) == ("12", 23687354, "SNV")
    assert row["P(present)"] == 1.0  # PHRED 0
    assert row["P(artifact)"] == 0.0  # PHRED inf
    assert row["P(denovo)"] == pytest.approx(10 ** (-2796.69 / 10))
    index = "case001subregion-N1-DNA1-WES1"
    assert (row[f"{index}: AF"], row[f"{index}: DP"], row[f"{index}: GT"]) == (0.5, 57, "0/1")
    assert row[f"{index}: observations"] == "26V2B|27V2B"
    assert row[f"{index}: AFD"] == "0=924.7,0.5=0.0,1=915.3"
    # The VCF fills neither ID, QUAL nor FILTER, so there are no such columns
    assert not {"id", "quality", "filter"} & set(frame.columns)
    assert table.columns["hints"].display == "available"

    assert table.columns["P(present)"].plot["heatmap"]["range"] == list(SEQUENTIAL)
    assert table.columns["P(artifact)"].plot["heatmap"]["range"] == list(AGAINST)
    sample = "case001subregionFather-N1-DNA1-WES1"
    assert table.columns[f"{sample}: AF"].plot["ticks"]["domain"] == [0, 1]
    assert table.columns[f"{sample}: DP"].plot["bars"]
    plot = table.columns[f"{sample}: observations"].custom_plot
    assert plot["spec"]["mark"]["type"] == "bar"
    assert plot["spec"]["$schema"] == VEGA_LITE and plot["spec"]["config"] == CELL_CONFIG
    assert 'value.split("|")' in plot["data"]
    # Every strength has its row and the cells of the column share the x scale
    encoding = plot["spec"]["encoding"]
    assert encoding["y"]["scale"]["domain"] == [
        "very strong",
        "strong",
        "positive",
        "barely",
        "equal",
    ]
    most = max(
        record["count"]
        for value in frame[f"{sample}: observations"].to_list()
        for record in observations(*value.split("|"))
    )
    assert encoding["x"]["scale"]["domain"] == [-most, most]
    assert table.columns[f"{sample}: observations"].display == "detail"  # no pedigree
    assert table.columns[f"{sample}: AFD"].custom_plot["spec"]["mark"]["type"] == "bar"
    assert table.columns[f"{sample}: AFD"].display == "detail"
    assert table.columns[f"{sample}: biases"].display == "available"

    cfg = _validate(summary)
    # The plots are in the event views; the all view leaves them out
    column = cfg["views"]["present"]["render-table"]["columns"][f"{sample}: observations"]
    assert f"{sample}: observations" not in cfg["views"]["all"]["render-table"]["columns"]
    assert json.loads(column["custom-plot"]["spec"])["width"] == 110
    assert cfg["views"]["all"]["max-in-memory-rows"] == 1000


def test_other_callers_get_genotype_frequency_and_depth_columns(tmp_path):
    path = tmp_path / "calls.vcf"
    path.write_text(GENERIC_VCF, encoding="utf-8")

    summary = summarize(str(path))
    frame = summary.table.frame

    assert frame.columns == [
        "chromosome",
        "position",
        "type",
        "ref",
        "alt",
        "S1: GT",
        "S1: AF",
        "S1: alt reads",
        "S1: DP",
        "S1: GQ",
        "S2: GT",
        "S2: AF",
        "S2: alt reads",
        "S2: DP",
        "S2: GQ",
        "id",
        "quality",
        "filter",
    ]
    assert frame["type"].to_list() == ["SNV", "SNV", "INS"]
    assert frame["S1: AF"].to_list() == pytest.approx([0.5, 0.25, 0.0])
    assert frame["S1: alt reads"].to_list() == [20, 10, 0]
    assert frame["S2: AF"].to_list()[:2] == pytest.approx([0.25, 0.0])
    assert frame["S2: AF"][2] is None and frame["S2: DP"][2] is None
    assert frame["id"].to_list() == ["rs1", "rs1", ""]
    assert frame["quality"].to_list() == [50.5, 50.5, None]
    assert frame["filter"].to_list() == ["PASS", "PASS", "low_qual"]
    assert summary.table.columns["filter"].display is None  # a filter besides PASS is shown
    assert summary.table.columns["id"].display == "available"
    assert summary.table.columns["S1: GQ"].display == "available"
    assert summary.table.columns["S1: AF"].plot["ticks"]["domain"] == [0, 1]
    _validate(summary)


def test_summary_of_an_empty_vcf_has_the_columns_and_no_rows(tmp_path):
    header = "\n".join(
        line for line in CALLS.read_text(encoding="utf-8").splitlines() if line.startswith("#")
    )
    path = tmp_path / "empty.vcf"
    path.write_text(header + "\n", encoding="utf-8")

    summary = summarize(str(path))

    assert summary.table.frame.height == 0
    assert "P(present)" in summary.table.frame.columns
    assert summary.table.frame.columns[:5] == ["chromosome", "position", "type", "ref", "alt"]
    _validate(summary)


ANNOTATED = FIXTURES / "calls.annotated.vcf"
MUTECT2 = FIXTURES / "calls.mutect2.vcf"
LINKOUTS = tuple(str(link) for link in Linkout)


def test_annotated_calls_get_gene_consequence_impact_hgvs_and_the_other_entries():
    summary = summarize(str(ANNOTATED))
    frame = summary.table.frame
    columns = summary.table.columns

    assert frame.columns[10:17] == [
        "gene",
        "consequence",
        "impact",
        "hgvs.c",
        "hgvs.p",
        "transcript",
        "biotype",
    ]
    row = frame.row(0, named=True)
    # The first entry of the highest impact for the allele, percent-decoded
    assert (row["gene"], row["consequence"], row["impact"]) == ("SOX5", "synonymous_variant", "LOW")
    assert (row["hgvs.c"], row["hgvs.p"]) == ("c.2052C>T", "p.(Ser684=)")
    assert row["transcript"] == "ENST00000309359.1" and row["exon"] == "15/15"
    assert row["gene_id"] == "ENSG00000134532" and row["hgnc_id"] == "11201"
    assert row["other annotations"].startswith("ENST00000367206.3: downstream_gene_variant; ")
    assert columns["impact"].plot["pills"]["domain"] == ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    assert columns["hgvs.p"].ellipsis == 28
    assert columns["transcript"].display == "available"
    assert columns["other annotations"].display == "available"
    assert list(columns)[-2:] == ["other annotations", "hints"]  # last in the row details
    assert columns["gene_id"].display == "hidden"
    _validate(summary)


def test_consequences_are_pills_split_at_ampersands_and_colored_by_their_own_impact():
    """A JS renderer, as datavzrd's pills stay in one piece: the terms of a call in the
    annotator's order, colored by their own impact (unknown terms as MODIFIERs), the
    ``_variant`` suffix dropped from the label with the term on hover, on one line with a ``+N``
    pill for the terms beyond the first two, and linked to the Sequence Ontology if selected."""
    summary = summarize(str(ANNOTATED))
    renderer = summary.table.columns["consequence"].custom

    assert renderer == annotation.consequence_column({}).custom
    assert renderer.startswith("function(value, row) {") and 'split("&")' in renderer
    colors = json.loads(re.search(r"const colors = (\{.*?\});", renderer).group(1))
    assert colors["stop_gained"] == STATUS["critical"]
    assert colors["missense_variant"] == STATUS["serious"]
    assert colors["synonymous_variant"] == STATUS["warning"]
    assert "intron_variant" not in colors and f'|| "{MUTED}"' in renderer
    assert 'term.endsWith("_variant") ? term.slice(0, -8) : term' in renderer
    assert 'title="${term}"' in renderer and "white-space: nowrap" in renderer
    assert "const shown = terms.length <= 3 ? terms.length : 2;" in renderer
    assert "+${rest.length}" in renderer and "const link = null;" in renderer
    assert summary.table.columns["consequence"].plot is None
    _validate(summary)

    linked = summarize(str(ANNOTATED), Options(linkouts=LINKOUTS)).table.columns["consequence"]
    assert (
        'const link = "https://www.sequenceontology.org/browser/current_release/search/?q=";'
        in linked.custom
    )
    assert '<a href="${link}${encodeURIComponent(term)}" target="_blank"' in linked.custom


def test_frequencies_and_clinvar_come_from_the_fields_mehari_adds():
    frame = summarize(str(ANNOTATED)).table.frame

    # AF = (2 hom + het + hemi) / AN; gnomAD AF is the larger of exomes and genomes
    assert frame["gnomAD exomes AF"].to_list()[:2] == pytest.approx([0.5, 1e-5])
    assert frame["gnomAD genomes AF"].to_list()[:2] == pytest.approx([0.2, 0.0])
    assert frame["gnomAD AF"].to_list()[:2] == pytest.approx([0.5, 1e-5])
    assert frame["gnomAD AF"][2] is None  # not observed
    assert frame["gnomAD exomes AF"][3] == pytest.approx((400 + 200 + 50) / 800)  # with hemi
    assert frame["gnomAD hom"].to_list()[:2] == [110, 0]
    assert frame["ClinVar"].to_list()[4] == "Pathogenic,Uncertain significance"
    assert frame["ClinVar VCV"][4] == "VCV000012345.6,VCV000067890.1"
    assert frame["clinvar_id"][4] == "12345" and frame["clinvar_id"][0] == ""
    assert frame["id"].to_list() == ["", "rs123", "", "", "rs456", ""]


def test_cadd_and_spliceai_come_from_the_fields_mehari_adds():
    summary = summarize(str(ANNOTATED))
    frame = summary.table.frame
    columns = summary.table.columns

    assert frame["CADD"].to_list() == [12.3, 3.5, None, None, 61.0, None]
    assert frame["CADD raw"].to_list()[:2] == [1.234, -0.21]
    # The largest delta score over the effects and genes; the predictions per gene on demand
    assert frame["SpliceAI"].to_list() == [0.35, 0.8, None, None, None, None]
    assert frame["SpliceAI predictions"][1] == (
        "SOX5: AG 0.00 (+1), AL 0.00 (+0), DG 0.00 (+0), DL 0.05 (-7); "
        "GENEB: AG 0.10 (-2), AL 0.00 (+0), DG 0.00 (+0), DL 0.80 (+15)"
    )
    assert frame["SpliceAI predictions"][4] == ""
    after_clinvar = list(columns).index("clinvar_id") + 1
    assert list(columns)[after_clinvar : after_clinvar + 4] == [
        "CADD",
        "CADD raw",
        "SpliceAI",
        "SpliceAI predictions",
    ]
    assert columns["CADD"].plot["ticks"]["domain"] == [0, 99]
    assert columns["CADD"].display is None
    assert columns["CADD raw"].display == "available"
    assert columns["SpliceAI predictions"].display == "available"
    assert "`CADD`: the PHRED-scaled CADD score" in summary.table.desc
    _validate(summary)


def test_the_tasks_display_modes_override_the_defaults(capsys):
    options = Options(columns={"CADD raw": "normal", "gene": "pinned", "nothing": "hidden"})
    summary = summarize(str(ANNOTATED), options)

    assert summary.table.columns["CADD raw"].display == "normal"
    assert summary.table.columns["gene"].display == "pinned"
    assert summary.table.columns["gene"].links is None  # the other settings stay
    assert "no column 'nothing'" in capsys.readouterr().err
    _validate(summary)


def test_linkouts_follow_the_assembly_and_the_selection():
    grch37 = summarize(str(ANNOTATED), Options(linkouts=LINKOUTS, assembly="GRCh37")).table
    grch38 = summarize(str(ANNOTATED), Options(linkouts=LINKOUTS, assembly="GRCh38")).table
    unknown = summarize(str(ANNOTATED), Options(linkouts=LINKOUTS)).table
    some = summarize(
        str(ANNOTATED), Options(linkouts=("genecards", "dbsnp"), assembly="GRCh38")
    ).table

    assert grch37.columns["position"].links["gnomAD"].endswith("?dataset=gnomad_r2_1")
    assert grch38.columns["position"].links["gnomAD"] == (
        "https://gnomad.broadinstitute.org/variant/{gnomad_chromosome}-{position}-{ref}-{alt}"
        "?dataset=gnomad_r4"
    )
    assert "[chrpos37]" in grch37.columns["position"].links["ClinVar"]
    assert "[Base Position Previous]" in grch37.columns["position"].links["dbSNP"]
    assert grch38.columns["gene"].links["Ensembl"].startswith("https://www.ensembl.org/")
    assert grch37.columns["gene"].links["Ensembl"].startswith("https://grch37.ensembl.org/")
    assert set(grch38.columns["gene"].links) == {"Ensembl", "GeneCards", "OMIM"}
    assert grch38.columns["id"].links == {"dbSNP": "https://www.ncbi.nlm.nih.gov/snp/{value}"}
    assert grch38.columns["ClinVar"].links == {
        "ClinVar": "https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/"
    }
    assert grch38.frame["gnomad_chromosome"][0] == "12"
    assert grch38.columns["gnomad_chromosome"].display == "hidden"
    # Without a build the position and Ensembl links are left out and the description says so
    assert unknown.columns["position"].links is None
    assert set(unknown.columns["gene"].links) == {"GeneCards", "OMIM"}
    assert "gnomad_chromosome" not in unknown.frame.columns
    assert "no annotation task in the chain" in unknown.desc
    assert set(some.columns["gene"].links) == {"GeneCards"}
    assert set(some.columns["position"].links) == {"dbSNP"}
    assert "inspect" not in some.frame.columns
    assert summarize(str(ANNOTATED), Options(linkouts=())).table.columns["gene"].links is None


def test_the_inspect_link_carries_the_observations_without_names():
    table = summarize(str(ANNOTATED), Options(linkouts=("varlociraptor_inspect",))).table

    url = table.frame["inspect_url"][0]
    assert url.startswith(
        "https://varlociraptor.github.io/varlociraptor-inspect/?PROB_DENOVO=2796.69&"
    )
    assert "&PROB_ARTIFACT=inf&" in url
    assert "&OBS_s1=" in url and "&AFD_s3=" in url and "case001" not in url
    assert "%" in url  # the observation string is URL-encoded
    assert table.columns["inspect"].links == {"varlociraptor-inspect": "{inspect_url}"}
    assert table.columns["inspect_url"].display == "hidden"
    assert table.frame["inspect"][0] == "inspect"
    assert "numbered in the order of the VCF" in table.desc


def test_annotated_mutect2_calls_mix_the_generic_sample_columns_with_the_annotation():
    summary = summarize(str(MUTECT2), Options(linkouts=LINKOUTS, assembly="GRCh38"))
    frame = summary.table.frame

    assert frame["gene"].to_list() == ["SOX5", "KRAS", "ITPR2"]
    assert frame["impact"].to_list() == ["MODIFIER", "MODERATE", "MODIFIER"]
    assert frame["hgvs.p"][1] == "p.(Gly12Asp)"
    assert "P(present)" not in frame.columns and "gnomAD AF" not in frame.columns
    assert "case001subregion-T1-DNA1-WES1: AF" in frame.columns
    assert frame["filter"].to_list() == ["PASS", "weak_evidence", "PASS"]
    assert "quality" not in frame.columns and "inspect" not in frame.columns
    _validate(summary)


def test_varlociraptor_calls_get_a_view_per_event_ranked_by_posterior_rarity_and_impact():
    summary = summarize(str(ANNOTATED), Options(linkouts=LINKOUTS, assembly="GRCh37"))
    table, views = summary.table, {view.name: view for view in summary.views}

    # One view per event that is not artifact or absent; the site opens with the first that has calls
    assert list(views) == ["denovo", "present"] and summary.default_view == "present"
    assert views["denovo"].frame.height == 0 and views["present"].frame.height == 6
    assert table.frame["event"].to_list() == ["present"] * 6
    assert table.columns["event"].plot["pills"]["domain"] == [
        "denovo",
        "present",
        "artifact",
        "absent",
    ]
    assert table.columns["event"].plot["pills"]["range"] == [
        "#eb6834",
        "#2a78d6",
        "#c3c2b7",
        "#c3c2b7",
    ]
    # Equal posteriors: unobserved variants first, then by gnomAD frequency, then impact
    present = views["present"].frame
    assert present["gnomAD AF"].to_list()[:2] == [None, None]
    assert present["gnomAD AF"].to_list()[2:] == pytest.approx([1e-5, 0.00015, 0.5, 0.8125])
    assert present["position"].to_list()[:2] == [23737566, 24982726]  # genomic order among equals
    assert "by that posterior, then gnomAD frequency (unobserved first), then impact" in (
        views["present"].desc
    )
    # Cross-links on the position, the heavy columns only in the event views and the result TSV
    assert views["present"].links == {"all calls": ("position", "all/position", False)}
    assert table.links == {
        "denovo": ("position", "denovo/position", True),
        "present": ("position", "present/position", True),
    }
    assert "inspect_url" in table.exclude and "other annotations" in table.exclude
    assert "case001subregion-N1-DNA1-WES1: observations" in table.exclude
    assert "case001subregion-N1-DNA1-WES1: AF" not in table.exclude
    assert table.frame["position"].to_list() == sorted(table.frame["position"].to_list())
    assert (
        "- The observation and density plots, the other annotations and the inspect link are "
        "in the other views" in table.desc
    )
    assert table.desc.startswith("All calls in genomic order.\n\n**Columns**\n- `event`")
    assert (
        "\n\n**Links**\n- `position`: gnomAD, ClinVar, dbSNP\n- `gene`: Ensembl, GeneCards, OMIM"
        in table.desc
    )
    assert "\n\n**Navigation**\n- " in table.desc
    cfg = _validate(summary)
    assert cfg["default-view"] == "present"
    assert "inspect" not in cfg["views"]["all"]["render-table"]["columns"]
    assert "inspect" in cfg["views"]["present"]["render-table"]["columns"]
    link = cfg["datasets"]["all"]["links"]["present"]
    assert (link["column"], link["table-row"], link["optional"]) == (
        "position",
        "present/position",
        True,
    )


def test_the_all_view_names_the_in_memory_limit_when_the_calls_exceed_it():
    summary = summarize(str(ANNOTATED), Options(max_in_memory_rows=5))

    assert "With 6 rows, above the in-memory limit of 5" in summary.table.desc
    assert "in-memory limit" not in summarize(str(ANNOTATED)).table.desc


def test_annotated_calls_without_events_get_a_ranked_view():
    summary = summarize(str(MUTECT2))
    (view,) = summary.views

    assert view.name == "ranked" and summary.default_view == "ranked"
    assert view.frame["impact"].to_list() == ["MODERATE", "MODIFIER", "MODIFIER"]
    assert view.frame["position"].to_list() == [25398284, 23893753, 26634226]
    assert "by impact, then position" in view.desc
    assert summary.table.exclude == ("other annotations",)
    assert "- The other annotations are in the other views" in summary.table.desc


def test_event_of_takes_the_most_probable_event_and_contigs_sort_by_the_header(tmp_path):
    from snappy_wrappers.tools.calls_report.views import contig_ranks, event_of, genomic_order

    assert event_of({"P(a)": 0.2, "P(b)": 0.7, "P(c)": None}, ["A", "B", "C"]) == "b"
    assert event_of({"P(a)": 0.5, "P(b)": 0.5}, ["A", "B"]) == "a"  # ties: header order
    assert event_of({}, ["A"]) == ""

    frame = pl.DataFrame({"chromosome": ["2", "X", "1", "10", "MT"], "position": [5, 1, 9, 2, 3]})
    assert contig_ranks(frame, ["1", "2", "X"]) == {"1": 0, "2": 1, "X": 2, "10": 3, "MT": 4}
    ordered = genomic_order(frame, ["1", "2", "X"])
    assert ordered["chromosome"].to_list() == ["1", "2", "X", "10", "MT"]
    natural = genomic_order(frame, [])
    assert natural["chromosome"].to_list() == ["1", "2", "10", "MT", "X"]


def test_a_vcf_without_calls_still_has_its_views():
    header = "\n".join(
        line for line in ANNOTATED.read_text(encoding="utf-8").splitlines() if line.startswith("#")
    )
    path = FIXTURES.parent / "empty_annotated.vcf"
    try:
        path.write_text(header + "\n", encoding="utf-8")
        summary = summarize(str(path), Options(linkouts=LINKOUTS, assembly="GRCh38"))
    finally:
        path.unlink()

    assert [view.name for view in summary.views] == ["denovo", "present"]
    assert summary.default_view is None  # no view has calls, so the site opens with all calls
    assert summary.table.frame.height == 0
    _validate(summary)


def test_structural_variants_get_their_length_and_a_hidden_end(tmp_path):
    path = tmp_path / "sv.vcf"
    path.write_text(SV_VCF)

    summary = summarize(str(path), Options())

    frame = summary.table.frame
    assert frame.columns[:7] == ["chromosome", "position", "type", "SV length", "end", "ref", "alt"]
    assert frame["type"].to_list() == ["DEL", "INS", "BND"]
    assert frame["SV length"].to_list() == [500, 50, None]
    assert frame["end"].to_list() == [1500, 3000, 5000]
    assert summary.table.columns["end"].display == "hidden"
    assert summary.table.columns["SV length"].precision == 0
    assert "`SV length`: the length of a structural variant" in summary.table.desc
    _validate(summary)
    # Small variants have neither column, although varlociraptor declares the SV fields
    assert "SV length" not in summarize(str(CALLS), Options()).table.frame.columns
