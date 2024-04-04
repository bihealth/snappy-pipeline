# -*- coding: utf-8 -*-

import re

import vcfpy
from snakemake import shell

inp = snakemake.input.vcf
out = snakemake.output.vcf
escape = snakemake.config["step_config"][snakemake.config["pileline_step"]["name"]]["vep"]["escape"]

# Nomenclature taken from ENSEMBL 110 (http://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
variant_classes_vep = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "__end_of_HIGH_IMPACT",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "__end_of_MODERATE_IMPACT",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "__end_of_LOW_IMPACT",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
    "__end_of_MODIFIER_IMPACT",
    "",
)

appris_codes = ("P1", "A1", "P2", "A2", "P3", "P4", "P5", "")
tsl_codes = ("1", "2", "3", "4", "5", "NA", "")
criteria = ("MANE", "Consequence", "APPRIS", "TSL")

reduced_escape_mapping = [
    (";", "%3B"),
    ("\r", "%0D"),
    ("\n", "%0A"),
    ("\t", "%09"),
]

sep = re.compile("\|")
prefix = re.compile("^Consequence annotations from Ensembl VEP. Format: ")
ampersand = re.compile("&")


def get_value(criterion, x):
    if criterion == "MANE":
        return _get_value_MANE(x)
    if criterion == "Consequence":
        return _get_value_Consequence(x)
    if criterion == "APPRIS":
        return _get_value_APPRIS(x)
    if criterion == "TSL":
        return _get_value_TSL(x)


def _get_value_MANE(x):
    if x:
        return 0
    else:
        return 1


def _get_value_Consequence(x):
    consequences = ampersand.split(x)
    worst = None
    for consequence in consequences:
        try:
            i = variant_classes_vep.index(consequence)
            if worst is None or i < worst:
                worst = i
        except ValueError:
            print("WARNING- unknown consequence {}".format(consequence))
    if worst is None:
        return len(variant_classes_vep) - 1
    else:
        return worst


def _get_value_APPRIS(x):
    try:
        return appris_codes.index(x)
    except ValueError:
        print("WARNING- unknown APPRIS code {}".format(x))
        return len(appris_codes) - 1


def _get_value_TSL(x):
    try:
        return tsl_codes.index(str(x))
    except ValueError:
        print("WARNING- unknown TSL code {}".format(x))
        return len(tsl_codes) - 1


reader = vcfpy.Reader.from_path(inp)
header = reader.header.copy()

if (
    not escape
    and header.lines[0].key == "fileformat"
    and header.lines[0].value in ("VCFv4.1", "VCFv4.2")
):
    vcfpy.record.ESCAPE_MAPPING = reduced_escape_mapping

csq = header.get_info_field_info("CSQ")
titles = sep.split(prefix.sub("", header.get_info_field_info("CSQ").description))
csq.mapping["Number"] = 1

writer = vcfpy.Writer.from_path(out, header)

for record in reader:
    codes = {
        "MANE": 1,
        "Consequence": len(variant_classes_vep) - 1,
        "APPRIS": len(appris_codes) - 1,
        "TSL": len(tsl_codes) - 1,
    }
    chosen = None
    first = None
    for oneAnnotation in record.INFO["CSQ"]:
        fields = sep.split(str(oneAnnotation))
        assert len(fields) == len(titles)
        fields = {titles[i]: fields[i] for i in range(len(fields))}

        if not first:
            first = fields

        selected = False
        for criterion in criteria:
            result = get_value(criterion, fields[criterion])
            if result is None or result == codes[criterion]:
                continue
            selected = result < codes[criterion]
            break

        if selected:
            for criterion in criteria:
                codes[criterion] = get_value(criterion, fields[criterion])
            chosen = fields.values()

    assert first is not None
    if chosen is None:
        chosen = first.values()

    record.INFO["CSQ"] = "|".join(chosen)

    writer.write_record(record)

writer.close()

shell(
    r"""
tabix {snakemake.output.vcf}

d=$(dirname {snakemake.output.vcf})
pushd
f=$(basename {snakemake.output.vcf})
md5sum $f > $f.md5
md5sum $f.tbi > $f.tbi.md5
popd
"""
)
