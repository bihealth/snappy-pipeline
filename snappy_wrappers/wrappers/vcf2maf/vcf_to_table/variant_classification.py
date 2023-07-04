# import exceptions
# from action import Action


# Mapping taken from vcftomaf.pl routine GetVariantClassification
# Github repo: https://github.com/mskcc/vcf2maf hash d13e404 2022-11-02
# No change since tag 1.6.21 hash 754d68a 2021-04-24
variant_classes_vcftomaf = (
    ("splice_acceptor_variant", "Splice_Site"),
    ("splice_donor_variant", "Splice_Site"),
    ("transcript_ablation", "Splice_Site"),
    ("exon_loss_variant", "Splice_Site"),
    ("stop_gained", "Nonsense_Mutation"),
    #
    ("stop_lost", "Nonstop_Mutation"),
    ("initiator_codon_variant", "Translation_Start_Site"),
    ("start_lost", "Translation_Start_Site"),
    #
    ("missense_variant", "Missense_Mutation"),
    ("coding_sequence_variant", "Missense_Mutation"),
    ("conservative_missense_variant", "Missense_Mutation"),
    ("rare_amino_acid_variant", "Missense_Mutation"),
    ("transcript_amplification", "Intron"),
    ("intron_variant", "Intron"),
    ("INTRAGENIC", "Intron"),
    ("intragenic_variant", "Intron"),
    ("splice_region_variant", "Splice_Region"),
    ("incomplete_terminal_codon_variant", "Silent"),
    ("synonymous_variant", "Silent"),
    ("stop_retained_variant", "Silent"),
    ("NMD_transcript_variant", "Silent"),
    ("mature_miRNA_variant", "RNA"),
    ("exon_variant", "RNA"),
    ("non_coding_exon_variant", "RNA"),
    ("non_coding_transcript_exon_variant", "RNA"),
    ("non_coding_transcript_variant", "RNA"),
    ("nc_transcript_variant", "RNA"),
    ("5_prime_UTR_variant", "5'UTR"),
    ("5_prime_UTR_premature_start_codon_gain_variant", "5'UTR"),
    ("3_prime_UTR_variant", "3'UTR"),
    ("TF_binding_site_variant", "IGR"),
    ("regulatory_region_variant", "IGR"),
    ("regulatory_region", "IGR"),
    ("intergenic_variant", "IGR"),
    ("intergenic_region", "IGR"),
    ("upstream_gene_variant", "5'Flank"),
    ("downstream_gene_variant", "3'Flank"),
)
# Missing from VEP 109 description
# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
variant_classes_vep109 = (
    ("splice_donor_5th_base_variant", "Splice_Region"),
    ("splice_donor_region_variant", "Splice_Region"),
    ("splice_polypyrimidine_tract_region", "Splice_Region"),
    ("start_retained_variant", "Silent"),
    ("TFBS_ablation", "IGR"),
    ("TFBS_amplification", "IGR"),
    ("regulatory_region_ablation", "IGR"),
    ("regulatory_region_amplification", "IGR"),
    ("feature_elongation", "Targeted_Region"),
    ("feature_truncation", "Targeted_Region"),
)
variant_classes_all = variant_classes_vcftomaf + variant_classes_vep109


def variant_classification(x, args=None):
    if x is None:
        return None
    assert isinstance(x, list) and len(x) == 4
    if x[0] is None:
        return None
    assert isinstance(x[0], list) and len(x[0]) > 0
    if x[1] is None:
        x[1] = [None] * len(x[0])
    assert isinstance(x[1], list) and len(x[1]) == len(x[0])
    assert isinstance(x[2], list) and len(x[2]) == len(x[0])
    assert isinstance(x[3], list) and len(x[3]) == len(x[0])
    tcga = list()
    for i in range(len(x[0])):
        ensembl = x[0][i]
        if ensembl is None:
            tcga.append(None)
            continue
        ensembl = ensembl.split("&")[0]
        # Convert to Sequence Ontology controlled vocabulary
        indel = x[1][i]
        if indel == "DEL":
            indel = "deletion"
        if indel == "INS":
            indel = "insertion"
        ref = x[2][i].replace("-", "")
        alt = x[3][i].replace("-", "")
        inFrame = (abs(len(ref) - len(alt)) % 3) == 0

        found = False
        for (k, v) in variant_classes_all:
            if ensembl == k:
                tcga.append(v)
                found = True
                break
        if found:
            continue

        if indel == "deletion" and (
            ensembl == "frameshift_variant"
            or (ensembl == "protein_altering_variant" and not inFrame)
        ):
            tcga.append("Frame_Shift_Del")
        elif indel == "insertion" and (
            ensembl == "frameshift_variant"
            or (ensembl == "protein_altering_variant" and not inFrame)
        ):
            tcga.append("Frame_Shift_Ins")
        elif ensembl == "inframe_deletion" or (
            ensembl == "protein_altering_variant" and inFrame and indel == "deletion"
        ):
            tcga.append("In_Frame_Del")
        elif ensembl == "inframe_insertion" or (
            ensembl == "protein_altering_variant" and inFrame and indel == "insertion"
        ):
            tcga.append("In_Frame_Ins")
        else:
            tcga.append("Targeted_Region")

    return tcga
