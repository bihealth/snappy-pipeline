# import packages

import json  # for writing out the json file
from cyvcf2 import VCF  # for reading vcf file
import tabix
from optparse import OptionParser


parser = OptionParser()
parser.add_option("--rawvcf", help="raw vcf file")
parser.add_option("--filtered-vcf", help="filtered vcf file")
parser.add_option(
    "--exom-bedfile",
    help="exom regions bed file",
)
parser.add_option("--padding", help="it's padding", default=0)
parser.add_option(
    "--repeated-region-bedfile",
    help="bed file contains repeated or difficult to map regions",
)
parser.add_option(
    "--output",
    help="json file contains information about vcf file",
)
(options, args) = parser.parse_args()

# from pybedtools import Interval,BedTool

# I will use for loop here
# It could be better by using multi threading.
# TODO:
# - Number of variants called (from the `*.full.vcf.gz`file) [X] Get all variants from full file
# - Number of variants that pass all filters [X] Get all variants from fitlered file
# - Number of variants passing the filters that are inside & outside of the exome bed file (with padding) [X]
# - Number of variants in repeated or other difficult to map regions [ ] Need a bed file containing repeated
# - Number of SNVs & indels [X]
# - Indel lengths statistics [X]
# - Number of variants with more than 1 ALT allele (I _think_ that this value should be 0 for all variants passing the filters) [ ]
# - Number of SNVs with minimal support (only one read supporting the variant), possibly separating inside & outside exome bed [X]
# - Number of SNVs with limited support (5 reads or less supporting the variant), possibly separating inside & outside exome bed [X]
# - Number of SNVs with strong support (at least 10 reads, VAF greater or equal to 10%) [X] VAF need to be joined
# - Number of SNVs in each mutation class (C>A, C>G, C>T, T>A, T>C, T>G), for minimal, limited, average & strong support [X] not for minimal,limited,average yet
# - Metrics on strand bias could also be collected (F1R2 & F2R1 vcf FORMAT) [ ] PAUSE
# - the VAF distribution. Perhaps we can split the variants (both SNVs & indels) into those with very little support (<1%), subclonal (<10%), the main part (<40%), affected by CNV (>40%)[ ]


# Adding A>C G>T grouping
# Table for VAF
# Genotype check > report different genotypes
# Checking number multi allelic
# Support alternative allel in the normal
# #alt in normal + average ALT + max ALT
# snappy pipeline for this

##### SNAPPY PIPELINE
## WRAPPER
# - Checking the configuration file
# - Checking input file
## PIPELINE
# - __init__ :


def check_variant_in_bed(v_chrom, v_start, v_end, bed_intervals, padding=0):
    # print(f'variant chrom: {v_chrom} - start: {v_start} - end:{v_end}')
    record = bed_intervals.query(v_chrom, v_start - padding, v_end + padding)
    if record:
        return True
    else:
        return False


# def calculate_exoms_length(bed_intervals):
#     contigs= ["chr" + str(i) for i in range(1, 23)]
#     contigs.append("chrX")
#     contigs.append("chrY")
#     length = 0
#     for contig in contigs:
#         records = bed_intervals.query(contig,0,250000000)
#         for record in records:
#             length = length + int(record[2]) - int(record[1])
#     return length


def get_variant_type(ref, alt):
    if (len(alt) == 1) and (len(ref) == len(alt)):
        return "snp"
    elif len(alt) != len(ref):
        return "indels"


def check_sp_read(variant):
    if get_variant_type(variant.REF, variant.ALT[0]) == "snp":
        dp = variant.INFO["DP"]
        if dp == 1:
            return "minimal"
        elif (dp > 1) and (dp <= 5):
            return "limited"
        if dp >= 10:
            return "strong"


def assign_class_snvs(variant, mt_mat):
    temp = str(variant.REF) + ">" + str(variant.ALT[0])
    # Following result from plot-vcfstats with bcftools
    match temp:
        case "A>C" | "C>A":
            mt_mat[0] += 1
        case "A>G" | "G>A":
            mt_mat[1] += 1
        case "A>T" | "T>A":
            mt_mat[2] += 1
        case "C>G" | "G>C":
            mt_mat[3] += 1
        case "C>T" | "T>C":
            mt_mat[4] += 1
        case "G>T" | "T>G":
            mt_mat[5] += 1
    return mt_mat


def process_vcf_file(vcf_file, bed_file="", filter_file=False, padding=0):
    total_v_count = 0
    if filter_file:
        # VAF
        vaf = []
        # numb of variants in and out of exoms
        v_inside_exom = 0
        v_outside_exom = 0
        # numb of snps and indels
        n_snps = 0
        n_indels = 0
        indels_length = []
        # support read informations
        minimal_rp_snvs_exom = 0
        minimal_rp_snvs_nexom = 0
        limited_rp_snvs_exom = 0
        limited_rp_snvs_nexom = 0
        strong_rp_snvs_exom = 0
        strong_rp_snvs_nexom = 0
        # mutation classes
        mt_classes = [0] * 6
        for variant in vcf_file:
            vaf.append(variant.format("AF")[1])
            total_v_count += 1
            # Counting variants in or out of exom and variants with different support reads
            if "chrUn" in variant.CHROM:
                v_outside_exom += 1
                if check_sp_read(variant) == "minimal":
                    minimal_rp_snvs_nexom += 1
                elif check_sp_read(variant) == "limited":
                    limited_rp_snvs_nexom += 1
                else:
                    strong_rp_snvs_nexom += 1
            elif check_variant_in_bed(
                variant.CHROM, variant.start, variant.end, bed_file, padding=0
            ):
                v_inside_exom += 1
                if check_sp_read(variant) == "minimal":
                    minimal_rp_snvs_exom += 1
                elif check_sp_read(variant) == "limited":
                    limited_rp_snvs_exom += 1
                else:
                    strong_rp_snvs_exom += 1

            else:
                v_outside_exom += 1
                if check_sp_read(variant) == "minimal":
                    minimal_rp_snvs_nexom += 1
                elif check_sp_read(variant) == "limited":
                    limited_rp_snvs_nexom += 1
                else:
                    strong_rp_snvs_nexom += 1
            # get number of snvs and indels
            # Need to check multi allelic. Users shouldn't input multi allelic vcf file.
            if get_variant_type(variant.REF, variant.ALT[0]) == "snp":
                n_snps += 1
                mt_classes = assign_class_snvs(variant, mt_classes)
            else:
                # More for indels
                n_indels += 1
                indels_length.append(abs(len(variant.REF) - len(variant.ALT[0])))

        return (
            total_v_count,
            v_inside_exom,
            v_outside_exom,
            n_snps,
            n_indels,
            indels_length,
            minimal_rp_snvs_exom,
            minimal_rp_snvs_nexom,
            limited_rp_snvs_exom,
            limited_rp_snvs_nexom,
            strong_rp_snvs_exom,
            strong_rp_snvs_nexom,
            mt_classes,
            vaf,
        )
    else:
        for variant in vcf_file:
            total_v_count += 1
        return total_v_count


def main():
    classes = ["A>C", "A>G", "A>T", "C>G", "C>T", "G>T"]
    path_full_vcf = options.rawvcf
    path_filtered_vcf = options.filtered_vcf
    path_bedfile = options.exom_bedfile
    outpath = options.output
    full_vcf = VCF(path_full_vcf)
    filter_vcf = VCF(path_filtered_vcf)
    # header = filter_vcf.raw_header
    sample = filter_vcf.samples[1]
    # read bed file
    bed_intervals = tabix.open(path_bedfile)
    full_v_count = process_vcf_file(full_vcf)
    (
        passed_v_count,
        v_inside_exom,
        v_outside_exom,
        n_snps,
        n_indels,
        indels_length,
        minimal_rp_snvs_exom,
        minimal_rp_snvs_nexom,
        limited_rp_snvs_exom,
        limited_rp_snvs_nexom,
        strong_rp_snvs_exom,
        strong_rp_snvs_nexom,
        mt_classes,
        vaf,
    ) = process_vcf_file(filter_vcf, bed_intervals, filter_file=True, padding=options.padding)

    # After one for loop through vcf file. The VCF will automatically be closed.
    # length_exoms = calculate_exoms_length(bed_intervals)
    # for variant in filter_vcf:
    #     print(variant.format("AF"))
    summary = {
        "sample": sample,
        "number_full_variants": full_v_count,
        "number_passed_variants": passed_v_count,
        "number_snvs_inside_exom": v_inside_exom,
        "number_snvs_outside_exom": v_outside_exom,
        "number_snps": n_snps,
        "number_indels": n_indels,
        "indels_length": indels_length,
        "minimal_rp_snvs_exom": minimal_rp_snvs_exom,
        "minimal_rp_snvs_nexom": minimal_rp_snvs_nexom,
        "limited_rp_snvs_exom": limited_rp_snvs_exom,
        "limited_rp_snvs_nexom": limited_rp_snvs_nexom,
        "strong_rp_snvs_exom": strong_rp_snvs_exom,
        "strong_rp_snvs_nexom": strong_rp_snvs_nexom,
        "mutation_classes": mt_classes,
        "classes": classes,
        "VAF": vaf,
    }
    json_object = json.dumps(summary, indent=4)

    with open(outpath, "w") as outfile:
        outfile.write(json_object)


if __name__ == "__main__":
    main()
