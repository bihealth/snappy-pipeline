# import packages
import argparse
import gzip
import json  # for writing out the json file

from cyvcf2 import VCF  # for reading vcf file
import tabix


def check_variant_in_bed(v_chrom, v_start, v_end, bed_intervals, padding=0):
    try:
        tmp = []
        records = bed_intervals.query(v_chrom, v_start - int(padding), v_end + int(padding))
        for record in records:
            tmp.append(record)
        if len(tmp) >= 1:
            return True
        else:
            return False
    except OSError:
        print("You should index the bed file with tabix first!")


def get_contigs_from_bed_file(bedfile):
    contigs = []
    try:
        with gzip.open(bedfile, "r") as intervals:
            for interval in intervals:
                line = interval.decode()
                tmp = line.split("\t")[0]
                if tmp not in contigs:
                    contigs.append(tmp)
        return contigs
    except OSError:
        print("Error reading the GZIP file.")


def get_variant_type(ref, alt):
    if len(ref) == len(alt):
        if len(alt) == 1:
            return "SNV"
        else:
            return "ONV"
    elif len(alt) < len(ref):
        return "indel"


def check_sp_read(variant, pos_sample, minimal, limited):
    dp = variant.format("AD")[1][pos_sample]
    if dp <= minimal:
        return "minimal"
    elif (dp > minimal) and (dp <= limited):
        return "limited"
    elif dp >= limited:
        return "strong"


def assign_class_snvs(variant, mt_mat):
    temp = str(variant.REF) + ">" + str(variant.ALT[0])
    # Following result from plot-vcfstats with bcftools
    if temp in ["A>C", "T>G"]:
        mt_mat[0] += 1
    elif temp in ["A>G", "T>C"]:
        mt_mat[1] += 1
    elif temp in ["A>T", "T>A"]:
        mt_mat[2] += 1
    elif temp in ["C>G", "G>C"]:
        mt_mat[3] += 1
    elif temp in ["C>T", "G>A"]:
        mt_mat[4] += 1
    elif temp in ["G>T", "C>A"]:
        mt_mat[5] += 1
    return mt_mat


def process_vcf_file(
    vcf_file,
    pos_sample=1,
    contigs=[],
    hard_contigs=[],
    minimal=1,
    limited=5,
    bed_file="",
    hard_file="",
    filter_file=False,
    padding=0,
    id_af="AF",
):
    infor = {
        "total_v_count": 0,
        "v_inside_exom": 0,
        "v_outside_exom": 0,
        "v_in_hard_regions": 0,
        "n_snps": 0,
        "n_onvs": 0,
        "n_indels": 0,
        "minimal_rp_snvs_exom": 0,
        "minimal_rp_snvs_nexom": 0,
        "limited_rp_snvs_exom": 0,
        "limited_rp_snvs_nexom": 0,
        "strong_rp_snvs_exom": 0,
        "strong_rp_snvs_nexom": 0,
        "indels_length": [],
        "mt_classes": [0] * 6,
        "vaf": [],
    }
    if filter_file:
        for variant in vcf_file:
            infor["total_v_count"] += 1
            # Gathering information of variants in comparison to interested regions
            if variant.CHROM in contigs:
                if len(variant.ALT) > 1:
                    raise ValueError(
                        "Multiallelic variants are not supported"
                    )  # checking for multiallelic
                infor["vaf"].append(float(variant.format(id_af)[1][0]))
                # Counting variants in or out of exom and variants with different support reads
                if check_variant_in_bed(
                    variant.CHROM, variant.start, variant.end, bed_file, padding
                ):
                    infor["v_inside_exom"] += 1
                    infor[
                        str(check_sp_read(variant, pos_sample, minimal, limited)) + "_rp_snvs_exom"
                    ] += 1
                else:
                    infor["v_outside_exom"] += 1
                    infor[
                        str(check_sp_read(variant, pos_sample, minimal, limited)) + "_rp_snvs_nexom"
                    ] += 1

                # Need to check multi allelic. Users shouldn't input multi allelic vcf file.
                if get_variant_type(variant.REF, variant.ALT[0]) == "SNV":
                    infor["n_snps"] += 1
                    infor["mt_classes"] = assign_class_snvs(variant, infor["mt_classes"])
                elif get_variant_type(variant.REF, variant.ALT[0]) == "indel":
                    # More for indels
                    infor["n_indels"] += 1
                    infor["indels_length"].append(abs(len(variant.REF) - len(variant.ALT[0])))
                elif get_variant_type(variant.REF, variant.ALT[0]) == "ONV":
                    infor["n_onvs"] += 1
            # Gathering information of variants in comparison to hard mapped regions
            else:
                infor["v_outside_exom"] += 1
            if variant.CHROM in hard_contigs:
                if check_variant_in_bed(
                    variant.CHROM, variant.start, variant.end, hard_file, padding
                ):
                    infor["v_in_hard_regions"] += 1
        return infor
    else:
        total_v_count = 0
        for variant in vcf_file:
            total_v_count += 1
        return total_v_count


parser = argparse.ArgumentParser(description="Gathering information of variants")
parser.add_argument(
    "--rawvcf",
    help="vcf file containing all called somatic variants, before filtration",
    required=True,
)
parser.add_argument(
    "--filtered-vcf", help="vcf file containing somatic variants passing all filters", required=True
)
parser.add_argument(
    "--sample",
    help="Name of tumor sample",
)
parser.add_argument("--exom-bedfile", help="exom regions bed file", required=True)
parser.add_argument("--padding", type=int, help="Size of padding around exom regions", default=0)
parser.add_argument(
    "--ignore-regions",
    help="bed file contains repeated or difficult to map regions",
)
parser.add_argument(
    "--variant-allele-frequency-id",
    default="AF",
    help="ID of allele frequency in the vcf file",
)
parser.add_argument(
    "-m",
    "--minimal",
    type=int,
    default=1,
    help="threshold for defining a variant that has minimal support reads",
)
parser.add_argument(
    "-l",
    "--limited",
    type=int,
    default=5,
    help="threshold for defining a variant that has limited support reads",
)
parser.add_argument(
    "--output",
    help="json file contains information about vcf file, if not given, the tool will print it out",
)
args = parser.parse_args()


def main():
    classes = ["T>G", "T>C", "T>A", "C>G", "C>T", "C>A"]
    # reading input
    full_vcf = VCF(args.rawvcf)
    filter_vcf = VCF(args.filtered_vcf)
    pos_sample = filter_vcf.samples.index(args.sample)
    ######################
    bed_intervals = tabix.open(args.exom_bedfile)
    contigs = get_contigs_from_bed_file(args.exom_bedfile)
    full_v_count = process_vcf_file(full_vcf)
    # read bed file
    if args.ignore_regions:
        path_hard_regions = args.ignore_regions
        hard_intervals = tabix.open(path_hard_regions)
        hard_contigs = get_contigs_from_bed_file(path_hard_regions)
        infor = process_vcf_file(
            filter_vcf,
            pos_sample,
            contigs,
            hard_contigs,
            args.minimal,
            args.limited,
            bed_intervals,
            hard_intervals,
            filter_file=True,
            padding=args.padding,
            id_af=args.variant_allele_frequency_id,
        )
    else:
        infor = process_vcf_file(
            filter_vcf,
            pos_sample,
            contigs,
            [],
            args.minimal,
            args.limited,
            bed_intervals,
            "",
            filter_file=True,
            padding=args.padding,
            id_af=args.variant_allele_frequency_id,
        )

    summary = {
        "sample": args.sample,
        "total_variants_number": full_v_count,
        "number_variants_passing_filters": infor["total_v_count"],
        "number_variants_in_enrichment_regions": infor["v_inside_exom"],
        "number_variants_outside_enrichment_regions": infor["v_outside_exom"],
        "number_variants_in_masked_regions": infor["v_in_hard_regions"],
        "number_snvs": infor["n_snps"],
        "number_onvs": infor["n_onvs"],
        "number_indels": infor["n_indels"],
        "number_variants_in_enriched_with_minimal_support": infor["minimal_rp_snvs_exom"],
        "number_variants_outside_enriched_with_minimal_support": infor["minimal_rp_snvs_nexom"],
        "number_variants_in_enriched_with_limited_support": infor["limited_rp_snvs_exom"],
        "number_variants_outside_enriched_with_limited_support": infor["limited_rp_snvs_nexom"],
        "number_variants_in_enriched_with_strong_support": infor["strong_rp_snvs_exom"],
        "number_variants_outside_enriched_with_strong_support": infor["strong_rp_snvs_nexom"],
        classes[0]: infor["mt_classes"][0],
        classes[1]: infor["mt_classes"][1],
        classes[2]: infor["mt_classes"][2],
        classes[3]: infor["mt_classes"][3],
        classes[4]: infor["mt_classes"][4],
        classes[5]: infor["mt_classes"][5],
        "indels_length": infor["indels_length"],
        "VAF": infor["vaf"],
    }
    if not args.ignore_regions:
        del summary["number_variants_in_masked_regions"]
    json_object = json.dumps(summary)
    if args.output:
        with open(args.output, "w") as outfile:
            outfile.write(json_object)
        outfile.close()
    else:
        print(json_object)


if __name__ == "__main__":
    main()
