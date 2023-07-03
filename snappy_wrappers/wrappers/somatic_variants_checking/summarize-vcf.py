# import packages
import argparse
import gzip
import json  # for writing out the json file

from cyvcf2 import VCF  # for reading vcf file
import tabix

parser = argparse.ArgumentParser(description="Gathering information of variants")
parser.add_argument(
    "--rawvcf",
    help="raw vcf file",
)
parser.add_argument(
    "--filtered-vcf",
    help="filtered vcf file",
)
parser.add_argument(
    "--exom-bedfile",
    help="exom regions bed file",
)
parser.add_argument("--padding", help="it's padding", default=0)
parser.add_argument(
    "--ignore-regions",
    help="bed file contains repeated or difficult to map regions",
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
    help="json file contains information about vcf file",
)
args = parser.parse_args()


def check_variant_in_bed(v_chrom, v_start, v_end, bed_intervals, padding=0):
    try:
        tmp = []
        records = bed_intervals.query(v_chrom, v_start - padding, v_end + padding)
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
    if (len(alt) == 1) and (len(ref) == len(alt)):
        return "snp"
    elif len(alt) != len(ref):
        return "indels"


def check_sp_read(variant, minimal, limited):
    dp = variant.INFO["DP"]
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
    contigs=[],
    hard_contigs=[],
    minimal=1,
    limited=5,
    bed_file="",
    hard_file="",
    filter_file=False,
    padding=0,
):
    infor = {
        "total_v_count": 0,
        "v_inside_exom": 0,
        "v_outside_exom": 0,
        "v_in_hard_regions": 0,
        "n_snps": 0,
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
                infor["vaf"].append(float(variant.format("AF")[1][0]))
                # Counting variants in or out of exom and variants with different support reads
                if check_variant_in_bed(
                    variant.CHROM, variant.start, variant.end, bed_file, padding=0
                ):
                    infor["v_inside_exom"] += 1
                    infor[str(check_sp_read(variant, minimal, limited)) + "_rp_snvs_exom"] += 1
                else:
                    infor["v_outside_exom"] += 1
                    infor[str(check_sp_read(variant, minimal, limited)) + "_rp_snvs_nexom"] += 1

                # Need to check multi allelic. Users shouldn't input multi allelic vcf file.
                if get_variant_type(variant.REF, variant.ALT[0]) == "snp":
                    infor["n_snps"] += 1
                    infor["mt_classes"] = assign_class_snvs(variant, infor["mt_classes"])
                else:
                    # More for indels
                    infor["n_indels"] += 1
                    infor["indels_length"].append(abs(len(variant.REF) - len(variant.ALT[0])))
            # Gathering information of variants in comparison to hard mapped regions

            if variant.CHROM in hard_contigs:
                if check_variant_in_bed(
                    variant.CHROM, variant.start, variant.end, hard_file, padding=0
                ):
                    infor["v_in_hard_regions"] += 1
        return infor
    else:
        total_v_count = 0
        for variant in vcf_file:
            total_v_count += 1
        return total_v_count


def main():
    classes = ["T>G", "T>C", "T>A", "C>G", "C>T", "C>A"]
    # reading input
    full_vcf = VCF(args.rawvcf)
    filter_vcf = VCF(args.filtered_vcf)
    sample = filter_vcf.samples[1]
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
            contigs,
            hard_contigs,
            args.minimal,
            args.limited,
            bed_intervals,
            hard_intervals,
            filter_file=True,
            padding=args.padding,
        )
    else:
        infor = process_vcf_file(
            filter_vcf,
            contigs,
            [],
            args.minimal,
            args.limited,
            bed_intervals,
            "",
            filter_file=True,
            padding=args.padding,
        )

    summary = {
        "sample": sample,
        "number_full_variants": full_v_count,
        "number_passed_variants": infor["total_v_count"],
        "number_snvs_inside_exom": infor["v_inside_exom"],
        "number_snvs_outside_exom": infor["v_outside_exom"],
        "number_snvs_in_hard_regions": infor["v_in_hard_regions"],
        "number_snps": infor["n_snps"],
        "number_indels": infor["n_indels"],
        "minimal_rp_snvs_exom": infor["minimal_rp_snvs_exom"],
        "minimal_rp_snvs_nexom": infor["minimal_rp_snvs_nexom"],
        "limited_rp_snvs_exom": infor["limited_rp_snvs_exom"],
        "limited_rp_snvs_nexom": infor["limited_rp_snvs_nexom"],
        "strong_rp_snvs_exom": infor["strong_rp_snvs_exom"],
        "strong_rp_snvs_nexom": infor["strong_rp_snvs_nexom"],
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
        del summary["number_snvs_in_hard_regions"]
    json_object = json.dumps(summary)

    with open(args.output, "w") as outfile:
        outfile.write(json_object)


if __name__ == "__main__":
    main()
