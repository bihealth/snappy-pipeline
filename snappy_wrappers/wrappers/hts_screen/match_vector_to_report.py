#!/usr/bin/env python3
"""Convert match_vector.tsv file from snappy-hts_screen.sh"""

import argparse
import sys

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

KEYS = ["unmapped", "one_one", "one_multi", "multi_one", "multi_multi"]


def init_counts(genomes):
    return [dict((key, 0) for key in KEYS) for g in genomes]


def run(args):
    counts = []
    num_reads = 0
    num_reads_no_hits = 0
    for lineno, line in enumerate(args.in_file):
        line = line.strip()
        fields = line.split("\t")
        if lineno == 0:  # initialize counts and genomes
            genomes = fields[1:-1]
            counts = init_counts(genomes)
            continue
        num_reads += 1
        # compute values for determining match class
        arr = list(map(int, fields[1:-1]))
        if not sum(arr):
            num_reads_no_hits += 1
        is_multi, n_genomes = False, 0
        for i, n_alis in enumerate(arr):
            if n_alis >= 1:
                n_genomes += 1
            if n_alis > 1:
                is_multi = True
        # increment counts
        for i, n_alis in enumerate(arr):
            if n_alis == 0:
                counts[i]["unmapped"] += 1
            elif not is_multi:
                assert n_genomes >= 1
                if genomes == 1:
                    counts[i]["one_one"] += 1
                else:
                    counts[i]["one_multi"] += 1
            else:
                if genomes == 1:
                    counts[i]["multi_one"] += 1
                else:
                    counts[i]["multi_multi"] += 1
    # print results
    print("#hts_screen version 0.1\t#Reads in subset: ${num_reads}", file=args.out_file)
    print(
        (
            "Genome\t#Reads_processed\t#Unmapped\t%Unmapped\t"
            "#One_hit_one_genome\t%One_hit_one_genome\t"
            "#Multiple_hits_one_genome\t%Multiple_hits_one_genome\t"
            "#One_hit_multiple_genomes\t%One_hit_multiple_genomes\t"
            "Multiple_hits_multiple_genomes\t"
            "%Multiple_hits_multiple_genomes"
        ),
        file=args.out_file,
    )
    for i, genome in enumerate(genomes):
        tpl = "{}\t{}\t{}\t{:0.2f}\t{}\t{:0.2f}\t{}\t{:0.2f}\t" "{}\t{:0.2f}\t{}\t{:0.2f}"
        print(
            tpl.format(
                genome,
                num_reads,
                counts[i]["unmapped"],
                100 * counts[i]["unmapped"] / num_reads,
                counts[i]["one_one"],
                100 * counts[i]["one_one"] / num_reads,
                counts[i]["multi_one"],
                100 * counts[i]["multi_one"] / num_reads,
                counts[i]["one_multi"],
                100 * counts[i]["one_multi"] / num_reads,
                counts[i]["multi_multi"],
                100 * counts[i]["multi_multi"] / num_reads,
            ),
            file=args.out_file,
        )

    perc_no_hits = 100 * num_reads_no_hits / num_reads
    print(file=args.out_file)
    print("%Hit_no_genomes {:.02f}".format(perc_no_hits), file=args.out_file)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument(
        "--in-file",
        default=sys.stdin,
        type=argparse.FileType("rt"),
        help="Input file, defaults to stdin",
    )
    parser.add_argument(
        "--out-file",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output file, defaults to stdout",
    )
    run(parser.parse_args())


def __main__():
    return main()


if __name__ == "__main__":
    main(sys.argv)
