import pathlib
import argparse


def main(args):
    try:
        from snakemake.shell import shell
    except ImportError:
        import subprocess

        def shell(cmd):
            return subprocess.run(
                cmd.split(" "), shell=True, check=True, text=True, capture_output=True
            )

    # We will first call the GATK tool.  However, there are issues with the reliability of the contig
    # ploidy calls.  This causes an exception in the  JointGermlineCNVSegmentation further downstream.
    # We thus rewrite the genotype calls based on the ones from the PED files as a workaround.
    #
    # cf. https://github.com/broadinstitute/gatk/issues/8164

    out_path = pathlib.Path(args.output_done).parent

    MALE = "male"
    FEMALE = "female"

    sex_map = {}
    for ped_path in args.pedigrees:
        with open(ped_path, "rt") as inputf:
            for line in inputf:
                arr = line.strip().split("\t")
                sex_map[arr[1]] = {"1": MALE, "2": FEMALE}[arr[4]]

    ploidy_x = {MALE: 1, FEMALE: 2}
    ploidy_y = {MALE: 1, FEMALE: 0}

    paths_tsv = " ".join(args.tsvs)

    # Add interval block list for PAR regions if configured.
    par_intervals = args.par_intervals
    if par_intervals:
        par_args = f"-XL {par_intervals}"
    else:
        par_args = ""

    threads = args.threads
    interval_list = args.interval_list
    output_done = args.output_done

    shell(
        rf"""
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" ERR EXIT

    export MKL_NUM_THREADS={threads}
    export OMP_NUM_THREADS={threads}
    export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"

    # Get contig name style
    egrep "^@SQ\s+SN:chr[0-9XY]{{1,2}}\s+" {interval_list} > /dev/null && true
    exit_status=$?
    if [[ $exit_status == 0 ]]; then
        STYLE="chr"
    else
        STYLE=""
    fi


    PRIORS=$TMPDIR/ploidy_priors.tsv
    echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" \
    > $PRIORS
    for i in {{1..22}}; do
        echo -e "$STYLE$i\t0\t0.01\t0.98\t0.01" >> $PRIORS
    done
    echo -e "${{STYLE}}X\t0.01\t0.49\t0.49\t0.01" >> $PRIORS
    echo -e "${{STYLE}}Y\t0.495\t0.495\t0.01\t0" >> $PRIORS

    set -x

    gatk DetermineGermlineContigPloidy \
        -L {interval_list} \
        {par_args} \
        --interval-merging-rule OVERLAPPING_ONLY \
        $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
        --contig-ploidy-priors $PRIORS \
        --output $(dirname {output_done}) \
        --output-prefix ploidy
    """
    )

    ploidy_calls = out_path / "ploidy-calls"
    for sample_dir in ploidy_calls.glob("SAMPLE_*"):
        path_name = sample_dir / "sample_name.txt"
        with path_name.open("rt") as inputf:
            sample_name = inputf.read().strip()
        sample_sex = sex_map[sample_name]
        path_call = sample_dir / "contig_ploidy.tsv"
        with path_call.open("rt") as inputf:
            call_lines = inputf.readlines()
        for i in range(len(call_lines)):
            line = call_lines[i]
            arr = line.split("\t")
            chrom = arr[0].lower()
            if "x" in chrom or "y" in chrom:
                ploidy = int(arr[1])
                if "x" in arr[0].lower():
                    if ploidy != ploidy_x[sample_sex]:
                        arr[1] = str(ploidy_x[sample_sex])
                        arr[2] = f"42.000{ploidy}"  # magic marker
                elif "y" in arr[0].lower():
                    if ploidy != ploidy_y[sample_sex]:
                        arr[1] = str(ploidy_y[sample_sex])
                        arr[2] = f"42.000{ploidy}"  # magic marker
            call_lines[i] = "\t".join(arr)
        with path_call.open("wt") as outputf:
            for line in call_lines:
                print(line.strip(), file=outputf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine Germline Contig Ploidy")
    parser.add_argument(
        "--pedigree",
        dest="pedigrees",
        required=True,
        help="Path(s) to pedigree file(s)",
        action="append",
    )
    parser.add_argument(
        "--tsv", dest="tsvs", required=True, action="append", help="Path(s) to TSV file(s)"
    )
    parser.add_argument(
        "--interval-list", required=True, help="Path to interval list"
    )
    parser.add_argument(
        "--par-intervals", required=False, help="Path to PAR intervals BED file"
    )
    parser.add_argument(
        "--threads", required=True, type=int, help="Number of threads", default=1
    )
    parser.add_argument(
        "--output-done", required=True, help="Path to .done output marker"
    )
    args = parser.parse_known_args()
    main(args)
