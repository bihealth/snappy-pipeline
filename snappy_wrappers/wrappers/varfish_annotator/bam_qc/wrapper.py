import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

DEF_HELPER_FUNCS = r"""
compute-md5()
{
    if [[ $# -ne 2 ]]; then
        >&2 echo "Invalid number of arguments: $#"
        exit 1
    fi
    md5sum $1 \
    | awk '{ gsub(/.*\//, "", $2); print; }' \
    > $2
}
"""

with tempfile.NamedTemporaryFile("wt") as json_tmpf:
    # Write library name to identifier that should be used in output file into a JSON file.
    # Dictionary expected structure. Key: library name; Value: identifier to be used in file.
    # For snappy-based analysis, both keys and values will be the library name. For analysis using
    # externally generated data, the values will be the sample name as provided by the external
    # source.
    # {
    #    "P001-N1-DNA1-WGS1": "P001",
    #    "P002-N1-DNA1-WGS1": "P002",
    #    "P003-N1-DNA1-WGS1": "P003",
    #  }
    for library_name, identifier in snakemake.params.args.items():
        print(f"{library_name} {identifier}", file=json_tmpf)
        json_tmpf.flush()

    # Actually run the script.
    shell(
        r"""
    set -x

    # Write files for reproducibility -------------------------------------------------------------

    {DEF_HELPER_FUNCS}

    # Write out information about conda and save a copy of the wrapper with picked variables
    # as well as the environment.yaml file.
    conda list >{snakemake.log.conda_list}
    conda info >{snakemake.log.conda_info}
    compute-md5 {snakemake.log.conda_list} {snakemake.log.conda_list_md5}
    compute-md5 {snakemake.log.conda_info} {snakemake.log.conda_info_md5}
    cp {__real_file__} {snakemake.log.wrapper}
    compute-md5 {snakemake.log.wrapper} {snakemake.log.wrapper_md5}
    cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
    compute-md5 {snakemake.log.env_yaml} {snakemake.log.env_yaml_md5}

    # Also pipe stderr to log file ----------------------------------------------------------------

    if [[ -n "{snakemake.log.log}" ]]; then
        if [[ "$(set +e; tty; set -e)" != "" ]]; then
            rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
            exec 2> >(tee -a "{snakemake.log.log}" >&2)
        else
            rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
            echo "No tty, logging disabled" >"{snakemake.log.log}"
        fi
    fi

    # Create auto-cleaned temporary directory
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    # Run actual tools ----------------------------------------------------------------------------

    # --------------------- #
    # Get output identifier #
    # --------------------- #
    get_output_identifier() {{
        # Parameters:
        # $1: file being processed
        # $2: path to temporary file with identifiers hash
        library_name=$(basename $1 | cut -d . -f 2)
        sample=$(grep $library_name $2 | cut -f2 -d' ')
        if [[ -z "$sample" ]]; then
            echo "Could not find '$1' in $2"
            cat $2
            exit 1
        fi
        echo $sample
        return 0
    }}

    # Extract "samtools stats" outputs into JSON
    for bamstats in {snakemake.input.bamstats}; do
        # Define identifier to be used in output
        sample=$(get_output_identifier $bamstats {json_tmpf.name})

        grep '^SN' $bamstats \
        | cut -f 2-3 \
        | sed -e 's/:\t/\t/g' \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ print "{{ \"" sample "\": {{ \"bamstats\": {{"; }}
            {{ if (NR > 1) print ","; print "\"" $1 "\":" $2 }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
        > $TMPDIR/bamstats.$sample
    done

    # Extract "samtools idxstats" into JSON
    for idxstats in {snakemake.input.idxstats}; do
        # Define identifier to be used in output
        sample=$(get_output_identifier $idxstats {json_tmpf.name})

        cat $idxstats \
        | cut -f 1,3,4 \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ print "{{ \"" sample "\": {{ \"idxstats\": {{"; }}
            {{ if (NR > 1) print ","; print "\"" $1 "\":{{\"mapped\":" $2 ",\"unmapped\":" $3 "}}" }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
        > $TMPDIR/idxstats.$sample
    done

    # Extract coverage statistics into JSON.
    for cov_qc in {snakemake.input.cov_qc}; do
        # Define identifier to be used in output
        sample=$(get_output_identifier $cov_qc {json_tmpf.name})

        grep '^SN' $cov_qc \
        | cut -f 2-4 \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ nout=1; print "{{ \"" sample "\": {{ \"summary\": {{" }}
            ($1 == "mean coverage") {{ if (nout > 1) print ","; print "\"" $1 "\":" $2; nout += 1; }}
            ($1 == "number of BED intervals") {{ if (nout > 1) print ","; print "\"target count\":" $2; nout += 1; }}
            ($1 == "total target size") {{ if (nout > 1) print ","; print "\"total target size\":" $2; nout += 1; }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
        > $TMPDIR/cov_qc.summary.$sample

        grep '^MBC' $cov_qc \
        | cut -f 2-4 \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ print "{{ \"" sample "\": {{ \"min_cov_base\": {{"; }}
            {{ if (NR > 1) print ","; print "\"" $1 "\":" $2 }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
        > $TMPDIR/cov_qc.mbc.$sample

        grep '^AMC' $cov_qc \
        | cut -f 2-4 \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ print "{{ \"" sample "\": {{ \"min_cov_target\": {{"; }}
            {{ if (NR > 1) print ","; print "\"" $1 "\":" $2 }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
        > $TMPDIR/cov_qc.amc.$sample
    done

    # Merge JSON files
    out_bam_qc={snakemake.output.bam_qc}

    echo -en "case_id\tset_id\tbam_stats\n.\t.\t" \
    >${{out_bam_qc%.gz}}

    jq -c --slurp 'reduce .[] as $item ({{}}; . * $item)' $TMPDIR/* \
    | sed -e 's/"/\"""/g' \
    >> ${{out_bam_qc%.gz}}

    # Compress output files and create MD5 sums.
    gzip ${{out_bam_qc%.gz}}

    # Compute MD5 sums on output files
    compute-md5 {snakemake.output.bam_qc} {snakemake.output.bam_qc_md5}

    # Create output links --------------------------------------------------------------------------

    for path in {snakemake.output.output_links}; do
      dst=$path
      src=work/${{dst#output/}}
      ln -sr $src $dst
    done
    """
    )

# Compute MD5 sums of logs.
shell(
    r"""
{DEF_HELPER_FUNCS}

sleep 1s  # try to wait for log file flush
compute-md5 {snakemake.log.log} {snakemake.log.log_md5}
"""
)
