# -*- coding: utf-8 -*-
"""Wrapper for running collecting coverage QC data."""
import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

this_file = __file__

# cov_qc bamstats flagstats idxstats

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
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    for library_name, identifier in snakemake.params.args.items():
        print(f"{library_name} {identifier}", file=json_tmpf)
        json_tmpf.flush()

    # Actually run the script.
    shell(
        r"""
    set -x

    # TODO: remove this again, is for fail early
    # Additional logging for transparency & reproducibility
    # Logging: Save a copy this wrapper (with the pickle details in the header)
    cp {this_file} {snakemake.log.wrapper}

    # Write out information about conda installation.
    conda list > {snakemake.log.conda_list}
    conda info > {snakemake.log.conda_info}
    md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
    md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}

    # Also pipe stderr to log file
    if [[ -n "{snakemake.log}" ]]; then
        if [[ "$(set +e; tty; set -e)" != "" ]]; then
            rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
            exec 2> >(tee -a "{snakemake.log}" >&2)
        else
            rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
            echo "No tty, logging disabled" >"{snakemake.log}"
        fi
    fi

    # Create auto-cleaned temporary directory
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

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

    pushd $(dirname {snakemake.output.bam_qc}) && \
        md5sum $(basename {snakemake.output.bam_qc}) > $(basename {snakemake.output.bam_qc}).md5
    """
    )
