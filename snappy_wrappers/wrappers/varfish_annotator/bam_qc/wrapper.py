# -*- coding: utf-8 -*-
"""Wrapper for running collecting coverage QC data."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

this_file = __file__

# cov_qc bamstats flagstats idxstats

shell(
    r"""
set -x

# TODO: remove this again, is for fail early
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} {snakemake.log.wrapper}

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

# Extract "samtools stats" outputs into JSON
for bamstats in {snakemake.input.bamstats}; do
    sample=$(basename $bamstats | cut -d . -f 2)

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
    sample=$(basename $idxstats | cut -d . -f 2)

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
    sample=$(basename $cov_qc | cut -d . -f 2)

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
