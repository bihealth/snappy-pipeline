# -*- coding: utf-8 -*-
"""Wrapper for running scalpel in somatic mode
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

this_file = __file__

shell(
    r"""
# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# TODO: Fix link problems of tabix. What the link?
export LD_LIBRARY_PATH=$(dirname $(which tabix))/../lib

export TMPDIR=$(mktemp -d)
#trap "rm -rf $TMPDIR" EXIT KILL TERM INT HUP

# Arguments:
# - output VCf
# - only results in target
# - parallel processing with 16 cores
# - improve sensitivity in low-coverage regions
# - improve results using two-pass
# - avoid oversampling in low-coverage regions
# - minimum reportable frequency set to 10 here
# - minimal coverage of 3
# - reference sequence and temporary output directory
# - BED file with regions to call indels for
# - PAIRED analysis only
scalpel-discovery \
    --format vcf \
    --intarget \
    --numprocs 16 \
    --covthr 3 \
    --lowcov 1 \
    --two-pass \
    --seed 1234567 \
    --pathlimit 10000 \
    --outratio 0.1 \
    --mincov 3 \
    --ref {snakemake.config[static_data_config][reference][path]} \
    --dir $TMPDIR/scalpel.tmp \
    --bed \
        {snakemake.config[step_config][somatic_variant_calling][scalpel][path_target_regions]} \
    --somatic \
    --normal {snakemake.input.normal_bam} \
    --tumor {snakemake.input.tumor_bam}

# Infer directory name to create output from
if [[ -f $TMPDIR/scalpel.tmp/twopass/common.indel.vcf ]]; then
    dirname=twopass
else
    dirname=main
fi

# Obtain fixed contig header lines
awk '{{ printf("##contig=<ID=%s,length=%d>\n", $1, $2); }}' \
    {snakemake.config[static_data_config][reference][path]}.fai \
> $TMPDIR/contig_headers.txt

# join and transform output file for tumor/normal pairs
bgzip $TMPDIR/scalpel.tmp/$dirname/common.indel.vcf
tabix $TMPDIR/scalpel.tmp/$dirname/common.indel.vcf.gz
bgzip $TMPDIR/scalpel.tmp/$dirname/somatic.indel.vcf
tabix $TMPDIR/scalpel.tmp/$dirname/somatic.indel.vcf.gz

bcftools concat \
    -a \
    $TMPDIR/scalpel.tmp/$dirname/common.indel.vcf.gz \
    $TMPDIR/scalpel.tmp/$dirname/somatic.indel.vcf.gz \
| awk -v fname=$TMPDIR/contig_headers.txt \
        'BEGIN {{ OFS="\t"; first=1; }}
        /^##contig/ {{
            if (first) {{
                file=fname;
                while ((getline < file) > 0) {{ print }}
                first=0;
            }}
        }}
        !/^##contig/ {{ print $0; }}' \
| sed 's/sample_name/{snakemake.wildcards.tumor_library}\t{snakemake.params.normal_lib_name}/' \
| awk 'BEGIN {{ OFS="\t" }}
        /^#/  {{ print $0 }}
        !/^#/ {{
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, "."
        }}' \
| bgzip -c \
> {snakemake.output.full_vcf}
tabix {snakemake.output.full_vcf}

# split out somatic variants
bcftools view \
    --include 'FILTER=="PASS" & INFO/INH=="no" & INFO/SOMATIC==1' \
    {snakemake.output.full_vcf} \
| bgzip -c \
> {snakemake.output.vcf}
tabix {snakemake.output.vcf}

tar -czf {snakemake.output.tar} $TMPDIR/scalpel.tmp

for f in {snakemake.output}; do
    if [[ -f $f ]]; then
        md5sum $f >$f.md5
    fi
done

# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper.py

# Logging: Save a permanent copy of the environment file used
cp $(dirname {this_file})/environment.yaml $(dirname {snakemake.log.log})/environment_wrapper.yaml
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
