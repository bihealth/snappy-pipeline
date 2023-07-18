# -*- coding: utf-8 -*-
"""Wrapper for running Minimap2."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

this_file = __file__

seq_platform = snakemake.params.args["extra_infos"]["seqPlatform"]
library_kit = snakemake.params.args["extra_infos"]["libraryKit"]

shell(
    r"""
set -x

# Write out information about conda and save a copy of the wrapper with picked variables ----------
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Create auto-cleaned temporary directory ---------------------------------------------------------
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Determine preset --------------------------------------------------------------------------------

if [[ "{library_kit}" == "PacBio HiFi" ]]; then
    platform=PACBIO
    preset=map-hifi
elif [[ "{library_kit}" == "PacBio CLR" ]]; then
    platform=PACBIO
    preset=map-pb
elif [[ "{library_kit}" == guppy/* ]] || [[ "{library_kit}" == dorado/* ]]; then
    platform=ONT
    preset=map-ont
else
    >&2 echo "Unknown library kit {library_kit}"
    exit 1
fi

# Function Definitions ----------------------------------------------------------------------------

# Under the assumption that we have single-ended reads, write out lines as FASTQ.
input_files()
{{
    for fname in $(find $(dirname {snakemake.input}) -name '*.bam' -or -name '*.fast?.gz'); do
        basename=$(basename $fname .bam)

        if [[ "$fname" == *.bam ]]; then
            samtools fastq -F 2048 $fname;
        else \
            zcat $fname;
        fi
    done
}}

# Actually run Minimap2.
run_minimap2()
{{
    minimap2 \
        -t {snakemake.config[step_config][ngs_mapping][minimap2][mapping_threads]} \
        -x $preset \
        -a {snakemake.config[step_config][ngs_mapping][minimap2][path_index]} \
        -Y \
        --MD \
        /dev/stdin \
    | samtools addreplacerg \
        -r "@RG\tID:{snakemake.wildcards.library_name}.$i\tSM:{snakemake.wildcards.library_name}\tPL:$PLATFORM" - \
}}

# Alignment postprocessing
#
# This function makes extensive use of tee-piping to prevent the BAM to be read from the disk
# again for statistics.
postproc_bam()
{{
    set -x

    samtools sort \
        -T $TMPDIR/sort_bam \
        -m {snakemake.config[step_config][ngs_mapping][bwa_mem2][memory_bam_sort]} \
        -@ {snakemake.config[step_config][ngs_mapping][bwa_mem2][num_threads_bam_sort]} \
        -O BAM \
        -o /dev/stdout \
        /dev/stdin \
    | tee >(samtools stats    /dev/stdin >{snakemake.output.report_bamstats_txt}) \
    | tee >(samtools flagstat /dev/stdin >{snakemake.output.report_flagstats_txt}) \
    | tee >(stdbuf -o0 samtools index /dev/stdin -o /dev/stdout >{snakemake.output.bam_bai}) \
    | tee >(compute-md5       /dev/stdin {snakemake.output.bam_md5}) \
    > {snakemake.output.bam}


    touch \
        {snakemake.output.bam} \
        {snakemake.output.bam_md5} \
        {snakemake.output.bam_bai} \
        {snakemake.output.bam_bai_md5}

    # We use stdbuf on samtools index above toether with stdout redirection.  We have seen
    # some interesting timing issues with clustered/distributed file systems and the small
    # index file not having been written by the time that samtools idxstats uses it.  In
    # addition we attempt to run it 5 times with incremental waits to make things more
    # robust.
    for i in {{1..5}}; do
        set +e
        samtools idxstats {snakemake.output.bam} >{snakemake.output.report_idxstats_txt}
        ret=$?
        set -e
        if [[ "$ret" -eq 0 ]]; then
            break
        elif [[ "$i" -lt 5 ]]; then
            sleep "${{i}}m"
        fi
    done
    if [[ "$ret" -ne 0 ]]; then
        >&2 echo "samtools idxstats failed 5 times - something wrong with the index?"
        exit 1
    fi
    compute-md5 {snakemake.output.bam_bai} {snakemake.output.bam_bai_md5}

    compute-md5 {snakemake.output.report_bamstats_txt}  {snakemake.output.report_bamstats_txt_md5}
    compute-md5 {snakemake.output.report_flagstats_txt} {snakemake.output.report_flagstats_txt_md5}
    compute-md5 {snakemake.output.report_idxstats_txt}  {snakemake.output.report_idxstats_txt_md5}
}}

# Create output links -----------------------------------------------------------------------------

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
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
