from snakemake import shell

shell(
    r"""
set -x

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

picard -Xmx6g -Djava.io.tmpdir=$TMPDIR \
    CollectHsMetrics \
    I={snakemake.input.bam} \
    O={snakemake.output.txt} \
    R={snakemake.config[static_data_config][reference][path]} \
    MINIMUM_MAPPING_QUALITY=0 \
    TARGET_INTERVALS={snakemake.config[step_config][ngs_mapping][picard_hs_metrics][path_targets_interval_list]} \
    BAIT_INTERVALS={snakemake.config[step_config][ngs_mapping][picard_hs_metrics][path_baits_interval_list]} \
    VALIDATION_STRINGENCY=SILENT

pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
popd
"""
)
