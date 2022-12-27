# -*- coding: utf-8 -*-
"""Wrapper for running Jannovar variant annotation
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

# Get shortcuts to static data and step configuration
static_config = snakemake.config["static_data_config"]
anno_config = snakemake.config["step_config"]["somatic_variant_annotation"]

# Build list of arguments to pass to Jannovar
annotation_args = []

# Add configuration for dbNSFP if configured
if static_config.get("dbnsfp", {}).get("path") and anno_config["dbnsfp"]["columns"]:
    annotation_args += [
        "--dbnsfp-tsv {}".format(repr(static_config["dbnsfp"]["path"])),
        "--dbnsfp-col-contig {}".format(repr(anno_config["dbnsfp"]["col_contig"])),
        "--dbnsfp-col-pos {}".format(repr(anno_config["dbnsfp"]["col_pos"])),
        "--dbnsfp-columns {}".format(repr(",".join(anno_config["dbnsfp"]["columns"]))),
    ]

# Add configuration for variant databases
if static_config.get("dbsnp", {}).get("path"):
    annotation_args.append("--dbsnp-vcf {}".format(static_config["dbsnp"]["path"]))
if static_config.get("exac", {}).get("path"):
    annotation_args.append("--exac-vcf {}".format(static_config["exac"]["path"]))
if static_config.get("gnomad_exomes", {}).get("path"):
    annotation_args.append("--gnomad-exomes-vcf {}".format(static_config["gnomad_exomes"]["path"]))
if static_config.get("gnomad_genomes", {}).get("path"):
    annotation_args.append(
        "--gnomad-genomes-vcf {}".format(static_config["gnomad_genomes"]["path"])
    )
if static_config.get("uk10k", {}).get("path"):
    annotation_args.append("--uk10k-vcf {}".format(static_config["uk10k"]["path"]))
if static_config.get("clinvar", {}).get("path"):
    annotation_args.append("--clinvar-vcf {}".format(static_config["clinvar"]["path"]))
if static_config.get("cosmic", {}).get("path"):
    annotation_args.append("--cosmic-vcf {}".format(static_config["cosmic"]["path"]))

# Add off target flag
if anno_config.get("flag_off_target"):
    annotation_args.append("--enable-off-target-filter")

# Add configuration for BED-based annotation
for conf in anno_config["annotation_tracks_bed"]:
    col_token = "" if not conf.get("label_conf") else ":{}".format(conf.get("label_conf"))
    annotation_args.append(
        "--bed-annotation "
        + repr("{path}:{info_field}:{description}{col_token}".format(col_token=col_token, **conf))
    )

# Add configuration for TSV-based annotation
for conf in anno_config["annotation_tracks_tsv"]:
    offset = 1 if conf["one_based"] else 0
    # pathToTsvFile:oneBasedOffset:colContig:colStart:colEnd:colRef(or=0):colAlt(or=0):colValue:fieldType:fieldName:fieldDescription: accumulationStrategy
    annotation_args.append(
        (
            "--tsv-annotation "
            + repr(
                "{path}:{offset}:{col_contig}:{col_start}:{col_end}:{col_ref}:"
                "{col_alt}:{annotated_alleles}:{col_value}:{field_type}:"
                "{field_name}:{description}:{accumulation_strategy}"
            ).format(offset=offset, **conf)
        )
    )

# Add configuration for VCF-based annotation
for conf in anno_config["annotation_tracks_vcf"]:
    annotation_args.append(
        "--vcf-annotation "
        + repr(
            "{path}:{prefix}:{joined_fields}".format(joined_fields=",".join(conf["fields"]), **conf)
        )
    )

# Build the snippet to put into Jannovar
annotation_snippet = " \\\n    ".join(annotation_args)

# Build intervals argument
arg_intervals = " ".join(
    ["--interval {}".format(interval) for interval in snakemake.params["args"]["intervals"]]
)

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

# See the following for the memory-related massaging
#
# http://bugs.java.com/view_bug.do?bug_id=8043516

# Call jannovar statistics
MALLOC_ARENA_MAX=4 \
jannovar \
    annotate-vcf \
    -XX:MaxHeapSize=5g \
    -XX:CompressedClassSpaceSize=512m \
    -XX:+UseG1GC \
    --input-vcf {snakemake.input.vcf} \
    --output-vcf {snakemake.output.vcf} \
    --database {snakemake.config[step_config][somatic_variant_annotation][path_jannovar_ser]} \
    --ref-fasta {snakemake.config[static_data_config][reference][path]} \
    {arg_intervals} \
    {annotation_snippet}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi}).md5
"""
)
