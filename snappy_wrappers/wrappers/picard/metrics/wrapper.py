# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for picard metrics collection: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = args["reference"]

collect_multiple_metrics_programs = {
    "CollectAlignmentSummaryMetrics",
    "CollectBaseDistributionByCycle",
    "CollectGcBiasMetrics",
    "CollectInsertSizeMetrics",
    "CollectQualityYieldMetrics",
    "CollectSequencingArtifactMetrics",
    "MeanQualityByCycle",
    "QualityScoreDistribution",
}
collect_multiple_metrics = " ".join(
    [
        f"-PROGRAM {pgm}"
        for pgm in collect_multiple_metrics_programs.intersection(set(args["programs"]))
    ]
)

prefix = args.get("prefix", "")

name = args.get("bait_name", "null")

ShellWrapper(snakemake).run(
    r"""
set -x

d=$(ls $CONDA_PREFIX/share | grep picard)
picard_jar="$CONDA_PREFIX/share/$d/picard.jar"
if [[ ! -r $picard_jar ]]
then
    echo "Can't find picard jar"
    exit -1
fi

d=$(dirname {snakemake.output[0]})

if [[ -n "{collect_multiple_metrics}" ]]
then
    # Setting METRIC_ACCUMULATION_LEVEL causes issues for some programs
    java -jar $picard_jar CollectMultipleMetrics \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectMultipleMetrics \
        -R {reference} \
        -FILE_EXTENSION .txt \
        -PROGRAM null \
        {collect_multiple_metrics}
fi

if [[ "{args[programs]}" == *"EstimateLibraryComplexity"* ]]
then
    java -jar $picard_jar EstimateLibraryComplexity \
        -I {snakemake.input.bam} \
        -O $d/{prefix}EstimateLibraryComplexity.txt
fi

if [[ "{args[programs]}" == *"CollectJumpingLibraryMetrics"* ]]
then
    java -jar $picard_jar CollectJumpingLibraryMetrics \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectJumpingLibraryMetrics.txt
fi

if [[ "{args[programs]}" == *"CollectOxoGMetrics"* ]]
then
    if [[ -n "{args[dbsnp]}" ]]
    then
        dbsnp="-DB_SNP {args[dbsnp]}"
    else
        dbsnp=""
    fi
    java -jar $picard_jar CollectOxoGMetrics \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectOxoGMetrics.txt \
        -R {reference} \
        $dbsnp
fi

if [[ "{args[programs]}" == *"CollectHsMetrics"* ]]
then
    java -jar $picard_jar CollectHsMetrics \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectHsMetrics.txt \
        -R {reference} \
        -BAIT_SET_NAME {name} \
        -BAIT_INTERVALS {snakemake.input.baits} \
        -TARGET_INTERVALS {snakemake.input.targets}
fi

if [[ "{args[programs]}" == *"CollectRawWgsMetrics"* ]]
then
    java -jar $picard_jar CollectRawWgsMetrics \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectRawWgsMetrics.txt \
        -R {reference}
fi

if [[ "{args[programs]}" == *"CollectWgsMetrics"* ]]
then
    java -jar $picard_jar CollectWgsMetrics \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectWgsMetrics.txt \
        -R {reference}
fi

if [[ "{args[programs]}" == *"CollectWgsMetricsWithNonZeroCoverage"* ]]
then
    java -jar $picard_jar CollectWgsMetricsWithNonZeroCoverage \
        -I {snakemake.input.bam} \
        -O $d/{prefix}CollectWgsMetricsWithNonZeroCoverage.txt \
        -CHART $d/{prefix}CollectWgsMetricsWithNonZeroCoverage.pdf \
        -R {reference}
fi
"""
)

