import math
import numbers
import re
import subprocess

from snakemake.shell import shell
import tabix
import vcfpy

NAME_PATTERN = re.compile("^seg_([0-9]+)_cn_([0-9]+)$")

LFC = {"ID": "LFC", "Number": 1, "Type": "Float", "Description": "Log2 fold change from coverage"}
CN = {
    "ID": "CN",
    "Number": 1,
    "Type": "Integer",
    "Description": "Number of allele called by the CNV caller",
}

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

bed = tabix.open(snakemake.input.cnv)

ps = []
nProc = 0
ps += [
    subprocess.Popen(
        ["bcftools", "merge", snakemake.input.normal, snakemake.input.tumor], stdout=subprocess.PIPE
    )
]
nProc += 1
ps += [
    subprocess.Popen(
        ["bcftools", "filter", "--include", "N_ALT=2 & FORMAT/AD[:2]=0"],
        stdin=ps[nProc - 1].stdout,
        stdout=subprocess.PIPE,
        text=True,
    )
]
ps[nProc - 1].stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
nProc += 1
if config["excluded_regions"]:
    ps += [
        subprocess.Popen(
            ["bcftools", "view", "--targets-file", "^" + config["excluded_regions"]],
            stdin=ps[nProc - 1].stdout,
            stdout=subprocess.PIPE,
            text=True,
        )
    ]
    ps[nProc - 1].stdout.close()
    nProc += 1
reader = vcfpy.Reader.from_stream(ps[nProc - 1].stdout)

sample = snakemake.wildcards["library_name"]
iSample = reader.header.samples.name_to_idx[sample]

reader.header.add_format_line(LFC)
reader.header.add_format_line(CN)

writer = vcfpy.Writer.from_path(snakemake.output.vcf, reader.header)

segments = {}
for variant in reader:
    logFoldChange = None
    cn = None

    # Add default (missing) values for all samples
    variant.add_format("LFC", logFoldChange)
    variant.add_format("CN", cn)

    try:
        locii = list(bed.query(variant.CHROM, variant.affected_start, variant.affected_end + 1))
    except tabix.TabixError as e:
        print(
            "Variant {}:{}{}>{} not covered by segmentation of sample {}".format(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value, sample
            )
        )
        continue

    if len(locii) == 0:
        print(
            "Variant {}:{}{}>{} not covered by segmentation of sample {}".format(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value, sample
            )
        )
    elif len(locii) == 1:
        for locus in locii:
            logFoldChange = float(locus[4])
            if locus[5] == "-":
                logFoldChange = -logFoldChange
            m = NAME_PATTERN.match(locus[3])
            assert m
            cn = int(m.groups()[1])
    else:
        print(
            "Variant {}:{}{}>{} crosses segments in sample {}".format(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value, sample
            )
        )
        for locus in locii:
            m = NAME_PATTERN.match(locus[3])
            assert m
            if cn is None:
                cn = int(m.groups()[1])
            else:
                if cn != int(m.groups()[1]):
                    cn = None
                    break

    variant.calls[iSample].data["LFC"] = logFoldChange
    variant.calls[iSample].data["CN"] = cn

    writer.write_record(variant)

    if len(locii) == 1:
        locus = locii[0]
        if locus[3] not in segments:
            segments[locus[3]] = {
                "logFoldChange": logFoldChange,
                "CHROM": locus[0],
                "start": locus[1],
                "stop": locus[2],
                "BAFs": [],
            }
        depth = variant.call_for_sample[sample].data["AD"]
        baf = min(depth[0], depth[1]) / (depth[0] + depth[1])
        segments[locus[3]]["BAFs"] += [baf]

writer.close()

# Implmentation of R's quantile function for continuous quantiles.
def quantile(x, probs, na_rm=False, method=7):
    if na_rm:
        x = list(
            filter(
                lambda y: y is not None
                and isinstance(y, Number)
                and not math.isnan(y)
                and not math.isinf(y),
                x,
            )
        )
    x.sort()
    n = len(x)
    if method == 4:
        m = [0 for p in probs]
    elif method == 5:
        m = [0.5 for p in probs]
    elif method == 6:
        m = [p for p in probs]
    elif method == 7:
        m = [1 - p for p in probs]
    elif method == 8:
        m = [(p + 1) / 3 for p in probs]
    elif method == 9:
        m = [p / 4 + 3 / 8 for p in probs]
    else:
        raise NotImplementedError("Only continuous methods 4 to 9 (R) are implmented")
    qs = [None] * len(probs)
    for i in range(len(probs)):
        j = math.floor(n * probs[i] + m[i])
        if j < 1:
            qs[i] = x[j]
        elif j >= n:
            qs[i] = x[n - 1]
        else:
            g = n * probs[i] + m[i] - j
            qs[i] = (1 - g) * x[j - 1] + g * x[j]
    return qs


with open(snakemake.output.tsv, "wt") as f:
    print(
        "Number\tCHROM\tstart\tstop\tCN\tLFC\tNvariants\tMin\t1st quartile\tMedian\t3rd quartile\tMax",
        file=f,
    )
    for segId, segment in segments.items():
        m = NAME_PATTERN.match(segId)
        n = int(m.groups()[0])
        cn = int(m.groups()[1])

        logFoldChange = segment["logFoldChange"]
        nVar = len(segment["BAFs"])
        bafs = sorted(segment["BAFs"])

        qs = quantile(segment["BAFs"], (0.25, 0.5, 0.75), method=8)
        values = [
            n,
            segment["CHROM"],
            segment["start"],
            segment["stop"],
            cn,
            segment["logFoldChange"],
            nVar,
            bafs[0],
            qs[0],
            qs[1],
            qs[2],
            bafs[nVar - 1],
        ]
        print("\t".join(map(str, values)), file=f)

shell(
    r"""
tabix {snakemake.output.vcf}

md5() {{
    filename=$1
    fn=$(basename $filename)
    pushd $(dirname $filename) 1> /dev/null 2>&1
    rslt=$(md5sum $fn)
    popd 1> /dev/null 2>&1
    echo "$rslt"
}}

md5 {snakemake.output.vcf} > {snakemake.output.vcf_md5}
md5 {snakemake.output.tbi} > {snakemake.output.tbi_md5}
md5 {snakemake.output.tsv} > {snakemake.output.tsv_md5}
"""
)
