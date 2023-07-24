import os
import re
import subprocess

import tabix
import vcfpy

# Can't add anything if there's no CNV results
if "cnv" not in snakemake.input:
    os.symlink(snakemake.input.vcf, snakemake.output.vcf)
    os.symlink(snakemake.input.vcf_tbi, snakemake.output.vcf_tbi)
    return

NAME_PATTERN = re.compile("^seg_([0-9]+)_cn_([0-9]+)$")

LFC = {"ID": "LFC", "Number": 1, "Type": "Float", "Description": "Log2 fold change from coverage"}
CN = {"ID": "CN", "Number": 1, "Type": "Integer", "Description": "Number of allele called by the CNV caller"}

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

bed = tabix.open(snakemake.input.cnv)

ps = []
nProc = 0
ps += [subprocess.Popen(["bcftools", "merge", snakemake.input.normal, snakemake.input.tumor], stdout=subprocess.PIPE)]
nProc += 1
ps += [subprocess.Popen(["bcftools", "filter", "--include", "N_ALT=3 & FORMAT/AD[:2]=0"], stdin=ps[nProc-1].stdout, stdout=subprocess.PIPE, text=True)]
ps[nProc-1].stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
nProc += 1
if config["excluded_regions"]:
    ps += [subprocess.Popen(["bcftools", "view", "--regions-file", "^" + config["excluded_regions"], stdin=ps[nProc-1].stdout, stdout=subprocess.PIPE, text=True)]
    ps[nProc-1].stdout.close()
    nProc += 1
reader = vcfpy.Reader.from_stream(ps[nProc-1].stdout)

sample = snakemake.params["args"]["library_name"]
iSample = reader.header.samples.name_to_idx[sample]

reader.header.add_format_line(LFC)
reader.header.add_format_line(CN)

writer = vcfpy.Writer.from_path(snakemake.output.vcf, reader.header)

for variant in reader:
    logFoldChange = None
    cn = None

    # Add default (missing) values for all samples
    variant.add_format("LFC", logFoldChange)
    variant.add_format("CN", cn)

    try:
        locii = list(bed.query(variant.CHROM, variant.affected_start, variant.affected_end+1))
    except tabix.TabixError as e:
        print("Variant {}:{}{}>{} not covered by segmentation of sample {}".format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value, sample))
        continue

    if len(locii) == 0:
        print("Variant {}:{}{}>{} not covered by segmentation of sample {}".format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value, sample))
    elif len(locii) == 1:
        for locus in locii:
            logFoldChange = float(locus[4])
            if locus[5] == "-":
                logFoldChange = -logFoldChange
            m = NAME_PATTERN.match(locus[3])
            assert m
            cn = int(m.groups()[1])
    else:
        print("Variant {}:{}{}>{} crosses segments in sample {}".format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value, sample))
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

writer.close()
