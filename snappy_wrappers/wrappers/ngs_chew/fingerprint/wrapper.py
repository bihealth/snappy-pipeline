# -*- coding: utf-8 -*-
"""Wrapper for running ``ngs-chew fingerprint``."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

path_ref = snakemake.config["static_data_config"]["reference"]["path"]
if "hg19" in path_ref or "37" in path_ref:
    genome_release = "GRCh37"
else:
    genome_release = "GRCh38"

shell(
    r"""
set -x

# Write out information about conda and save a copy of the wrapper with picked variables
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

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/{{out,sorted,sort.tmp}}

if [[ "{library_kit}" == "PacBio HiFi" ]]; then
    preset=map-hifi
elif [[ "{library_kit}" == "PacBio CLR" ]]; then
    preset=map-pb
elif [[ "{library_kit}" == ONT* ]]; then
    preset=map-ont
else
    >&2 echo "Unknown library kit {library_kit}"
    exit 1
fi

ngs-chew fingerprint \
    --min-coverage 5 \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --output-fingerprint {snakemake.ouput.npz} \
    --input-bam {snakemake.input.bam}
    --genome-release {genome_release}


# Compute MD5 sums
pushd $(dirname $out_bam)
md5sum $(basename $out_bam) >$(basename $out_bam).md5
md5sum $(basename $out_bam).bai >$(basename $out_bam).bai.md5
popd

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}

# call plot-bamstats
mkdir $TMPDIR/bamstats.d
plot-bamstats \
    -p $TMPDIR/bamstats.d/ \
    {snakemake.output.report_bamstats_txt} \
|| true  # ignore failure

# Patch inline-html if necessary.
cat >$TMPDIR/inline-html.diff <<EOF
diff --git a/inline_html/inline_html.py b/inline_html/inline_html.py
index 893086c..cbef6dd 100644
--- a/inline_html/inline_html.py
+++ b/inline_html/inline_html.py
@@ -20,7 +20,10 @@ def resource_to_data(path, in_file):
     mime, _ = mimetypes.guess_type(path)
     with open(path, 'rb') as fp:
         data = fp.read()
-        data64 = b''.join(base64.encodestring(data).splitlines())
+        try:
+            data64 = b''.join(base64.encodestring(data).splitlines())
+        except AttributeError:
+            data64 = b''.join(base64.encodebytes(data).splitlines())
         return 'data:%s;base64,%s' % (mime, data64.decode('ascii'))
EOF
pushd $(python3 -c 'import inline_html; print(inline_html.__path__[0])')
if ! grep encodebytes inline_html.py; then
    patch -p2 <$TMPDIR/inline-html.diff
fi
popd

# Convert HTML report into one file.
inline-html \
    --in-file $TMPDIR/bamstats.d/index.html \
    --out-file {snakemake.output.report_bamstats_html} \
|| touch {snakemake.output.report_bamstats_html}

# Build MD5 files for the reports
md5sum {snakemake.output.report_bamstats_html} > {snakemake.output.report_bamstats_html_md5}
md5sum {snakemake.output.report_bamstats_txt} > {snakemake.output.report_bamstats_txt_md5}
md5sum {snakemake.output.report_flagstats_txt} >{snakemake.output.report_flagstats_txt_md5}
md5sum {snakemake.output.report_idxstats_txt} > {snakemake.output.report_idxstats_txt_md5}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
