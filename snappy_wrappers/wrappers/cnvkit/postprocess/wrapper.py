# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

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

set -x

# -----------------------------------------------------------------------------

cat << _EOF | R --vanilla --slave
segments <- read.table("{snakemake.input.segment}", sep="\t", header=1, stringsAsFactors=FALSE)
stopifnot(all(c("chromosome", "start", "end") %in% colnames(segments)))
calls <- read.table("{snakemake.input.call}", sep="\t", header=1, stringsAsFactors=FALSE)
stopifnot(all(c("chromosome", "start", "end", "cn") %in% colnames(calls)))

r1 <- GenomicRanges::GRanges(
    seqnames=segments[["chromosome"]],
    ranges=IRanges::IRanges(start=segments[["start"]], end=segments[["end"]]),
    strand="*"
)
r2 <- GenomicRanges::GRanges(
    seqnames=calls[["chromosome"]],
    ranges=IRanges::IRanges(start=calls[["start"]], end=calls[["end"]]),
    strand="*"
)

i <- GenomicRanges::findOverlaps(r1, r2, ignore.strand=TRUE)
stopifnot(!any(duplicated(S4Vectors::queryHits(i))))

# Default value: diploid segment (number of copies = 2)
segments[,"cn"] <- 2
segments[S4Vectors::queryHits(i),"cn"] <- calls[S4Vectors::subjectHits(i),"cn"]

write.table(segments, file="{snakemake.output.final}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
_EOF

d=$(dirname "{snakemake.output.final}")
pushd $d
fn=$(basename "{snakemake.output.final}")
md5sum $fn > $fn.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
