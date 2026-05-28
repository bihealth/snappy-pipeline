# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

# WARNING- In rare cases, the parallel invokation of the prepare_panel
#     wrapper creates multiple instances of variants located in padding
#     regions. Most of them are removed when the distinct regions are
#     merged, but it may happen that the variant in padding region is
#     called with support with different ALT. In this case, the variant
#     is not merged, the same locus appears multiple times in the vcf
#     produced by the prepare_panel wrapper, and the creation of the
#     panel of normals fails.
#     These rare occurences can be fixed by manually removing the
#     multiple variants.

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

if java_options := args.get("java_options", ""):
    java_options = f"--java-options '{java_options}'"

extra_arguments = " ".join(args.get("extra_arguments", []))

ShellWrapper(snakemake).run(
    r"""
set -x

export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

mkdir -p $TMPDIR/out
mkdir -p $TMPDIR/vcfs

out_base=$TMPDIR/out/$(basename {snakemake.output.vcf} .vcf.gz)
mkdir -p $out_base

vcfs=$(echo "{snakemake.input.normals}" | tr ' ' '\n')

# Create a file with the list of contigs & vcf list for GenomicsDBImport
rm -f $TMPDIR/contigs.txt
cmd=""
for vcf in ${{vcfs}}
do
    bcftools view -h ${{vcf}} \
        | grep "^##contig=<" \
        | sed -re "s/.*ID=([^,>]+),length=([^,>]+).*/\1:1-\2/" \
        >> $TMPDIR/contigs_all.list
    cmd="$cmd -V ${{vcf}} "
done
sort $TMPDIR/contigs_all.list | uniq > $TMPDIR/contigs.list

# Create the genomicsdb
rm -rf $TMPDIR/pon_db
gatk {java_options} GenomicsDBImport \
    --tmp-dir $TMPDIR \
    --reference {snakemake.input.reference} \
    --genomicsdb-workspace-path $TMPDIR/pon_db \
    --intervals $TMPDIR/contigs.list \
    {extra_arguments} \
    $cmd

# Create the panel of normals vcf
gatk CreateSomaticPanelOfNormals \
    --tmp-dir $TMPDIR \
    --reference {snakemake.input.reference} \
    --germline-resource "{snakemake.input.germline_resource}" \
    --variant gendb://$TMPDIR/pon_db \
    --output ${{out_base}}.vcf

bgzip ${{out_base}}.vcf
tabix -f ${{out_base}}.vcf.gz

# Make a copy of the genomics database for PureCN
# NOTE: the sleep & true commands are required to work around
#       a tar error triggered by a cephfs bug/feature
#       (https://ceph-users.ceph.narkive.com/th0JxsKR/cephfs-tar-archiving-immediately-after-writing)
#       The bug is probably triggered because GATK genomicsdb is large is size & can contain 100000s files
sleep 10
tar -zcvf {snakemake.output.db} -C $TMPDIR pon_db || true

# Copy the results to destination & compute checksums
cp ${{out_base}}.vcf.gz {snakemake.output.vcf}
cp ${{out_base}}.vcf.gz.tbi {snakemake.output.vcf}.tbi

pushd $(dirname {snakemake.output.vcf})
f=$(basename {snakemake.output.vcf})
md5sum $f > $f.md5
md5sum $f.tbi > $f.tbi.md5
popd

pushd $(dirname {snakemake.output.db})
f=$(basename {snakemake.output.db})
md5sum $f > $f.md5
popd
"""
)

