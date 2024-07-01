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

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf ${{tmpdir}}" EXIT
mkdir -p ${{tmpdir}}/out
mkdir -p ${{tmpdir}}/vcfs

out_base=${{tmpdir}}/out/$(basename {snakemake.output.vcf} .vcf.gz)
mkdir -p $out_base

vcfs=$(echo "{snakemake.input.normals}" | tr ' ' '\n')

# Create a file with the list of contigs & vcf list for GenomicsDBImport
rm -f ${{tmpdir}}/contigs.txt
cmd=""
for vcf in ${{vcfs}}
do
    bcftools view -h ${{vcf}} \
        | grep "^##contig=<" \
        | sed -re "s/.*ID=([^,>]+),length=([^,>]+).*/\1:1-\2/" \
        >> ${{tmpdir}}/contigs_all.list
    cmd="$cmd -V ${{vcf}} "
done
sort ${{tmpdir}}/contigs_all.list | uniq > ${{tmpdir}}/contigs.list

# Create the genomicsdb
rm -rf ${{tmpdir}}/pon_db
gatk --java-options '-Xms10000m -Xmx20000m' GenomicsDBImport \
    --tmp-dir ${{tmpdir}} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --genomicsdb-workspace-path ${{tmpdir}}/pon_db \
    --intervals ${{tmpdir}}/contigs.list \
    $cmd

# Create the panel of normals vcf
gatk CreateSomaticPanelOfNormals \
    --tmp-dir ${{tmpdir}} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --germline-resource "{snakemake.config[step_config][panel_of_normals][mutect2][germline_resource]}" \
    --variant gendb://${{tmpdir}}/pon_db \
    --output ${{out_base}}.vcf

bgzip ${{out_base}}.vcf
tabix -f ${{out_base}}.vcf.gz

# Make a copy of the genomics database for PureCN
# NOTE: the sleep & true commands are required to work around 
#       a tar error triggered by a cephfs bug/feature
#       (https://ceph-users.ceph.narkive.com/th0JxsKR/cephfs-tar-archiving-immediately-after-writing)
#       The bug is probably triggered because GATK genomicsdb is large is size & can contain 100000s files 
sleep 10
tar -zcvf {snakemake.output.db} -C ${{tmpdir}} pon_db || true

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

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
