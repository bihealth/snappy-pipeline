# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py
"""

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

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# TODO: add through shell.prefix
export TMPDIR=/fast/users/$USER/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf ${{TMPDIR}}" EXIT
mkdir -p ${{TMPDIR}}/out
mkdir -p ${{TMPDIR}}/vcfs

out_base=${{TMPDIR}}/out/$(basename {snakemake.output.vcf} .vcf.gz)

vcfs=$(echo "{snakemake.input.vcf}" | tr ' ' '\n')

# Create a file with the list of contigs & vcf list for GenomicsDBImport
rm -f ${{TMPDIR}}/contigs.txt
cmd=""
for vcf in ${{vcfs}}
do
    bcftools view -h ${{vcf}} \
        | grep "^##contig=<" \
        | sed -re "s/.*ID=([^,>]+),length=([^,>]+).*/\1:1-\2/" \
        >> ${{TMPDIR}}/contigs_all.list
    cmd="$cmd -V ${{vcf}} "

    rm -f ${{vcf}}.tbi ${{vcf}}.tbi.md5
    gatk IndexFeatureFile -I ${{vcf}}
done
sort ${{TMPDIR}}/contigs_all.list | uniq > ${{TMPDIR}}/contigs.list

# Create the genomicsdb
rm -rf pon_db
gatk --java-options '-Xms10000m -Xmx20000m' GenomicsDBImport \
    --tmp-dir ${{TMPDIR}} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --genomicsdb-workspace-path pon_db \
    --intervals ${{TMPDIR}}/contigs.list \
    $cmd

# Create the panel of normals vcf
gatk CreateSomaticPanelOfNormals \
    --tmp-dir ${{TMPDIR}} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --germline-resource "{snakemake.config[step_config][panel_of_normals][mutect2][germline_resource]}" \
    --variant gendb://pon_db \
    --output ${{out_base}}.vcf

bgzip ${{out_base}}.vcf
tabix -f ${{out_base}}.vcf.gz

pushd $TMPDIR && \
    for f in ${{out_base}}.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv ${{out_base}}.* $(dirname {snakemake.output.vcf})
"""
)
