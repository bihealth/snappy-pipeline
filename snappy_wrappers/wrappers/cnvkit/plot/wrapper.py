# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py plot
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

mapper = snakemake.wildcards.mapper
if "library_name" in snakemake.wildcards.keys():
    libname = snakemake.wildcards.library_name
    tool_and_mapper = mapper + ".cnvkit.plot"
    make_per_chr_plots = "true"
    opt_filetype = "pdf"
elif "cancer_library" in snakemake.wildcards.keys():
    libname = snakemake.wildcards.cancer_library
    tool_and_mapper = mapper + ".control_freec"
    make_per_chr_plots = "false"
    opt_filetype = "png"
else:
    raise Exception("Unsupported naming")

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

# -----------------------------------------------------------------------------

unset DISPLAY

d="work/{tool_and_mapper}.{libname}/out"

fn="{tool_and_mapper}.{libname}.diagram.pdf"
cnvkit.py diagram -s {snakemake.input.cns} {snakemake.input.cnr} -o $d/$fn
pushd $d ; md5sum $fn > $fn.md5 ; popd

fn="{tool_and_mapper}.{libname}.scatter.{opt_filetype}"
cnvkit.py scatter -s {snakemake.input.cns} {snakemake.input.cnr} -o $d/$fn
pushd $d ; md5sum $fn > $fn.md5 ; popd

fn="{tool_and_mapper}.{libname}.heatmap.{opt_filetype}"
cnvkit.py heatmap {snakemake.input.cnr} -o $d/$fn
pushd $d ; md5sum $fn > $fn.md5 ; popd

if [ "{make_per_chr_plots}" = true ] ; then
    chroms=$(tail -n +2 "{snakemake.input.cnr}" | cut -f 1 | sort | uniq | grep -E "^(chr)?([0-9]+|X|Y)")

    for chr in $chroms ; do
        c=$(echo "$chr" | sed -e "s/^chr//")
        fn="{tool_and_mapper}.{libname}.scatter.chr${{c}}.{opt_filetype}"
        cnvkit.py scatter -s {snakemake.input.cns} --chromosome $chr {snakemake.input.cnr} -o $d/$fn
        pushd $d ; md5sum $fn > $fn.md5 ; popd
    done

    for chr in $chroms ; do
        c=$(echo "$chr" | sed -e "s/^chr//")
        fn="{tool_and_mapper}.{libname}.heatmap.chr${{c}}.{opt_filetype}"
        cnvkit.py heatmap --chromosome $chr {snakemake.input.cnr} -o $d/$fn
        pushd $d ; md5sum $fn > $fn.md5 ; popd
    done
else
    echo "Not creating per-chromosome heatmap or scatter plots."
fi
"""
)
