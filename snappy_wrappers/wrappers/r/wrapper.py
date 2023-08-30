"""CUBI+Snakemake wrapper code for non-conda package installation
"""

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

if "packages" in snakemake.params.keys():
    packages = snakemake.params["packages"]
elif "packages" in config.keys():
    packages = config["packages"]
else:
    packages = None

assert packages is not None
assert isinstance(packages, list)
assert len(packages) > 0

to_install = {"github": [], "bitbucket": [], "local": [], "bioconductor": [], "cran": []}
for package in packages:
    package = dict(package)
    if "repo" in package.keys() and package["repo"] in to_install.keys():
        to_install[package["repo"]].append(package["name"])
    else:
        to_install["cran"].append(package["name"])
to_install = {k: "c({})".format(", ".join(['"{}"'.format(vv) for vv in v])) for k, v in to_install.items()}

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

R --vanilla --slave << __EOF
if (length({to_install[cran]}) > 0) {{
    install.packages({to_install[cran]}, lib=dirname("{snakemake.output.done}"), update=FALSE, ask=FALSE)
}}
if (length({to_install[bioconductor]}) > 0) {{
    BiocManager::install({to_install[bioconductor]}, lib=dirname("{snakemake.output.done}"), update=FALSE, ask=FALSE)
}}
if (length({to_install[github]}) > 0) {{
    remotes::install_github({to_install[github]}, lib=dirname("{snakemake.output.done}"), upgrade="never")
}}
if (length({to_install[bitbucket]}) > 0) {{
    remotes::install_bitbucket({to_install[bitbucket]}, lib=dirname("{snakemake.output.done}"), upgrade="never")
}}
if (length({to_install[local]}) > 0) {{
    remotes::install_local({to_install[local]}, lib=dirname("{snakemake.output.done}"), upgrade="never")
}}
__EOF
touch {snakemake.output.done}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
