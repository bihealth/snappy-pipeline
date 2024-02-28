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

to_install = []
for package in packages:
    package = dict(package)
    to_install.append('list(name="{}", repo="{}")'.format(package["name"], package["repo"]))
to_install = "list({})".format(", ".join(to_install))

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
for (pkg in {to_install}) {{
    if (pkg[["repo"]] == "cran") install.packages(pkg[["name"]], lib=dirname("{snakemake.output.done}"), update=FALSE, ask=FALSE)
    if (pkg[["repo"]] == "bioconductor") BiocManager::install(pkg[["name"]], lib=dirname("{snakemake.output.done}"), update=FALSE, ask=FALSE)
    if (pkg[["repo"]] == "github") remotes::install_github(pkg[["name"]], lib=dirname("{snakemake.output.done}"), upgrade="never")
    if (pkg[["repo"]] == "bitbucket") remotes::install_bitbucket(pkg[["name"]], lib=dirname("{snakemake.output.done}"), upgrade="never")
    if (pkg[["repo"]] == "local") remotes::install_local(pkg[["name"]], lib=dirname("{snakemake.output.done}"), upgrade="never")
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
