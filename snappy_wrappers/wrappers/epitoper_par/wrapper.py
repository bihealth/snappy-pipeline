import json
import os
import sys
import tempfile
import textwrap
from math import ceil

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ResourceUsage,
    SgeResourceUsageConverter,
    gib,
    hours,
    in_working_dir,
)

# TODO: call on overlapping windows, on merge make unique

# Naming clash limbo...
snake_job = snakemake
del snakemake
from snakemake import shell  # noqa: C0411, E402
from snakemake import snakemake as run_snakemake  # noqa: C0411, E402

snakemake = snake_job


# Create Temp Work Directory ----------------------------------------------------------------------

tmpdir = tempfile.mkdtemp("snake_par")

# Perform Splitting -------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -ex;")

# Figure out the number of chunks
chunk_size = snakemake.config["step_config"]["epitope_prediction"]["epitoper"]["chunk_size"]
record_count = int(shell('zgrep -v "^#" {snakemake.input} | wc -l', read=True))
num_chunks = int(ceil(record_count / chunk_size))

# Split input into chunks of the configured size
with in_working_dir(tmpdir, print_chdir=True):
    shell(
        r"""
    # Hack: get back bin directory of base/root environment.
    export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

    mkdir -p splitting

    SUFFIX=.vcf \
    NUM_LINES={snakemake.config[step_config][epitope_prediction][epitoper][chunk_size]} \
        snappy-vcf_split \
            splitting/ \
            {snakemake.input}

    i=0
    for fname in $(ls splitting/*.vcf | sort); do
        mkdir -p job_in.$i.d
        mv $fname job_in.$i.d/input.vcf
        bgzip -f job_in.$i.d/input.vcf
        let "++i"
    done
    """
    )

# Generate Snakefile chunk-by-chunk ---------------------------------------------------------------
result_files = ["job_out.{jobno}.d/.done".format(jobno=jobno) for jobno in range(num_chunks)]
chunks = [
    textwrap.dedent(
        r"""
    shell.executable("/bin/bash")
    shell.prefix("set -ex;")

    configfile: 'config.json'

    localrules: all

    rule all:
        input: {all_results}
    """
    )
    .lstrip()
    .format(all_results=", ".join(map(repr, result_files)))
]


def esc(s):
    return s.replace("{", "{{").replace("}", "}}")


#: Extensions to generate
key_ext = {"txt_gz": "txt.gz", "txt_gz_md5": "txt.gz.md5"}
resources = ResourceUsage(cores=8, memory=gib(30.0), duration=hours(6.0))
for jobno in range(num_chunks):
    output = {
        key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
        for key, ext in key_ext.items()
    }
    chunks.append(
        textwrap.dedent(
            r"""
    rule chunk_{jobno}:
        input:
            vcf='job_in.{jobno}.d/input.vcf.gz',
        output:
            touch("job_out.{jobno}.d/.done"),
            **{output}
        params:
            args={args}
        wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/epitoper'


    cluster_config['chunk_{jobno}'] = {resources}
    """
        ).format(
            jobno=jobno,
            output=repr(output),
            args=repr(snakemake.params["args"]),
            wrapper_prefix="file://" + base_dir,
            resources=repr(SgeResourceUsageConverter(resources).to_res_dict()),
        )
    )

with in_working_dir(tmpdir, print_chdir=True):
    with open("Snakefile", "wt") as snakefile:
        print("\n\n".join(chunks), file=sys.stderr)
        print("\n\n".join(chunks), file=snakefile)

# Write out config file
with in_working_dir(tmpdir, print_chdir=True):
    with open("config.json", "wt") as configfile:
        json.dump(snakemake.config, configfile)
    with open("config.json", "rt") as configfile:
        print(configfile.read(), file=sys.stderr)


with in_working_dir(tmpdir, print_chdir=True):
    # Launch execution on Snakefile
    # snakemake.config['step_config']['epitope_prediction']['epitoper']['use_drmaa']
    run_snakemake(
        "Snakefile",
        use_conda=True,
        cores=snakemake.config["step_config"]["epitope_prediction"]["epitoper"]["num_threads"],
        restart_times=snakemake.config["step_config"]["epitope_prediction"]["epitoper"][
            "restart_times"
        ],
        max_jobs_per_second=snakemake.config["step_config"]["epitope_prediction"]["epitoper"][
            "max_jobs_per_second"
        ],
        max_status_checks_per_second=snakemake.config["step_config"]["epitope_prediction"][
            "epitoper"
        ]["max_status_checks_per_second"],
    )

# Join results

max_jobno = num_chunks - 1
with in_working_dir(tmpdir, print_chdir=True):
    shell(
        textwrap.dedent(
            r"""
    set -euo pipefail

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

    out_base=$(dirname {snakemake.output.txt_gz})/$(basename {snakemake.output.txt_gz} .txt.gz)
    mkdir -p $(dirname {snakemake.output.txt_gz})

    # take first header -------------------------------------------------------
    zcat job_out.0.d/out/tmp_0.txt.gz | head -n 2 >$out_base.txt || true

    # append body contents ----------------------------------------------------
    for jobno in {{0..{max_jobno}}}; do
        zcat job_out.$jobno.d/out/tmp_$jobno.txt.gz | tail -n +3 \
        >>$out_base.txt
    done

    # bgzip output ------------------------------------------------------------
    bgzip -f $out_base.txt
    pushd $(dirname $out_base) && md5sum $(basename $out_base.txt).gz >$out_base.txt.gz.md5 && popd
    """
        ).lstrip()
    )
