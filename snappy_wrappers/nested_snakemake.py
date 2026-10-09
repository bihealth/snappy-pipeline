# -*- coding: utf-8 -*-
"""Run a nested Snakemake workflow from within a wrapper, e.g. for the mbcs meta-wrapper."""

import os
import shlex
from pathlib import Path

from snakemake.api import ResourceSettings, SnakemakeApi


class SnakemakeExecutionFailed(Exception):
    """Raised when nested snakemake execution failed"""


def run_snakemake(
    config,
    snakefile="Snakefile",
    cores=1,
    num_jobs=0,
    max_jobs_per_second=0,
    max_status_checks_per_second=0,
    job_name_token="",
    partition=None,
    profile=None,
):
    """Given a pipeline step's configuration, launch sequential or parallel Snakemake"""
    snakefile = Path(snakefile)
    if config["use_profile"]:
        workdir = Path(os.getcwd())
        print(
            f"Running with Snakemake profile on {num_jobs or config['num_jobs']} "
            f"cores in directory {workdir}"
        )
        os.mkdir(os.path.join(workdir, "slurm_log"))
        if partition:
            os.environ["SNAPPY_PIPELINE_DEFAULT_PARTITION"] = partition

        # Write Snakemake file: debug helper
        write_snakemake_debug_helper(
            profile=profile,
            jobs=str(num_jobs or config["num_jobs"]),
            restart_times=str(config["restart_times"]),
            job_name_token=job_name_token,
            max_jobs_per_second=str(max_jobs_per_second or config["max_jobs_per_second"]),
            max_status_checks_per_second=str(
                max_status_checks_per_second or config["max_status_checks_per_second"]
            ),
        )

        with SnakemakeApi() as api:
            result = (
                api.workflow(
                    snakefile=snakefile,
                    workdir=workdir,
                    resource_settings=ResourceSettings(
                        cores=cores, nodes=num_jobs or config["num_jobs"]
                    ),
                    # TODO properly choose remaining *_settings, if needed
                    # config_settings=None,
                    # storage_settings=None,
                    # workflow_settings=None,
                    # deployment_settings=None,
                    # storage_provider_settings=None,
                )
                .dag()
                .execute_workflow()
            )
    else:
        print(
            "Running locally with {num_jobs} jobs in directory {cwd}".format(
                num_jobs=config["num_jobs"], cwd=os.getcwd()
            )
        )
        with SnakemakeApi() as api:
            result = (
                api.workflow(
                    snakefile=snakefile,
                    resource_settings=ResourceSettings(cores=config["num_jobs"]),
                    # TODO properly choose remaining *_settings, if needed
                    # config_settings=None,
                    # storage_settings=None,
                    # workflow_settings=None,
                    # deployment_settings=None,
                    # storage_provider_settings=None,
                )
                .dag()
                .execute_workflow()
            )
    if result is False:
        raise SnakemakeExecutionFailed("Could not perform nested Snakemake call")


def write_snakemake_debug_helper(
    profile, jobs, restart_times, job_name_token, max_jobs_per_second, max_status_checks_per_second
):
    """Write Snakemake debug helper file

    When the temporary directory is kept, a failed execution can be restarted by calling snakemake
    in the temporary directory with the command line written to the file ``snakemake_call.sh``.

    :param profile: Snakemake profile name.
    :type profile: str

    :param jobs: Number of jobs argument,  ``--jobs``.
    :type jobs: str

    :param restart_times: Number of restarts argument, ``--restart-times``.
    :type restart_times: str

    :param job_name_token: Token included in job name, ``--jobname``.
    :type job_name_token: str

    :param max_jobs_per_second: Max number of jobs per second argument, ``--max-jobs-per-second``.
    :type max_jobs_per_second: str

    :param max_status_checks_per_second: Max status checks per second argument,
    ``--max-status-checks-per-second``.
    :type max_status_checks_per_second: str
    """
    with open(os.path.join(os.getcwd(), "snakemake_call.sh"), "wt") as f_call:
        print("/bin/bash")
        print("#SBATCH --output {}/slurm_log/%x-%J.log".format(os.getcwd()))
        print(
            " ".join(
                map(
                    str,
                    [
                        "snakemake",
                        "--directory",
                        os.getcwd(),
                        "--cores",
                        "--printshellcmds",
                        "--verbose",
                        "--software-deployment-method",  # sic! <- ?
                        "conda",
                        "--profile",
                        shlex.quote(profile),
                        "--jobs",
                        jobs,
                        "--restart-times",
                        restart_times,
                        "--jobname",
                        shlex.quote(
                            "snakejob{token}.{{rulename}}.{{jobid}}.sh".format(
                                token="." + job_name_token
                            )
                        ),
                        "--max-jobs-per-second",
                        max_jobs_per_second,
                        "--max-status-checks-per-second",
                        max_status_checks_per_second,
                    ],
                )
            ),
            file=f_call,
        )
