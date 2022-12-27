# -*- coding: utf-8 -*-
"""Code for running a parallel CUBI+Snakemake wrappers

Parallel running is performed in a fork/join parallelism manner.  The user can specify a "split",
a "work", and a "merge" task.
"""

from collections.abc import MutableMapping, MutableSequence
import contextlib
import datetime
import functools
import itertools
import json
import logging
import math
import os
import shlex
import shutil
import sys
import tempfile
import textwrap

from snakemake import get_profile_file, snakemake

from snappy_wrappers.tools.genome_windows import yield_regions

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


@contextlib.contextmanager
def in_working_dir(path, print_chdir=False):
    """Context manager to perform work in a given working directory"""
    cwd = os.getcwd()
    if print_chdir:
        print("cd {}".format(path), file=sys.stderr)
    os.chdir(path)
    yield
    if print_chdir:
        print("cd {}".format(cwd), file=sys.stderr)
    os.chdir(cwd)


def augment_parser(parser):
    """Augment parser from sequential wrapper with parallel wrapper"""

    group = parser.add_argument_group("Parallel Options", "Configuration of parallel execution")
    group.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for local processing/jobs to spawn concurrently",
    )
    group.add_argument(
        "--use-profile",
        default=None,
        help="Enables running the parallelization on cluster with the given snakemake profile",
    )
    group.add_argument(
        "--restart-times", type=int, default=5, help="Number of times to restart jobs automatically"
    )
    group.add_argument(
        "--max-jobs-per-second",
        type=int,
        default=10,
        help="Maximal number of jobs to launch per second in cluster mode",
    )
    group.add_argument(
        "--max-status-checks-per-second",
        type=int,
        default=10,
        help="Maximal number of cluster status checks for perform for second",
    )

    return parser


class WrapperInfoMixin:
    """Mixin for your parallel ``WrapperInfo`` classes"""

    def get_config_chunk(self):
        """Return chunk to pass into step configuration for the parallelism-specific configuration"""
        return {
            "num_jobs": self.args.num_jobs if hasattr(self.args, "num_jobs") else 1,
            "num_thread": self.args.num_threads if hasattr(self.args, "num_threads") else 1,
            "use_profile": self.args.use_profile,
            "restart_times": self.args.restart_times,
            "max_jobs_per_second": self.args.max_jobs_per_second,
            "max_status_checks_per_second": self.args.max_status_checks_per_second,
        }


class JobDescription:
    """Information that is required for executing a shell-based job

    The information can be used for generating a Snakefile for execution.  Note that the
    memory requirement here is given in total and not per core as for SGE.
    """

    def __init__(self, input_files, output_files, resource_usage, script, log_dir=None):
        #: list of input files; if relative, relative to temporary working
        #: directory
        self.input_files = input_files
        #: list of output files; if relative, relative to temporary working
        #: directory
        self.output_files = output_files
        #: ``ResourceUsage`` describing the resource usage
        self.resource_usage = resource_usage
        #: Bash script for executing the job, automatically dedented
        self.script = textwrap.dedent(script)
        #: Path to log directory, if any
        self.log_dir = log_dir


def kib(num):
    """Convert kilobytes to bytes.

    :param num: Number of kilobytes.
    :type num: int

    :return: Returns number of bytes for number of kilobytes.
    """
    return int(num) * 1024


def mib(num):
    """Convert megabytes to bytes.

    :param num: Number of megabytes.
    :type num: int

    :return: Returns number of bytes for number of megabytes.
    """
    return int(num) * 1024 * 1024


def gib(num):
    """Convert gigabytes to bytes.

    :param num: Number of gigabytes.
    :type num: int

    :return: Returns number of bytes for number of gigabytes.
    """
    return int(num) * 1024 * 1024 * 1024


def kib_to_string(num):
    """
    :param num: Number of kilobytes.
    :type num: int

    :return: Returns kilobytes string for memory usage with profile. Example: '200k'.
    """
    return f"{math.ceil(num)}k"


def mib_to_string(num):
    """
    :param num: Number of megabytes.
    :type num: int

    :return: Returns megabytes string for memory usage with profile. Example: '8M'.
    """
    return f"{math.ceil(num)}M"


def gib_to_string(num):
    """
    :param num: Number of gigabytes.
    :type num: int

    :return: Returns gigabytes string for memory usage with profile. Example: '12G'.
    """
    return f"{math.ceil(num)}G"


def minutes(num):
    """Get Timedelta for minutes.

    :param num: Number of minutes.
    :type num: int
    :return: Returns ``datetime.timedelta`` as string for ``num`` minutes.
    """
    return str(datetime.timedelta(minutes=num))


def hours(num):
    """Get Timedelta for hours.

    :param num: Number of hours.
    :type num: int
    :return: Returns ``datetime.timedelta`` as string for ``num`` hours.
    """
    output = str(datetime.timedelta(hours=num))
    # If necessary converts datetime to resource expected format.
    # Examples:
    # '1 day, 1:00:00' -> '1-01:00:00'
    # '2 days, 0:00:00' -> '2-00:00:00'
    if "day" in output:
        output = output.replace(" days, ", "-").replace(" day, ", "-")
        tmp_arr = output.split("-")
        output = tmp_arr[0] + "-" + ":".join([str(i).zfill(2) for i in tmp_arr[1].split(":")])
    return output


def days(num):
    """Get Timedelta for days.

    :param num: Number of days.
    :type num: int
    :return: Returns ``datetime.timedelta`` as string for ``num`` days.
    """
    return str(datetime.timedelta(days=num))


class SnakemakeExecutionFailed(Exception):
    """Raised when nested snakemake execution failed"""


def run_snakemake(
    config,
    snakefile="Snakefile",
    num_jobs=0,
    max_jobs_per_second=0,
    max_status_checks_per_second=0,
    job_name_token="",
    partition=None,
    profile=None,
):
    """Given a pipeline step's configuration, launch sequential or parallel Snakemake"""
    if config["use_profile"]:
        print(
            f"Running with Snakemake profile on {num_jobs or config['num_jobs']} "
            f"cores in directory {os.getcwd()}"
        )
        os.mkdir(os.path.join(os.getcwd(), "slurm_log"))
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

        result = snakemake(
            snakefile,
            workdir=os.getcwd(),
            jobname="snakejob{token}.{{rulename}}.{{jobid}}.sh".format(token="." + job_name_token),
            cores=1,
            nodes=num_jobs or config["num_jobs"],
            max_jobs_per_second=max_jobs_per_second or config["max_jobs_per_second"],
            max_status_checks_per_second=max_status_checks_per_second
            or config["max_status_checks_per_second"],
            restart_times=config["restart_times"],
            verbose=True,
            use_conda=False,  # has to be done externally (no locking if True here) and is!
            jobscript=get_profile_file(profile, "slurm-jobscript.sh"),
            cluster=get_profile_file(profile, "slurm-submit.py"),
            cluster_status=get_profile_file(profile, "slurm-status.py"),
            cluster_sidecar=get_profile_file(profile, "slurm-sidecar.py"),
            cluster_cancel="scancel",
        )
    else:
        print(
            "Running locally with {num_jobs} jobs in directory {cwd}".format(
                num_jobs=config["num_jobs"], cwd=os.getcwd()
            )
        )
        result = snakemake(
            snakefile,
            cores=config["num_jobs"],
            max_jobs_per_second=config["max_jobs_per_second"],
            max_status_checks_per_second=config["max_status_checks_per_second"],
            restart_times=config["restart_times"],
            verbose=True,
            use_conda=False,  # has to be done externally (no locking if True here) and is!
        )
    if not result:
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
                        "--use-conda",  # sic!
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


def to_plain_python(obj):
    """Convert nested structure mixing dict/list objects with dict/lits (such as the output)
    of loading via ``ruamel.yaml`` to using builtin dict/list only.
    """
    if isinstance(obj, (list, MutableSequence)):
        return [to_plain_python(elem) for elem in obj]
    elif isinstance(obj, (dict, MutableMapping)):
        return {key: to_plain_python(elem) for key, elem in obj.items()}
    else:
        return obj


class ParallelBaseWrapper:
    """Base class for parallel wrapper classes.

    Serves mainly as a template class with methods to override to obtain the desired behaviour
    in the sub classes.  Also, it takes care of the dirty details such as creating the temporary
    directory and cleaning up when necessary.
    """

    #: The keys in ``snakemake.output`` that ``os.path.realpath`` should be applied to.
    realpath_output_keys = None
    realpath_log_keys = (
        "log",
        "log_md5",
        "conda_info",
        "conda_info_md5",
        "conda_list",
        "conda_list_md5",
    )
    #: The name of the step, used for selecting configuration settings.
    step_name = None
    #: The name of the tool, used for selecting configuration settings (if ``None`` then
    #: don't go further into config dict).
    tool_name = None
    #: The token to use for job names and temporary directories.  Constructed from ``step_name``
    #: and ``tool_name`` if not specified explicitly here.
    job_name_token = None
    #: The number of bases to pad the parallelization windows with.
    window_padding = 0
    job_resources = None
    merge_resources = None

    def __init__(self, snakemake):
        """Constructor.

        :param snakemake: Reference to ``snakemake`` object from ``wrapper.py``.
        :type snakemake: snakemake.script.Snakemake
        """
        #: Reference to ``snakemake`` object from ``wrapper.py``.
        self.snakemake = snakemake
        #: Base directory to wrappers
        self.wrapper_base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
        # Kick-off initialization
        self._apply_realpath_to_output()
        # Setup logging (will run in wrapper and its own process)
        if hasattr(snakemake.log, "log"):
            log_filename = self.snakemake.log.log
        else:
            log_filename = os.path.realpath(str(snakemake.log))
        logging.basicConfig(
            filename=log_filename,
            format="%(asctime)s %(name)s %(levelname)-7s %(message)s",
            datefmt="%Y-%m-%d %H:%M",
            level=logging.INFO,
        )  # TODO: make configurable?
        self.logger = logging.getLogger(self._job_name_token())
        #: Path to the main pipeline step workflow directory
        self.main_cwd = os.getcwd()

    def _apply_realpath_to_output(self):
        """Update output and log file paths to be realpaths."""
        for key in self.realpath_output_keys or []:
            if hasattr(self.snakemake.output, key):
                setattr(
                    self.snakemake.output,
                    key,
                    os.path.realpath(getattr(self.snakemake.output, key)),
                )
        for key in self.realpath_log_keys or []:
            if hasattr(self.snakemake.log, key):
                setattr(self.snakemake.log, key, os.path.realpath(getattr(self.snakemake.log, key)))

    def get_fai_path(self):
        """Return path to FAI file for reference to use.

        The base implementation uses the value in ``/static_data_config/reference/path`` to
        construct the path by appending ``".fai"``.
        """
        return self.snakemake.config["static_data_config"]["reference"]["path"] + ".fai"

    def get_all_log_files(self):
        """Return dict with realpaths of log files."""
        key_ext = {
            # Explicitly do not pass the log file, as conda list and info get written after merging
            # but the log file contains the overall logs from all steps. It gets passed to
            # logging.basicConfig() above.
            # 'log': 'log',
            "conda_info": "conda_info.txt",
            "conda_info_md5": "conda_info.txt.md5",
            "conda_list": "conda_list.txt",
            "conda_list_md5": "conda_list.txt.md5",
        }
        return {key: getattr(self.snakemake.log, key) for key in key_ext.keys()}

    def _get_config(self):
        if self.tool_name:
            return self._get_step_config()[self.tool_name]
        else:
            return self._get_step_config()

    def _get_step_config(self):
        return self.snakemake.config["step_config"][self.step_name]

    def get_job_mult_memory(self):
        """Get job memory multiplier.

        :return: Returns job memory multiplier, as defined in the provided config. If not defined,
        the default is 1.
        """
        return self._get_config().get("job_mult_memory", 1)

    def get_job_mult_time(self):
        """Get job running time multiplier.

        :return: Returns job running time multiplier, as defined in the provided config.
        If not defined, the default is 1.
        """
        return self._get_config().get("job_mult_time", 1)

    def get_merge_mult_memory(self):
        """Get job running memory for merger.

        :return: Returns job memory multiplier for merger, as defined in the provided config.
        If not defined, the default is 1.
        """
        return self._get_config().get("merge_mult_memory", 1)

    def get_merge_mult_time(self):
        """Get job running time multiplier for merger.

        :return: Returns job running time multiplier for merger, as defined in the provided config.
        If not defined, the default is 1.
        """
        return self._get_config().get("merge_mult_time", 1)

    def get_window_length(self):
        """Get window length.

        :return: Returns window length, as defined in the configuration.
        :raises KeyError: if window_length not defined in the configuration.
        """
        try:
            return self._get_config()["window_length"]
        except KeyError:
            error_msg = (
                "No default available for 'window_length', it must be defined in configuration."
            )
            self.logger.error(error_msg)
            raise

    def get_ignore_chroms(self):
        """Get list of ignored chromosomes.

        :return: Returns list of chromosomes names to be ignored, as defined in the configuration.
        If not defined, the default is an empty list.
        """
        return self._get_step_config().get("ignore_chroms", [])

    @functools.lru_cache(maxsize=16)
    def get_regions(self):
        """Return regions to process.

        This is constructed from the FAI (returned by ``get_fai_path()``), the window length (as
        returned by ``get_window_length()``, and the chromosomes to ignore
        (``get_ignore_chroms()``).
        """
        self.logger.info("Loading regions from FAI at %s", self.get_fai_path())
        with open(self.get_fai_path(), "rt") as fai_file:
            result = list(
                yield_regions(
                    fai_file,
                    self.get_window_length(),
                    ignore_chroms=self.get_ignore_chroms(),
                    padding=self.window_padding,
                )
            )
        # Users can truncate to a number of tokens by specifying in the configuration.
        self.logger.debug("Raw region list is %s", result)
        if self._get_config().get("debug_trunc_tokens", 0):
            self.logger.info(
                "debug_trunc_tokens has been set to %d, will limit number of created files",
                self._get_config()["debug_trunc_tokens"],
            )
            result = result[: self._get_config()["debug_trunc_tokens"]]
        self.logger.info("Region list is %s", result)
        return result

    def construct_parallel_result_files(self):
        """Construct list of parallel result files.

        This is based on the number of regions as returned by ``get_regions()``.
        """
        return ["job_out.{}.d/.done".format(jobno) for jobno in range(len(self.get_regions()))]

    def construct_preamble(self):
        """Return preamble for Snakefile.

        The default preamble configure shell executable to bash, sets prefix to ``"set -ex;"``,
        loads configuration, and starts with local rule ``all``.
        """
        return (
            textwrap.dedent(
                r"""
            shell.executable("/bin/bash")
            shell.prefix("set -ex;")

            configfile: 'config.json'

            localrules: all

            def resource_chunk_threads(wildcards):
                '''Return the number of threads to use for running one chunk.'''
                return {chunk_resources_threads}

            def resource_chunk_memory(wildcards):
                '''Return the memory to use for running one chunk.'''
                return {chunk_resources_memory}

            def resource_chunk_time(wildcards):
                '''Return the time to use for running one chunk.'''
                return {chunk_resources_time}

            def resource_chunk_partition(wildcards):
                '''Return the partition to use for running one chunk.'''
                return {chunk_resources_partition}

            def resource_merge_threads(wildcards):
                '''Return the number of threads to use for running merging.'''
                return {merge_resources_threads}

            def resource_merge_memory(wildcards):
                '''Return the memory to use for running merging.'''
                return {merge_resources_memory}

            def resource_merge_time(wildcards):
                '''Return the time to use for running merging.'''
                return {merge_resources_time}

            def resource_merge_partition(wildcards):
                '''Return the partition to use for running merging.'''
                return {merge_resources_partition}

            rule all:
                input: **{all_output}
        """
            )
            .lstrip()
            .format(
                all_output=repr(self.get_all_output()),
                chunk_resources_threads=repr(self.job_resources.threads),
                chunk_resources_time=repr(self.job_resources.time),
                chunk_resources_memory=repr(self.job_resources.memory),
                chunk_resources_partition=repr(self.job_resources.partition),
                merge_resources_threads=repr(self.merge_resources.threads),
                merge_resources_time=repr(self.merge_resources.time),
                merge_resources_memory=repr(self.merge_resources.memory),
                merge_resources_partition=repr(self.merge_resources.partition),
            )
        )

    def get_all_output(self):
        """Return overall output."""
        raise NotImplementedError("Override me!")

    def construct_parallel_rules(self):
        """Construct parallel rules."""
        raise NotImplementedError("Override me!")

    def construct_merge_rule(self):
        """Construct rule for merging chunks."""
        raise NotImplementedError("Override me!")

    def construct_epilogue(self):
        """Return epilogue for Snakefile, default is to return empty string"""
        return ""

    def joint_chunks(self):
        return "\n\n".join(
            itertools.chain(
                ["# PREAMBLE", self.construct_preamble(), "# PARALLEL WORK"],
                self.construct_parallel_rules(),
                [
                    "# JOIN PARALLEL RESULTS",
                    self.construct_merge_rule(),
                    "# EPILOGUE",
                    self.construct_epilogue(),
                ],
            )
        )

    def _job_name_token(self):
        """Construct job name token or take from class."""
        if self.job_name_token:
            return self.job_name_token
        else:
            if self.tool_name:
                return self.step_name + "_" + self.tool_name
            else:
                return self.step_name

    def run(self):
        # The setup of the temporary directory depends on whether it is to be kept (for debugging
        # purposes) or not.
        keep_tmpdir = self._get_config().get("keep_tmpdir", "never")
        # Either run with TemporaryDirectory as context manager and auto-cleanup on exit or
        # create temporary directory using mkdtemp().
        if keep_tmpdir == "never":
            self.logger.info("Running temporary directory with automated clean-up")
            with tempfile.TemporaryDirectory(self._job_name_token()) as tmpdir:
                self._do_run(tmpdir)
            return self  # short-circuit
        else:
            tmpdir = tempfile.mkdtemp(self._job_name_token())
        # Handle case of always keeping or cleanup on error only.
        if keep_tmpdir == "always":
            self.logger.info("Running temporary directory and WILL NOT CLEAN UP")
            self._do_run(tmpdir)
        else:  # keep_tmpdir == 'onerror'
            self.logger.info("Running in temporary directory, will cleanup in case of success")
            try:
                self._do_run(tmpdir)
            except SnakemakeExecutionFailed as e:
                self.logger.info(
                    "Caught error %s: %s, WILL NOT CLEAN UP temporary directory %s",
                    type(e),
                    e,
                    tmpdir,
                )
                raise  # re-raise e
            else:
                self.logger.info("Ran through successfully, cleaning up %s...", tmpdir)
                shutil.rmtree(tmpdir)
                self.logger.info("Done cleaning up.")
        return self

    def shutdown_logging(self):
        logging.shutdown()

    def _do_run(self, tmpdir):
        """Actual processing in ``tmpdir``.

        Separated from ``run()`` which handles the creation and cleanup of tmpdir.
        """
        with in_working_dir(tmpdir, print_chdir=True):
            self._write_config_file()
            self._write_snakefile()
            self._launch_execution()

    def _write_config_file(self):
        """Write out the configuration file (dump of ``snakemake.config``)."""
        with open("config.json", "wt") as configfile:
            self.logger.info(
                "Writing config.json with content >>>%s<<<",
                json.dumps(to_plain_python(self.snakemake.config), indent="  ", sort_keys=True),
            )
            json.dump(to_plain_python(self.snakemake.config), configfile)

    def _write_snakefile(self):
        with open("Snakefile", "wt") as snakefile:
            self.logger.info("Writing Snakefile with content >>>%s<<<", self.joint_chunks())
            print(self.joint_chunks(), file=snakefile)

    def _launch_execution(self):
        kwargs = {
            "num_jobs": self._get_config()["num_jobs"],
            "max_jobs_per_second": self._get_config()["max_jobs_per_second"],
            "max_status_checks_per_second": self._get_config()["max_status_checks_per_second"],
            "job_name_token": self._job_name_token(),
            "profile": os.getenv("SNAPPY_PIPELINE_SNAKEMAKE_PROFILE"),
        }
        self.logger.info("Launching excecution with args: %s", repr(kwargs))
        run_snakemake(self._get_config(), **kwargs)


class ParallelVcfOutputBaseWrapper(ParallelBaseWrapper):
    """Base class for wrappers generating VCF, splitting genome into windows."""

    #: The maximal number of files to merge in one go.  This must be lower than the maximal number
    #: of files that can be open on the system.  Thus, 1000 is a good default value.  We probably
    #: do not need to have a per-workflow or per-tool tuning algorithm in the step configuration.
    merge_block_size = 1000

    #: Relative path to wrapper to use.
    inner_wrapper = None

    #: Wrappers generating VCF output want to call ``realpath`` on these keys.
    realpath_output_keys = ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")

    #: Extensions to generate
    key_ext = {
        "vcf": "vcf.gz",
        "vcf_md5": "vcf.gz.md5",
        "vcf_tbi": "vcf.gz.tbi",
        "vcf_tbi_md5": "vcf.gz.tbi.md5",
    }

    def get_all_output(self):
        """Return dict with overall output."""
        return {key: getattr(self.snakemake.output, key) for key in self.key_ext.keys()}

    def construct_merge_rule(self):
        """Join the overall result files"""
        # Get list of parallel result ``.done`` files and transform into ``*.vcf.gz``.
        merge_input = [
            os.path.join(os.path.dirname(p), "out", "tmp_{}.vcf.gz".format(i))
            for i, p in enumerate(self.construct_parallel_result_files())
        ]
        # Distinguish cases of more than two levels, two levels, and one levels for merging.
        if len(merge_input) > self.merge_block_size * self.merge_block_size:
            # Too many files (>1M with merge_block_size of 1k)
            raise Exception("Number of output file requires more than two mergin steps!")
        elif len(merge_input) > self.merge_block_size:
            # We need two merge passes.
            merge_outputs = [
                "merge_out.{}.d/out/out.vcf.gz".format(i)
                for i in range(  # ceiled integer division below
                    (len(merge_input) + self.merge_block_size - 1) // self.merge_block_size
                )
            ]
            # Compute breaks (sequence of ``[start, end)`` for splitting merge_input into blocks).
            breaks = list(range(0, len(merge_input), self.merge_block_size)) + [len(merge_input)]
            return "\n\n".join(
                itertools.chain(
                    [
                        self._construct_level_one_merge_rule(chunk_no, merge_input[start:end])
                        for chunk_no, (start, end) in enumerate(zip(breaks, breaks[1:]))
                    ],
                    [self._construct_final_merge_rule(merge_outputs)],
                )
            )
        else:
            # We can do with one merge pass.
            return self._construct_final_merge_rule(merge_input)

    def _construct_level_one_merge_rule(self, chunk_no, merge_input):
        return (
            textwrap.dedent(
                r"""
            rule merge_chunk_{chunk_no}:
                input: {chunk_input}
                output:
                    vcf='merge_out.{chunk_no}.d/out/out.vcf.gz',
                    tbi='merge_out.{chunk_no}.d/out/out.vcf.gz.tbi',
                threads: resource_merge_threads
                resources:
                    time=resource_merge_time,
                    memory=resource_merge_memory,
                    partition=resource_merge_partition,
                shell:
                    r'''
                    set -euo pipefail  # inofficial Bash strict mode

                    bcftools concat \
                        --allow-overlaps \
                        -d none \
                        -o {{output.vcf}} \
                        -O z \
                        {{input}}

                    tabix -f {{output.vcf}}
                    '''
        """
            )
            .lstrip()
            .format(
                chunk_no=chunk_no,
                chunk_input=repr(merge_input),
            )
        )

    def _construct_final_merge_rule(self, merge_input):
        return (
            textwrap.dedent(
                r"""
            rule merge_all:
                input: {all_input}
                output: **{all_output}
                threads: resource_merge_threads
                resources:
                    time=resource_merge_time,
                    memory=resource_merge_memory,
                    partition=resource_merge_partition,
                log: **{all_log}
                shell:
                    r'''
                    set -euo pipefail  # Unofficial Bash strict mode

                    # Initialize output directory -----------------------------------------

                    outdir=$(basename {{output.vcf}})

                    mkdir -p output

                    # Concatenate files ---------------------------------------------------
                    bcftools concat \
                        --allow-overlaps \
                        -d none \
                        -o output/out.vcf.gz \
                        -O z \
                        {{input}}

                    tabix -f output/out.vcf.gz

                    pushd output
                    for f in *; do md5sum $f >$f.md5; done
                    popd

                    # Move to output directory --------------------------------------------
                    mkdir -p $(dirname {{output.vcf}})
                    mv output/out.vcf.gz {{output.vcf}}
                    mv output/out.vcf.gz.md5 {{output.vcf_md5}}
                    mv output/out.vcf.gz.tbi {{output.tbi}}
                    mv output/out.vcf.gz.tbi.md5 {{output.tbi_md5}}

                    # Write out information about conda installation.
                    conda list >{{log.conda_list}}
                    conda info >{{log.conda_info}}

                    pushd $(dirname {{log.conda_list}})
                    md5sum $(basename {{log.conda_list}}) >$(basename {{log.conda_list}}).md5
                    md5sum $(basename {{log.conda_info}}) >$(basename {{log.conda_info}}).md5
                    popd
                    '''
        """
            )
            .lstrip()
            .format(
                all_input=repr(merge_input),
                all_output=repr(self.get_all_output()),
                all_log=repr(self.get_all_log_files()),
            )
        )


class ParallelVariantCallingBaseWrapper(ParallelVcfOutputBaseWrapper):
    """Base class for wrappers performing variant calling, split genome into windows."""

    window_padding = 10000  # >=py36: 10_000

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
            params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "input_bam": repr(self.snakemake.input),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        {input_bam},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_chunk_threads
                    resources:
                        time=resource_chunk_time,
                        memory=resource_chunk_memory,
                        partition=resource_chunk_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'
            """
            ).format(**vals).lstrip()


class ParallelVariantAnnotationBaseWrapper(ParallelVcfOutputBaseWrapper):
    """Base class for wrappers performing variant annotation, split genome into windows."""

    #: Sensible padding when splitting along genome.
    window_padding = 1000  # >=py36: 1_000

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
            params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "input_": repr(
                    {
                        key: os.path.realpath(
                            os.path.join(self.main_cwd, getattr(self.snakemake.input, key))
                        )
                        for key in ("vcf", "tbi", "ped")
                    }
                ),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        **{input_},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_chunk_threads
                    resources:
                        time=resource_chunk_time,
                        memory=resource_chunk_memory,
                        partition=resource_chunk_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'
            """
            ).format(**vals).lstrip()


class ParallelSomaticVariantCallingBaseWrapper(ParallelVcfOutputBaseWrapper):
    """Base class for wrappers performing somatic variant calling, split genome into windows."""

    #: Sensible padding when splitting along genome.
    window_padding = 10000  # >=py36: 10_000

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
            params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
            # The parameters "normal_lib_name" and "tumor_lib_name" are only available for
            # Mutect 2, thus the conditional assignment.
            if hasattr(self.snakemake.params, "normal_lib_name"):
                params["normal_lib_name"] = self.snakemake.params.normal_lib_name
            if hasattr(self.snakemake.wildcards, "tumor_library"):
                params["tumor_lib_name"] = self.snakemake.wildcards.tumor_library
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "normal_bam": repr(self.snakemake.input.normal_bam),
                "tumor_bam": repr(self.snakemake.input.tumor_bam),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        tumor_bam={tumor_bam},
                        normal_bam={normal_bam},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_chunk_threads
                    resources:
                        time=resource_chunk_time,
                        memory=resource_chunk_memory,
                        partition=resource_chunk_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'
            """
            ).format(**vals).lstrip()


class ParallelSomaticVariantAnnotationBaseWrapper(ParallelVcfOutputBaseWrapper):
    """Base class for wrappers performing somatic variant annotation, split genome into windows."""

    #: Sensible padding when splitting along genome.
    window_padding = 1000  # >=py36: 1_000

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
            params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "input_": repr(
                    {
                        key: os.path.realpath(
                            os.path.join(self.main_cwd, getattr(self.snakemake.input, key))
                        )
                        for key in ("vcf", "tbi")
                    }
                ),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        **{input_},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_chunk_threads
                    resources:
                        time=resource_chunk_time,
                        memory=resource_chunk_memory,
                        partition=resource_chunk_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'
            """
            ).format(**vals).lstrip()
