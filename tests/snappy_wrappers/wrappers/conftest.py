"""Support code for the snappy pipeline wrapper tests."""
from collections import namedtuple
from itertools import chain
import os
import shutil
import subprocess
import textwrap
from unittest.mock import MagicMock

from pyfakefs import fake_filesystem
import pytest
from ruamel import yaml as yaml_ruamel
from snakemake.io import InputFiles, Log, OutputFiles, Params, Resources, Wildcards
from snakemake.script import Snakemake
import yaml

FORCE_RUN = os.environ.get("FORCE_RUN", "false") == "true"
DIFF_MASTER = os.environ.get("DIFF_MASTER", "false") == "true"
DIFF_LAST_COMMIT = os.environ.get("DIFF_LAST_COMMIT", "false") == "true"

if FORCE_RUN or DIFF_MASTER or DIFF_LAST_COMMIT:
    compare = "HEAD^" if DIFF_LAST_COMMIT else "origin/master"

    # check if wrapper is modified compared to master
    DIFF_FILES = set(
        subprocess.check_output(["git", "diff", compare, "--name-only"]).decode().split("\n")
    )

#: Allow running containerized (local only).
CONTAINERIZED = os.environ.get("CONTAINERIZED", "false") == "true"


class Skipped(Exception):
    pass


skip_if_not_modified = pytest.mark.xfail(raises=Skipped)


def run_workflow(wrapper, test_dir, cmd, tmpdir, check_log=None):
    d = str(tmpdir)
    origdir = os.getcwd()
    dst = os.path.join(d, "master")
    os.makedirs(dst, exist_ok=True)

    def copy(pth, src):
        shutil.copy(os.path.join(pth, src), os.path.join(dst, pth))

    used_wrappers = []
    wrapper_file = "used_wrappers.yaml"
    if os.path.exists(os.path.join(wrapper, wrapper_file)):
        # is meta wrapper
        with open(os.path.join(wrapper, wrapper_file), "r") as wf:
            wf = yaml.safe_load(wf)
            used_wrappers = wf["wrappers"]
    else:
        used_wrappers.append(wrapper)

    for w in used_wrappers:
        success = False
        for ext in ("py", "R", "Rmd"):
            script = "wrapper." + ext
            if os.path.exists(os.path.join(w, script)):
                os.makedirs(os.path.join(dst, w), exist_ok=True)
                copy(w, script)
                success = True
                break
        assert success, "No wrapper script found for {}".format(w)
        copy(w, "environment.yaml")

    if (
        not FORCE_RUN
        and (DIFF_MASTER or DIFF_LAST_COMMIT)
        and not any(
            any(f.startswith(w) for f in DIFF_FILES) for w in chain(used_wrappers, [wrapper])
        )
    ):
        raise Skipped("wrappers not modified")

    testdir = os.path.join(d, "test")
    # pkgdir = os.path.join(d, "pkgs")
    shutil.copytree(os.path.join(os.path.dirname(__file__), test_dir), testdir)
    # prepare conda package dir
    # os.makedirs(pkgdir)
    # switch to test directory
    os.chdir(testdir)
    if os.path.exists(".snakemake"):
        shutil.rmtree(".snakemake")
    cmd = cmd + ["--wrapper-prefix", "file://{}/master/".format(d), "--conda-cleanup-pkgs"]

    if CONTAINERIZED:
        # run snakemake in container
        cmd = [
            "sudo",
            "docker",
            "run",
            "-it",
            "-v",
            "{}:{}".format(os.getcwd(), "/workdir/test"),
            "-w",
            "/workdir",
            "snakemake/snakemake",
            " ".join(cmd),
        ]

    # env = dict(os.environ)
    # env["CONDA_PKGS_DIRS"] = pkgdir
    try:
        subprocess.check_call(cmd)
    except Exception as e:
        # go back to original directory
        os.chdir(origdir)
        logfiles = [
            os.path.join(d, f)
            for d, _, files in os.walk(os.path.join(testdir, "logs"))
            for f in files
        ]
        for path in logfiles:
            with open(path) as f:
                msg = "###### Logfile: " + path + " ######"
                print(msg, "\n")
                print(f.read())
                print("#" * len(msg))
        if check_log is not None:
            for f in logfiles:
                check_log(open(f).read())
        else:
            raise e
    finally:
        # cleanup environments to save disk space
        subprocess.check_call(
            r"for env in `conda env list | grep -P '\.snakemake/conda' | "
            r"cut -f1 | tr -d ' '`; do conda env remove --prefix $env; done",
            shell=True,  # nosec
        )
        # go back to original directory
        os.chdir(origdir)


@pytest.fixture
def fai_file_content():
    """Returns FAI file content based on hs37d5 (chromosome 1 only)."""
    return "1\t249250621\t52\t60\t61"


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config_mutect2():
    """Return YAML parsing result for (germline) configuration"""
    return yaml_ruamel.round_trip_load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          cosmic:
            path: /path/to/cosmic.vcf.gz
          dbsnp:
            path: /path/to/dbsnp.vcf.gz
        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa
          somatic_variant_calling:
            tools:
            - mutect2
            mutect2:
              panel_of_normals: ''      # Set path to panel of normals vcf if required
              germline_resource: REQUIRED # Germline variants resource (same as panel of normals)
              common_variants: REQUIRED # Common germline variants for contamination estimation
              # Parallelization configuration
              drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
              num_cores: 2              # number of cores to use locally
              window_length: 50000000   # split input into windows of this size, each triggers a job
              num_jobs: 500             # number of windows to process in parallel
              use_drmaa: true           # use DRMAA for parallel processing
              restart_times: 5          # number of times to re-launch jobs in case of failure
              max_jobs_per_second: 2    # throttling of job creation
              max_status_checks_per_second: 10   # throttling of status checks
              debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
              keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
              job_mult_memory: 2        # memory multiplier
              job_mult_time: 3          # running time multiplier
              merge_mult_memory: 4      # memory multiplier for merging
              merge_mult_time: 5        # running time multiplier for merging
              ignore_chroms:            # patterns of chromosome names to ignore
              - NC_007605    # herpes virus
              - hs37d5       # GRCh37 decoy
              - chrEBV       # Eppstein-Barr Virus
              - '*_decoy'    # decoy contig
              - 'HLA-*'      # HLA genes
              - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def snakemake_output_dict():
    """Returns dictionary that defined snakemake.output"""
    output_base_name = "work/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1"
    return {
        "raw": output_base_name + ".raw.vcf.gz",
        "raw_md5": output_base_name + ".raw.vcf.gz.md5",
        "raw_tbi": output_base_name + ".raw.vcf.gz.tbi",
        "raw_tbi_md5": output_base_name + ".raw.vcf.gz.tbi.md5",
        "stats": output_base_name + ".raw.vcf.stats",
        "stats_md5": output_base_name + ".raw.vcf.stats.md5",
        "f1r2": output_base_name + ".raw.f1r2_tar.tar.gz",
        "f1r2_md5": output_base_name + ".raw.f1r2_tar.tar.gz.md5",
    }


@pytest.fixture
def snakemake_obj(minimal_config_mutect2, snakemake_output_dict):
    """Returns Snakemake object."""
    # Define helper variables
    rule_name = "somatic_variant_calling_mutect2_run"
    threads = 2
    bench_iteration = 2
    scriptdir = "/work"
    input_dict = {
        "tumor_bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "tumor_bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }

    log_base_name = "work/bwa.mutect2.P001-T1-DNA1-WGS1/log/bwa.mutect2.P001-T1-DNA1-WGS1"
    log_dict = {
        "conda_info": log_base_name + ".conda_info.txt",
        "conda_info_md5": log_base_name + ".conda_info.txt.md5",
        "conda_list": log_base_name + ".conda_list.txt",
        "conda_list_md5": log_base_name + ".conda_list.txt.md5",
        "log": log_base_name + ".log",
        "log_md5": log_base_name + ".log.md5",
    }
    wildcards_dict = {"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"}
    params_dict = {"normal_lib_name": "P001-N1-DNA1-WGS1"}

    # Define Snakemake class input
    input_ = InputFiles(fromdict=input_dict)
    output_ = OutputFiles(fromdict=snakemake_output_dict)
    params_ = Params(fromdict=params_dict)
    log_ = Log(fromdict=log_dict)
    wildcards_ = Wildcards(fromdict=wildcards_dict)
    resources_ = Resources(fromdict={})

    return Snakemake(
        rulename=rule_name,
        threads=threads,
        bench_iteration=bench_iteration,
        input_=input_,
        output=output_,
        log=log_,
        params=params_,
        wildcards=wildcards_,
        config=minimal_config_mutect2,
        scriptdir=scriptdir,
        resources=resources_,
    )


@pytest.fixture
def fake_fs():
    """Return ``namedtuple`` with fake file system objects."""
    klass = namedtuple("FakeFsBundle", "fs os open inter_process_lock")
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    fake_lock = MagicMock()
    return klass(fs=fake_fs, os=fake_os, open=fake_open, inter_process_lock=fake_lock)


@pytest.fixture
def somatic_variant_fake_fs(fake_fs, fai_file_content):
    """Return fake file system setup with files for the somatic variant calling workflow."""
    # Create work directory
    fake_fs.fs.makedirs("/work", exist_ok=True)
    # Create static files
    fake_fs.fs.create_file("/path/to/ref.fa", create_missing_dirs=True)
    fake_fs.fs.create_file("/path/to/ref.fa.fai", contents=fai_file_content)
    return fake_fs


def patch_module_fs(module_name, fake_fs, mocker):
    """Helper function to mock out the file-system related things in the module with the given
    name using the given fake_fs and pytest-mock mocker
    """
    mocker.patch("{}.os".format(module_name), fake_fs.os)
    mocker.patch("{}.open".format(module_name), fake_fs.open, create=True)
    mocker.patch("{}.os".format(module_name), fake_fs.os)
