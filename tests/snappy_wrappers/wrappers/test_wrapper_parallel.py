import textwrap

import pytest

from snappy_wrappers.genome_regions import GenomeRegion
from snappy_wrappers.wrapper_parallel import (
    ParallelSomaticVariantCallingBaseWrapper,
    ResourceUsage,
    SgeResourceUsageConverter,
    days,
    gib,
    hours,
    kib,
    mib,
    minutes,
)

from .conftest import patch_module_fs


class FakeParallelMutect2Wrapper(ParallelSomaticVariantCallingBaseWrapper):
    """Fake parallel execution of Mutect2.

    It is used among other things as a way to test the class ``ParallelBaseWrapper``.
    """

    inner_wrapper = "mutect2/run"
    step_name = "somatic_variant_calling"
    tool_name = "mutect2"
    realpath_output_keys = (
        "raw",
        "raw_md5",
        "raw_tbi",
        "raw_tbi_md5",
        "stats",
        "stats_md5",
        "f1r2",
        "f1r2_md5",
    )
    key_ext = {
        "raw": "vcf.gz",
        "raw_md5": "vcf.gz.md5",
        "raw_tbi": "vcf.gz.tbi",
        "raw_tbi_md5": "vcf.gz.tbi.md5",
        "stats": "vcf.stats",
        "stats_md5": "vcf.stats.md5",
        "f1r2": "f1r2_tar.tar.gz",
        "f1r2_md5": "f1r2_tar.tar.gz.md5",
    }


@pytest.fixture
def resource_usage():
    """Returns ResourceUsage object."""
    return ResourceUsage(
        cores=2,
        memory=gib(1.0),
        duration=hours(3),
    )


@pytest.fixture
def fake_parallel_wrapper(snakemake_obj, somatic_variant_fake_fs, mocker):
    """Returns FakeParallelMutect2Wrapper object."""
    # Patch out file-system
    patch_module_fs("snappy_wrappers.wrapper_parallel", somatic_variant_fake_fs, mocker)
    return FakeParallelMutect2Wrapper(snakemake=snakemake_obj)


# Test isolated methods ----------------------------------------------------------------------------


def test_kib():
    """Tests wrapper_parallel.kib() call."""
    # Define expected dictionary
    expected_dict = {1: 1024, 2: 2048, 3: 3072, 1.5: 1024}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = kib(input_)
        assert actual == expected


def test_kib_exception():
    """Tests wrapper_parallel.kib() exception."""
    # Test raise ValueError
    with pytest.raises(ValueError):
        kib("One")
    # Test raise TypeError
    bad_input_list = [{1}, [1, 2]]
    for input_ in bad_input_list:
        with pytest.raises(TypeError):
            kib(input_)


def test_mib():
    """Tests wrapper_parallel.mib() call."""
    # Define expected dictionary
    expected_dict = {1: 1048576, 2: 2097152, 3: 3145728, 1.5: 1048576}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = mib(input_)
        assert actual == expected


def test_mib_exception():
    """Tests wrapper_parallel.mib() exception."""
    # Test raise ValueError
    with pytest.raises(ValueError):
        mib("One")
    # Test raise TypeError
    bad_input_list = [{1}, [1, 2]]
    for input_ in bad_input_list:
        with pytest.raises(TypeError):
            mib(input_)


def test_gib():
    """Tests wrapper_parallel.gib() call."""
    # Define expected dictionary
    expected_dict = {1: 1073741824, 2: 2147483648, 3: 3221225472, 1.5: 1073741824}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = gib(input_)
        assert actual == expected


def test_gib_exception():
    """Tests wrapper_parallel.gib() exception."""
    # Test raise ValueError
    with pytest.raises(ValueError):
        gib("One")
    # Test raise TypeError
    bad_input_list = [{1}, [1, 2]]
    for input_ in bad_input_list:
        with pytest.raises(TypeError):
            gib(input_)


def test_minutes():
    """Tests wrapper_parallel.minutes()"""
    # Define expected values
    expected_dict = {1: "0:01:00", 0.5: "0:00:30"}
    # Get actual values and assert
    for input_, expected in expected_dict.items():
        actual = str(minutes(input_))
        assert actual == expected


def test_hours():
    """Tests wrapper_parallel.hours()"""
    # Define expected values
    expected_dict = {1: "1:00:00", 0.5: "0:30:00"}
    # Get actual values and assert
    for input_, expected in expected_dict.items():
        actual = str(hours(input_))
        assert actual == expected


def test_days():
    """Tests wrapper_parallel.days()"""
    # Define expected values
    expected_dict = {1: "1 day, 0:00:00", 0.5: "12:00:00"}
    # Get actual values and assert
    for input_, expected in expected_dict.items():
        actual = str(days(input_))
        assert actual == expected


# Test SgeResourceUsageConverter  ------------------------------------------------------------------


def test_sge_resource_usage_converter_to_qsub_args(resource_usage):
    """Tests SgeResourceUsageConverter.to_qsub_args()"""
    sge_converter = SgeResourceUsageConverter(res_usage=resource_usage)
    # Define expected
    expected = ["--ntasks=2", "--time=03:00", "--mem=1024"]
    # Get actual and assert
    actual = sge_converter.to_qsub_args()
    assert actual == expected


def test_sge_resource_usage_converter_to_res_dict(resource_usage):
    """Tests SgeResourceUsageConverter.to_res_dict()"""
    sge_converter = SgeResourceUsageConverter(res_usage=resource_usage)
    # Define expected
    expected = {"ntasks": 2, "time": "03:00", "mem": 1024}
    # Get actual and assert
    actual = sge_converter.to_res_dict()
    assert actual == expected


# Test ParallelBaseWrapper   -----------------------------------------------------------------------


def test_parallel_base_wrapper_get_fai_path(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_fai_path()"""
    # Define expected - defined in `minimal_config`: static_data_config/reference/path
    expected = "/path/to/ref.fa.fai"
    # Get actual and assert
    actual = fake_parallel_wrapper.get_fai_path()
    assert actual == expected


def test_parallel_base_wrapper_get_all_log_files(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_all_log_files()"""
    # Define expected - defined in snakemake.log
    base_path = "/work/bwa.mutect2.P001-T1-DNA1-WGS1/log/bwa.mutect2.P001-T1-DNA1-WGS1"
    expected = {
        "conda_info": base_path + ".conda_info.txt",
        "conda_info_md5": base_path + ".conda_info.txt.md5",
        "conda_list": base_path + ".conda_list.txt",
        "conda_list_md5": base_path + ".conda_list.txt.md5",
    }
    # Get actual and assert
    actual = fake_parallel_wrapper.get_all_log_files()
    assert actual == expected


def test_parallel_base_wrapper_get_all_output(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_all_output()"""
    # Define expected - defined in snakemake.output
    base_path = "/work/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1"
    expected = {
        "raw": base_path + ".raw.vcf.gz",
        "raw_md5": base_path + ".raw.vcf.gz.md5",
        "raw_tbi": base_path + ".raw.vcf.gz.tbi",
        "raw_tbi_md5": base_path + ".raw.vcf.gz.tbi.md5",
        "stats": base_path + ".raw.vcf.stats",
        "stats_md5": base_path + ".raw.vcf.stats.md5",
        "f1r2": base_path + ".raw.f1r2_tar.tar.gz",
        "f1r2_md5": base_path + ".raw.f1r2_tar.tar.gz.md5",
    }
    # Get actual and assert
    actual = fake_parallel_wrapper.get_all_output()
    assert actual == expected


def test_parallel_base_wrapper_get_regions(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_regions()"""
    # Define expected
    genome_region_chr1_list = [
        (0, 50010000),
        (49990000, 100010000),
        (99990000, 150010000),
        (149990000, 200010000),
        (199990000, 249250621),
    ]
    expected = [
        GenomeRegion(chrom="1", begin=region[0], end=region[1])
        for region in genome_region_chr1_list
    ]
    # Get actual and assert
    actual = fake_parallel_wrapper.get_regions()
    assert actual == expected


def test_parallel_base_wrapper_get_job_mult_memory(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_job_mult_memory()"""
    # Define expected
    # from: `minimal_config`: step_config/somatic_variant_calling/mutect2/job_mult_memory
    expected = 2
    # Get actual and assert
    actual = fake_parallel_wrapper.get_job_mult_memory()
    assert actual == expected


def test_parallel_base_wrapper_get_job_mult_time(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_job_mult_time()"""
    # Define expected
    # from: `minimal_config`: step_config/somatic_variant_calling/mutect2/job_mult_time
    expected = 3
    # Get actual and assert
    actual = fake_parallel_wrapper.get_job_mult_time()
    assert actual == expected


def test_parallel_base_wrapper_get_merge_mult_memory(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_merge_mult_memory()"""
    # Define expected
    # from: `minimal_config`: step_config/somatic_variant_calling/mutect2/merge_mult_memory
    expected = 4
    # Get actual and assert
    actual = fake_parallel_wrapper.get_merge_mult_memory()
    assert actual == expected


def test_parallel_base_wrapper_get_merge_mult_time(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_merge_mult_time()"""
    # Define expected
    # from: `minimal_config`: step_config/somatic_variant_calling/mutect2/merge_mult_time
    expected = 5
    # Get actual and assert
    actual = fake_parallel_wrapper.get_merge_mult_time()
    assert actual == expected


def test_parallel_base_wrapper_get_window_length(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_window_length()"""
    # Define expected
    # from: `minimal_config`: step_config/somatic_variant_calling/mutect2/window_length
    expected = 50000000
    # Get actual and assert
    actual = fake_parallel_wrapper.get_window_length()
    assert actual == expected


def test_parallel_base_wrapper_get_ignore_chroms(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.get_ignore_chroms()"""
    # Define expected
    # from: `minimal_config`: step_config/somatic_variant_calling/mutect2/ignore_chroms
    expected = ["NC_007605", "hs37d5", "chrEBV", "*_decoy", "HLA-*", "GL000220.*"]
    # Get actual and assert
    actual = fake_parallel_wrapper.get_ignore_chroms()
    assert actual == expected


def test_parallel_base_wrapper_construct_parallel_result_files(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.construct_parallel_result_files()"""
    # Define expected - based on number of genomic regions in chr1 with window length of 50MB
    expected = [
        "job_out.0.d/.done",
        "job_out.1.d/.done",
        "job_out.2.d/.done",
        "job_out.3.d/.done",
        "job_out.4.d/.done",
    ]
    # Get actual and assert
    actual = fake_parallel_wrapper.construct_parallel_result_files()
    assert actual == expected


def test_parallel_base_wrapper_construct_preamble(fake_parallel_wrapper, snakemake_output_dict):
    """Tests ParallelBaseWrapper.construct_preamble()"""
    # Simplistic way to emulate `ParallelBaseWrapper._apply_realpath_to_output()`
    realpath_applied = {}
    for key, value in snakemake_output_dict.items():
        realpath_applied[key] = "/" + value
    # Define expected
    expected = (
        textwrap.dedent(
            """
        shell.executable("/bin/bash")
        shell.prefix("set -ex;")

        configfile: 'config.json'

        localrules: all

        rule all:
            input: **{out_dict}
        """
        )
        .lstrip()
        .format(out_dict=realpath_applied)
    )
    # Get actual and assert
    actual = fake_parallel_wrapper.construct_preamble()
    assert actual == expected


def test_parallel_base_wrapper_construct_epilogue(fake_parallel_wrapper):
    """Tests ParallelBaseWrapper.construct_epilogue()"""
    # Define expected
    expected = ""
    # Get actual and assert
    actual = fake_parallel_wrapper.construct_epilogue()
    assert actual == expected
