# -*- coding: utf-8 -*-
"""Tests for find_file module code"""

import copy
import json
from unittest.mock import MagicMock, patch

from pyfakefs import fake_filesystem
import pytest

from snappy_pipeline.find_file import FileSystemCrawler, FileSystemCrawlerResult, PatternSet


# TODO: test the invalidation_paths parameter/feature


# Test PatternSet ---------------------------------------------------------------------------------


def test_pattern_set_no_names():
    obj = PatternSet(("*_R1.fastq.gz", "*_R2.fastq.gz"))
    assert obj.patterns == ("*_R1.fastq.gz", "*_R2.fastq.gz")
    assert obj.names == tuple()
    assert obj.named_patterns is None


def test_pattern_set_with_names():
    obj = PatternSet(("*_R1.fastq.gz", "*_R2.fastq.gz"), ("first", "second"))
    assert obj.patterns == ("*_R1.fastq.gz", "*_R2.fastq.gz")
    assert obj.names == ("first", "second")
    assert obj.named_patterns == {"first": "*_R1.fastq.gz", "second": "*_R2.fastq.gz"}


# Test FileSystemCrawlerResult --------------------------------------------------------------------


def test_file_system_crawler_result_no_names():
    obj = FileSystemCrawlerResult("/base", ("foo_R1.fastq.gz", "foo_R2.fastq.gz"))
    assert obj.base_folder == "/base"
    assert obj.names == tuple()
    assert obj.named_files is None
    with pytest.raises(ValueError):
        obj.to_dict()
    assert str(obj) == (
        """FileSystemCrawlerResult('/base', ('foo_R1.fastq.gz', 'foo_R2.fastq.gz'), (), None)"""
    )


def test_file_system_crawler_result_with_names():
    obj = FileSystemCrawlerResult(
        "/base", ("foo_R1.fastq.gz", "foo_R2.fastq.gz"), ("first", "second")
    )
    assert obj.base_folder == "/base"
    assert obj.names == ("first", "second")
    assert obj.named_files == {"first": "foo_R1.fastq.gz", "second": "foo_R2.fastq.gz"}
    assert obj.to_dict() == {"first": "foo_R1.fastq.gz", "second": "foo_R2.fastq.gz"}
    # OrderedDict str representation changes depending on python implementation
    assert str(obj) == (
        "FileSystemCrawlerResult('/base', ('foo_R1.fastq.gz', 'foo_R2.fastq.gz'), "
        "('first', 'second'), OrderedDict([('first', 'foo_R1.fastq.gz'), "
        "('second', 'foo_R2.fastq.gz')]))"
    ) or str(obj) == (
               "FileSystemCrawlerResult('/base', ('foo_R1.fastq.gz', 'foo_R2.fastq.gz'), "
               "('first', 'second'), OrderedDict({'first': 'foo_R1.fastq.gz', "
               "'second': 'foo_R2.fastq.gz'}))"
           )


# Test FileSystemCrawler --------------------------------------------------------------------------

# Path to test file
CACHE_PATH = "/tmp/cache_file"


@pytest.fixture
def sample_cache_dict():
    return {
        "cache_version": 1,
        "root_dirs": {
            "/path": [
                ["P001", "flowcell", "lane", "P001_R1.fastq.gz"],
                ["P001", "flowcell", "lane", "P001_R2.fastq.gz"],
            ]
        },
    }


def test_file_system_crawler_invalidate_cache(sample_cache_dict):
    """Crawl file system. It starting with an existing CACHE_PATH that is older
    than the invalidation paths. Expected result: the original cache will be
    set as 'invalidated' and replaced in the first iteration; the file won't be
    replaced in the second iteration.

       The files must be present in the fake file system
    """
    # Initialise variables
    new_file = "/path/new_file.tsv"
    invalidation_paths = ["/path/P001", new_file]
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))

    # Prepare fake system
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict), create_missing_dirs=True
    )
    fake_fs.create_file("/path/P001/flowcell/lane/P001_R1.fastq.gz")
    fake_fs.create_file("/path/P001/flowcell/lane/P001_R2.fastq.gz")
    fake_fs.create_file(new_file, contents="original")
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        # Get the original modification time
        original_cache_file_time = fake_os.path.getmtime(CACHE_PATH)

        # Run the crawler
        # First iteration: should be invalidated, older than invalidation paths
        crawler_first = FileSystemCrawler(CACHE_PATH, invalidation_paths)
        _ = list(crawler_first.run("/path", "P001", (pattern_set,), False))
        crawler_first.save_cache()
        first_iteration_cache_file_time = fake_os.path.getmtime(CACHE_PATH)

        # Run the crawler
        # Second iteration: shouldn't be invalidated, newer than invalidation paths
        crawler_second = FileSystemCrawler(CACHE_PATH, invalidation_paths)
        _ = list(crawler_second.run("/path", "P001", (pattern_set,), False))
        crawler_second.save_cache()
        second_iteration_cache_file_time = fake_os.path.getmtime(CACHE_PATH)

    # Assert first iteration
    assert first_iteration_cache_file_time > original_cache_file_time
    assert crawler_first.cache_invalidated
    # Assert second iteration
    assert first_iteration_cache_file_time == second_iteration_cache_file_time
    assert not crawler_second.cache_invalidated


def test_file_system_crawler_construct_existing_cache(sample_cache_dict):
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        crawler.save_cache()
    assert crawler.cache_path == CACHE_PATH
    assert crawler.cache == sample_cache_dict
    assert crawler.invalidation_paths == []
    assert not crawler.cache_dirty
    assert crawler.lock_timeout == 60


def test_file_system_crawler_construct_no_existing_cache():
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_fs.create_dir(fake_os.path.dirname(CACHE_PATH))
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        crawler.save_cache()
    assert crawler.cache_path == CACHE_PATH
    empty_cache = {"cache_version": 1, "root_dirs": {}}
    assert crawler.cache == empty_cache
    assert crawler.invalidation_paths == []
    assert crawler.cache_dirty
    assert crawler.lock_timeout == 60
    assert fake_os.path.exists(CACHE_PATH)
    assert json.loads(fake_open(CACHE_PATH).read()) == empty_cache


def test_file_system_crawler_crawl_existing_cache(sample_cache_dict):
    """Crawl file system, starting off existing cache"""
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        res = list(crawler.run("/path", "P001", (pattern_set,), False))
    assert len(res) == 1
    assert res[0].base_folder == "/path/P001"
    assert res[0].files == (
        "/path/P001/flowcell/lane/P001_R1.fastq.gz",
        "/path/P001/flowcell/lane/P001_R2.fastq.gz",
    )
    assert res[0].names == ("left", "right")
    assert res[0].named_files == {
        "left": "/path/P001/flowcell/lane/P001_R1.fastq.gz",
        "right": "/path/P001/flowcell/lane/P001_R2.fastq.gz",
    }


def test_file_system_crawler_crawl_mismatching_file_count(sample_cache_dict):
    """Crawl file system (with existing cache) with mismatching files"""
    # Update cache to include the bogus file leading to imbalance
    sample_cache_dict = copy.deepcopy(sample_cache_dict)
    sample_cache_dict["root_dirs"]["/path"].append(
        ["P001", "flowcell", "lane2", "P001_R2.fastq.gz"]
    )
    # Go on as in the other tests
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        with pytest.raises(ValueError) as e:
            list(crawler.run("/path", "P001", (pattern_set,), False))
    assert e.match("Must have the same number.*")


def test_file_system_crawler_crawl_few_matches(sample_cache_dict, capsys):
    """Crawl file system (with existing cache), too few matches"""
    # Update cache to include the bogus file leading to imbalance
    sample_cache_dict = copy.deepcopy(sample_cache_dict)
    sample_cache_dict["root_dirs"]["/path"].pop()
    # Go on as in the other tests
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        res = list(crawler.run("/path", "P001", (pattern_set,), False))
    assert not res
    assert "WARNING: Skipping matches" in "\n".join(capsys.readouterr())


def test_file_system_crawler_crawl_no_existing_cache():
    """Crawl file system, starting without existing CACHE_PATH

    The files must be present in the fake file system
    """
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file("/path/P001/flowcell/lane/P001_R1.fastq.gz")
    fake_fs.create_file("/path/P001/flowcell/lane/P001_R2.fastq.gz")
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        res = list(crawler.run("/path", "P001", (pattern_set,), False))
    assert len(res) == 1
    assert res[0].base_folder == "/path/P001"
    assert res[0].files == (
        "/path/P001/flowcell/lane/P001_R1.fastq.gz",
        "/path/P001/flowcell/lane/P001_R2.fastq.gz",
    )
    assert res[0].names == ("left", "right")
    assert res[0].named_files == {
        "left": "/path/P001/flowcell/lane/P001_R1.fastq.gz",
        "right": "/path/P001/flowcell/lane/P001_R2.fastq.gz",
    }


# Test scenarious for finding SE files with PE pattern.


@pytest.fixture
def sample_cache_dict_se_only():  # only SE data => OK
    return {
        "cache_version": 1,
        "root_dirs": {
            "/path": [
                ["P001", "flowcell", "lane", "P001_R1.fastq.gz"],
                ["P001", "flowcell2", "lane", "P001_R1.fastq.gz"],
            ]
        },
    }


@pytest.fixture
def sample_cache_dict_se_pe_mixed():  # first PE then SE => error
    return {
        "cache_version": 1,
        "root_dirs": {
            "/path": [
                ["P001", "flowcell", "lane", "P001_R1.fastq.gz"],
                ["P001", "flowcell", "lane", "P001_R2.fastq.gz"],
                ["P001", "flowcell2", "lane", "P001_R1.fastq.gz"],
            ]
        },
    }


@pytest.fixture
def sample_cache_dict_pe_se_mixed():  # first SE then PE => error
    return {
        "cache_version": 1,
        "root_dirs": {
            "/path": [
                ["P001", "flowcell", "lane", "P001_R1.fastq.gz"],
                ["P001", "flowcell2", "lane", "P001_R1.fastq.gz"],
                ["P001", "flowcell2", "lane", "P001_R2.fastq.gz"],
            ]
        },
    }


def test_file_system_crawler_se_data_pe_pattern_good(sample_cache_dict_se_only):
    """Crawl file system, starting off existing cache"""
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict_se_only), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        res = list(crawler.run("/path", "P001", (pattern_set,), True))
    assert len(res) == 2
    assert res[0].base_folder == "/path/P001"
    assert res[0].files == ("/path/P001/flowcell/lane/P001_R1.fastq.gz",)
    assert res[0].names == ("left",)
    assert res[0].named_files == {"left": "/path/P001/flowcell/lane/P001_R1.fastq.gz"}
    assert res[1].base_folder == "/path/P001"
    assert res[1].files == ("/path/P001/flowcell2/lane/P001_R1.fastq.gz",)
    assert res[1].names == ("left",)
    assert res[1].named_files == {"left": "/path/P001/flowcell2/lane/P001_R1.fastq.gz"}


def test_file_system_crawler_se_data_pe_pattern_bad_se_pe(sample_cache_dict_se_pe_mixed):
    """Crawl file system, starting off existing cache"""
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict_se_pe_mixed), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        with pytest.raises(ValueError) as excinfo:
            list(crawler.run("/path", "P001", (pattern_set,), True))
        assert str(excinfo.value).startswith("Must have the same number of matches per pattern,")


def test_file_system_crawler_se_data_pe_pattern_bad_pe_se(sample_cache_dict_pe_se_mixed):
    """Crawl file system, starting off existing cache"""
    pattern_set = PatternSet(("*/*/*_R1.fastq.gz", "*/*/*_R2.fastq.gz"), ("left", "right"))
    fake_fs = fake_filesystem.FakeFilesystem()
    fake_fs.create_file(
        CACHE_PATH, contents=json.dumps(sample_cache_dict_pe_se_mixed), create_missing_dirs=True
    )
    fake_os = fake_filesystem.FakeOsModule(fake_fs)
    fake_open = fake_filesystem.FakeFileOpen(fake_fs)
    mock_lock = MagicMock()
    with patch("snappy_pipeline.find_file.os", fake_os), patch(
            "snappy_pipeline.find_file.InterProcessLock", mock_lock
    ), patch("snappy_pipeline.find_file.open", fake_open, create=True):
        crawler = FileSystemCrawler(CACHE_PATH, [])
        with pytest.raises(ValueError) as excinfo:
            list(crawler.run("/path", "P001", (pattern_set,), True))
        assert str(excinfo.value).startswith("Must have the same number of matches per pattern,")
