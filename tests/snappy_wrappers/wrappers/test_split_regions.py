"""The window mode of the ``split_regions`` wrapper (popdel): fixed-length padded windows, shared by the parts in order."""

import importlib.util
from pathlib import Path

import pytest

WRAPPER = Path(__file__).parents[3] / "snappy_wrappers/wrappers/split_regions/wrapper.py"


@pytest.fixture(scope="module")
def split():
    spec = importlib.util.spec_from_file_location("split_regions_wrapper", WRAPPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_contigs_are_cut_from_their_start(split):
    windows = split.cut_windows({"1": 25, "2": 8}, 10)

    assert windows == [
        split.Interval("1", 0, 10),
        split.Interval("1", 10, 20),
        split.Interval("1", 20, 25),
        split.Interval("2", 0, 8),
    ]


def test_parts_take_windows_in_order_with_about_equal_bases(split):
    windows = split.cut_windows({"1": 40, "2": 10, "3": 10}, 10)

    parts = split.partition_windows(windows, 3)

    assert [[(w.contig, w.start) for w in part] for part in parts] == [
        [("1", 0), ("1", 10)],
        [("1", 20), ("1", 30)],
        [("2", 0), ("3", 0)],
    ]


def test_every_part_gets_a_window(split):
    windows = [split.Interval("1", 0, 1), split.Interval("2", 0, 1), split.Interval("3", 0, 100)]

    parts = split.partition_windows(windows, 3)

    assert parts == [[window] for window in windows]


def test_fewer_windows_than_parts_is_an_error(split):
    with pytest.raises(ValueError, match="2 windows for 3 parts"):
        split.partition_windows(split.cut_windows({"1": 20}, 10), 3)


def test_the_windows_are_padded_within_their_contig(split, tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("1\t25\t3\t60\t61\n2\t8\t30\t60\t61\nhs37d5\t5\t50\t60\t61\n")
    outputs = [tmp_path / f"{i}-of-2.region.bed" for i in (1, 2)]

    split.main(str(fai), ["^hs37d5$"], 2, [str(o) for o in outputs], window_size=10)

    assert [o.read_text() for o in outputs] == ["1\t0\t12\n1\t8\t22\n", "1\t18\t25\n2\t0\t8\n"]
