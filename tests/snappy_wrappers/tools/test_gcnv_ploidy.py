"""Tests for ``snappy_wrappers.tools.gcnv_ploidy``"""

import pytest

from snappy_wrappers.tools import gcnv_ploidy

MALE = gcnv_ploidy.PLOIDY["1"]
HEADER = ["@RG\tID:GATKCopyNumber\tSM:index\n", "CONTIG\tPLOIDY\tPLOIDY_GQ\n"]


@pytest.mark.parametrize(
    "contig, kind",
    [("chrX", "x"), ("X", "x"), ("chry", "y"), ("Y", "y"), ("chr1", None), ("hs37d5", None)]
    + [("chrUn_decoy", None), ("chrX_random", None), ("chrEBV", None)],
)
def test_sex_chromosomes_match_exactly(contig, kind):
    assert gcnv_ploidy.sex_chromosome(contig) == kind


def test_x_and_y_follow_the_sex_and_keep_the_call_as_marker():
    lines = HEADER + ["chr1\t3\t50.1\n", "chrX\t2\t50.1\n", "chrY\t1\t50.1\n", "chr1_decoy\t0\t9\n"]
    assert gcnv_ploidy.fixed_calls(lines, MALE, autosomes_diploid=False) == [
        "@RG\tID:GATKCopyNumber\tSM:index",
        "CONTIG\tPLOIDY\tPLOIDY_GQ",
        "chr1\t3\t50.1",
        "chrX\t1\t42.0002",
        "chrY\t1\t50.1",
        "chr1_decoy\t0\t9",
    ]


def test_case_mode_sets_the_autosomes_to_two():
    lines = HEADER + ["1\t3\t50.1\n", "2\t2\t50.1\n", "X\t1\t50.1\n"]
    assert gcnv_ploidy.fixed_calls(lines, MALE, autosomes_diploid=True)[2:] == [
        "1\t2\t42.0003",
        "2\t2\t50.1",
        "X\t1\t50.1",
    ]


def _ploidy_dir(tmp_path, samples):
    for i, sample in enumerate(samples):
        sample_dir = tmp_path / "ploidy" / "ploidy-calls" / f"SAMPLE_{i}"
        sample_dir.mkdir(parents=True)
        (sample_dir / "sample_name.txt").write_text(f"{sample}\n")
        (sample_dir / "contig_ploidy.tsv").write_text("".join(HEADER) + "chrX\t1\t60\n")
    ped = tmp_path / "fam.ped"
    ped.write_text("FAM\tindex\t0\t0\t2\t2\n")
    return tmp_path / "ploidy", ped


def test_main_rewrites_the_calls_of_every_sample(tmp_path):
    ploidy_dir, ped = _ploidy_dir(tmp_path, ["index"])
    gcnv_ploidy.main([str(ploidy_dir), str(ped)])
    calls = (ploidy_dir / "ploidy-calls" / "SAMPLE_0" / "contig_ploidy.tsv").read_text()
    assert calls.endswith("chrX\t2\t42.0001\n")


def test_a_sample_without_sex_fails(tmp_path):
    ploidy_dir, ped = _ploidy_dir(tmp_path, ["index", "stranger"])
    with pytest.raises(ValueError, match="sample stranger is in no PED file"):
        gcnv_ploidy.main(["--autosomes-diploid", str(ploidy_dir), str(ped)])


def test_an_unknown_sex_fails(tmp_path):
    ped = tmp_path / "fam.ped"
    ped.write_text("FAM\tindex\t0\t0\t0\t2\n")
    with pytest.raises(ValueError, match="sample index has sex '0'"):
        gcnv_ploidy.read_sexes([str(ped)])
