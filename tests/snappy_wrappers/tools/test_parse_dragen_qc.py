from pathlib import Path

import pytest

from snappy_wrappers.tools.parse_dragen_qc import create_varfish_json, find_coverage_files

OVERALL_MEAN_COV_TEXT = "Average alignment coverage over ./../test.bed, 42.134"
FINE_HIST_TEXT = """Depth,Overall
0,1345062
1,154347
2,95633
3,83481
4,77150
5,71957
6,70995
7,74788
8,83818
9,95822
10,119129"""

COV_REPORT_TEXT = """#chrom	start	end	total_cvg	mean_cvg	Q1_cvg	median_cvg	Q3_cvg	min_cvg	max_cvg	pct_above_5	pct_above_15	pct_above_20	pct_above_30	pct_above_50	pct_above_100	pct_above_200	pct_above_300	pct_above_400	pct_above_500	pct_above_1000
1	100	200	35292	59.41	49.00	62.00	71.00	27	87	100.00	100.00	100.00	96.63	72.05	0.00	0.00	0.00	0.00	0.00	0.00
1	99999	199999	17515	40.17	18.00	28.00	64.00	8	108	100.00	86.70	68.35	47.02	30.50	4.59	0.00	0.00	0.00	0.00	0.00
1	10	11	99531	49.10	11.00	48.00	77.00	0	148	83.82	68.38	64.78	63.64	47.46	11.15	0.00	0.00	0.00	0.00	0.00
"""

MAPPING_METRICS_TEXT = """MAPPING/ALIGNING SUMMARY,,Total input reads,90000,100.00
MAPPING/ALIGNING SUMMARY,,Number of duplicate marked reads,900,11.77
MAPPING/ALIGNING SUMMARY,,Number of duplicate marked and mate reads removed,NA
MAPPING/ALIGNING SUMMARY,,Number of unique reads (excl. duplicate marked reads),900,88.23
MAPPING/ALIGNING SUMMARY,,Reads with mate sequenced,9999,100.00
MAPPING/ALIGNING SUMMARY,,Reads without mate sequenced,0,0.00
MAPPING/ALIGNING SUMMARY,,QC-failed reads,0,0.00
MAPPING/ALIGNING SUMMARY,,Mapped reads,9999,99.93
MAPPING/ALIGNING SUMMARY,,Mapped reads adjusted for filtered mapping,999,99.93
MAPPING/ALIGNING SUMMARY,,Mapped reads R1,999,99.99
MAPPING/ALIGNING SUMMARY,,Mapped reads R2,999,99.87
MAPPING/ALIGNING SUMMARY,,Number of unique & mapped reads (excl. duplicate marked reads),999,88.16
MAPPING/ALIGNING SUMMARY,,Unmapped reads,999,0.07
MAPPING/ALIGNING SUMMARY,,Unmapped reads adjusted for filtered mapping,999,0.07
MAPPING/ALIGNING SUMMARY,,Adjustment of reads matching non-reference decoys,0,0.00
MAPPING/ALIGNING SUMMARY,,Singleton reads (itself mapped; mate unmapped),999,0.06
MAPPING/ALIGNING SUMMARY,,Paired reads (itself & mate mapped),999,99.87
MAPPING/ALIGNING SUMMARY,,Properly paired reads,999,98.06
MAPPING/ALIGNING SUMMARY,,Not properly paired reads (discordant),999,1.81
MAPPING/ALIGNING SUMMARY,,Paired reads mapped to different chromosomes,999,1.30
MAPPING/ALIGNING SUMMARY,,Paired reads mapped to different chromosomes (MAPQ>=10),999,1.01
MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [40:inf),999,84.00
MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [30:40),999,1.00
MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [20:30),999,5.00
MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [10:20),999,5.00
MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [ 0:10),999,5.00
MAPPING/ALIGNING SUMMARY,,Reads with MAPQ NA (Unmapped reads),999,0.07
MAPPING/ALIGNING SUMMARY,,Reads with indel R1,999,3.28
MAPPING/ALIGNING SUMMARY,,Reads with indel R2,999,3.48
MAPPING/ALIGNING SUMMARY,,Total bases,10
MAPPING/ALIGNING SUMMARY,,Total bases R1,5
MAPPING/ALIGNING SUMMARY,,Total bases R2,5
MAPPING/ALIGNING SUMMARY,,Mapped bases,1
MAPPING/ALIGNING SUMMARY,,Mapped bases R1,1
MAPPING/ALIGNING SUMMARY,,Mapped bases R2,2
MAPPING/ALIGNING SUMMARY,,Soft-clipped bases,8,0.84
MAPPING/ALIGNING SUMMARY,,Soft-clipped bases R1,4,0.61
MAPPING/ALIGNING SUMMARY,,Soft-clipped bases R2,4,1.06
MAPPING/ALIGNING SUMMARY,,Hard-clipped bases,0,0.00
MAPPING/ALIGNING SUMMARY,,Hard-clipped bases R1,0,0.00
MAPPING/ALIGNING SUMMARY,,Hard-clipped bases R2,0,0.00
MAPPING/ALIGNING SUMMARY,,Mismatched bases R1,1,0.40
MAPPING/ALIGNING SUMMARY,,Mismatched bases R2,1,0.57
MAPPING/ALIGNING SUMMARY,,Mismatched bases R1 (excl. indels),1,0.34
MAPPING/ALIGNING SUMMARY,,Mismatched bases R2 (excl. indels),1,0.51
MAPPING/ALIGNING SUMMARY,,Q30 bases,1,92.45
MAPPING/ALIGNING SUMMARY,,Q30 bases R1,1,94.08
MAPPING/ALIGNING SUMMARY,,Q30 bases R2,1,90.81
MAPPING/ALIGNING SUMMARY,,Q30 bases (excl. dups & clipped bases),10
MAPPING/ALIGNING SUMMARY,,Total alignments,8
MAPPING/ALIGNING SUMMARY,,Secondary alignments,0
MAPPING/ALIGNING SUMMARY,,Supplementary (chimeric) alignments,3
MAPPING/ALIGNING SUMMARY,,Estimated read length,150.03
MAPPING/ALIGNING SUMMARY,,Bases in reference genome,31
MAPPING/ALIGNING SUMMARY,,Bases in target bed [% of genome],NA
MAPPING/ALIGNING SUMMARY,,Insert length: mean,400.00
MAPPING/ALIGNING SUMMARY,,Insert length: median,534.00
MAPPING/ALIGNING SUMMARY,,Insert length: standard deviation,9232.33
MAPPING/ALIGNING SUMMARY,,Provided sex chromosome ploidy,NA
MAPPING/ALIGNING SUMMARY,,Estimated sample contamination,0.0010
MAPPING/ALIGNING SUMMARY,,DRAGEN mapping rate [mil. reads/second],100.06"""

IDXSTATS_TEXT = """X\t155270560\t10\t2
Y\t59373566\t13\t1
MT\t41234123\t20\t99"""

FASTQC_METRICS_TEXT = """READ MEAN QUALITY,Read1,Q2 Reads,2736
READ MEAN QUALITY,Read1,Q13 Reads,2"""


@pytest.fixture
def dragen_case_complete(fs):
    """Create an artificial DRAGEN case."""
    sample_name = "XX00_1234"
    date = "2001-01-01"
    region_id = "1"
    region_name = f"qc-coverage-region-{region_id}"

    misc_path = f"/tmp/{sample_name}/misc/{date}"

    overall_mean_cov_path = f"{misc_path}/{sample_name}_dragen.{region_name}_overall_mean_cov.csv"
    fs.create_file(overall_mean_cov_path, create_missing_dirs=True, contents=OVERALL_MEAN_COV_TEXT)

    fine_hist_path = f"{misc_path}/{sample_name}_dragen.{region_name}_fine_hist.csv"
    fs.create_file(fine_hist_path, create_missing_dirs=True, contents=FINE_HIST_TEXT)

    cov_report_path = f"{misc_path}/{sample_name}_dragen.{region_name}_cov_report.bed"
    fs.create_file(cov_report_path, create_missing_dirs=True, contents=COV_REPORT_TEXT)

    mapping_metrics_path = f"{misc_path}/{sample_name}_dragen.mapping_metrics.csv"
    fs.create_file(mapping_metrics_path, create_missing_dirs=True, contents=MAPPING_METRICS_TEXT)

    idxstats_path = f"{misc_path}/{sample_name}_dragen.bam.idxstats.txt"
    fs.create_file(idxstats_path, create_missing_dirs=True, contents=IDXSTATS_TEXT)

    fastqc_path = f"{misc_path}/{sample_name}_dragen.fastqc_metrics.csv"
    fs.create_file(fastqc_path, create_missing_dirs=True, contents=FASTQC_METRICS_TEXT)

    wrong_cov_report_path = f"{misc_path}/{sample_name}_dragen.{region_name}_read_cov_report.bed"
    fs.create_file(wrong_cov_report_path, create_missing_dirs=True, contents="")

    yield fs


@pytest.mark.usefixtures("dragen_case_complete")
def test_create_varfish_json():

    region_id = "1"
    input_path = Path("/tmp/XX00_1234")
    all_paths = find_coverage_files(input_path)
    found_paths = all_paths[region_id]

    result_json = create_varfish_json("TEST", found_paths)
    expected_json = {
        "TEST": {
            "summary": {"mean coverage": 42.134, "target_count": 3, "total_target_size": 100101},
            "min_cov_base": {
                "0": 1.0,
                "10": 0.05242933884697617,
                "20": 0.0,
                "30": 0.0,
                "40": 0.0,
                "50": 0.0,
                "60": 0.0,
                "70": 0.0,
                "80": 0.0,
                "90": 0.0,
                "100": 0.0,
                "110": 0.0,
                "120": 0.0,
                "130": 0.0,
                "140": 0.0,
                "150": 0.0,
                "160": 0.0,
                "170": 0.0,
                "180": 0.0,
                "190": 0.0,
            },
            "min_cov_target": {
                "0": 1.0,
                "10": 0.05242933884697617,
                "20": 0.0,
                "30": 0.0,
                "40": 0.0,
                "50": 0.0,
                "60": 0.0,
                "70": 0.0,
                "80": 0.0,
                "90": 0.0,
                "100": 0.0,
                "110": 0.0,
                "120": 0.0,
                "130": 0.0,
                "140": 0.0,
                "150": 0.0,
                "160": 0.0,
                "170": 0.0,
                "180": 0.0,
                "190": 0.0,
            },
            "idxstats": {
                "X": {"mapped": 10, "unmapped": 2},
                "Y": {"mapped": 13, "unmapped": 1},
                "MT": {"mapped": 20, "unmapped": 99},
            },
        }
    }
    assert result_json["TEST"]["min_cov_base"] == expected_json["TEST"]["min_cov_base"]
    assert result_json["TEST"]["min_cov_target"] == expected_json["TEST"]["min_cov_target"]
    assert result_json["TEST"]["idxstats"] == expected_json["TEST"]["idxstats"]
    assert result_json["TEST"]["summary"] == expected_json["TEST"]["summary"]
