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
        }
    }
    assert result_json == expected_json
