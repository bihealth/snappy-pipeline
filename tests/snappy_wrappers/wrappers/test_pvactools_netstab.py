"""netMHCstabpan runs from the environment of the job (the pVACtools container), not through
apptainer."""

import importlib.util
import stat
from pathlib import Path

SCRIPT = Path(__file__).parents[3] / "snappy_wrappers/wrappers/pvactools/netstab/netMHCstabpan.py"

#: Stands in for the netMHCstabpan binary: one row per sequence, and the allele it used
FAKE_NETMHCSTABPAN = r"""#!/usr/bin/env bash
while [[ $# -gt 0 ]]; do
    case $1 in
        -a) allele=$2; shift 2 ;;
        -l) shift 2 ;;
        -xls) shift ;;
        -xlsfile) xls=$2; shift 2 ;;
        *) fasta=$1; shift ;;
    esac
done
echo "$TMPDIR" > "$NETMHCstabpan/tmpdir"
printf '\t\t\t%s\t\t\nPos\tPeptide\tID\tPred\tThalf(h)\tRank\tAve\tNB\n' "$allele" > "$xls"
paste - - < "$fasta" | while read -r id seq; do
    printf '0\t%s\t%s\t0.5\t1.0\t4.5\t0.5\t0\n' "$seq" "${id#>}" >> "$xls"
done
echo "$allele : Distance to traning data  0.000 (using nearest neighbor $allele)"
"""


def _script():
    spec = importlib.util.spec_from_file_location("netMHCstabpan", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_netmhcstabpan_runs_the_binary_of_the_job_environment(tmp_path):
    home = tmp_path / "netMHCstabpan-1.0" / "Linux_x86_64"
    (home / "bin").mkdir(parents=True)
    (home / "data").mkdir()
    (home / "data" / "MHC_pseudo.dat").write_text("HLA-A02:01 YFAMYGEKVAHTHVDTLYVRYHYYTW\n")
    binary = home / "bin" / "netMHCstabpan"
    binary.write_text(FAKE_NETMHCSTABPAN)
    binary.chmod(binary.stat().st_mode | stat.S_IXUSR)
    epitopes = [
        {"HLA Allele": "HLA-A*02:01", "MT Epitope Seq": "YLEPGPVTA"},
        {"HLA Allele": "HLA-B*07:02", "MT Epitope Seq": "KLVALGINA"},  # unknown to netMHCstabpan
    ]

    result = _script().run_netMHCstabpan(
        epitopes,
        hla_types={"HLA-A*02:01", "HLA-B*07:02"},
        lengths={9},
        base_path=str(tmp_path / "base"),
        netMHCstabpan_path=str(binary),
        netMHCpan_path=str(binary),
    )

    assert result[0]["Half Life"] == "1.0"
    assert result[0]["NetMHCstab allele"] == "HLA-A02:01 (distance: 0.0)"
    assert result[1]["Half Life"] == "NA"
    # Each run has a temporary directory of its own below the base path
    workers = tmp_path / "base" / "tmp" / "netMHCstabpan" / "workers"
    assert (home / "tmpdir").read_text().startswith(str(workers))
