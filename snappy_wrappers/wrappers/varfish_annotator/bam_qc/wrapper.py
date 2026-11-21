import gzip
import json
import shlex
import tempfile

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

COV_THRESHOLDS = [5, 15, 20, 30, 50, 100, 200, 300, 400, 500, 1000]

# The name of each library in the table, in the order of the samtools reports
sample_names = " ".join(shlex.quote(str(name)) for name in snakemake.params.args["sample_names"])

with tempfile.TemporaryDirectory() as tmpdir:
    # Extract the coverage data from the alfred qc files.
    result = {}
    for alfred_qc in snakemake.input.alfred_qc:
        with gzip.open(alfred_qc, "rt") as inputf:
            data = json.load(inputf)

        for sample in data["samples"]:
            sample_name = sample["id"]

            cols = sample["summary"]["data"]["columns"]
            rows = sample["summary"]["data"]["rows"]
            data0 = dict(zip(cols, rows[0]))
            total_aligned_in_target = data0["#AlignedBasesInBed"]
            total_target_size = data0["#TotalBedBp"]
            for metric in sample["readGroups"][0]["metrics"]:
                if metric["id"] == "targetCoverage":
                    target_cov_thresholds = {cov: 0 for cov in COV_THRESHOLDS}
                    idx_to_cov = {}
                    for idx, cov in enumerate(metric["x"]["data"][0]["values"]):
                        if cov in COV_THRESHOLDS:
                            idx_to_cov[idx] = cov
                    cov_to_frac = {}
                    for idx, frac in enumerate(metric["y"]["data"][0]["values"]):
                        if idx in idx_to_cov:
                            cov_to_frac[idx_to_cov[idx]] = frac * 100
                    break
            else:
                raise RuntimeError("Could not find targetCoverage metric")

            with open(f"{tmpdir}/alfred_qc.{sample_name}", "wt") as outputf:
                json.dump(
                    {
                        sample_name: {
                            "summary": {
                                "mean coverage": total_aligned_in_target / total_target_size,
                                "total target size": total_target_size,
                            },
                            "min_cov_target": cov_to_frac,
                        }
                    },
                    outputf,
                )

    # Actually run the script.
    ShellWrapper(snakemake).run(
        r"""
samples=({sample_names})
bamstats=({snakemake.input.bamstats})
idxstats=({snakemake.input.idxstats})

# Extract "samtools stats" outputs into JSON
for i in "${{!samples[@]}}"; do
    sample=${{samples[$i]}}

    grep '^SN' ${{bamstats[$i]}} \
        | cut -f 2-3 \
        | sed -e 's/:\t/\t/g' \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ print "{{ \"" sample "\": {{ \"bamstats\": {{"; }}
            {{ if (NR > 1) print ","; print "\"" $1 "\":" $2 }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
            > $TMPDIR/bamstats.$sample
done

# Extract "samtools idxstats" into JSON
for i in "${{!samples[@]}}"; do
    sample=${{samples[$i]}}

    cat ${{idxstats[$i]}} \
        | cut -f 1,3,4 \
        | awk -F $'\t' -v sample=$sample 'BEGIN {{ print "{{ \"" sample "\": {{ \"idxstats\": {{"; }}
            {{ if (NR > 1) print ","; print "\"" $1 "\":{{\"mapped\":" $2 ",\"unmapped\":" $3 "}}" }}
            END {{ print "}} }} }}" }}' \
        | jq '.' \
            > $TMPDIR/idxstats.$sample
done

# Merge JSON files
out_bam_qc={snakemake.output.bam_qc}

echo -en "case_id\tset_id\tbam_stats\n.\t.\t" \
    > ${{out_bam_qc%.gz}}

jq -c --slurp 'reduce .[] as $item ({{}}; . * $item)' {tmpdir}/alfred_qc.* $TMPDIR/* \
    | sed -e 's/"/\"""/g' \
        >> ${{out_bam_qc%.gz}}

# Compress output files and create MD5 sums.
gzip ${{out_bam_qc%.gz}}"""
    )
