"""Tests for the shell texts of the fgumi wrappers: ``fgumi_grouping`` and ``fgumi_consensus``"""

from types import SimpleNamespace

import pytest
from snakemake.utils import format as snakemake_format

from snappy_wrappers.tools import fgumi_consensus, fgumi_grouping

INDEX = [f"idx/reference{ext}" for ext in (".amb", ".ann", ".bwt.2bit.64", ".pac")]
REFERENCE = {
    "reference": "refs/genome.fa",
    "reference_fai": "refs/genome.fa.fai",
    "reference_dict": "refs/genome.dict",
}
ALIGNER = {"num_threads": 16, "bases_per_batch": 150000000, "extra_arguments": ["--fast"]}
GROUP = {"edits": 1, "min_map_q": 1, "allow_unmapped": False, "extra_arguments": []}
FILTER = {
    "max_read_error_rate": 0.025,
    "max_base_error_rate": 0.1,
    "min_base_quality": 40,
    "max_no_call_fraction": 0.2,
    "extra_arguments": ["--min-mean-base-quality", "20"],
}
DUPLEX = {
    "min_reads": [1],
    "min_input_base_quality": 20,
    "filter_min_reads": [3, 2, 2],
    "require_single_strand_agreement": True,
}
SIMPLEX = {"min_reads": 1, "min_input_base_quality": 20, "filter_min_reads": 3}


def _snakemake(args, inputs, outputs):
    """A stand-in for the ``snakemake`` object of a wrapper, as ``ShellWrapper`` formats it."""
    return SimpleNamespace(
        params=SimpleNamespace(args=args),
        input=SimpleNamespace(get=lambda key, default=None: inputs.get(key, default), **inputs),
        output=SimpleNamespace(**outputs),
        log=SimpleNamespace(log="work/L1/log/L1.log"),
        threads=16,
    )


def _script(text, snakemake, extra_arguments=""):
    return snakemake_format(text, snakemake=snakemake, extra_arguments=extra_arguments)


def _grouping(correct, read_groups=("L1",), reads=None):
    reads = reads or {"reads_left": ["L1.R1.fastq.gz"], "reads_right": ["L1.R2.fastq.gz"]}
    args = {
        "sample_name": "L1",
        "read_groups": list(read_groups),
        "read_structures": ["8M+T", "8M+T"],
        "strategy": "paired",
        "num_threads": 4,
        "platform": "ILLUMINA",
        "aligner": ALIGNER,
        "correct": correct,
        "group": GROUP,
    }
    inputs = {**reads, "index": INDEX, **REFERENCE}
    outputs = {
        "grouped_bam": "work/L1/out/L1.grouped.bam",
        "group_family_sizes": "work/L1/report/umi_qc/L1.group.family_sizes.txt",
    }
    if correct is not None:
        inputs["umis"] = "kit/umis.txt"
        outputs["correct_metrics"] = "work/L1/report/umi_qc/L1.correct.txt"
    snakemake = _snakemake(args, inputs, outputs)
    return _script(fgumi_grouping.command(snakemake), snakemake)


def test_grouping_extracts_aligns_sorts_and_groups_in_one_pipeline():
    script = _grouping(correct=None)
    assert "declare -a reads_left=(L1.R1.fastq.gz)" in script
    assert "--read-structures 8M+T 8M+T" in script
    assert "--read-group-id ${read_groups[$i]}" in script and "--sample L1" in script
    assert "mv $TMPDIR/extracted.bam $TMPDIR/unmapped.bam" in script
    assert "fgumi correct" not in script
    assert (
        "bwa-mem3 mem \\\n        --bam=0 \\\n        -p \\\n        -Y \\\n        -t 16" in script
    )
    assert "-K 150000000 \\\n        --fast \\\n        idx/reference \\\n        - \\" in script
    assert "--unmapped $TMPDIR/unmapped.bam \\\n        --reference refs/genome.fa" in script
    assert "--order template-coordinate" in script and "--tmp-dir $TMPDIR" in script
    assert "--strategy paired \\\n        --edits 1 \\\n        --min-map-q 1" in script
    assert "--allow-unmapped false \\\n        --realign-window true" in script
    assert "--metrics work/L1/report/umi_qc/L1.group \\" in script
    assert script.endswith("--output work/L1/out/L1.grouped.bam")


def test_grouping_corrects_the_umis_against_the_whitelist():
    script = _grouping(correct={"max_mismatches": 1, "min_distance": 2})
    assert (
        "fgumi correct \\\n    --input $TMPDIR/extracted.bam \\\n    --umi-files kit/umis.txt"
        in script
    )
    assert "--max-mismatches 1 \\\n    --min-distance 2" in script
    assert "--metrics work/L1/report/umi_qc/L1.correct.txt" in script
    assert "--output $TMPDIR/unmapped.bam" in script
    assert "mv $TMPDIR/extracted.bam" not in script


def test_grouping_merges_the_read_groups_under_one_header():
    script = _grouping(
        correct=None,
        read_groups=("L1_L001", "L1_L002"),
        reads={
            "reads_left": ["a.R1.fq.gz", "b.R1.fq.gz"],
            "reads_right": ["a.R2.fq.gz", "b.R2.fq.gz"],
        },
    )
    assert "declare -a read_groups=(L1_L001 L1_L002)" in script
    assert "samtools cat -h $TMPDIR/merged.hdr.bam -o $TMPDIR/extracted.bam" in script
    # Temp files are named by the index of the read group, never by its name
    assert "$TMPDIR/$i.extracted.bam" in script


def test_grouping_reads_the_files_another_task_wrote():
    reads = {"reads_left": ["trimmed/L1.R1.fastq.gz"], "reads_right": ["trimmed/L1.R2.fastq.gz"]}
    args = {
        "sample_name": "L1",
        "read_groups": ["L1"],
        "read_structures": ["8M+T", "8M+T"],
        "strategy": "adjacency",
        "num_threads": 4,
        "platform": "ILLUMINA",
        "aligner": ALIGNER,
        "correct": None,
        "group": GROUP,
        "input": reads,
    }
    snakemake = _snakemake(
        args,
        {"reads": ["tasks/trimming/work/L1/out"], "index": INDEX, **REFERENCE},
        {"grouped_bam": "g.bam", "group_family_sizes": "x.group.family_sizes.txt"},
    )
    script = _script(fgumi_grouping.command(snakemake), snakemake)
    assert "declare -a reads_left=(trimmed/L1.R1.fastq.gz)" in script
    assert "--strategy adjacency" in script


def _consensus(mode, caller, clip=None, extra_arguments=""):
    args = {
        "consensus": mode,
        "num_threads": 4,
        "compression_level": 6,
        "aligner": ALIGNER,
        "caller": caller,
        "filter": FILTER,
        "realign": {"chunk_size": 10000000},
        "clip": clip,
    }
    outputs = {
        "bam": "work/L1/out/L1.bam",
        "consensus_family_sizes": "work/L1/report/umi_qc/L1.consensus.family_sizes.txt",
        "consensus_stats": "work/L1/report/umi_qc/L1.consensus.txt",
        "filter_stats": "work/L1/report/umi_qc/L1.filter.txt",
        "report_bamstats_txt": "work/L1/report/bam_qc/L1.bam.bamstats.txt",
        "report_flagstats_txt": "work/L1/report/bam_qc/L1.bam.flagstats.txt",
        "report_idxstats_txt": "work/L1/report/bam_qc/L1.bam.idxstats.txt",
    }
    if clip is not None:
        outputs["clip_metrics"] = "work/L1/report/umi_qc/L1.clip.txt"
    snakemake = _snakemake(
        args, {"bam": "work/L1/out/L1.grouped.bam", "index": INDEX, **REFERENCE}, outputs
    )
    return _script(fgumi_consensus.command(snakemake), snakemake, extra_arguments)


def test_duplex_consensus_calls_realigns_filters_and_sorts_in_one_pipeline():
    script = _consensus("duplex", DUPLEX, extra_arguments="--trim true")
    assert (
        "fgumi duplex \\\n    --input work/L1/out/L1.grouped.bam \\\n    --min-reads 1 \\" in script
    )
    assert "--min-input-base-quality 20 \\\n    --output-per-base-tags true" in script
    assert "--metrics work/L1/report/umi_qc/L1.consensus \\" in script
    assert "--stats work/L1/report/umi_qc/L1.consensus.txt" in script
    assert "--trim true \\\n    --output /dev/stdout" in script
    # The realigner's bwa-mem3 reads the index next to the linked reference, with -Y and the
    # chunk size as -K; its thread count is literal, since --aligner::threads is rejected with a command
    assert "ln -s $(realpath refs/genome.fa) $TMPDIR/ref/genome.fa\n" in script
    assert (
        "for f in idx/reference.amb idx/reference.ann idx/reference.bwt.2bit.64 idx/reference.pac; do"
        in script
    )
    assert "$TMPDIR/ref/genome.fa${f#idx/reference}" in script
    assert (
        '--aligner::command "bwa-mem3 mem --bam=0 -p -Y -t 16 -K 10000000 --fast {ref} /dev/stdin"'
        in script
    )
    assert "--aligner::chunk-size 10000000" in script
    assert "--aligner::threads" not in script
    commands = [line for line in script.splitlines() if not line.strip().startswith("#")]
    assert not [line for line in commands if "--reverse-per-base-tags" in line]
    assert "--min-reads 3,2,2 \\\n        --max-read-error-rate 0.025" in script
    assert "--max-no-call-fraction 0.2 \\\n        --require-single-strand-agreement true" in script
    assert "--min-mean-base-quality 20 \\\n        --output /dev/stdout" in script
    assert "clip_reads() {\n    cat\n}" in script
    assert "--order coordinate \\\n        --write-index true" in script
    assert "--compression-level 6 \\\n        --output work/L1/out/L1.bam" in script
    assert script.endswith(
        "samtools idxstats work/L1/out/L1.bam > work/L1/report/bam_qc/L1.bam.idxstats.txt"
    )


def test_simplex_consensus_takes_single_min_reads_and_clips():
    clip = {
        "read_one_five_prime": 15,
        "read_one_three_prime": 8,
        "read_two_five_prime": 15,
        "read_two_three_prime": 8,
        "overlapping_reads": False,
        "extra_arguments": [],
    }
    script = _consensus("simplex", SIMPLEX, clip=clip)
    assert (
        "fgumi simplex \\\n    --input work/L1/out/L1.grouped.bam \\\n    --min-reads 1 \\"
        in script
    )
    assert "--min-reads 3 \\\n        --max-read-error-rate 0.025" in script
    assert "--require-single-strand-agreement" not in script
    assert (
        "fgumi clip \\\n        --input /dev/stdin \\\n        --reference refs/genome.fa" in script
    )
    assert "--clipping-mode hard \\\n        --auto-clip-attributes true" in script
    assert "--read-one-five-prime 15 \\\n        --read-one-three-prime 8" in script
    assert "--clip-overlapping-reads false" in script
    assert "--metrics work/L1/report/umi_qc/L1.clip.txt" in script


@pytest.mark.parametrize("module", [fgumi_grouping, fgumi_consensus])
def test_the_shell_texts_remove_only_files_of_the_temp_directory(module):
    texts = [getattr(module, name) for name in dir(module) if name.isupper()]
    for text in texts:
        for line in text.splitlines():
            if line.strip().startswith("rm "):
                assert all(arg.startswith("$TMPDIR/") for arg in line.split()[1:]), line
