"""Tests for the helpers of the wrappers: ``shards``, ``java``, ``bwa_mem``"""

import ast
import pathlib
from types import SimpleNamespace

import pytest
from snakemake.utils import format as snakemake_format

from snappy_wrappers.tools import bwa_mem, gcnv_ploidy, java, shards


def test_intervals_exclude_each_ignored_contig_and_pad():
    text = shards.intervals("ref.fa", "ref.fa.fai", 1000, 50, ["^chrUn", "_decoy$"])
    assert "--bin-length 1000" in text
    assert "for ignore_contig in ^chrUn _decoy$; do" in text
    assert 'awk "(\\$1 ~ /$ignore_contig/) { printf(' in text  # braces are single
    assert "            ref.fa.fai\n" in text
    assert "-v PADDING=50" in text
    assert text.endswith("mkdir -p $TMPDIR/shards-output")


def test_run_shards_passes_job_number_and_interval():
    text = shards.run_shards(8)
    assert "num_threads=8\n" in text
    assert text.endswith("| parallel --plain -j $num_threads 'run-shard {#} {}'")


@pytest.mark.parametrize("check_ref, line", [("", None), ("ws", "        -c ws \\\n")])
def test_gather_normalizes_against_the_reference(check_ref, line):
    text = shards.gather("ref.fa", "out.vcf.gz", check_ref=check_ref)
    assert ("-c " in text) == bool(line)
    assert (line or "") + "        -f ref.fa \\\n" in text
    assert text.endswith("        -o out.vcf.gz\ntabix out.vcf.gz")


@pytest.mark.parametrize(
    "mem_mb, shards_, heap", [(4000, 1, 3000), (12345, 1, 9259), (22000, 4, 4125), (12345, 7, 1323)]
)
def test_java_heap_is_three_quarters_of_the_memory_per_shard(mem_mb, shards_, heap):
    snakemake = SimpleNamespace(resources=SimpleNamespace(mem_mb=mem_mb))
    assert java.java_heap_mb(snakemake, shards=shards_) == heap


def test_bwa_mem_formats_lists_as_snakemake_does():
    text = bwa_mem.read_functions(["l1.fq.gz", "l2.fq.gz"], "", "rg1 rg2", True, 4, False)
    assert "declare -a reads_left=(l1.fq.gz l2.fq.gz)\ndeclare -a reads_right=()\n" in text
    assert 'if [[ "True" == "True" ]]; then\n        trimadap-mt -p 4' in text
    assert "${#reads_right[@]}" in text


def test_bwa_mem_pipeline_calls_the_mapper():
    output = SimpleNamespace(
        bam="x.bam",
        bam_bai="x.bam.bai",
        report_bamstats_txt="s",
        report_flagstats_txt="f",
        report_idxstats_txt="i",
    )
    text = bwa_mem.mapping_pipeline("run_bwa_mem3", output, 2, "4G", 3)
    assert text.count("| run_bwa_mem3 ") == 2
    assert "samtools idxstats x.bam > i\n" in text
    assert text.endswith("check_input_line_counts")


def test_bwa_mem_command_renders_through_the_shell_wrapper():
    """The final shell texts are pasted into the template, which ShellWrapper formats once more:
    their braces must survive that pass, and the placeholders it fills must be the only ones left."""
    args = {
        "read_groups": ["rg1"],
        "trim_adapters": True,
        "num_threads_trimming": 2,
        "mask_duplicates": True,
        "split_as_secondary": False,
        "num_threads_align": 4,
        "bases_per_batch": 100,
        "num_threads_bam_view": 1,
        "memory_bam_sort": "1G",
        "num_threads_bam_sort": 1,
        "sample_name": "S1",
        "platform": "ILLUMINA",
    }
    reads = {"reads_left": ["l1.fq.gz"], "reads_right": ["r1.fq.gz"]}
    snakemake = SimpleNamespace(
        params=SimpleNamespace(args=args),
        input=SimpleNamespace(index=["ref.fa.amb"], get=lambda key, default=None: reads.get(key)),
        output=SimpleNamespace(
            bam="x.bam",
            bam_bai="x.bam.bai",
            report_bamstats_txt="s",
            report_flagstats_txt="f",
            report_idxstats_txt="i",
        ),
    )

    script = snakemake_format(
        bwa_mem.command(snakemake, "bwa-mem3"), snakemake=snakemake, extra_arguments="-Y"
    )

    assert "input_files() {\n    set -x\n    i=$1\n" in script
    assert "${#reads_right[@]}" in script and "{{" not in script and "}}" not in script
    assert "run_bwa_mem3() {\n" in script and "| run_bwa_mem3 0 \\\n" in script
    assert "    bwa-mem3 mem \\\n        ref.fa \\\n" in script
    assert "\\tSM:S1\\tPL:ILLUMINA" in script and "        -Y \\\n" in script


@pytest.mark.parametrize("module", [shards, java, bwa_mem, gcnv_ploidy])
def test_modules_run_on_the_oldest_python_that_snakemake_uses_for_wrappers(module):
    """The wrappers that import them run with the Python of their environment."""
    from snakemake.script import MIN_PY_VERSION

    ast.parse(pathlib.Path(module.__file__).read_text(), feature_version=MIN_PY_VERSION)
