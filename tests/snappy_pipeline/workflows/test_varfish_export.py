"""Tests for the VarFish export: the files of each pedigree, the settings of the mehari
commands, and the BAM QC table."""

import gzip
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from snakemake.iocontainers import Wildcards

from snappy_pipeline.workflows.varfish_export import MehariStepPart, VarfishExportWorkflow
from snappy_pipeline.workflows.varfish_export.model import VarfishExport
from snappy_wrappers.utils import mehari_options

WRAPPERS = Path(__file__).resolve().parents[3] / "snappy_wrappers" / "wrappers"

DEPENDS_ON = {
    "alignments": "mapping",
    "variants": "calling",
    "reference": "genome",
    "database": "db",
}


def _upstream(field, **kwargs):
    name = kwargs.get("library_name", "")
    return {
        "variants": SimpleNamespace(vcf=f"calls/{name}.vcf.gz"),
        "database": SimpleNamespace(
            hgnc_tsv="db/hgnc.tsv",
            transcript_db=["db/a.bin.zst", "db/b.bin.zst"],
            clinvar_db="db/clinvar",
            frequency_db="db/freqs",
        ),
        "alignments": SimpleNamespace(
            bamstats=f"{name}.bam.bamstats.txt", idxstats=f"{name}.bam.idxstats.txt"
        ),
        "qc": SimpleNamespace(alfred_qc=f"{name}.target_coverage.json.gz"),
        "exons": SimpleNamespace(bed="exons.bed"),
    }[field]


def _part(**settings):
    config = VarfishExport(depends_on={**DEPENDS_ON, **settings.pop("depends_on", {})}, **settings)
    part = object.__new__(MehariStepPart)
    part.config = config
    part.parent = SimpleNamespace(
        config=config,
        get_upstream_paths=_upstream,
        reference_inputs=lambda: {"reference": "ref.fa"},
        pedigree_members=lambda index, include_background: pd.DataFrame(
            {"library_name": ["child-DNA1", "father-DNA1", "mother-DNA1"]}
        ),
        coverage_qc_libraries={"child-DNA1", "mother-DNA1"},
    )
    return part


def test_the_export_writes_one_file_set_per_pedigree_index():
    part = _part()

    assert part._get_output_files_annotate_seqvars() == {
        "ped": "work/{library_name}/out/{library_name}.ped",
        "ped_md5": "work/{library_name}/out/{library_name}.ped.md5",
        "gts": "work/{library_name}/out/{library_name}.seqvars.gts.tsv.gz",
        "gts_md5": "work/{library_name}/out/{library_name}.seqvars.gts.tsv.gz.md5",
    }
    assert (
        part._get_output_files_bam_qc()["bam_qc"]
        == "work/{library_name}/out/{library_name}.bam_qc.tsv.gz"
    )
    assert part.result_keys["annotate_strucvars"] == ("gts", "gts_md5")
    assert VarfishExportWorkflow.get_output_paths(None, library_name="L") == {
        "ped": "work/L/out/L.ped",
        "gts": "work/L/out/L.seqvars.gts.tsv.gz",
    }


def test_the_seqvars_annotation_reads_the_databases_and_passes_the_varfish_defaults():
    part = _part(depends_on={"exons": "exons"})
    wildcards = Wildcards(fromdict={"library_name": "child-DNA1"})

    inputs = part._get_input_files_annotate_seqvars(wildcards)
    params = part._get_params_annotate_seqvars(wildcards)

    assert inputs["ped"] == "work/write_pedigree.child-DNA1/out/child-DNA1.ped"
    assert inputs["transcript_db"] == ["db/a.bin.zst", "db/b.bin.zst"]
    assert inputs["exon_bed"] == "exons.bed"
    config = params["config"]
    assert (config["pick"], config["pick_order"], config["keep_intergenic"]) == (
        "gene",
        ["mane-select", "mane-plus-clinical", "length"],
        True,
    )
    assert "extra_arguments" not in config
    assert mehari_options(config) == [
        "--threads 1",
        "--assembly GRCh38",
        "--pick gene",
        "--pick-order mane-select,mane-plus-clinical,length",
        "--keep-intergenic",
    ]


def test_the_strucvars_export_gets_the_assembly_the_seed_and_its_arguments():
    part = _part(
        mehari={"assembly": "GRCh37"},
        strucvars={"rng_seed": 7, "extra_arguments": ["--slack-bnd 100"]},
    )

    assert part._get_params_annotate_strucvars(SimpleNamespace()) == {
        "assembly": "GRCh37",
        "rng_seed": 7,
        "extra_arguments": ["--slack-bnd 100"],
    }


def test_the_bam_qc_reads_the_samtools_reports_of_the_pedigree_and_the_coverage_of_the_kits():
    inputs = _part()._get_input_files_bam_qc(Wildcards(fromdict={"library_name": "child-DNA1"}))

    assert inputs["bamstats"] == [
        f"{lib}.bam.bamstats.txt" for lib in ("child-DNA1", "father-DNA1", "mother-DNA1")
    ]
    assert inputs["alfred_qc"] == [
        "child-DNA1.target_coverage.json.gz",
        "mother-DNA1.target_coverage.json.gz",
    ]
    assert "flagstats" not in inputs


def test_mehari_options_leave_unset_fields_to_mehari_and_repeat_plugins():
    config = {
        "pick": None,
        "pick_order": [],
        "flag_pick": False,
        "plugins": ["nmd", "aa-sequence-alt"],
        "hgvs_parentheses": "never",
    }

    assert mehari_options(config) == [
        "--plugin nmd",
        "--plugin aa-sequence-alt",
        "--hgvs-parentheses never",
    ]


def _bam_qc_module():
    spec = importlib.util.spec_from_file_location(
        "bam_qc_wrapper", WRAPPERS / "mehari" / "bam_qc" / "wrapper.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_bam_qc_table_holds_the_stats_of_every_library_as_json(tmp_path):
    bam_qc = _bam_qc_module()
    (tmp_path / "L1.bam.bamstats.txt").write_text(
        "# comment\nSN\traw total sequences:\t1000\nSN\taverage length:\t150.5\nSN\tis sorted:\t1\t# yes\nFFQ\t1\t2\n"
    )
    (tmp_path / "L1.bam.idxstats.txt").write_text("1\t249250621\t900\t10\n*\t0\t0\t90\n")
    alfred = {
        "samples": [
            {
                "id": "L1",
                "summary": {
                    "data": {"columns": ["#AlignedBasesInBed", "#TotalBedBp"], "rows": [[300, 100]]}
                },
                "readGroups": [
                    {
                        "metrics": [
                            {
                                "id": "targetCoverage",
                                "x": {"data": [{"values": [0, 5, 10, 15]}]},
                                "y": {"data": [{"values": [1.0, 0.9, 0.5, 0.25]}]},
                            }
                        ]
                    }
                ],
            }
        ]
    }
    with gzip.open(tmp_path / "L1.target_coverage.json.gz", "wt") as f:
        json.dump(alfred, f)

    table = bam_qc.bam_qc_table(
        [str(tmp_path / "L1.bam.bamstats.txt")],
        [str(tmp_path / "L1.bam.idxstats.txt")],
        [str(tmp_path / "L1.target_coverage.json.gz")],
    )

    header, row = table.rstrip("\n").split("\n")
    assert header == "case_id\tset_id\tbam_stats"
    case_id, set_id, bam_stats = row.split("\t")
    assert (case_id, set_id) == (".", ".")
    assert json.loads(bam_stats.replace('"""', '"')) == {
        "L1": {
            "bamstats": {"raw total sequences": 1000, "average length": 150.5, "is sorted": 1},
            "idxstats": {"1": {"mapped": 900, "unmapped": 10}, "*": {"mapped": 0, "unmapped": 90}},
            "summary": {"mean coverage": 3.0, "total target size": 100},
            "min_cov_target": {"5": 90.0, "15": 25.0},
        }
    }


def test_the_bam_qc_without_coverage_has_only_the_samtools_stats(tmp_path):
    bam_qc = _bam_qc_module()
    (tmp_path / "L2.bam.bamstats.txt").write_text("SN\treads mapped:\t5\n")
    (tmp_path / "L2.bam.idxstats.txt").write_text("1\t10\t5\t0\n")

    table = bam_qc.bam_qc_table(
        [str(tmp_path / "L2.bam.bamstats.txt")], [str(tmp_path / "L2.bam.idxstats.txt")], []
    )

    assert json.loads(table.split("\n")[1].split("\t")[2].replace('"""', '"')) == {
        "L2": {"bamstats": {"reads mapped": 5}, "idxstats": {"1": {"mapped": 5, "unmapped": 0}}}
    }
