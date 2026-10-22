"""Tests for the params and the files of the MELT step part."""

from types import SimpleNamespace

from snakemake.io import is_flagged

from snappy_pipeline.workflows.common.melt import MeltStepPart


def _params(action, step="sv_calling_targeted"):
    # The params only read the MELT config section and the step name.
    part = object.__new__(MeltStepPart)
    part.config = SimpleNamespace(
        melt=SimpleNamespace(
            me_refs_infix="1KGP_Hg19",
            exclude_contigs=["hs37d5", "NC_007605"],
            read_length=150,
            max_no_call_percent=100,
        )
    )
    part.parent = SimpleNamespace(name=step)
    return getattr(part, f"_get_params_{action}")(SimpleNamespace(me_type="ALU"))


def test_indiv_analysis_gets_exome_for_targeted_data():
    assert _params("indiv_analysis") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "read_length": 150,
        "exclude_contigs": ["hs37d5", "NC_007605"],
        "exome": True,
    }


def test_indiv_analysis_gets_no_exome_for_wgs():
    assert "exome" not in _params("indiv_analysis", step="sv_calling_wgs")


def test_genotype_takes_neither_read_length_nor_no_call_percent():
    assert _params("genotype") == {"me_refs_infix": "1KGP_Hg19", "me_type": "ALU"}


def test_group_analysis_gets_the_read_length():
    assert _params("group_analysis") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "read_length": 150,
    }


def test_make_vcf_gets_the_no_call_percent():
    assert _params("make_vcf") == {
        "me_refs_infix": "1KGP_Hg19",
        "me_type": "ALU",
        "max_no_call_percent": 100,
    }


def _part():
    # The input functions only read the MELT files, the pedigrees and the upstream paths.
    part = object.__new__(MeltStepPart)
    part.config = SimpleNamespace(melt=SimpleNamespace(jar_file="melt.jar", me_refs_path="refs"))
    part.parent = SimpleNamespace(
        pedigrees=lambda: {"P1": ["P1", "F1"]},
        reference_inputs=lambda: {},
        get_upstream_paths=lambda field: SimpleNamespace(bed="genes.bed"),
    )
    return part


def test_indiv_analysis_writes_one_directory_per_library():
    outputs = _part()._get_output_files_indiv_analysis()
    assert outputs == {"out": "work/indiv_analysis.{library_name}.{me_type}/out"}
    assert is_flagged(outputs["out"], "directory")


def test_group_analysis_reads_the_directories_of_all_pedigree_members():
    wildcards = SimpleNamespace(index_library_name="P1", me_type="ALU")
    inputs = _part()._get_input_files_group_analysis(wildcards)
    assert inputs["indiv_analysis"] == [
        "work/indiv_analysis.P1.ALU/out",
        "work/indiv_analysis.F1.ALU/out",
    ]


def test_genotype_and_make_vcf_read_the_sites_of_group_analysis():
    part = _part()
    part.parent.pedigrees = lambda: {"P1": ["P1"]}
    wildcards = SimpleNamespace(index_library_name="P1", library_name="P1", me_type="ALU")
    pre_geno = "work/group_analysis.P1.ALU/out/ALU.pre_geno.tsv"
    assert part._get_output_files_group_analysis()["pre_geno"] == (
        "work/group_analysis.{index_library_name}.{me_type}/out/{me_type}.pre_geno.tsv"
    )
    assert part._get_input_files_genotype(wildcards)["pre_geno"] == pre_geno
    assert part._get_input_files_make_vcf(wildcards)["pre_geno"] == pre_geno


def test_group_analysis_and_make_vcf_write_no_marker():
    part = _part()
    assert "done" not in part._get_output_files_group_analysis()
    assert "done" not in part._get_output_files_make_vcf()
