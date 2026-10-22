# -*- coding: utf-8 -*-
"""Tests for the step parts of ``somatic_wgs_cnv_calling``."""

from types import SimpleNamespace

import pandas as pd
import pytest
from snakemake.iocontainers import Wildcards

from snappy_pipeline.workflows.somatic_wgs_cnv_calling import (
    CnvettiSomaticWgsStepPart,
    CnvkitSomaticWgsStepPart,
    ControlFreecSomaticWgsStepPart,
    SomaticWgsCnvCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_wgs_cnv_calling.model import (
    CnvkitWgs,
    ControlFreec,
    SomaticWgsCnvCalling,
    Tool,
)

WILDCARDS = Wildcards(fromdict={"library_name": "T1-DNA1-WGS1"})
TUMOR_WILDCARDS = Wildcards(fromdict={"tumor_library": "T1-DNA1-WGS1"})


class _Parent:
    """The upstream files of a task with a cnvkit panel of normals."""

    def get_upstream_paths(self, field, library_name=None):
        if field == "panel_of_normals":
            prefix = "pon/work/panel_of_normals/out/panel_of_normals"
            return SimpleNamespace(
                panel_of_normals=prefix + ".cnn",
                target=prefix + ".target.bed",
                antitarget=prefix + ".antitarget.bed",
            )
        assert field == "alignments"
        return SimpleNamespace(bam=f"{library_name}.bam", bai=f"{library_name}.bam.bai")

    def reference_inputs(self):
        return {"reference": "genome.fa", "reference_fai": "genome.fa.fai"}


def _cnvkit_part():
    # The input functions only read the upstream paths, so skip the full step setup.
    part = object.__new__(CnvkitSomaticWgsStepPart)
    part.parent = _Parent()
    return part


def test_cnvkit_postprocess_gets_the_library_name():
    params = _cnvkit_part()._get_params_postprocess(WILDCARDS)

    assert params == {"library_name": "T1-DNA1-WGS1"}


def test_cnvkit_coverage_reads_the_regions_of_the_panel_of_normals():
    inputs = _cnvkit_part()._get_input_files_coverage(WILDCARDS)

    assert inputs == {
        "bam": "T1-DNA1-WGS1.bam",
        "bai": "T1-DNA1-WGS1.bam.bai",
        "reference": "genome.fa",
        "reference_fai": "genome.fa.fai",
        "target": "pon/work/panel_of_normals/out/panel_of_normals.target.bed",
        "antitarget": "pon/work/panel_of_normals/out/panel_of_normals.antitarget.bed",
    }


def test_cnvkit_fix_reads_the_panel_of_normals():
    inputs = _cnvkit_part()._get_input_files_fix(WILDCARDS)

    assert inputs["ref"] == "pon/work/panel_of_normals/out/panel_of_normals.cnn"


def test_cnvkit_postprocess_writes_the_dnacopy_segments_from_segments_and_calls():
    part = _cnvkit_part()

    assert set(part._get_input_files_postprocess(WILDCARDS)) == {"segment", "call"}
    assert part._get_output_files_postprocess()["final"].endswith(".dnacopy.seg")


def test_cnvkit_declares_the_files_its_wrappers_write():
    part = _cnvkit_part()

    assert not [key for key in part._get_output_files_plot() if key.startswith("heatmap")]
    assert set(part._get_output_files_export()) == {
        key + md5 for key in ("bed", "bed_tbi", "seg", "vcf", "vcf_tbi") for md5 in ("", "_md5")
    }


def test_cnvkit_fix_has_no_edge_correction_for_whole_genomes():
    assert CnvkitWgs().fix.edge_correction is False


def test_cnvkit_needs_a_panel_of_normals():
    with pytest.raises(ValueError, match="depends_on.panel_of_normals"):
        SomaticWgsCnvCalling(
            tool="cnvkit",
            cnvkit={},
            depends_on={"alignments": "mapping", "reference": "genome"},
        )


def test_cnvetti_segment_declares_the_files_its_wrapper_writes():
    outputs = object.__new__(CnvettiSomaticWgsStepPart)._get_output_files_segment()

    assert outputs["segments_bcf"] == "work/{tumor_library}/out/{tumor_library}.bcf"
    assert outputs["windows_bcf"] == "work/{tumor_library}/out/{tumor_library}.windows.bcf"
    # Downstream tasks read the segments
    paths = SomaticWgsCnvCallingWorkflow.get_output_paths(
        SimpleNamespace(tool=Tool.cnvetti), library_name="{tumor_library}"
    )
    assert (paths["bcf"], paths["bcf_csi"]) == (
        outputs["segments_bcf"],
        outputs["segments_bcf_csi"],
    )


def test_control_freec_plot_declares_the_plots_of_cnvkit_plot():
    # The plot rule runs the cnvkit/plot wrapper, which writes a diagram and a scatter plot
    outputs = object.__new__(ControlFreecSomaticWgsStepPart)._get_output_files_plot()

    assert set(outputs) == {"diagram", "diagram_md5", "scatter", "scatter_md5"}


def _control_freec_part(sex, **cfg):
    # The params only read the Control-FREEC config section and the library table.
    part = object.__new__(ControlFreecSomaticWgsStepPart)
    part.config = SimpleNamespace(control_freec=ControlFreec(convert={}, **cfg))
    table = pd.DataFrame({"library_name": ["T1-DNA1-WGS1"], "sex": [sex]})
    part.parent = SimpleNamespace(
        get_upstream_paths=lambda field: SimpleNamespace(mappability="mappability.gem"),
        build_library_dataframe=lambda: table,
    )
    return part


def test_control_freec_run_passes_the_settings_of_the_ini_file():
    part = _control_freec_part("male", ploidy=3)

    assert part._get_params_run(TUMOR_WILDCARDS) == {
        "path_mappability_enabled": True,
        "window_size": -1,
        "ploidy": 3,
        "break_point_threshold": 0.8,
        "coefficient_of_variation": 0.05,
        "contamination_adjustment": True,
        "contamination": 0.4,
        "minimal_subclone_presence": 0.2,
        "unique_match": True,
        "force_gc_content_normalization": 0,
        "min_cna_length": 2,
        "sex": "XY",
    }


@pytest.mark.parametrize(
    "sex, unknown_sex, freec_sex",
    [
        ("male", "XX", "XY"),
        ("female", "XY", "XX"),
        ("unknown", "XY", "XY"),
        ("unknown", "XX", "XX"),
    ],
)
def test_control_freec_takes_the_sex_from_the_sample_sheet(sex, unknown_sex, freec_sex):
    part = _control_freec_part(sex, unknown_sex=unknown_sex)

    assert part._get_params_run(TUMOR_WILDCARDS)["sex"] == freec_sex


def test_control_freec_transform_passes_the_amplification_threshold():
    part = object.__new__(ControlFreecSomaticWgsStepPart)
    part.config = SimpleNamespace(control_freec=ControlFreec(convert={"amplification": 8}))
    wildcards = Wildcards(fromdict={"tumor_library": "T1-DNA1-WGS1"})

    assert part._get_params_transform(wildcards)["amplification"] == 8
    assert ControlFreec(convert={}).convert.amplification == 6
