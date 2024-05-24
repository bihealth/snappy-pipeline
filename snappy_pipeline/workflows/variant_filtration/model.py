import re
from typing import Annotated, Self

from pydantic import Field, AfterValidator, model_validator

from snappy_pipeline.models import SnappyStepModel, SnappyModel


class Threshold(SnappyModel):
    min_gq: int
    min_dp_het: int
    min_dp_hom: int
    include_expressions: list[str]


class Frequencies(SnappyModel):
    af_dominant: float = 0.001
    """AF (allele frequency) values"""

    af_recessive: float = 0.01
    """AF (allele frequency) values"""

    ac_dominant: int = 3
    """AC (allele count in gnomAD) values"""


class ScoreThreshold(SnappyModel):
    require_coding: bool
    require_gerpp_gt2: bool
    min_cadd: int | None


def check_combination(s: str) -> str:
    """
    A very simple validator that checks if the string has 5 dots.
    The actual validation is done in the model validator, because the sets of valid pattern strings
    can only be known at that point.
    """
    assert s.count(".") == 5
    return s


FilterCombination = Annotated[str, AfterValidator(check_combination)]


class VariantFiltration(SnappyStepModel):
    path_variant_annotation: Annotated[str, Field(examples=["../variant_annotation"])]

    tools_ngs_mapping: list[str] = []
    """defaults to ngs_mapping tool"""

    tools_variant_calling: list[str] = []
    """defaults to variant_annotation tool"""

    thresholds: dict[str, Threshold] = {
        "conservative": Threshold(
            **dict(
                min_gq=40,
                min_dp_het=10,
                min_dp_hom=5,
                include_expressions=["'MEDGEN_COHORT_INCONSISTENT_AC=0'"],
            )
        ),
        "relaxed": Threshold(
            **dict(
                min_gq=20,
                min_dp_het=6,
                min_dp_hom=3,
                include_expressions=["'MEDGEN_COHORT_INCONSISTENT_AC=0'"],
            )
        ),
    }
    """quality filter sets, "keep_all" implicitly defined"""

    frequencies: Frequencies = Frequencies()

    region_beds: Annotated[
        dict[str, str],
        Field(
            examples=[
                {
                    "all_tads": "/fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed",
                    "all_genes": "/fast/projects/medgen_genomes/static_data/GRCh37/gene_bed/ENSEMBL_v75.bed.gz",
                    "limb_tads": "/fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed",
                    "lifted_enhancers": "/fast/projects/medgen_genomes/static_data/GRCh37/all_but_onlyMB.bed",
                    "vista_enhancers": "/fast/projects/medgen_genomes/static_data/GRCh37/vista_limb_enhancers.bed",
                }
            ]
        ),
    ] = {}
    """regions to filter to, "whole_genome" implicitly defined"""

    score_thresholds: dict[str, ScoreThreshold] = {
        "coding": ScoreThreshold(
            **dict(require_coding=True, require_gerpp_gt2=False, min_cadd=None)
        ),
        "conservative": ScoreThreshold(
            **dict(require_coding=False, require_gerpp_gt2=False, min_cadd=0)
        ),
        "conserved": ScoreThreshold(
            **dict(require_coding=False, require_gerpp_gt2=True, min_cadd=None)
        ),
    }

    filter_combinations: list[FilterCombination] = [
        "conservative.de_novo.dominant_freq.lifted_enhancers.all_scores.passthrough",
        "conservative.de_novo.dominant_freq.lifted_enhancers.conserved.passthrough",
        "conservative.de_novo.dominant_freq.limb_tads.all_scores.passthrough",
        "conservative.de_novo.dominant_freq.limb_tads.coding.passthrough",
        "conservative.de_novo.dominant_freq.limb_tads.conserved.passthrough",
        "conservative.de_novo.dominant_freq.vista_enhancers.all_scores.passthrough",
        "conservative.de_novo.dominant_freq.vista_enhancers.conserved.passthrough",
        "conservative.de_novo.dominant_freq.whole_genome.all_scores.passthrough",
        "conservative.de_novo.dominant_freq.whole_genome.coding.passthrough",
        "conservative.de_novo.dominant_freq.whole_genome.conserved.passthrough",
        "conservative.dominant.dominant_freq.lifted_enhancers.all_scores.passthrough",
        "conservative.dominant.dominant_freq.lifted_enhancers.conserved.passthrough",
        "conservative.dominant.dominant_freq.limb_tads.all_scores.passthrough",
        "conservative.dominant.dominant_freq.limb_tads.coding.passthrough",
        "conservative.dominant.dominant_freq.limb_tads.conserved.passthrough",
        "conservative.dominant.dominant_freq.vista_enhancers.all_scores.passthrough",
        "conservative.dominant.dominant_freq.vista_enhancers.conserved.passthrough",
        "conservative.dominant.dominant_freq.whole_genome.all_scores.passthrough",
        "conservative.dominant.dominant_freq.whole_genome.coding.passthrough",
        "conservative.dominant.dominant_freq.whole_genome.conserved.passthrough",
        "conservative.dominant.recessive_freq.lifted_enhancers.all_scores.intervals500",
        "conservative.dominant.recessive_freq.lifted_enhancers.conserved.intervals500",
        "conservative.dominant.recessive_freq.lifted_enhancers.conserved.tads",
        "conservative.dominant.recessive_freq.limb_tads.all_scores.intervals500",
        "conservative.dominant.recessive_freq.limb_tads.coding.gene",
        "conservative.dominant.recessive_freq.limb_tads.conserved.intervals500",
        "conservative.dominant.recessive_freq.limb_tads.conserved.tads",
        "conservative.dominant.recessive_freq.vista_enhancers.all_scores.intervals500",
        "conservative.dominant.recessive_freq.vista_enhancers.conserved.intervals500",
        "conservative.dominant.recessive_freq.vista_enhancers.conserved.tads",
        "conservative.dominant.recessive_freq.whole_genome.all_scores.intervals500",
        "conservative.dominant.recessive_freq.whole_genome.coding.gene",
        "conservative.dominant.recessive_freq.whole_genome.conserved.intervals500",
        "conservative.dominant.recessive_freq.whole_genome.conserved.tads",
        "conservative.recessive_hom.recessive_freq.lifted_enhancers.all_scores.passthrough",
        "conservative.recessive_hom.recessive_freq.lifted_enhancers.conserved.passthrough",
        "conservative.recessive_hom.recessive_freq.limb_tads.all_scores.passthrough",
        "conservative.recessive_hom.recessive_freq.limb_tads.coding.passthrough",
        "conservative.recessive_hom.recessive_freq.limb_tads.conserved.passthrough",
        "conservative.recessive_hom.recessive_freq.vista_enhancers.all_scores.passthrough",
        "conservative.recessive_hom.recessive_freq.vista_enhancers.conserved.passthrough",
        "conservative.recessive_hom.recessive_freq.whole_genome.all_scores.passthrough",
        "conservative.recessive_hom.recessive_freq.whole_genome.coding.passthrough",
        "conservative.recessive_hom.recessive_freq.whole_genome.conserved.passthrough",
        # The following are for input to variant_combination.
        "conservative.dominant.recessive_freq.whole_genome.coding.passthrough",
        "conservative.dominant.recessive_freq.whole_genome.conserved.passthrough",
    ]
    """dot-separated {thresholds}.{inherit}.{freq}.{region}.{score}.{het_comp}"""

    @model_validator(mode="after")
    def ensure_filter_combinations_are_valid(self: Self) -> Self:
        thresholds: set[str] = set(self.thresholds.keys())
        inherit: set[str] = {"de_novo", "dominant", "recessive_hom"}
        freq: set[str] = {"dominant_freq", "recessive_freq"}
        region: set[str] = set(self.region_beds.keys()) | {"whole_genome"}
        score: set[str] = set(self.score_thresholds.keys()) | {"all_scores"}
        het_comp: set[str] = {"passthrough", "intervals500", "tads", "gene"}
        pattern: str = r".".join(
            f'({"|".join(p)})' for p in [thresholds, inherit, freq, region, score, het_comp]
        )
        pattern: re.Pattern[str] = re.compile(pattern)
        for combination in self.filter_combinations:
            assert pattern.fullmatch(combination) is not None

        return self
