import enum
import logging
from typing import Annotated, Self

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import (
    ConfigPath,
    EnumField,
    ReadBy,
    RequiredBy,
    SnappyModel,
    SnappyStepModel,
    ToggleModel,
    ToolModel,
)
from snappy_pipeline.models.calls_report import CallsReport
from snappy_pipeline.models.gatk import GATK
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    CommonVariants,
    Dbsnp,
    DnaAlignments,
    ExpectedPathSchema,
    ExpectedVcf,
    GermlineResource,
    Kits,
    Reference,
)

logger = logging.getLogger(__name__)


class ExpectedCandidates(BaseModel):
    """The candidate variants of a varlociraptor task: a VCF or BCF, read once in order."""

    vcf: str


class BafFileGeneration(ToggleModel, ToolModel):
    min_dp: Annotated[int, Field(ge=1)] = 10
    """minimal DP of variant, must be >=1"""


class BcftoolsStats(ToggleModel, ToolModel):
    pass


class BcftoolsRoh(ToggleModel, ToolModel):
    path_af_file: ConfigPath

    ignore_homref: bool = False

    skip_indels: bool = False

    rec_rate: float = 1e-8


class Tool(enum.StrEnum):
    bcftools_call = "bcftools_call"
    gatk3_hc = "gatk3_hc"
    gatk3_ug = "gatk3_ug"
    gatk4_hc_joint = "gatk4_hc_joint"
    gatk4_hc_gvcf = "gatk4_hc_gvcf"
    mutect2 = "mutect2"
    freebayes = "freebayes"
    varlociraptor = "varlociraptor"


class IgnoreContigs(SnappyModel):
    """The contigs that a caller skips"""

    ignore_contigs: list[str] = ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-"]
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""


class IntervalPadding(SnappyModel):
    """The padding around the windows that a caller processes in parallel"""

    padding: Annotated[int, Field(ge=0)] = 1000
    """Bases added on both sides of each window (the windows then overlap; duplicate calls are merged)"""


class WindowedCaller(ToolModel, IgnoreContigs, IntervalPadding):
    """A caller that processes the windows of the genome in parallel"""

    window_length: int = 10000000
    num_threads: int = 16


class GatkCaller(WindowedCaller):
    """A GATK caller"""

    allow_seq_dict_incompatibility: bool = False


class Assembly(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class BcftoolsCall(WindowedCaller):
    assembly: Annotated[Assembly, EnumField(Assembly)]
    """The genome build: the ploidy of chrX and chrY by sex (``bcftools call --ploidy``)"""

    max_depth: int = 250
    max_indel_depth: int = 250

    annotate: str = "FORMAT/AD,FORMAT/DP"
    """Tags that ``bcftools mpileup`` adds (``--annotate``)"""

    per_sample_mf: bool = True
    """Apply the minimum number of gapped reads per sample (``bcftools mpileup --per-sample-mF``)"""

    redo_baq: bool = True
    """Recalculate BAQ on the fly (``bcftools mpileup --redo-BAQ``)"""

    variants_only: bool = True
    """Output variant sites only (``bcftools call --variants-only``)"""


class Gatk3Hc(GatkCaller):
    sample_ploidy: Annotated[int, Field(ge=1)] = 2
    """Ploidy of the samples (``--sample_ploidy``)"""


class Gatk3Ug(GatkCaller):
    downsample_to_coverage: int = 250

    sample_ploidy: Annotated[int, Field(ge=1)] = 2
    """Ploidy of the samples (``--sample_ploidy``)"""


class Gatk4HcJoint(GatkCaller):
    annotation_groups: list[str] = ["StandardAnnotation", "StandardHCAnnotation"]
    """Annotation groups of ``HaplotypeCaller`` (``-G``)"""


class Gatk4HcGvcf(GatkCaller):
    discover_annotation_groups: list[str] = [
        "StandardAnnotation",
        "AS_StandardAnnotation",
        "StandardHCAnnotation",
        "AlleleSpecificAnnotation",
    ]
    """Annotation groups of ``HaplotypeCaller`` (``-G``)"""

    combine_gvcfs_annotation_groups: list[str] = [
        "StandardAnnotation",
        "AS_StandardAnnotation",
        "StandardHCAnnotation",
        "AlleleSpecificAnnotation",
    ]
    """Annotation groups of ``CombineGVCFs`` (``-G``)"""

    genotype_annotation_groups: list[str] = [
        "StandardAnnotation",
        "AS_StandardAnnotation",
        "StandardHCAnnotation",
    ]
    """Annotation groups of ``GenotypeGVCFs`` (``-G``)"""


class TumorNormalMode(enum.StrEnum):
    AUTOMATIC = "automatic"
    PAIRED = "paired"
    TUMOR_ONLY = "tumor_only"
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class Contamination(ToggleModel, GATK):
    java_options: str = "-Xms4000m -Xmx8000m"
    """Options for the java run-time of ``CalculateContamination``"""

    pileup: GATK = GATK()
    """Parameters for GetPileupSummaries used on the tumor bam (& normalk if present)"""


class Mutect2(GATK, IgnoreContigs):
    num_jobs: Annotated[int, Field(ge=1)] = 24
    """The number of regions the genome is split into for scattered calling"""

    num_cores: Annotated[int, Field(ge=1)] = 2
    """The threads of each scattered calling job"""

    contamination: Contamination
    """Estimation of contamination using GetPileupSummaries & CalculateContamination"""

    filtration: GATK = GATK()
    """Additional arguments for filtration"""

    padding: int = 5000
    """Padding around intervals for scatter/gather"""

    tumor_normal_mode: TumorNormalMode = TumorNormalMode.AUTOMATIC
    """Whether to call variants in paired, tumor_only, or automatic mode."""


class Freebayes(ToolModel, IgnoreContigs):
    """Candidate calling with ``freebayes --pooled-continuous``, jointly over the libraries of an
    entity, in parallel over genome windows"""

    min_alternate_fraction: Annotated[float, Field(gt=0, le=1)] = 0.05
    """Minimum fraction of the observations of a position that support an alternate allele
    (``--min-alternate-fraction``)"""

    min_alternate_count: Annotated[int, Field(ge=1)] = 2
    """Minimum number of observations that support an alternate allele (``--min-alternate-count``)"""

    atomize: bool = True
    """Split complex alleles into their SNVs and indels (``bcftools norm --atomize``).
    Varlociraptor then evaluates each atom instead of realigning a long complex allele,
    which is much cheaper on the alleles freebayes builds in polymorphic regions.
    Annotation and filters see atomic variants.
    Off, adjacent substitutions stay one record and one codon change"""

    num_jobs: Annotated[int, Field(ge=1)] = 24
    """Number of genome windows that are called in parallel"""


class ScenarioKind(enum.StrEnum):
    germline = "germline"
    somatic = "somatic"


class Species(SnappyModel):
    """The priors of the ``germline`` scenario (its ``species`` section)"""

    heterozygosity: Annotated[float, Field(gt=0, lt=1)] = 0.001
    """Expected heterozygosity of the population"""

    germline_mutation_rate: Annotated[float, Field(gt=0, lt=1)] = 1.2e-8
    """Rate of de novo germline mutations per base and generation (``germline-mutation-rate``).
    With the full prior, it weighs ``denovo`` against inheritance from an unobserved parent"""

    ploidy: dict[str, dict[str, int]] = {
        "male": {"all": 2, "X": 1, "Y": 1},
        "female": {"all": 2, "X": 2, "Y": 0},
    }
    """Ploidy per sex and contig, ``all`` for the other contigs, in the contig names of the
    reference (``chrX`` and ``chrY`` on hg38); a library of unknown sex is diploid everywhere"""


class Preprocess(ToolModel):
    """``varlociraptor preprocess variants``"""

    max_depth: Annotated[int, Field(ge=1)] = 200
    """Observations per library and position that varlociraptor keeps (``--max-depth``); more reads are downsampled.
    Raise it where the expected allele frequencies are low, such as subclonal somatic variants or mosaics:
    at 1 percent, two of 200 observations support the variant, and a single one counts as ambiguous.
    Preprocessing and calling take longer with more observations"""


class Varlociraptor(ToolModel):
    """Calling with varlociraptor against a scenario; ``extra_arguments`` go to
    ``varlociraptor call variants``"""

    scenario: Annotated[ScenarioKind, EnumField(ScenarioKind)]
    """The scenario preset: ``germline`` calls the pedigree with mendelian priors, ``somatic``
    each tumor with its matched normal; it sets the tags of the calls, the grouping and the
    libraries taken"""

    template: ConfigPath = ""
    """A YTE template that replaces the YAML of the preset, rendered with ``samples`` (the
    library table rows of the entity's libraries), ``library_name`` (the entity) and
    ``settings`` (``species``, ``purity``)"""

    tumor_normal_mode: TumorNormalMode = TumorNormalMode.AUTOMATIC
    """``somatic`` only: ``paired`` leaves out tumors without a matched normal, ``tumor_only``
    calls every tumor alone, ``automatic`` calls each tumor with its matched normal if it has one"""

    species: Species = Species()
    """``germline`` only: the priors"""

    full_prior: bool = True
    """Pass ``--full-prior`` to ``varlociraptor call variants``.
    The prior then weighs every event, with mendelian inheritance and the mutation rate.
    Without it, the prior only tells absent from present,
    and a child with unobserved parents gets a fixed share of ``denovo``"""

    preprocess: Preprocess = Preprocess()

    num_jobs: Annotated[int, Field(ge=1)] = 24
    """Number of chunks of the candidates that are called in parallel"""

    @model_validator(mode="after")
    def check_scenario_settings(self) -> Self:
        """Reject a changed setting that the chosen scenario does not use (a written-out default
        passes, as ``snappy task add`` writes every field)."""
        unused = {"germline": ("tumor_normal_mode",), "somatic": ("species",)}[self.scenario]
        fields = type(self).model_fields
        changed = [name for name in unused if getattr(self, name) != fields[name].get_default()]
        if changed:
            raise ValueError(f"{', '.join(changed)}: not used by the {self.scenario} scenario")
        return self


class VariantCallingDependsOn(SnappyModel):
    alignments: DnaAlignments

    reference: Reference

    dbsnp: Annotated[
        Dbsnp, RequiredBy(Tool.gatk3_hc, Tool.gatk3_ug, Tool.gatk4_hc_joint, Tool.gatk4_hc_gvcf)
    ] = ""

    panel_of_normals: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"pon", "mutect2"})),
        ExpectedPathSchema(ExpectedVcf),
        ReadBy(Tool.mutect2),
    ] = ""
    """
    The panel VCF that Mutect2 uses to filter recurrent artefacts: a ``panel_of_normals`` task
    with ``tool: mutect2``, or an ``external_data`` task with ``{type: models, tags: [pon,
    mutect2]}`` and the key ``vcf`` (GATK resource bundle: ``1000g_pon.hg38.vcf.gz``)
    """

    germline_resource: Annotated[GermlineResource, ReadBy(Tool.mutect2)] = ""
    """Population allele frequencies (af-only gnomAD VCF)"""

    common_variants: Annotated[
        CommonVariants, RequiredBy(Tool.mutect2, when="mutect2.contamination")
    ] = ""
    """Common germline variants, for the contamination estimate"""

    candidates: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedCandidates),
        RequiredBy(Tool.varlociraptor),
    ] = ""
    """The candidate variants that ``varlociraptor`` evaluates: a ``variant_calling`` task with
    ``tool: freebayes``, an SV calling task (the calls are then tagged ``sv``), or another task
    that produces ``variants`` with the key ``vcf``; with the same ``group_by`` and
    ``library_selection``"""

    kits: Annotated[Kits, ReadBy(when="bcftools_roh"), ReadBy(Tool.freebayes)] = ""
    """Restricts ``bcftools_roh`` and the candidates of ``freebayes`` to the targets of the
    (targeted) libraries' kit"""


class VariantCalling(SnappyStepModel):
    depends_on: VariantCallingDependsOn

    calls_report: CallsReport | None = None
    """The calls of each VCF as a table and a datavzrd site for the report; unset: no report"""

    tool: Annotated[Tool, EnumField(Tool)]

    baf_file_generation: BafFileGeneration = BafFileGeneration()

    bcftools_stats: BcftoolsStats | None = None

    bcftools_roh: BcftoolsRoh | None = None

    bcftools_call: BcftoolsCall | None = None

    gatk3_hc: Gatk3Hc | None = None

    gatk3_ug: Gatk3Ug | None = None

    gatk4_hc_joint: Gatk4HcJoint | None = None

    gatk4_hc_gvcf: Gatk4HcGvcf | None = None

    mutect2: Mutect2 | None = None
    """Configuration for MuTect 2"""

    freebayes: Freebayes | None = None
    """Candidate calling with freebayes, for a ``varlociraptor`` task"""

    varlociraptor: Varlociraptor | None = None
    """Calling with varlociraptor on the candidates of ``depends_on.candidates``"""

    @model_validator(mode="after")
    def check_varlociraptor_settings(self) -> Self:
        """The varlociraptor VCF has no PL or AD, which the ROH and BAF reports need, and the
        somatic scenario calls each tumor on its own."""
        if self.tool != Tool.varlociraptor or self.varlociraptor is None:
            return self
        for name in ("bcftools_roh", "baf_file_generation"):
            section = getattr(self, name)
            if section is not None and section.enabled:
                raise ValueError(f"{name}: not available with tool: varlociraptor")
        if (
            self.varlociraptor.scenario == ScenarioKind.somatic
            and self.group_by == "related"
            and not self.varlociraptor.template
        ):
            raise ValueError(
                "group_by: related with the somatic scenario needs a template, the joint "
                "scenario of a donor's tumors; the preset calls each tumor with its matched normal"
            )
        return self
