import dataclasses
import enum
import re
from typing import Annotated

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToggleModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.hla_typing.model import (
    MHCIClassDnaTool,
    MHCIClassRnaTool,
    MHCIIClassDnaTool,
    MHCIIClassRnaTool,
)
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation.model import ExpectedVariantVcf


@dataclasses.dataclass
class MHC_CLASS:
    name: str
    genes: tuple[str]
    dirname: str
    filename: str
    prefix: str
    pattern: re.Pattern = dataclasses.field(init=False)

    def __post_init__(self):
        self.pattern = re.compile(
            r"^(?P<valid>({})\*[0-9]+:[0-9]+N?)(?P<suppl>.*)$".format("|".join(self.genes))
        )


MHC_CLASS_I = MHC_CLASS(
    name="class_i",
    genes=("A", "B", "C"),
    dirname="MHC_Class_I",
    filename="MHC_I",
    prefix="HLA-",
)

MHC_CLASS_II = MHC_CLASS(
    name="class_ii",
    genes=("DPA1", "DPB1", "DPA2", "DPB2", "DQA1", "DQB1"),
    dirname="MHC_Class_II",
    filename="MHC_II",
    prefix="",
)


class SupportedPredictionTool(enum.StrEnum):
    PVACSEQ = "pvacseq"
    PVACFUSE = "pvacfuse"
    PVACSPLICE = "pvacsplice"


class SupportedPileupTool(enum.StrEnum):
    STAR = "star"


class SupportedExpressionTool(enum.StrEnum):
    SALMON = "salmon"
    OTHER = "other"


class SupportedGeneFusionTool(enum.StrEnum):
    ARRIBA = "arriba"


class SupportedGermlineVariantCallingTool(enum.StrEnum):
    GATK4_HC = "gatk4_hc"


class SupportedGermlineVariantAnnotationTool(enum.StrEnum):
    VEP = "vep"


class HlaTypingDnaTool(SnappyModel):
    class_i: MHCIClassDnaTool | None = None
    class_ii: MHCIIClassDnaTool | None = None


class HlaTypingRnaTool(SnappyModel):
    class_i: MHCIClassRnaTool | None = None
    class_ii: MHCIIClassRnaTool | None = None


class HlaTypingTool(SnappyModel):
    dna: HlaTypingDnaTool = HlaTypingDnaTool()
    rna: HlaTypingRnaTool = HlaTypingRnaTool()


class InputVariantType(enum.StrEnum):
    CALLING = "calling"
    ANNOTATION = "annotation"
    FILTRATION = "filtration"


class Algorithm(enum.StrEnum):
    BigMHC_EL = "BigMHC_EL"
    BigMHC_IM = "BigMHC_IM"
    DeepImmuno = "DeepImmuno"
    ImmunoScope_IM = "ImmunoScope_IM"
    MHCflurry = "MHCflurry"
    MHCflurryEL = "MHCflurryEL"
    MHCnuggetsI = "MHCnuggetsI"
    MHCnuggetsII = "MHCnuggetsII"
    MixMHC2pred = "MixMHC2pred"
    MixMHCpred = "MixMHCpred"
    NNalign = "NNalign"
    NetMHC = "NetMHC"
    NetMHCIIpan = "NetMHCIIpan"
    NetMHCIIpanEL = "NetMHCIIpanEL"
    NetMHCcons = "NetMHCcons"
    NetMHCpan = "NetMHCpan"
    NetMHCpanEL = "NetMHCpanEL"
    PRIME = "PRIME"
    PickPocket = "PickPocket"
    SMM = "SMM"
    SMMPMBEC = "SMMPMBEC"
    SMMalign = "SMMalign"


class AllAlgorithms(enum.StrEnum):
    AllClassI = "all_class_i"
    AllClassII = "all_class_ii"
    All = "all"


class PercentageThresholdStrategy(enum.StrEnum):
    CONSERVATIVE = "conservative"
    EXPLORATORY = "exploratory"


class TopScoreMetric(enum.StrEnum):
    MEDIAN = "median"
    LOWEST = "lowest"


class TopScoreMetric2(enum.StrEnum):
    COMBINED_PERCENTILE = "combined_percentile"
    BINDING_PERCENTILE = "binding_percentile"
    IMMUNOGENICITY_PERCENTILE = "immunogenicity_percentile"
    PRESENTATION_PERCENTILE = "presentation_percentile"
    IC50 = "ic50"


class TranscriptPrioritizationStrategy(enum.StrEnum):
    MANE_SELECT = "mane_select"
    CANONICAL = "canonical"
    TSL = "tsl"


class NetChopMethod(enum.StrEnum):
    CTERM = "cterm"
    TWENTY_S = "20s"


class NetMHCIIpanVersion(enum.StrEnum):
    FourZero = "4.0"
    FourOne = "4.1"
    FourTwo = "4.2"
    FourThree = "4.3"


class AnchorType(enum.StrEnum):
    A = "A"
    D = "D"
    NDA = "NDA"


class ProteomeSimilarityByFile(SnappyModel):
    path: str


class ProteomeSimilarityByBlast(SnappyModel):
    todo: str


class ProteomeSimilarityByCreateProteome(SnappyModel):
    path: str = "../create_proteome"


class ProteomeSimilarity(ToggleModel):
    file: ProteomeSimilarityByFile | None = None
    blast: ProteomeSimilarityByBlast | None = None
    create_proteome: ProteomeSimilarityByCreateProteome | None = None

    @model_validator(mode="after")
    def ensureOneAndOnlyOneModeDefined(self):
        if self.enabled:
            if (
                (self.file and (self.blast or self.create_proteome))
                or (self.blast and (self.create_proteome or self.file))
                or (self.create_proteome and (self.file or self.blast))
            ):
                raise ValueError("Only one proteome similarity mode can be defined")
            if not (self.file or self.blast or self.create_proteome):
                raise ValueError("One proteome similarity mode must be defined")
        return self


class NetChop(ToggleModel):
    path_netchop: str | None = None
    threshold: float = 0.5
    method: NetChopMethod = NetChopMethod.CTERM

    @model_validator(mode="after")
    def ensurePathIsDefined(self):
        if self.enabled and not self.path_netchop:
            raise ValueError("Path to netchop binary must be defined")
        return self


class NetMHCStab(ToggleModel):
    path_netmhc_stab: str | None = None
    path_netmhc_pan: str | None = None

    @model_validator(mode="after")
    def ensurePathIsDefined(self):
        if self.enabled and not (self.path_netmhc_stab and self.path_netmhc_pan):
            raise ValueError("Path to netMHCstanpan & netMHCpan binaries must be defined")
        return self


class PVACtools(SnappyModel):
    path_container: str | None = None
    n_threads: int = 1

    algorithms: list[Algorithm] | AllAlgorithms = AllAlgorithms.AllClassI
    netmhciipan_version: NetMHCIIpanVersion = NetMHCIIpanVersion.FourOne

    class_i_epitope_length: list[int] = [8, 9, 10, 11]
    class_ii_epitope_length: list[int] = []

    use_normalized_percentiles: bool = False

    binding_threshold: int = 500
    binding_percentile_threshold: float = 2.0
    presentation_percentile_threshold: float = 2.0
    immunogenicity_percentile_threshold: float = 2.0
    percentile_threshold_strategy: PercentageThresholdStrategy = (
        PercentageThresholdStrategy.CONSERVATIVE
    )
    allele_specific_binding_thresholds: bool = False

    top_score_metric: TopScoreMetric = TopScoreMetric.MEDIAN
    top_score_metric2: list[TopScoreMetric2] = [
        TopScoreMetric2.IC50,
        TopScoreMetric2.COMBINED_PERCENTILE,
    ]

    net_chop: NetChop = NetChop()
    netmhc_stab: NetMHCStab = NetMHCStab()

    expn_val: float = 1.0

    problematic_amino_acids: list[str] = []

    genes_of_interest_file: str | None = None

    fasta_size: int = 200

    exclude_NAs: bool = False

    aggregate_inclusion_binding_threshold: int = 5000
    aggregate_inclusion_count_limit: int = 25

    @model_validator(mode="after")
    def ensure_one_class_length_defined(self):
        if not (self.class_i_epitope_length or self.class_ii_epitope_length):
            raise ValueError("Epitope lengths must be defined for at least one MHC class")
        return self


class PVACseq(PVACtools):
    use_all_transcripts: bool = False

    normal_cov: int = 25
    tdna_cov: int = 25
    trna_cov: int = 2
    normal_vaf: float = 0.02
    tdna_vaf: float = 0.1
    trna_vaf: float = 0.25
    minimum_fold_change: float = 0.0

    transcript_prioritization_strategy: list[TranscriptPrioritizationStrategy] = [
        TranscriptPrioritizationStrategy.CANONICAL,
        TranscriptPrioritizationStrategy.MANE_SELECT,
        TranscriptPrioritizationStrategy.TSL,
    ]
    maximum_transcript_support_level: int = 1
    biotypes: list[str] = ["protein_coding"]
    allow_incomplete_transcript: bool = False

    allele_specific_anchors: bool = False
    anchor_contribution_threshold: float = 0.8

    downstream_sequence_length: int = 1000

    run_ml_predictions: bool = False
    ml_threshold_accept: float = 0.55
    ml_threshold_reject: float = 0.3


class PVACfuse(PVACtools):
    tool_somatic_gene_fusion_calling: SupportedGeneFusionTool = SupportedGeneFusionTool.ARRIBA

    downstream_sequence_length: int = 1000

    read_support: int = 5


class PVACsplice(PVACtools):
    use_all_transcripts: bool = False

    normal_cov: int = 25
    tdna_cov: int = 25
    trna_cov: int = 2
    normal_vaf: float = 0.02
    tdna_vaf: float = 0.1
    trna_vaf: float = 0.25

    transcript_prioritization_strategy: list[TranscriptPrioritizationStrategy] = [
        TranscriptPrioritizationStrategy.CANONICAL,
        TranscriptPrioritizationStrategy.MANE_SELECT,
        TranscriptPrioritizationStrategy.TSL,
    ]
    maximum_transcript_support_level: int = 1
    biotypes: list[str] = ["protein_coding"]
    allow_incomplete_transcript: bool = False

    junction_score: int = 10
    variant_distance: int = 100
    anchor_types: list[AnchorType] = [AnchorType.A, AnchorType.D, AnchorType.NDA]


class BAQ(enum.StrEnum):
    NO = "no"
    FULL = "full"
    REDO = "redo"


class FORMATTAGS(enum.StrEnum):
    AD = "FORMAT/AD"
    ADF = "FORMAT/ADF"
    ADR = "FORMAT/ADR"
    DP = "FORMAT/DP"
    NMBZ = "FORMAT/NMBZ"
    QS = "FORMAT/QS"
    SP = "FORMAT/SP"
    SCR = "FORMAT/SCR"


class INFOTAGS(enum.StrEnum):
    AD = "INFO/AD"
    ADF = "INFO/ADF"
    ADR = "INFO/ADR"
    BQBZ = "INFO/BQBZ"
    FS = "INFO/FS"
    IDV = "INFO/IDV"
    IMF = "INFO/IMF"
    MIN_PL_SUM = "INFO/MIN_PL_SUM"
    MQ0F = "INFO/MQ0F"
    MQBZ = "INFO/MQBZ"
    MQSBZ = "INFO/MQSBZ"
    NM = "INFO/NM"
    NMBZ = "INFO/NMBZ"
    RPBZ = "INFO/RPBZ"
    SCBZ = "INFO/SCBZ"
    SCR = "INFO/SCR"
    SGB = "INFO/SGB"
    VDB = "INFO/VDB"


class EnsemblVersion(enum.StrEnum):
    NONE = "none"
    GENE = "gene"
    TRANSCRIPT = "transcript"
    BOTH = "both"


class RnaMapping(ToggleModel):
    tool_rna_mapping: SupportedPileupTool = SupportedPileupTool.STAR

    baq: BAQ | None = None
    max_depth: int = 250
    min_MQ: int = 0
    adjust_MQ: int = 0
    min_BQ: int = 1
    max_BQ: int = 60
    delta_BQ: int = 30
    annotate: list[INFOTAGS | FORMATTAGS] = []
    snp_indel_genotype: list[str] = []


class RnaQuantification(ToggleModel):
    tool_gene_expression_quantification: SupportedExpressionTool = SupportedExpressionTool.SALMON

    duplicate_transcripts_table: str | None = None

    ensembl_id: bool = True
    use_ensembl_version: EnsemblVersion = EnsemblVersion.NONE

    annotation: str = "CSQ"
    annotation_description_regex: str = (
        r"^Consequence annotations from Ensembl VEP. Format: (?P<titles>.+)$"
    )
    annotation_separator: str = r"\|"
    annotation_gene_id: str | int = "Gene"
    annotation_transcript_id: str | int = "Feature"

    @model_validator(mode="after")
    def ensure_ensemb_id_enabled_for_version(self):
        if self.use_ensembl_version != EnsemblVersion.NONE and not self.ensembl_id:
            raise ValueError("'use_ensembl_version' not allowed unless 'ensembl_id' is enabled")
        return self


class Phasing(ToggleModel):
    tool_ngs_mapping: str = "bwa"


class Proteome(ToggleModel):
    add_unmutated: bool = True
    external_proteome: str | None = None

    @model_validator(mode="after")
    def ensure_at_least_one_proteome_source(self):
        if self.enabled and not self.add_unmutated and not self.external_proteome:
            raise ValueError(
                "One proteome source must be configured when proteome similarity is enabled"
            )
        return self


class SomaticNeoepitopePredictionDependsOn(SnappyModel):
    hla_typing: Annotated[str, DataSignature(DataType.TABULAR, frozenset({"hla"}))]
    somatic_variant_annotation: Annotated[
        str,
        DataSignature(DataType.VARIANTS),
        ExpectedPathSchema(ExpectedVariantVcf),
    ] = ""
    ngs_mapping: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS),
        ExpectedPathSchema(ExpectedAlignments),
    ] = ""
    gene_expression_quantification: Annotated[
        str, DataSignature(DataType.EXPRESSION, frozenset({"rna"}))
    ] = ""
    combine_variants: Annotated[str, DataSignature(DataType.VARIANTS)] = ""
    somatic_gene_fusion_calling: Annotated[
        str, DataSignature(DataType.VARIANTS, frozenset({"somatic", "fusion", "rna"}))
    ] = ""
    germline_variant: Annotated[str, DataSignature(DataType.VARIANTS)] = ""


class SomaticNeoepitopePrediction(SnappyStepModel):
    depends_on: SomaticNeoepitopePredictionDependsOn

    tool: SupportedPredictionTool = SupportedPredictionTool.PVACSEQ

    tool_hla_typing: HlaTypingTool = HlaTypingTool()

    pileup: RnaMapping = RnaMapping()
    quantification: RnaQuantification = RnaQuantification()
    phasing: Phasing = Phasing()
    proteome: Proteome = Proteome()

    pvacseq: PVACseq = PVACseq()
    pvacfuse: PVACfuse = PVACfuse()
    pvacsplice: PVACsplice = PVACsplice()

    @model_validator(mode="after")
    def ensure_at_least_one_tool_configured(self):
        if self.tool_hla_typing.dna.class_i is None and self.tool_hla_typing.dna.class_ii is None:
            raise ValueError("No HLA typing tool has been defined for DNA data")
        return self

    @model_validator(mode="after")
    def ensure_tool_dependencies_satisfied(self):
        """Validate tool-specific depends_on requirements."""
        deps = self.depends_on

        # HLA typing is mandatory across all prediction tools
        if not deps.hla_typing:
            raise ValueError(
                f"depends_on.hla_typing is required for neoepitope prediction (tool: {self.tool!r})"
            )

        match self.tool:
            case SupportedPredictionTool.PVACSEQ:
                if not deps.somatic_variant_annotation:
                    raise ValueError(
                        "depends_on.somatic_variant_annotation is required when tool is 'pvacseq'"
                    )

            case SupportedPredictionTool.PVACSPLICE:
                if not deps.somatic_variant_annotation:
                    raise ValueError(
                        "depends_on.somatic_variant_annotation is required when tool is 'pvacsplice'"
                    )
                if not deps.ngs_mapping:
                    raise ValueError("depends_on.ngs_mapping is required when tool is 'pvacsplice'")

            case SupportedPredictionTool.PVACFUSE:
                if not deps.somatic_gene_fusion_calling:
                    raise ValueError(
                        "depends_on.somatic_gene_fusion_calling is required when tool is 'pvacfuse'"
                    )

        return self
