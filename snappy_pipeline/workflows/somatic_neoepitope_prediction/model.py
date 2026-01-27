import enum

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.hla_typing.model import Tool as HlaTypingTool


class SupportedPredictionTool(enum.StrEnum):
    PVACSEQ = "pVACseq"


class SupportedPileupTool(enum.StrEnum):
    STAR = "star"


class SupportedExpressionTool(enum.StrEnum):
    SALMON = "salmon"
    OTHER = "other"


class Algorithm(enum.StrEnum):
    BigMHC_EL = "BigMHC_EL"
    BigMHC_IM = "BigMHC_IM"
    DeepImmuno = "DeepImmuno"
    MHCflurry = "MHCflurry"
    MHCflurryEL = "MHCflurryEL"
    MHCnuggetsI = "MHCnuggetsI"
    MHCnuggetsII = "MHCnuggetsII"
    NNalign = "NNalign"
    NetMHC = "NetMHC"
    NetMHCIIpan = "NetMHCIIpan"
    NetMHCIIpanEL = "NetMHCIIpanEL"
    NetMHCcons = "NetMHCcons"
    NetMHCpan = "NetMHCpan"
    NetMHCpanEL = "NetMHCpanEL"
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
    PERCENTILE = "percentile"
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


class PVACseq(SnappyModel):
    path_container: str | None = None
    n_threads: int = 1

    algorithms: list[Algorithm] | AllAlgorithms = AllAlgorithms.AllClassI
    netmhciipan_version: NetMHCIIpanVersion = NetMHCIIpanVersion.FourOne

    class_i_epitope_length: list[int] = [8, 9, 10, 11]
    class_ii_epitope_length: list[int] = []

    binding_threshold: int = 500
    percentile_threshold: float | None = 1.0
    percentile_threshold_strategy: PercentageThresholdStrategy = (
        PercentageThresholdStrategy.CONSERVATIVE
    )
    allele_specific_binding_thresholds: bool = False

    top_score_metric: TopScoreMetric = TopScoreMetric.MEDIAN
    top_score_metric2: TopScoreMetric2 = TopScoreMetric2.PERCENTILE

    normal_cov: int = 25
    tdna_cov: int = 25
    trna_cov: int = 2
    normal_vaf: float = 0.02
    tdna_vaf: float = 0.1
    trna_vaf: float = 0.25
    minimum_fold_change: float = 0.0
    expn_val: float = 1.0

    transcript_prioritization_strategy: TranscriptPrioritizationStrategy = (
        TranscriptPrioritizationStrategy.MANE_SELECT
    )
    maximum_transcript_support_level: int | None = None
    biotypes: list[str] = []
    allow_incomplete_transcript: bool = False

    net_chop_method: NetChopMethod | None = NetChopMethod.CTERM
    netmhc_stab: bool = True
    allele_specific_anchors: bool = False
    anchor_contribution_threshold: float = 0.8

    problematic_amino_acids: list[str] = []

    run_reference_proteome_similarity: bool = False
    peptide_fasta: str | None = None

    genes_of_interest_file: str | None = None

    exclude_NAs: bool = False

    downstream_sequence_length: int = 1000
    aggregate_inclusion_binding_threshold: int = 5000
    aggregate_inclusion_count_limit: int = 25

    @model_validator(mode="after")
    def ensure_peptide_fasta_exists(self):
        if self.run_reference_proteome_similarity and not self.peptide_fasta:
            raise ValueError(
                "Missing peptide fasta file required when enabling 'run_reference_proteome_similarity'"
            )
        return self


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


class RnaMapping(SnappyModel):
    enabled: bool = False

    path_ngs_mapping: str = "../ngs_mapping"
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


class RnaQuantification(SnappyModel):
    enabled: bool = False

    path_gene_expression_quantification: str = "../gene_expression_quantification"
    tool_gene_expression_quantification: SupportedExpressionTool = SupportedExpressionTool.SALMON

    ensembl_id: bool = True
    use_ensembl_version: bool = False

    annotation: str = "CSQ"
    annotation_description_regex: str = (
        r"^Consequence annotations from Ensembl VEP. Format: (?P<titles>.+)$"
    )
    annotation_separator: str = r"\|"
    annotation_gene_id: str | int = "Gene"
    annotation_transcript_id: str | int = "Feature"


class SomaticNeoepitopePrediction(SnappyStepModel):
    tools: list[SupportedPredictionTool] = [SupportedPredictionTool.PVACSEQ]

    path_somatic_variant_annotation: str = "../somatic_variant_annotation"
    is_filtered: bool = True

    path_hla_typing: str = "../hla_typing"
    tool_hla_typing: HlaTypingTool = HlaTypingTool.optitype

    pileup: RnaMapping = RnaMapping()
    quantification: RnaQuantification = RnaQuantification()

    pVACseq: PVACseq = PVACseq()
