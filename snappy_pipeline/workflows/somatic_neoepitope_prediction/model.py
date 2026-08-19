import dataclasses
import enum
import re

from typing import Union

from pydantic import model_validator

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToggleModel
from snappy_pipeline.workflows.hla_typing.model import (
    MHCIClassDnaTool,
    MHCIClassRnaTool,
    MHCIIClassDnaTool,
    MHCIIClassRnaTool,
)


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


class HlaTypingDnaTools(SnappyModel):
    class_i: list[MHCIClassDnaTool] = []
    class_ii: list[MHCIIClassDnaTool] = []


class HlaTypingRnaTools(SnappyModel):
    class_i: list[MHCIClassRnaTool] = []
    class_ii: list[MHCIIClassRnaTool] = []


class HlaTypingTools(SnappyModel):
    dna: HlaTypingDnaTools = HlaTypingDnaTools()
    rna: HlaTypingRnaTools = HlaTypingRnaTools()


class InputVariantType(enum.StrEnum):
    CALLING = "calling"
    ANNOTATION = "annotation"
    FILTRATION = "filtration"


class ClassIAlgorithm(enum.StrEnum):
    AllClassI = "all_class_i"

    DeepImmuno = "DeepImmuno"
    BigMHC_EL = "BigMHC_EL"
    BigMHC_IM = "BigMHC_IM"
    MHCflurry = "MHCflurry"
    MHCflurryEL = "MHCflurryEL"
    MixMHCpred = "MixMHCpred"
    PRIME = "PRIME"
    MHCnuggetsI = "MHCnuggetsI"
    NetMHC = "NetMHC"
    NetMHCpan = "NetMHCpan"
    NetMHCpanEL = "NetMHCpanEL"
    SMMPMBEC = "SMMPMBEC"
    SMM = "SMM"
    NetMHCcons = "NetMHCcons"
    PickPocket = "PickPocket"
    TLBind = "TLBind"
    TLImm = "TLImm"


class ClassIIAlgorithm(enum.StrEnum):
    AllClassII = "all_class_ii"

    MHCnuggetsII = "MHCnuggetsII"
    ImmunoScope_IM = "ImmunoScope_IM"
    NetMHCIIpan = "NetMHCIIpan"
    NetMHCIIpanEL = "NetMHCIIpanEL"
    NNalign = "NNalign"
    SMMalign = "SMMalign"
    MixMHC2pred = "MixMHC2pred"


class MLPredictions(ToggleModel):
    accept: float = 0.55
    reject: float = 0.3


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


class NetChopMethod(enum.StrEnum):
    CTERM = "cterm"
    TWENTY_S = "20s"


class NetChop(ToggleModel):
    path_netchop: str | None = None
    threshold: float = 0.5
    method: NetChopMethod = NetChopMethod.CTERM
    all_epitopes: bool = False

    @model_validator(mode="after")
    def ensurePathIsDefined(self):
        if self.enabled and not self.path_netchop:
            raise ValueError("Path to netchop binary must be defined")
        return self


class NetMHCStab(ToggleModel):
    all_epitopes: bool = False


class PVACtools(SnappyModel):
    path_container: str | None = None
    n_threads: int = 1
    use_all_transcripts: bool = False

    algorithms: list[Union[ClassIAlgorithm, ClassIIAlgorithm]] = [
        ClassIAlgorithm.AllClassI,
        ClassIIAlgorithm.AllClassII,
    ]

    class_i_epitope_length: list[int] = [8, 9, 10, 11]
    class_ii_epitope_length: list[int] = [12, 13, 14, 15, 16, 17, 18]

    use_normalized_percentiles: bool = False
    reference_score_path: str | None = None

    net_chop: NetChop = NetChop()
    netmhc_stab: NetMHCStab = NetMHCStab()

    genes_of_interest_file: str | None = None

    extra_args: list[str] = []

    @model_validator(mode="after")
    def ensure_valid_algorithm_and_length(self):
        if any(map(lambda a: a in ClassIAlgorithm, self.algorithms)):
            if ClassIAlgorithm.AllClassI in self.algorithms:
                for a in self.algorithms:
                    if a != ClassIAlgorithm.AllClassI and a in ClassIAlgorithm:
                        raise ValueError(
                            f"Algorithm {a} not compatible with all_class_i (already class I)"
                        )
            if len(self.class_i_epitope_length) == 0:
                raise ValueError("No epitope lengths supplied for class I algorithms")
        if any(map(lambda a: a in ClassIIAlgorithm, self.algorithms)):
            if ClassIIAlgorithm.AllClassII in self.algorithms:
                for a in self.algorithms:
                    if a != ClassIIAlgorithm.AllClassII and a in ClassIIAlgorithm:
                        raise ValueError(
                            f"Algorithm {a} not compatible with all_class_ii (already class II)"
                        )
            if len(self.class_ii_epitope_length) == 0:
                raise ValueError("No epitope lengths supplied for class II algorithms")
        return self

    @model_validator(mode="after")
    def ensure_reference_scores(self):
        if self.use_normalized_percentiles and not self.reference_score_path:
            raise ValueError(
                "A path for the reference score must be set when using normalized percentiles"
            )
        return self


class PVACseq(PVACtools):
    ml_predictions: MLPredictions = MLPredictions()


class PVACfuse(PVACtools):
    path_somatic_gene_fusion_calling: str = "../somatic_gene_fusion_calling"
    tool_somatic_gene_fusion_calling: SupportedGeneFusionTool = SupportedGeneFusionTool.ARRIBA


class PVACsplice(PVACtools):
    pass


class EnsemblVersion(enum.StrEnum):
    NONE = "none"
    GENE = "gene"
    TRANSCRIPT = "transcript"
    BOTH = "both"


class RnaMapping(ToggleModel):
    path_ngs_mapping: str = "../ngs_mapping"
    tool_rna_mapping: SupportedPileupTool = SupportedPileupTool.STAR

    extra_args: list[str] = []


class RnaQuantification(ToggleModel):
    path_gene_expression_quantification: str = "../gene_expression_quantification"
    tool_gene_expression_quantification: SupportedExpressionTool = SupportedExpressionTool.SALMON

    duplicate_transcripts_table: str | None = None

    ensembl_id: bool = True
    use_ensembl_version: EnsemblVersion = EnsemblVersion.NONE

    extra_args: list[str] = []

    @model_validator(mode="after")
    def ensure_ensemb_id_enabled_for_version(self):
        if self.use_ensembl_version != EnsemblVersion.NONE and not self.ensembl_id:
            raise ValueError("'use_ensembl_version' not allowed unless 'ensembl_id' is enabled")
        return self


class Phasing(ToggleModel):
    tool_ngs_mapping: str = "bwa"
    path_combine_variants: str = "../combine_variants"


class GermlineVariantStep(enum.StrEnum):
    CALL = "germline_variant_calling"
    ANNOTATION = "germline_variant_annotation"
    FILTER = "germline_variant_filtration"


class Proteome(ToggleModel):
    path_germline_variants: str | None = None
    tool_ngs_mapping: str = "bwa"
    tool_germline_variant_calling: str = "gatk4_hc"
    tool_germline_variant_annotation: str | None = None
    is_filtered: bool = True
    germline_variant_step: GermlineVariantStep = GermlineVariantStep.FILTER

    add_unmutated: bool = True
    external_proteome: str | None = None

    @model_validator(mode="after")
    def ensure_valid_variant_configuration(self):
        if self.enabled and self.path_germline_variants:
            match self.germline_variant_step:
                case GermlineVariantStep.CALL:
                    if self.is_filtered | self.tool_germline_variant_annotation:
                        raise ValueError(
                            "Filtration & annotation tool must be unset in calling mode"
                        )
                case GermlineVariantStep.FILTER:
                    if not self.is_filtered:
                        raise ValueError("Filtration must be set in filtration mode")
                case GermlineVariantStep.ANNOTATION:
                    if not self.tool_germline_variant_annotation:
                        raise ValueError("Annotation tool must be set in annotation mode")
        return self

    @model_validator(mode="after")
    def ensure_at_least_one_proteome_source(self):
        if (
            self.enabled
            and not self.path_germline_variants
            and not self.add_unmutated
            and not self.external_proteome
        ):
            raise ValueError(
                "One proteome source must be configured when proteome similarity is enabled"
            )
        return self


class SomaticNeoepitopePrediction(SnappyStepModel):
    tools: list[SupportedPredictionTool] = [SupportedPredictionTool.PVACSEQ]

    path_somatic_variant_annotation: str = "../somatic_variant_annotation"
    is_filtered: bool = True

    path_hla_typing: str = "../hla_typing"
    tools_hla_typing: HlaTypingTools = HlaTypingTools()

    pileup: RnaMapping = RnaMapping()
    quantification: RnaQuantification = RnaQuantification()
    phasing: Phasing = Phasing()
    proteome: Proteome = Proteome()

    pvacseq: PVACseq = PVACseq()
    pvacfuse: PVACfuse = PVACfuse()
    pvacsplice: PVACsplice = PVACsplice()

    @model_validator(mode="after")
    def ensure_at_least_one_tool_configured(self):
        if self.tools_hla_typing.dna.class_i is None and self.tools_hla_typing.dna.class_ii is None:
            raise ValueError("No HLA typing tools has been defined for DNA data")
        return self
