import enum
from typing import Annotated, Union

from pydantic import BaseModel, Field, model_validator

from snappy_pipeline.models import ConfigPath, SnappyModel, SnappyStepModel, ToggleModel
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedVcf,
    Features,
    Reference,
)
from snappy_pipeline.workflows.hla_typing.model import ExpectedHlaTyping
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments, ExpectedRnaAlignments


class Tool(enum.StrEnum):
    pvactools = "pvactools"


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
    """pVACseq's machine learning predictions (``--run-ml-predictions``)"""

    accept: float = 0.55
    reject: float = 0.3


class NetChopMethod(enum.StrEnum):
    CTERM = "cterm"
    TWENTY_S = "20s"


class NetChop(ToggleModel):
    """Proteasomal cleavage prediction with a local NetChop installation"""

    path: ConfigPath | None = None
    """The ``netchop`` executable"""

    threshold: float = 0.5
    method: NetChopMethod = NetChopMethod.CTERM
    all_epitopes: bool = False
    """Predict for all epitopes, not only the filtered ones"""

    @model_validator(mode="after")
    def ensure_path_when_enabled(self):
        if self.enabled and not self.path:
            raise ValueError("net_chop needs the path of the netchop executable")
        return self


class NetMHCStab(ToggleModel):
    """Stability prediction with netMHCstabpan in the container (MHC class I)"""

    all_epitopes: bool = False
    """Predict for all epitopes, not only the filtered ones (or those of NetChop)"""


class PvacTool(ToggleModel):
    """Options shared by pVACseq, pVACfuse and pVACsplice"""

    n_threads: int = 1
    algorithms: list[Union[ClassIAlgorithm, ClassIIAlgorithm]] = [
        ClassIAlgorithm.AllClassI,
        ClassIIAlgorithm.AllClassII,
    ]
    class_i_epitope_length: list[int] = [8, 9, 10, 11]
    class_ii_epitope_length: list[int] = [12, 13, 14, 15, 16, 17, 18]
    net_chop: NetChop = NetChop()
    netmhc_stab: NetMHCStab = NetMHCStab()
    genes_of_interest: ConfigPath | None = None
    """A file with the genes of interest (``--genes-of-interest-file``)"""

    extra_args: list[str] = []

    @model_validator(mode="after")
    def ensure_valid_algorithm_and_length(self):
        for algorithms, all_name, lengths, name in (
            (ClassIAlgorithm, ClassIAlgorithm.AllClassI, self.class_i_epitope_length, "I"),
            (ClassIIAlgorithm, ClassIIAlgorithm.AllClassII, self.class_ii_epitope_length, "II"),
        ):
            chosen = [a for a in self.algorithms if a in algorithms]
            if all_name in chosen and len(chosen) > 1:
                raise ValueError(f"{all_name} already runs all class {name} algorithms")
            if chosen and not lengths:
                raise ValueError(f"No epitope lengths for the class {name} algorithms")
        return self


class Pileup(ToggleModel):
    """Read counts of the tumor RNA at the somatic variants (``bcftools mpileup``)"""

    extra_args: list[str] = []


class Expression(ToggleModel):
    """Gene and transcript TPMs of the tumor RNA (salmon), added to the VCF"""

    duplicate_transcripts_table: ConfigPath | None = None
    extra_args: list[str] = []


class PVACseq(PvacTool):
    enabled: bool = True
    """pVACseq is enabled by default"""

    ml_predictions: MLPredictions = MLPredictions()
    pileup: Pileup = Pileup()
    """Needs ``depends_on.alignments_rna``"""

    expression: Expression = Expression()
    """Needs ``depends_on.expression``"""

    phasing: ToggleModel = ToggleModel()
    """Phase proximal germline variants (GATK 3 ReadBackedPhasing); needs
    ``depends_on.combined_variants`` and ``depends_on.alignments_dna``"""

    def __init__(self, enabled: bool = True, **kwargs):
        super().__init__(enabled=enabled, **kwargs)


class PVACfuse(PvacTool):
    pass


class PVACsplice(PvacTool):
    pass


class PvacTools(SnappyModel):
    container: ConfigPath | None = None
    """The pVACtools container image (``.sif``); without it, the task builds one from
    ``docker://griffithlab/pvactools:7.1.3``"""

    pvacseq: PVACseq = PVACseq()
    pvacfuse: PVACfuse = PVACfuse()
    pvacsplice: PVACsplice = PVACsplice()

    @model_validator(mode="after")
    def ensure_one_prediction(self):
        if not (self.pvacseq.enabled or self.pvacfuse.enabled or self.pvacsplice.enabled):
            raise ValueError("Enable at least one of pvacseq, pvacfuse and pvacsplice")
        return self


class ExpectedSalmonExpression(BaseModel):
    """Consumer-driven contract: the salmon TPMs of an RNA library."""

    gene_sf: str
    transcript_sf: str


class ExpectedFusions(BaseModel):
    """Consumer-driven contract: the fusions table of an RNA library (arriba)."""

    fusions: str


class ExpectedProteome(BaseModel):
    """Consumer-driven contract: a reference proteome."""

    fasta: str


class SomaticNeoepitopePredictionDependsOn(SnappyModel):
    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "annotated", "vep"})),
        ExpectedPathSchema(ExpectedVcf),
    ] = ""
    """Somatic variants annotated by VEP with the ``Frameshift`` and ``Wildtype`` plugins
    (pvacseq, pvacsplice)"""

    hla_types: Annotated[
        list[str],
        DataSignature(DataType.TABULAR, frozenset({"hla"})),
        ExpectedPathSchema(ExpectedHlaTyping),
    ] = Field(min_length=1)
    """HLA typing tasks; the HLA types of all of them are merged"""

    reference: Reference

    features: Features = ""
    """Gene annotation (pvacsplice)"""

    # Two inputs of one kind, so qualified by molecule; a step that needs only one uses plain
    # ``alignments``, narrowed by the tags of its signature
    alignments_rna: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedRnaAlignments),
    ] = ""
    """STAR alignments of the tumor RNA, with their strandedness (pvacseq pileup, pvacsplice)"""

    alignments_dna: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignments),
    ] = ""
    """Alignments of the tumor DNA (pvacseq phasing)"""

    combined_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"combined", "annotated", "vep"})),
        ExpectedPathSchema(ExpectedVcf),
    ] = ""
    """Somatic and germline variants combined and annotated by VEP (pvacseq phasing)"""

    expression: Annotated[
        str,
        DataSignature(DataType.EXPRESSION, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedSalmonExpression),
    ] = ""
    """Salmon TPMs of the tumor RNA (pvacseq expression)"""

    fusions: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "fusion", "rna"})),
        ExpectedPathSchema(ExpectedFusions),
    ] = ""
    """Arriba fusions of the tumor RNA (pvacfuse)"""

    proteome: Annotated[
        str,
        DataSignature(DataType.RAW, frozenset({"proteome"})),
        ExpectedPathSchema(ExpectedProteome),
    ] = ""
    """A reference proteome for the similarity check of the epitopes (optional)"""


class SomaticNeoepitopePrediction(SnappyStepModel):
    depends_on: SomaticNeoepitopePredictionDependsOn

    tool: Tool

    pvactools: PvacTools = PvacTools()

    @model_validator(mode="after")
    def ensure_inputs_of_enabled_predictions(self):
        """Each enabled prediction and option has the ``depends_on`` fields it reads."""
        deps, tools = self.depends_on, self.pvactools
        needed = {}
        if tools.pvacseq.enabled:
            needed["pvacseq"] = ["variants"]
            if tools.pvacseq.pileup.enabled:
                needed["pvacseq.pileup"] = ["alignments_rna"]
            if tools.pvacseq.expression.enabled:
                needed["pvacseq.expression"] = ["expression"]
            if tools.pvacseq.phasing.enabled:
                needed["pvacseq.phasing"] = ["combined_variants", "alignments_dna"]
        if tools.pvacfuse.enabled:
            needed["pvacfuse"] = ["fusions"]
        if tools.pvacsplice.enabled:
            needed["pvacsplice"] = ["variants", "alignments_rna", "features"]
        missing = [
            f"{name} needs depends_on.{field}"
            for name, fields in needed.items()
            for field in fields
            if not getattr(deps, field)
        ]
        if missing:
            raise ValueError("; ".join(missing))
        return self
