import enum
from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel


class Format(enum.StrEnum):
    stringtie = "stringtie"
    snappy_custom = "snappy_custom"
    cufflinks = "cufflinks"
    kallisto = "kallisto"
    custom = "custom"


class Preparation(SnappyModel):
    format: Annotated[Format, EnumField(Format)] = Format.snappy_custom
    "The file format of the expression file to process. (stringtie, kallisto, cufflinks, snappy_custom, custom)"
    "Use `custom` to process file formats not explicitly supported."
    "The `custom` option requires the use of the --id-column and --expression-column arguments."
    path_features: str = ""
    """Gencode file path, required for star and snappy format"""
    mode: Annotated[str, Field(examples=["gene", "transcript"])] = "gene"
    """Determine whether the expression file contains gene or transcript TPM values."""
    full_vep_annotation: bool = True
    """The somatic_variant_annotation has been done in fully annotation mode"""
    id_column: str = ""
    """Gene/Transcript id column name. Need when using the `custom` format."""
    expression_column: str = ""
    """Expression column name. Need when using the `custom` format."""
    ignore_ensembl_id_version: bool = True
    """Ignore the ensemble id version"""
    max_depth: int = 4000
    "max_depth option of bcftools mpileup"


class Algorithms(enum.StrEnum):
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
    all_class_i = "all_class_i"


class Prediction(SnappyModel):
    CLASS_I_EPITOPE_LENGTH: list[int] = [8, 9, 10, 11]
    "Length of MHC Class I subpeptides (neoepitopes) to predict."
    "Multiple epitope lengths can be specified using a comma-separated list."
    "Typical epitope lengths vary between 8-15."
    algorithms: Annotated[list[Algorithms], EnumField(Algorithms, [])]
    "The epitope prediction algorithms to use."
    "For running with all prediction assign to ['all_class_i']"
    # Following parameter please follow the documentation of pvactools to know its function
    BINDING_THRESHOLD: int = 500
    percentile_threshold: float | None = Field(default=None)
    allele_specific_binding_thresholds: bool = False
    aggregate_inclusion_binding_threshold: int = 5000
    netmhc_stab: bool = False
    NET_CHOP_THRESHOLD: float = 0.5
    PROBLEMATIC_AMINO_ACIDS: list[str] | None = Field(default=None)
    # top-score-metric
    # net-chop-method
    # net-chop-threshold

    # E.g amino_acid:position - G:2 would check for a Glycine at
    # the second position of the epitope
    # run_reference_proteome_similarity: bool = False
    # blastp_path
    # blastp-db
    # peptide-fasta
    FAST_SIZE: int = 200
    exclude_NAs: bool = False
    NORMAL_COV: int = 5
    TDNA_COV: int = 10
    TRNA_COV: int = 10
    NORMAL_VAF: float = 0.02
    maximum_transcript_support_level: int = Field(gt=0, lt=6, default=1)
    # allele-specific-anchors
    # anchor_contribution_threshold: float=0.8
    pass_only: bool = False
    tumor_purity: float | None = Field(default=None)


class SomaticNeoepitopePrediction(SnappyStepModel):
    path_somatic_variant_annotation: Annotated[
        str, Field(examples=["../somatic_variant_annotation"])
    ]
    path_container: Annotated[
        str, Field(examples=["../somatic_neoepitope_prediction/work/containers/out/pvactools.simg"])
    ]
    """
    Running somatic neoepitope prediction with pvactools 
    is required,with the container
    """
    path_hla_typing: Annotated[str, Field(examples=["../hla_typing"])]
    path_rna_ngs_mapping: Annotated[str, Field(examples=["../ngs_mapping"])]
    tools_somatic_variant_annotation: list[str] = ["vep"]
    tools_ngs_mapping: list[str] = []
    tools_somatic_variant_calling: list[str] = []
    tools_rna_mapping: list[str] = []
    """Deafult to those configured for ngs_mapping"""
    tools_ngs_mapping: list[str] = []
    """Deafult to those configured for ngs_mapping"""
    tools_somatic_variant_calling: list[str] = []
    """Deafult to those configured for somatic_variant_calling"""
    preparation: Preparation | None = None
    prediction: Prediction | None = None
