from typing import Annotated, Literal

from snappy_pipeline.models import SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.models.annotation import Mehari, MehariPickCriterion
from snappy_pipeline.workflows.abstract.protocol import (
    RELATED_LEVEL,
    DataSignature,
    DataType,
    Exons,
    ExpectedPathSchema,
    ExpectedVariants,
    Reference,
)
from snappy_pipeline.workflows.ngs_data_qc.model import ExpectedAlfredQc
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignmentQc


class ExpectedMehariDatabase(SnappyModel):
    """The mehari databases of the export: the HGNC table of the TSV, the transcript database
    (one or several), the ClinVar and the gnomAD frequency database"""

    hgnc_tsv: str
    transcript_db: str | list[str]
    clinvar_db: str
    frequency_db: str


class VarfishMehari(Mehari):
    """``mehari annotate seqvars`` for the VarFish table, which holds one row per variant and
    gene: the MANE transcript of each gene, else the longest, and intergenic variants too"""

    pick: Literal["gene", "transcript", "allele"] | None = "gene"

    pick_order: list[MehariPickCriterion] = ["mane-select", "mane-plus-clinical", "length"]

    keep_intergenic: bool | None = True


class VarfishStrucvars(ToolModel):
    """``mehari annotate strucvars``; its ``extra_arguments`` go to that command"""

    rng_seed: int = 42
    """Seed of the UUIDs of the structural variants, so that a rerun writes the same table"""


class VarfishExportDependsOn(SnappyModel):
    alignments: Annotated[
        str,
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        ExpectedPathSchema(ExpectedAlignmentQc),
    ]
    """The ``ngs_mapping`` task whose ``samtools stats`` and ``idxstats`` reports the BAM QC exports"""

    qc: Annotated[
        str,
        DataSignature(DataType.QC, frozenset({"alfred"})),
        ExpectedPathSchema(ExpectedAlfredQc),
    ] = ""
    """Optional ``ngs_data_qc`` task with ``tool: alfred``: the coverage of the targets of the
    libraries with a kit; the BAM QC export needs it. ``ngs_data_qc`` checks no libraries of
    background data sets, so pedigree members from them have no target coverage in the export."""

    variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline"})),
        RELATED_LEVEL,
        ExpectedPathSchema(ExpectedVariants),
    ]
    """The germline calls of each pedigree"""

    structural_variants: Annotated[
        str,
        DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"})),
        ExpectedPathSchema(ExpectedVariants),
    ] = ""
    """Optional ``sv_calling_targeted`` or ``sv_calling_wgs`` task whose calls are exported"""

    reference: Reference

    database: Annotated[
        str,
        DataSignature(DataType.DATABASE, frozenset({"mehari"})),
        ExpectedPathSchema(ExpectedMehariDatabase),
    ]
    """``external_data`` task of type ``database`` tagged ``mehari``: the HGNC table and the
    mehari databases, built by mehari 0.47 or newer"""

    exons: Exons = ""
    """Optional task with a BED file of exons; the export keeps only the small variants near them"""


class VarfishExport(SnappyStepModel):
    """Export for VarFish: the tables that ``varfish-cli`` imports, one set per pedigree."""

    depends_on: VarfishExportDependsOn

    mehari: VarfishMehari = VarfishMehari()
    """The annotation of the small variants; ``assembly`` must match the genome build that
    ``varfish-cli`` imports into"""

    strucvars: VarfishStrucvars = VarfishStrucvars()
    """The export of the structural variants"""
