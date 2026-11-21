import enum
from typing import Annotated, Literal

from pydantic import BaseModel, Field

from snappy_pipeline.models import ConfigPath, ContainerImage, ReadBy, RequiredBy
from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, ToolModel
from snappy_pipeline.models.cnvkit import PanelOfNormals as CnvKit
from snappy_pipeline.models.gatk import GATK
from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    DnaAlignments,
    ExpectedPathSchema,
    GermlineResource,
    Kits,
    Reference,
)


class ExpectedGenomicsDb(BaseModel):
    """Consumer-driven contract: the GenomicsDB of a mutect2 panel of normals task."""

    genomicsdb: str


class ExpectedAccess(BaseModel):
    """Consumer-driven contract: the accessible regions (BED)."""

    access: str


class Tool(enum.StrEnum):
    mutect2 = "mutect2"
    cnvkit = "cnvkit"
    purecn = "purecn"
    access = "access"


class GenomicsDb(GATK):
    java_options: str = "-Xms16g -Xmx32g"
    """
    GenomicsDBImport needs lots of memory.
    (mind the caveats in https://gatk.broadinstitute.org/hc/en-us/articles/35967560824859-GenomicsDBImport)
    """


class Mutect2(GATK):
    num_jobs: Annotated[int, Field(ge=1)] = 24
    """The number of regions the genome is split into for scattered calling"""

    path_normals_list: ConfigPath = ""
    """Optional file listing libraries to include in panel"""

    padding: int = 5000
    """Padding around intervals for scatter/gather"""

    genomicsdb: GenomicsDb = GenomicsDb()
    """
    Parameters for the creation of the (transient?) genomics database.
    The 2nd part of panel of normals creation (CreateSomaticPanelOfNormals) is not under user control
    """


class Access(ToolModel):
    """Creates access file for cnvkit, based on genomic sequence & excluded regions (optionally)"""

    exclude: list[str] = []
    """[access] Bed file of regions to exclude (mappability, blacklisted, ...)"""

    min_gap_size: int = 0
    """[access] Minimum gap size between accessible sequence regions (0: use default value)"""


class GenomeName(enum.StrEnum):
    hg18 = "hg18"
    hg19 = "hg19"
    hg38 = "hg38"
    mm9 = "mm9"
    mm10 = "mm10"
    rn4 = "rn4"
    rn5 = "rn5"
    rn6 = "rn6"
    canFam3 = "canFam3"


class PureCn(ToolModel):
    container: ContainerImage = (
        "docker://markusriester/purecn"
        "@sha256:5d676079969281d343199e291cfbc529bf1b2c69c70342311998f90386acddce"
    )
    """
    The image that Snakemake pulls to run the PureCN jobs, also those of the calling tasks that read
    this panel; the default is ``markusriester/purecn:latest`` of 2024-05-28, pinned by its digest.
    The path of a local ``.sif`` file works too.
    """

    path_normals_list: ConfigPath = ""
    """Optional file listing libraries to include in panel"""

    genome_name: Annotated[
        GenomeName | Literal["unknown"],
        EnumField(GenomeName, json_schema_extra={"options": {"unknown"}}),
    ] = "unknown"
    """The genome of PureCN's annotations, for the calling tasks that read this panel, too"""

    enrichment_kit_name: str = "unknown"
    """For filename only..."""

    mappability: str = ""
    """
    GRCh38:
    ``mappability.bw`` from the PureCN extdata (hg38)
    """

    reptiming: str = ""
    """Nothing for GRCh38"""

    seed: int = 1234567
    """The seed of PureCN's random numbers, for the calling tasks that read this panel, too"""


class PanelOfNormalsDependsOn(SnappyModel):
    alignments: Annotated[DnaAlignments, RequiredBy(Tool.mutect2, Tool.cnvkit, Tool.purecn)] = ""

    panel_of_normals: Annotated[
        str,
        DataSignature(DataType.MODELS, frozenset({"pon", "mutect2"})),
        ExpectedPathSchema(ExpectedGenomicsDb),
        RequiredBy(Tool.purecn),
    ] = ""
    """
    A ``panel_of_normals`` task with ``tool: mutect2``, whose genomicsDB PureCN's NormalDB.R
    reads
    """

    access: Annotated[
        str,
        DataSignature(DataType.RAW, frozenset({"access"})),
        ExpectedPathSchema(ExpectedAccess),
        ReadBy(Tool.cnvkit),
    ] = ""
    """
    The accessible regions: a ``panel_of_normals`` task with ``tool: access``, or an
    ``external_data`` task that produces ``{type: raw, tags: [access]}`` with the key ``access``
    """

    reference: Reference

    germline_resource: Annotated[GermlineResource, RequiredBy(Tool.mutect2)] = ""
    """Population allele frequencies (af-only gnomAD VCF)"""

    kits: Annotated[Kits, RequiredBy(Tool.purecn), ReadBy(Tool.cnvkit)] = ""
    """
    For ``purecn``, the baits of the kit of the normals are the bait regions (MergedProbes for
    Agilent SureSelect, recommended by the PureCN author).
    For ``cnvkit``, the baits are the targets of the exome panel; without kits, the panel is
    built for whole genome data.
    All normals must be targeted libraries (WES or Panel_seq) of one kit.
    """


class PanelOfNormals(SnappyStepModel):
    depends_on: PanelOfNormalsDependsOn

    tool: Annotated[Tool, EnumField(Tool)]

    ignore_contigs: Annotated[
        list[str],
        Field(
            examples=[
                ["^NC_007605$", "^hs37d5$", "^chrEBV$", "_decoy$", "^HLA-", "^GL000220"],
                [
                    "^chrEBV$",
                    "^HPV",
                    "^CMV$",
                    "^HBV$",
                    "^HCV-",
                    "^HIV-",
                    "^KSHV$",
                    "^HTLV-1$",
                    "^MCV$",
                    "_decoy$",
                    "^chrUn_GL00220",
                    "^SV40$",
                ],
            ]
        ),
    ] = []
    """Regular expressions matched with ``re.search`` against contig names; matching contigs are ignored."""

    mutect2: Mutect2 | None = None

    cnvkit: CnvKit | None = None

    access: Access = Access()

    purecn: PureCn | None = None
