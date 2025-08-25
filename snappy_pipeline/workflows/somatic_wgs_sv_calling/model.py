import enum
from typing import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel, validators


class Tool(enum.StrEnum):
    manta = "manta"
    delly2 = "delly2"


class MantaExtraConfig(SnappyModel):
    minCandidateVariantSize: int = 8
    rnaMinCandidateVariantSize: int = 1000
    """
    Run discovery and candidate reporting for all SVs/indels at or above this size
    Separate option (to provide different default) used for runs in RNA-mode
    """

    minEdgeObservations: int = 3
    """
    Remove all edges from the graph unless they're supported by this many 'observations'.
    Note that one supporting read pair or split read usually equals one observation, but evidence is sometimes downweighted.
    """

    graphNodeMaxEdgeCount: int = 10
    """If both nodes of an edge have an edge count higher than this, then skip evaluation of the edge. Set to 0 to turn this filtration off"""

    minCandidateSpanningCount: int = 3
    """Run discovery and candidate reporting for all SVs/indels with at least this many spanning support observations"""

    minScoredVariantSize: int = 50
    """After candidate identification, only score and report SVs/indels at or above this size"""

    minDiploidVariantScore: int = 10
    """Minimum VCF "QUAL" score for a variant to be included in the diploid vcf"""

    minPassDiploidVariantScore: int = 20
    """VCF "QUAL" score below which a variant is marked as filtered in the diploid vcf"""

    minPassDiploidGTScore: int = 15
    """Minimum genotype quality score below which single samples are filtered for a variant in the diploid vcf"""

    minSomaticScore: int = 10
    """Somatic quality scores below this level are not included in the somatic vcf"""

    minPassSomaticScore: int = 30
    """Somatic quality scores below this level are filtered in the somatic vcf"""

    enableRemoteReadRetrievalForInsertionsInGermlineCallingModes: int = 1
    enableRemoteReadRetrievalForInsertionsInCancerCallingModes: int = 0
    """
    Remote read retrieval is used ot improve the assembly of putative insertions by retrieving any mate reads in remote
    locations with poor mapping quality, which pair to confidently mapping reads near the insertion locus. These reads
    can help to fully assemble longer insertions, under certain circumstances this feature can add a very large runtime
    burden. For instance, given the very high chimeric pair rates found in degraded FFPE samples, the runtime of the read
    retrieval process can be unpredicable. For this reason the feature is disabled by default for somatic variant calling.
    This feature can be enabled/disabled separately for germline and cancer calling below.

    Here "CancerCallingModes" includes tumor-normal subtraction and tumor-only calling. "GermlineCallingModes" includes all other calling modes.
    """

    useOverlapPairEvidence: int = 0
    """Set if an overlapping read pair will be considered as evidence. Set to 0 to skip overlapping read pairs"""

    enableEvidenceSignalFilter: int = 1
    """Set the filter on candidates of insignificant evidence signal. This is forced to 0 for runs in RNA-mode"""


class MantaExtendedOptions(SnappyModel):
    existingAlignStatsFile: str | None = None
    """Pre-calculated alignment statistics file. Skips alignment stats calculation."""

    useExistingChromDepths: bool = False
    """Use pre-calculated chromosome depths."""

    generateEvidenceBam: bool = False
    """Generate a bam of supporting reads for all SVs."""

    outputContig: bool = False
    """Output assembled contig sequences in VCF file."""

    region: list[str] = []
    """
    Limit the analysis to a region of the genome for debugging purposes.

    If this argument is provided multiple times all specified regions will be analyzed together.
    All regions must be non-overlapping to get a meaningful result. Examples:
    - '--region chr20' (whole chromosome),
    - '--region chr2:100-2000 --region chr3:2500-3000' (two regions)

    If this option is specified (one or more times) together with the --callRegions BED file,
    then all region arguments will be intersected with the callRegions BED track.
    """

    retainTempFiles: bool = False
    """Keep all temporary files (for workflow debugging)"""

    scanSizeMb: int = 12
    """Maximum sequence region size (in megabases) scanned by each task during SV Locus graph generation."""

    callMemMb: int | None = None
    """
    Set default task memory requirement (in megabytes) for common tasks.

    This may benefit an analysis of unusual depth, chimera rate, etc..

    'Common' tasks refers to most compute intensive scatter-phase tasks of graph creation and candidate generation.
    """


class Manta(SnappyModel):
    ignore_chroms: list[str] = []

    unstrandedRNA: bool = False
    """Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand"""

    extra_config: MantaExtraConfig = MantaExtraConfig()
    extended_options: MantaExtendedOptions = MantaExtendedOptions()


class Delly2(SnappyModel):
    path_exclude_tsv: str | None = None
    max_threads: int = 16


class SomaticWgsSvCalling(SnappyStepModel, validators.ToolsMixin):
    path_ngs_mapping: str = "../ngs_mapping"
    tools: Annotated[list[Tool], EnumField(Tool, [Tool.manta], min_length=1)]

    manta: Manta | None = None

    delly2: Delly2 | None = None
