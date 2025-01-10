import enum
import re
from typing import Annotated
from pydantic import Field, model_validator

from snappy_pipeline.models import SnappyModel


class FormatAnnotation(enum.StrEnum):
    AD = "FORMAT/AD"
    ADF = "FORMAT/ADF"
    ADR = "FORMAT/ADR"
    DP = "FORMAT/DP"
    NMBZ = "FORMAT/NMBZ"
    QS = "FORMAT/QS"
    SP = "FORMAT/SP"
    SCR = "FORMAT/SCR"


class InfoAnnotation(enum.StrEnum):
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


class Config(enum.strEnum):
    BGI = "bgi"
    ILLUMINA = "illumina"
    ILLUMINA_1_18 = "illumina-1.18"
    ONT = "ont"
    ONT_SUP = "ont-sup"
    PACBIO_CCS = "pacbio-ccs"
    PACBIO_CCS_1_18 = "pacbio-css-1.18"
    ULTIMA = "ultima"
    ONE_TWELVE = "1.12"


class BAQAction(enum.StrEnum):
    DEFAULT = "default"
    DISABLE = "disable"
    REDO = "redo"


class BAQ(SnappyModel):
    action: BAQAction = BAQAction.DEFAULT
    max_read_len: int = 500
    """Maximum length of read to pass to BAQ algorithm"""


class IndelAmbiguousReadAction(enum.StrEnum):
    DROP = "drop"
    INCAD = "incAD"
    INCAD0 = "incAD0"


class IndelCallingModel(enum.StrEnum):
    DISABLED = "no_indels_cns"
    CNS = "indels_cns"
    TWO_ZERO = "indels-2.0"


class Indel(SnappyModel):
    enable: bool = True
    """Perform indel calling"""

    ext_prob: int = 20
    """Phred-scaled gap extension seq error probability"""
    gap_frac: float = 0.05
    """Minimum fraction of gapped reads"""
    tandem_qual: int = 500
    """Coefficient for homopolymer errors"""
    max_idepth: int = 250
    """Maximum per-file depth for INDEL calling"""
    min_ireads: int = 2
    """Minimum number gapped reads for indel candidates"""
    open_prob: int = 40
    """Phred-scaled gap open seq error probability"""
    ambig_reads: IndelAmbiguousReadAction = IndelAmbiguousReadAction.DROP
    """What to do with ambiguous indel reads"""
    indel_bias: float = 1.0
    """Raise to favour recall over precision"""
    del_bias: float = 0.0
    """Relative likelihood of insertion to deletion"""
    score_vs_ref: float = 0.0
    """Ratio of score vs ref (1) or 2nd-best allele (0)"""
    indel_size: int = 110
    """Approximate maximum indel size considered"""
    seqq_offset: int = 120
    """Indel-cns tuning for indel seq-qual scores"""
    poly_mqual: bool = False
    """(Edlib mode) Use minimum quality within homopolymers"""


class GermlineVariantCalling(SnappyModel):
    min_variant_depth: int = 20
    """Minimum read depth for a SNV to be displayed in the b-allele frequency plot."""
    zygosity_freq: float | None = None
    """Ignore VCF's genotypes (GT field) and instead infer zygosity from allele frequencies."""

    illumina_1_3: bool = False
    """Quality is in the Illumina-1.3+ encoding"""
    count_orphans: bool = False
    """Include anomalous read pairs, with flag PAIRED but not PROPER_PAIR set"""
    exclude_flags_any: int | None = None
    """Exclude reads with any of those flags set"""
    include_flags_any: int | None = None
    """Include reads with any of those flags unset"""
    exclude_flags_all: int | None = None
    """Exclude reads with all of those flags set"""
    include_flags_all: int | None = None
    """Include reads with all of those flags unset"""
    max_coverage: int = 250
    """Maximum coverage"""
    min_map_qual: int = 0
    """Minimum map quality"""
    min_base_qual: int = 1
    """Minimum base quality"""
    max_base_qual: int = 60
    """Maximum base quality"""
    delta_base_qual: int = 30
    """Maximum base quality"""

    config: Config | None = None
    """Specify platform profile"""
    baq: BAQ = BAQ()
    indels: Indel = Indel()

    include_annotations: list[FormatAnnotation | InfoAnnotation] = []
    exclude_annotations: list[FormatAnnotation | InfoAnnotation] = []

    min_genotype_qual: int = 20
    """Minimum genotype quality to consider germline variant"""
    min_map_qual: int = 35
    """Minimum mapping quality allowed during pileup"""

    path_regions: Annotated[
        str | None,
        Field(
            examples=[
                "/data/cephfs-1/work/groups/cubi/projects/biotools/GATK_Best_Practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/GATK_Best_Practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz",
                "/data/cephfs-1/work/groups/cubi/projects/biotools/DRAGEN/hg38/hg38_1000G_phase1.snps.high_confidence.vcf.gz",
            ]
        ),
    ] = None
    """
    Path to the regions where to perform germline variant calling.
    """

    path_read_groups: str | None = None
    """Select or exclude read groups listed in the file"""


class BcfToolsOutputFormat(enum.StrEnum):
    UNCOMPRESSED_BCF = "u"
    UNCOMPRESSED_VCF = "v"
    COMPRESSED_BCF = "b"
    COMPRESSED_VCF = "z"
    VCF_ZERO = "z0"
    VCF_ONE = "z1"
    VCF_TWO = "z2"
    VCF_THREE = "z3"
    VCF_FOUR = "z4"
    VCF_FIVE = "z5"
    VCF_SIX = "z6"
    VCF_SEVEN = "z7"
    VCF_EIGHT = "z8"
    VCF_NINE = "z9"
    BCF_ZERO = "b0"
    BCF_ONE = "b1"
    BCF_TWO = "b2"
    BCF_THREE = "b3"
    BCF_FOUR = "b4"
    BCF_FIVE = "b5"
    BCF_SIX = "b6"
    BCF_SEVEN = "b7"
    BCF_EIGHT = "b8"
    BCF_NINE = "b9"


class BcfToolsIndexFormat(enum.StrEnum):
    TBI = "tbi"
    CSI = "csi"


class BcfToolsOverlap(enum.StrEnum):
    POS = "0"
    RECORD = "1"
    VARIANT = "2"


class BcfToolsVariantType(enum.StrEnum):
    SNPS = "snps"
    INDELS = "indels"
    BOTH = "both"
    ALL = "all"

    SOME = "some"
    EXACT = "exact"
    NONE = "none"
    """Used in merge"""
    SNP_INS_DEL = "snp-ins-del"
    """Used in merge"""
    ID = "id"
    """Used in merge"""
    STAR = "*"
    """Used in merge"""
    STAR_STAR = "**"
    """Used in merge"""
    MNPS = "mnps"
    """Used in view"""
    REF = "ref"
    """Used in view"""
    BND = "bnd"
    """Used in view"""
    OTHER = "other"
    """Used in view"""


class BcfToolsMergeLogic(enum.StrEnum):
    PASS = "x"
    APPLY = "+"


class BcfToolsMergegVCFBlocks(enum.StrEnum):
    VOID = "-"
    REF = "REF.FA"


class BcfToolsNormCheckRef(enum.StrEnum):
    EXIT = "e"
    WARN = "w"
    EXCLUDE = "x"
    SET = "s"


class BcfToolsNormAtomOverlap(enum.StrEnum):
    STAR = "*"
    MISSING = "."


class BcfToolsNormMultiallelics(enum.StrEnum):
    SPLIT_SNPS = "+snps"
    SPLIT_INDELS = "+indels"
    SPLIT_BOTH = "+both"
    SPLIT_ANY = "+any"
    JOIN_SNPS = "-snps"
    JOIN_INDELS = "-indels"
    JOIN_BOTH = "-both"
    JOIN_ANY = "-any"


class BcfToolsNormMultiOverlaps(enum.StrEnum):
    REF = "0"
    MISSING = "."


class BcfToolsNormSortOrder(enum.StrEnum):
    CHR_POS = "chr_pos"
    LEX = "lex"


class BcfToolsViewGenotype(enum.StrEnum):
    HOM = "hom"
    HET = "het"
    MISSING = "miss"


class BcfToolsCallPloidy(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    X = "X"
    Y = "Y"
    HAPLOID = "1"
    DIPLOID = "2"


class BcfToolsCallAnnotation(enum.StrEnum):
    STRAND_BIAS_P_VALUES = "INFO/PV4"
    GENOTYPE_QUALITY = "FORMAT/GQ"
    GENOTYPE_POSTERIOR_PROBABILITIES = "FORMAT/GP"


class BcfToolsCallConstrain(enum.StrEnum):
    ALLELES = "alleles"
    TRIO = "trio"


class BcfToolsCsqPhase(enum.StrEnum):
    ASIS = "a"
    MERGE = "m"
    REQUIRE = "r"
    CREATE_NON_REF = "R"
    SKIP = "s"


class BcfToolsCsqChrNames(enum.StrEnum):
    PASS = "0"
    UNIFY = "1"


class BcfToolsFilterMode(enum.StrEnum):
    ADD = "+"
    RESET = "x"


class BcfToolsMpileupAmbiguousIndel(enum.StrEnum):
    DROP = "drop"
    INC_AD = "incAD"
    INC_AD0 = "incAD0"


class BcfToolsMixinOutput:
    output_type: BcfToolsOutputFormat = BcfToolsOutputFormat.UNCOMPRESSED_VCF
    """u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]"""
    write_index: BcfToolsIndexFormat | None = None
    """Automatically index the output files [off]"""
    # output: str | None = None  # FILE
    # """Write output to a file [standard output]"""


class BcfToolsMixinLocus:
    regions: list[str] = []
    """Restrict to comma-separated list of regions"""
    regions_overlap: BcfToolsOverlap = BcfToolsOverlap.RECORD
    """Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]"""
    targets: list[str] = []
    """Similar to -r but streams rather than index-jumps"""
    targets_overlap: BcfToolsOverlap = BcfToolsOverlap.POS
    """Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]"""
    # regions_file: str | None = None  # FILE
    # """Restrict to regions listed in FILE"""
    # targets_file: str | None = None
    # """Similar to -R but streams rather than index-jumps"""


class BcfToolsMixinFilter:
    exclude: str | None = None
    """EXPR: Exclude sites for which the expression is true (see man page for details)"""
    include: str | None = None
    """EXPR: Select sites for which the expression is true (see man page for details)"""


class BcfToolsMixinSample:
    samples: list[str] = []
    """Comma separated list of samples to annotate (or exclude with "^" prefix)"""
    # samples_file: [^]FILE
    # """File of samples to annotate (or exclude with "^" prefix)"""


BcfToolsAnnotateMinOverlap = enum.StrEnum(
    "BcfToolsAnnotateMinOverlap",
    {"ANN": "ANN", "VCF": "VCF", "RECIPROCAL": "ANN:VCF"},
)


class BcfToolsAnnotate(
    SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample
):
    columns: list[str] = []
    """List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details"""
    force: bool = False
    """Continue despite parsing error (at your own risk!)"""
    header_line: str | None = None
    """Header line which should be appended to the VCF header, can be given multiple times"""
    set_id: str | None = None
    """[+]FORMAT: Set ID column using a `bcftools query`-like expression, see man page for details"""
    keep_sites: bool = False
    """Leave -i/-e sites unchanged instead of discarding them"""
    merge_logic: str | None = None
    """TAG:TYPE: Merge logic for multiple overlapping regions (see man page for details), EXPERIMENTAL"""
    mark_sites: str | None = None
    """[+-]TAG: Add INFO/TAG flag to sites which are ("+") or are not ("-") listed in the -a file"""
    min_overlap: BcfToolsAnnotateMinOverlap = BcfToolsAnnotateMinOverlap.RECIPROCAL
    """Required overlap as a fraction of variant in the -a file (ANN), the VCF (:VCF), or reciprocal (ANN:VCF)"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    pair_logic: BcfToolsVariantType = BcfToolsVariantType.SOME
    """Matching records by <snps|indels|both|all|some|exact>, see man page for details [some]"""
    single_overlaps: bool = False
    """Keep memory low by avoiding complexities arising from handling multiple overlapping intervals"""
    remove: list[str] = []
    """List of annotations (e.g. ID,INFO/DP,FORMAT/DP,FILTER) to remove (or keep with "^" prefix). See man page for details"""

    @model_validator(mode="after")
    def ensure_valid_variant_type(self):
        if self.pair_logic not in (
            BcfToolsVariantType.SNPS,
            BcfToolsVariantType.INDELS,
            BcfToolsVariantType.BOTH,
            BcfToolsVariantType.ALL,
            BcfToolsVariantType.SOME,
            BcfToolsVariantType.EXACT,
        ):
            raise ValueError(f"Illegal variant type '{self.pair_logic}'")
        return self

    # annotations: str | None = None
    # """VCF file or tabix-indexed FILE with annotations: CHR\tPOS[\tVALUE]+"""
    # columns_file: str | None = None
    # """Read -c columns from FILE, one name per row, with optional --merge-logic TYPE: NAME[ TYPE]"""
    # header_lines: str | None = None
    # """Lines which should be appended to the VCF header"""
    # rename_annots: str | None = None  # FILE
    # """Rename annotations: TYPE/old\tnew, where TYPE is one of FILTER,INFO,FORMAT"""
    # rename_chrs: str | None = None  # FILE
    # """Rename sequences according to the mapping: old\tnew"""


class BcfToolsConcat(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinLocus):
    allow_overlaps: bool = False
    """First coordinate of the next file can precede last record of the current file."""
    compact_PS: bool = False
    """Do not output PS tag at each site, only at the start of a new phase set block."""
    rm_dups: BcfToolsVariantType | None = None
    """Output duplicate records present in multiple files only once: <snps|indels|both|all|exact>"""
    drop_genotypes: bool = False
    """Drop individual genotype information."""
    ligate: bool = False
    """Ligate phased VCFs by matching phase at overlapping haplotypes"""
    ligate_force: bool = False
    """Ligate even non-overlapping chunks, keep all sites"""
    ligate_warn: bool = False
    """Drop sites in imperfect overlaps"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    naive: bool = False
    """Concatenate files without recompression, a header check compatibility is performed"""
    naive_force: bool = False
    """Same as --naive, but header compatibility is not checked. Dangerous, use with caution."""
    min_PQ: int = 30
    """Break phase set if phasing quality is lower than <int> [30]"""
    verbose: 0 | 1
    """Set verbosity level [1]"""

    @model_validator(mode="after")
    def ensure_valid_variant_type(self):
        if self.rm_dups not in (
            BcfToolsVariantType.SNPS,
            BcfToolsVariantType.INDELS,
            BcfToolsVariantType.BOTH,
            BcfToolsVariantType.ALL,
            BcfToolsVariantType.EXACT,
        ):
            raise ValueError(f"Illegal variant type '{self.rm_dups}'")
        return self

    # remove_duplicates: bool = False
    # """Alias for -d exact"""
    # file_list: str | None = None
    # """Read the list of files from a file."""


class BcfToolsConvert(
    SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample
):
    no_version: bool = False
    """Do not append version and command line to the header"""

    three_n_plus_6: bool = False
    """Use 3*N+6 column format instead of the old 3*N+5 column format"""
    tag: str = "GT"
    """Tag to take values for .gen file: GT,PL,GL,GP [GT]"""
    chrom: bool = False
    """Output chromosome in first column instead of CHROM:POS_REF_ALT"""
    keep_duplicates: bool = False
    """Keep duplicate positions"""
    vcf_ids: bool = False
    """Output VCF IDs in second column instead of CHROM:POS_REF_ALT"""
    gvcf2vcf: bool = False
    """Expand gVCF reference blocks"""

    haploid2diploid: bool = False
    """Convert haploid genotypes to diploid homozygotes"""
    vcf_ids: bool = False
    """Output VCF IDs instead of CHROM:POS_REF_ALT"""

    columns: str = "ID,CHROM,POS,AA"
    """Columns of the input tsv file, see man page for details [ID,CHROM,POS,AA]"""

    # gensample2vcf: <PREFIX>|<GEN-FILE>,<SAMPLE-FILE>
    # gensample: <PREFIX>|<GEN-FILE>,<SAMPLE-FILE>
    # sex: str | None = None
    # """Output sex column in the sample-file, input format is: Sample\t[MF]"""
    # fasta_ref: str | None = None
    # """Reference sequence in fasta format"""
    # hapsample2vcf: <PREFIX>|<HAP-FILE>,<SAMPLE-FILE>
    # hapsample: <PREFIX>|<HAP-FILE>,<SAMPLE-FILE>
    # haplegendsample2vcf: <PREFIX>|<HAP-FILE>,<LEGEND-FILE>,<SAMPLE-FILE>
    # haplegendsample: <PREFIX>|<HAP-FILE>,<LEGEND-FILE>,<SAMPLE-FILE>
    # tsv2vcf: FILE


class BcfToolsHead(SnappyModel, BcfToolsMixinSample):
    headers: int | None = None
    """Display INT header lines [all]"""
    records: int | None = None
    """Display INT variant record lines [none]"""


class BcfToolsIsec(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus):
    collapse: BcfToolsVariantType = BcfToolsVariantType.NONE
    """Treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]"""
    complement: bool = False
    """Output positions present only in the first file but missing in the others"""
    apply_filters: list[str] = []
    """Require at least one of the listed FILTER strings (e.g. "PASS,.")"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    nfiles: str | None = None
    """Output positions present in this many (=), this many or more (+), this many or fewer (-), the exact (~) files"""

    @model_validator(mode="after")
    def ensure_valid_variant_type(self):
        if self.collapse not in (
            BcfToolsVariantType.SNPS,
            BcfToolsVariantType.INDELS,
            BcfToolsVariantType.BOTH,
            BcfToolsVariantType.ALL,
            BcfToolsVariantType.SOME,
            BcfToolsVariantType.NONE,
        ):
            raise ValueError(f"Illegal variant type '{self.collapse}'")
        return self

    # file_list: str | None = None
    # """Read the input file names from the file"""
    # prefix: DIR
    # """If given, subset each of the input files accordingly, see also -w"""
    # threads: INT
    # """Use multithreading with <int> worker threads [0]"""
    # write: list = []
    # """List of files to write with -p given as 1-based indexes. By default, all files are written"""


BcfToolsMergeFilterLogic = enum.StrEnum("BcfToolsMergeFilterLogic", {"PASS": "x", "APPLY": "+"})


class BcfToolsMerge(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinLocus):
    force_no_index: bool = False
    """Merge unindexed files, synonymous to --no-index"""
    force_samples: bool = False
    """Resolve duplicate sample names"""
    force_single: bool = False
    """Run even if there is only one file on input"""
    print_header: bool = False
    """Print only the merged header and exit"""
    missing_to_ref: bool = False
    """Assume genotypes at missing sites are 0/0"""
    apply_filters: list[str] = []
    """Require at least one of the listed FILTER strings (e.g. "PASS,.")"""
    filter_logic: BcfToolsMergeFilterLogic = BcfToolsMergeFilterLogic.APPLY
    """Remove filters if some input is PASS ("x"), or apply all filters ("+") [+]"""
    gvcf: BcfToolsMergegVCFBlocks | None = None
    """Merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,MIN_DP:min,I16:sum,IDV:max,IMF:max -M PL:max,AD:0"""
    info_rules: list[str] = ["DP:sum", "DP4:sum"]
    """TAG:METHOD,...: Rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [DP:sum,DP4:sum]"""
    local_alleles: int = 0
    """If more than INT alt alleles are encountered, drop FMT/PL and output LAA+LPL instead; 0=unlimited [0]"""
    merge: BcfToolsVariantType = BcfToolsVariantType.BOTH
    """Allow multiallelic records for snps,indels,both,snp-ins-del,all,none,id,*,**; see man page for details [both]"""
    missing_rules: list[str] = ["."]
    """Rules for replacing missing values in numeric vectors (.,0,max) when unknown allele <*> is not present [.]"""
    no_index: bool = False
    """Merge unindexed files, the same chromosomal order is required and -r/-R are not allowed"""
    no_version: bool = False
    """Do not append version and command line to the header"""

    @model_validator(mode="after")
    def ensure_valid_variant_type(self):
        if self.collapse not in (
            BcfToolsVariantType.SNPS,
            BcfToolsVariantType.INDELS,
            BcfToolsVariantType.BOTH,
            BcfToolsVariantType.ALL,
            BcfToolsVariantType.NONE,
            BcfToolsVariantType.SNP_INS_DEL,
            BcfToolsVariantType.ID,
            BcfToolsVariantType.STAR,
            BcfToolsVariantType.STAR_STAR,
        ):
            raise ValueError(f"Illegal variant type '{self.collapse}'")
        return self

    @model_validator(mode="after")
    def ensure_valid_info_rules(self):
        INFO_RULE_PATTERN = re.compile("^([A-Za-z_][0-9A-Za-z_.]*|1000G):(sum|avg|min|max|join)$")
        if len(self.info_rules) == 0:
            raise ValueError("Missing info-rules parameter")
        if len(self.info_rules) != 1 or self.info_rules[0] != "-":
            for rule in self.info_rules:
                if not INFO_RULE_PATTERN.match(rule):
                    raise ValueError(f"Illegal info rule '{rule}'")
        return self

    @model_validator(mode="after")
    def ensure_valid_missing_rules(self):
        MISSING_RULE_PATTERN = re.compile("^([A-Za-z ][0-9A-Za-z .]*|1000G):(.|0|max)$")
        if len(self.missing_rules) == 0:
            raise ValueError("Missing missing-rules parameter")
        for rule in self.missing_rules:
            if not MISSING_RULE_PATTERN.match(rule):
                raise ValueError(f"Illegal missing rule '{rule}'")
        return self

    # use_header: str | None = None
    # """Use the provided header"""
    # file_list: str | None = None
    # """Read file names from the file"""


class BcfToolsNorm(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus):
    atomize: bool = False
    """Decompose complex variants (e.g. MNVs become consecutive SNVs)"""
    atom_overlaps: BcfToolsNormAtomOverlap = BcfToolsNormAtomOverlap.STAR
    """Use the star allele (*) for overlapping alleles or set to missing (.) [*]"""
    check_ref: BcfToolsNormCheckRef = BcfToolsNormCheckRef.EXIT
    """Check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [e]"""
    remove_duplicates: bool = False
    """Remove duplicate lines of the same type."""
    rm_dup: BcfToolsVariantType | None = None
    """Remove duplicate snps|indels|both|all|exact"""
    force: bool = False
    """Try to proceed even if malformed tags are encountered. Experimental, use at your own risk"""
    keep_sum: list[str] = []
    """Keep vector sum constant when splitting multiallelics (see github issue #360)"""
    multiallelics: BcfToolsNormMultiallelics | None = None
    """-|+TYPE: Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]"""
    multi_overlaps: BcfToolsNormMultiOverlaps = BcfToolsNormMultiOverlaps.REF
    """Fill in the reference (0) or missing (.) allele when splitting multiallelics [0]"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    do_not_normalize: bool = False
    """Do not normalize indels (with -m or -c s)"""
    old_rec_tag: str | None = None
    """Annotate modified records with INFO/STR indicating the original variant"""
    strict_filter: bool = False
    """When merging (-m+), merged site is PASS only if all sites being merged PASS"""
    sort: BcfToolsNormSortOrder = BcfToolsNormSortOrder.CHR_POS
    """Sort order: chr_pos,lex [chr_pos]"""
    verbose: int = 1
    """Verbosity level (0-2) of GFF parsing [1]"""
    site_win: int = 1000
    """Buffer for sorting lines which changed position during realignment [1000]"""

    @model_validator(mode="after")
    def ensure_valid_remove_duplicates(self):
        if self.rm_dup is not None and self.rm_dup not in (
            BcfToolsVariantType.SNPS,
            BcfToolsVariantType.INDELS,
            BcfToolsVariantType.BOTH,
            BcfToolsVariantType.ALL,
            BcfToolsVariantType.EXACT,
        ):
            raise ValueError(f"Illegal variant type '{self.rm_dup}'")
        return self

    # fasta_ref: str | None = None
    # """Reference sequence"""
    # gff_annot: str | None = None
    # """Follow HGVS 3'rule and right-align variants in transcripts on the forward strand"""


class BcfToolsQuery(SnappyModel, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample):
    force_samples: bool = False
    """Only warn about unknown subset samples"""
    print_filtered: str | None = None
    """Output STR for samples failing the -i/-e filtering expression"""
    format: str
    """See man page for details"""
    print_header: bool = False
    """Print header, -HH to omit column indices"""
    list_samples: bool = False
    """Print the list of samples and exit"""
    disable_automatic_newline: bool = False
    """Disable automatic addition of newline character when not present"""
    allow_undef_tags: bool = False
    """Print "." for undefined tags"""

    # vcf_list: str | None = None
    # """Process multiple VCFs listed in the file"""


# BcfToolsReheader(SnappyModel):
# fai: str | None = None
# """update sequences and their lengths from the .fai file"""
# header: str | None = None
# """new header"""
# output: str | None = None
# """write output to a file [standard output]"""
# samples: str | None = None
# """new sample names"""
# temp_prefix: PATH
# """ignored; was template for temporary file name"""
# threads: INT
# """use multithreading with <int> worker threads (BCF only) [0]"""


class BcfToolsSort(SnappyModel, BcfToolsMixinOutput):
    pass


class BcfToolsView(
    SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample
):
    drop_genotypes: bool = False
    """Drop individual genotype information (after subsetting if -s option set)"""
    header_only: bool = False
    """Print only the header in VCF output (equivalent to bcftools head)"""
    no_header: bool = False
    """Suppress the header in VCF output"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    trim_unseen_allele: bool = False
    """Remove '<*>' or '<NON_REF>' at variant (-A) or at all (-AA) sites"""
    trim_alt_alleles: bool = False
    """Trim ALT alleles not seen in the genotype fields (or their subset with -s/-S)"""
    no_update: bool = False
    """Do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)"""
    force_samples: bool = False
    """Only warn about unknown subset samples"""
    min_ac: str | None = None
    """Minimum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]"""
    max_ac: str | None = None
    """Maximum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]"""
    apply_filters: list[str] = []
    """Require at least one of the listed FILTER strings (e.g. "PASS,.")"""
    genotype: BcfToolsViewGenotype | None = None
    """Require one or more hom/het/missing genotype or, if prefixed with "^", exclude such sites"""
    known: bool = False
    """Select known sites only (ID is not '.')"""
    novel: bool = False
    """Select novel sites only (ID is is '.')"""
    min_allele: int | None = None
    """Minimum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)"""
    max_allele: int | None = None
    """Maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)"""
    phased: bool = False
    """Select sites where all samples are phased"""
    exclude_phased: bool = False
    """Exclude sites where all samples are phased"""
    min_af: str | None = None
    """Minimum frequency for non-reference (nref), 1st alternate (alt1), least frequen (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]"""
    max_af: str | None = None
    """Maximum frequency for non-reference (nref), 1st alternate (alt1), least frequen (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]"""
    uncalled: bool = False
    """Select sites without a called genotype"""
    exclude_uncalled: bool = False
    """Exclude sites without a called genotype"""
    types: list[BcfToolsVariantType] | None = None
    """Select comma-separated list of variant types: snps,indels,mnps,ref,bnd,other [null]"""
    exclude_types: list[BcfToolsVariantType] | None = None
    """Exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other [null]"""
    private: bool = False
    """Select sites sites where the non-reference alleles are exclusive (private) to the subset samples"""
    exclude_private: bool = False
    """Exclude sites sites where the non-reference alleles are exclusive (private) to the subset samples"""
    write_index: BcfToolsIndexFormat | None = None
    """Automatically index the output files [off]"""

    @model_validator(mode="after")
    def ensure_valid_types(self):
        for types in (self.types, self.exclude_types):
            if types is not None:
                for one_type in types:
                    if one_type is not None and one_type not in (
                        BcfToolsVariantType.SNPS,
                        BcfToolsVariantType.INDELS,
                        BcfToolsVariantType.MNPS,
                        BcfToolsVariantType.REF,
                        BcfToolsVariantType.BND,
                        BcfToolsVariantType.OTHER,
                    ):
                        raise ValueError(f"Illegal variant type '{one_type}'")
        return self

    @model_validator(mode="after")
    def ensure_valid_counts(self):
        COUNT_PATTERN = re.compile("^([0-9]+)(:(nref|alt1|minor|major|nonmajor))?$")
        for counts in (self.min_ac, self.max_ac):
            if counts is not None and not COUNT_PATTERN.match(counts):
                raise ValueError(f"Illegal min/max counts '{counts}'")
        return self

    @model_validator(mode="after")
    def ensure_valid_frequency(self):
        FREQUENCY_PATTERN = re.compile(
            r"^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([EeDd][+-]?[0-9]+)?(:nref|alt1|minor|major|nonmajor))?$"
        )
        for freq in (self.min_af, self.max_af):
            if freq is not None and not FREQUENCY_PATTERN.match(freq):
                raise ValueError(f"Illegal min/max frequency '{freq}'")
        return self


class BcfToolsCall(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinLocus, BcfToolsMixinSample):
    no_version: bool = False
    """Do not append version and command line to the header"""
    ploidy: BcfToolsCallPloidy | None = None
    """Predefined ploidy, 'list' to print available settings, append '?' for details [2]"""

    keep_alts: bool = False
    """Keep all possible alternate alleles at variant sites"""
    keep_unseen_allele: bool = False
    """Keep the unobserved allele <*> or <NON_REF>"""
    annotate: list[BcfToolsCallAnnotation] = []
    """Optional tags to output (lowercase allowed); '?' to list available tags"""
    prior_freqs: tuple[str, str] | None = None
    """Use prior allele frequencies AN & AC, determined from these pre-filled tags"""
    group_samples_tag: str | None = None
    """The tag to use with -G, by default FORMAT/QS and FORMAT/AD are checked automatically"""
    gvcf: list[int] = []
    """Group non-variant sites into gVCF blocks by minimum per-sample DP"""
    insert_missed: bool = False
    """Output also sites missed by mpileup but present in -T"""
    keep_masked_ref: bool = False
    """Keep sites with masked reference allele (REF=N)"""
    skip_variants: BcfToolsVariantType | None = None
    """Skip indels/snps"""
    variants_only: bool = False
    """Output variant sites only"""
    write_index: BcfToolsIndexFormat | None = None
    """Automatically index the output files [off]"""

    consensus_caller: bool = False
    """The original calling method (conflicts with -m)"""
    constrain: BcfToolsCallConstrain | None = None
    """One of: alleles, trio (see manual)"""
    multiallelic_caller: bool = False
    """Alternative model for multiallelic and rare-variant calling (conflicts with -c)"""
    novel_rate: tuple[float, float, float] = (1e-8, 1e-9, 1e-9)
    """Likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]"""
    pval_threshold: float = 0.5
    """Variant if P(ref|D)<FLOAT with -c [0.5]"""
    prior: float = 1.1e-3
    """Mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]"""

    @model_validator(mode="after")
    def ensure_valid_types(self):
        if self.skip_variants is not None and self.skip_variants not in (
            BcfToolsVariantType.SNPS,
            BcfToolsVariantType.INDELS,
        ):
            raise ValueError(f"Illegal variant type '{self.skip_variants}'")
        return self

    @model_validator(mode="after")
    def avoid_caller_conflict(self):
        if self.consensus_caller and self.multiallelic_caller:
            raise ValueError("Consensus & multiallelic callers cannot be simultaneously selected")

    # ploidy_file: str | None = None
    # """Space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY"""
    # group_samples: FILE|-
    # """Group samples by population (file with "sample\tgroup") or "-" for single-sample calling."""
    #                 This requires FORMAT/QS or other Number=R,Type=Integer tag such as FORMAT/AD


class BcfToolsConsensus(SnappyModel, BcfToolsMixinFilter, BcfToolsMixinSample):
    absent: str | None = None
    """Replace positions absent from VCF with CHAR"""
    haplotype: str | None
    """
    Choose which allele to use from the FORMAT/GT field, note the codes are case-insensitive:
        N: N={1,2,3,..} is the index of the allele from GT, regardless of phasing (e.g. "2")
        R: REF allele in het genotypes
        A: ALT allele
        I: IUPAC code for all genotypes
        LR,LA: longer allele and REF/ALT if equal length
        SR,SA: shorter allele and REF/ALT if equal length
        NpIu: index of the allele for phased and IUPAC code for unphased GTs (e.g. "2pIu")
    """
    iupac_codes: bool = False
    """Output IUPAC codes based on FORMAT/GT, use -s/-S to subset samples"""
    mark_del: str | None = None
    """Instead of removing sequence, insert character CHAR for deletions"""
    mark_ins: str | None = None
    """Highlight insertions in uppercase (uc), lowercase (lc), or use CHAR, leaving the rest as is"""
    mark_snv: str | None = None
    """Highlight substitutions in uppercase (uc), lowercase (lc), or use CHAR, leaving the rest as is"""
    mask_with: str | None = None
    """Replace with CHAR (skips overlapping variants); change to uppercase (uc) or lowercase (lc)"""
    missing: str | None = None
    """Output CHAR instead of skipping a missing genotype './.'"""

    @model_validator(mode="after")
    def ensure_one_char(self):
        for s in (self.absent, self.make_del, self.missing):
            if s is not None and len(s) != 1:
                raise ValueError(f"String {s} should be made of a single character")
        return self

    @model_validator(mode="after")
    def ensure_one_char_or_case(self):
        for s in (self.mark_ins, self.mark_snv, self.mask_with):
            if s is not None and len(s) != 1 and s not in ("uc", "lc"):
                raise ValueError(f"String {s} should be a single charaacter, or 'lc' or 'uc'")
        return self

    @model_validator(mode="after")
    def ensure_valid_haplotype(self):
        HAPLOTYPE_PATTERN = re.compile("^([0-9]+(pIu)?|[RAI]|[LS][RA])$")
        if self.haplotype is not None and not HAPLOTYPE_PATTERN.matches(self.haplotype):
            raise ValueError(f"Illegal haplotype definition {self.haplotype}")
        return self

    # chain: str | None = None
    # """Write a chain file for liftover"""
    # fasta_ref: str | None = None
    # """Reference sequence in fasta format"""
    # mask: str | None = None
    # """Replace regions according to the next --mask-with option. The default is --mask-with N"""


class BcfToolsCnv(SnappyModel, BcfToolsMixinLocus):
    control_sample: str | None = None
    """Optional control sample name to highlight differences"""
    plot_threshold: float | None = None
    """Plot aberrant chromosomes with quality at least FLOAT"""
    query_sample: str | None = None
    """Query samply name"""
    aberrant: tuple[float, float] = (1.0, 1.0)
    """Fraction of aberrant cells in query and control [1.0,1.0]"""
    BAF_weight: float = 1.0
    """Relative contribution from BAF [1]"""
    BAF_dev: tuple[float, float] = (0.04, 0.04)
    """Expected BAF deviation in query and control [0.04,0.04]"""
    err_prob: float = 1e-4
    """Uniform error probability [1e-4]"""
    LRR_dev: tuple[float, float] = (0.2, 0.2)
    """Expected LRR deviation [0.2,0.2]"""
    LRR_weight: float = 0.2
    """Relative contribution from LRR [0.2]"""
    LRR_smooth_win: int = 10
    """Window of LRR moving average smoothing [10]"""
    optimize: float = 1.0
    """Estimate fraction of aberrant cells down to FLOAT [1.0]"""
    same_prob: float = 0.5
    """Prior probability of -s/-c being the same [0.5]"""
    xy_prob: float = 1e-9
    """P(x|y) transition probability [1e-9]"""

    # AF_file: str | None = None
    # """Read allele frequencies from file (CHR\tPOS\tREF,ALT\tAF)"""
    # output_dir: PATH


class BcfToolsCsq(
    SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample
):
    trim_protein_seq: int | None = None
    """Abbreviate protein-changing predictions to max INT aminoacids"""
    custom_tag: str | None = None
    """Use this tag instead of the default BCSQ"""
    local_csq: bool = False
    """Localized predictions, consider only one VCF record at a time"""
    ncsq: int = 15
    """Maximum number of per-haplotype consequences to consider for each site [15]"""
    phase: BcfToolsCsqPhase = BcfToolsCsqPhase.REQUIRE
    """
    How to handle unphased heterozygous genotypes: [r]
        a: take GTs as is, create haplotypes regardless of phase (0/1 -> 0|1)
        m: merge *all* GTs into a single haplotype (0/1 -> 1, 1/2 -> 1)
        r: require phased GTs, throw an error on unphased het GTs
        R: create non-reference haplotypes if possible (0/1 -> 1|1, 1/2 -> 1|2)
        s: skip unphased hets
    """
    force: bool = False
    """Run even if some sanity checks fail"""
    unify_chr_names: BcfToolsCsqChrNames = BcfToolsCsqChrNames.UNIFY
    """Automatically unify chromosome naming (e.g. chrX vs X) in GFF, fasta, and VCF [1]"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    verbose: int = 1
    """Verbosity level 0-2 [1]"""

    # fasta_ref: str | None = None
    # """Reference file in fasta format"""
    # gff_annot: str | None = None
    # """GFF3 annotation file"""
    # dump_gff: FILE.gz
    # """Dump the parsed GFF file (for debugging purposes)"""


BcfToolsFilterFailedGenotype = enum.StrEnum(
    "BcfToolsFilterFailedGenotype", {"REF": "0", "MISSING": "."}
)


class BcfToolsFilter(
    SnappyModel, BcfToolsMixinOutput, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample
):
    SnpGap: str | None = None
    """Filter SNPs within <int> base pairs of an indel (the default) or any combination of indel,mnp,bnd,other,overlap"""
    IndelGap: int | None = None
    """Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass"""
    mask: list[str] = []
    """Soft filter regions, "^" to negate"""
    mask_overlap: BcfToolsOverlap = BcfToolsOverlap.RECORD
    """Mask if POS in the region (0), record overlaps (1), variant overlaps (2) [1]"""
    mode: BcfToolsFilterMode | None = None
    """"+": do not replace but add to existing FILTER; "x": reset filters at sites which pass"""
    no_version: bool = False
    """Do not append version and command line to the header"""
    soft_filter: str | None = None
    """Annotate FILTER column with <string> or unique filter name ("Filter%d") made up by the program ("+")"""
    set_GTs: BcfToolsFilterFailedGenotype = BcfToolsFilterFailedGenotype.REF
    """Set genotypes of failed samples to missing (.) or ref (0)"""

    @model_validator(mode="after")
    def ensure_valid_snp_gap(self):
        SNP_GAP_PATTERN = re.compile("^([0-9]+)(:(indel|mnp|bnd|other|overlap))?$")
        if self.SnpGap is not None and not SNP_GAP_PATTERN.matches(self.SnpGap):
            raise ValueError(f"Illegal SNP gap '{self.SnpGap}'")
        return self

    # mask_file: [^]FILE
    # """Soft filter regions listed in a file, "^" to negate"""


class BcfToolsGtcheck(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinLocus, BcfToolsMixinSample):
    distinctive_sites: list | None = None
    """
    Find sites that can distinguish between at least NUM sample pairs.
        NUM[,MEM[,TMP]]          If the number is smaller or equal to 1, it is interpreted as the fraction of pairs.
        The optional MEM string sets the maximum memory used for in-memory sorting [500M]
        and TMP is a prefix of temporary files used by external sorting [/tmp/bcftools.XXXXXX]
    """
    error_probability: int = 40
    """Phred-scaled probability of genotyping error, 0 for faster but less accurate results [40]"""
    homs_only: bool = False
    """Homozygous genotypes only, useful with low coverage data (requires -g)"""
    n_matches: int = 0
    """
    Print only top INT matches for each sample (sorted by average score), 0 for unlimited.
    Use negative value to sort by HWE probability rather than by discordance [0]
    """
    no_HWE_prob: bool = False
    """Disable calculation of HWE probability"""
    pairs: list[str] = []
    """Comma-separated sample pairs to compare (qry,gt[,qry,gt..] with -g or qry,qry[,qry,qry..] w/o)"""
    use: list[str] = ["PL", "GT"]
    """Which tag to use in the query file (TAG1) and the -g file (TAG2) [PL,GT]"""

    # genotypes: str | None = None  # FILE
    # """Genotypes to compare against"""
    # pairs_file: str | None = None
    # """File with tab-delimited sample pairs to compare (qry,gt with -g or qry,qry w/o)"""


class BcfToolsMpileup(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinLocus, BcfToolsMixinSample):
    illumina1_3: bool = False
    """Quality is in the Illumina-1.3+ encoding"""
    count_orphans: bool = False
    """Include anomalous read pairs, with flag PAIRED but not PROPER_PAIR set"""
    no_BAQ: bool = False
    """Disable BAQ (per-Base Alignment Quality)"""
    adjust_MQ: int = 0
    """Adjust mapping quality [0]"""
    full_BAQ: bool = False
    """Apply BAQ everywhere, not just in problematic regions"""
    max_depth: int = 250
    """Max raw per-file depth; avoids excessive memory usage [250]"""
    redo_BAQ: bool = False
    """Recalculate BAQ on the fly, ignore existing BQs"""
    no_reference: bool = False
    """Do not require fasta reference file"""
    min_MQ: int = 0
    """Skip alignments with mapQ smaller than INT [0]"""
    min_BQ: int = 1
    """Skip bases with baseQ/BAQ smaller than INT [1]"""
    max_BQ: int = 60
    """Limit baseQ/BAQ to no more than INT [60]"""
    delta_BQ: int = 30
    """Use neighbour_qual + INT if less than qual [30]"""
    ignore_RG: bool = False
    """Ignore RG tags (one BAM = one sample)"""
    skip_all_set: int | None = None
    """Skip reads with all of the bits set"""
    skip_any_set: int | None = 1796  # (4 + 256 + 512 + 1024)
    """Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP]"""
    skip_all_unset: int | None = None
    """Skip reads with all of the bits unset"""
    skip_any_unset: int | None = None
    """Skip reads with any of the bits unset"""
    ignore_overlaps: bool = False
    """Disable read-pair overlap detection"""
    seed: int = 0
    """Random number seed used for sampling deep regions [0]"""
    annotate: list = []
    r"""Optional tags to output; '\?' to list available tags []"""
    gvcf: list[int] = []
    """Group non-variant sites into gVCF blocks according to minimum per-sample DP"""
    no_version: bool = False
    """Do not append version and command line to the header"""

    config: str | None = None
    """Specify platform profile (use "-X list" for details)"""
    ext_prob: int = 20
    """Phred-scaled gap extension seq error probability [20]"""
    gap_frac: float = 0.05
    """Minimum fraction of gapped reads [0.05]"""
    tandem_qual: int = 500
    """Coefficient for homopolymer errors [500]"""
    skip_indels: bool = False
    """Do not perform indel calling"""
    max_idepth: int = 250
    """Maximum per-file depth for INDEL calling [250]"""
    min_ireads: int = 2
    """Minimum number gapped reads for indel candidates [2]"""
    max_read_len: int = 500
    """Maximum length of read to pass to BAQ algorithm [500]"""
    open_prob: int = 40
    """Phred-scaled gap open seq error probability [40]"""
    per_sample_mF: bool = False
    """Apply -m and -F per-sample for increased sensitivity"""
    platforms: list[str] = ["all"]
    """Comma separated list of platforms for indels [all]"""
    ambig_reads: BcfToolsMpileupAmbiguousIndel = BcfToolsMpileupAmbiguousIndel.DROP
    """What to do with ambiguous indel reads: drop,incAD,incAD0 [drop]"""
    indel_bias: float = 1.0
    """Raise to favour recall over precision [1.00]"""
    del_bias: float = 0.0
    """Relative likelihood of insertion to deletion [0.00]"""
    score_vs_ref: float = 0.0
    """Ratio of score vs ref (1) or 2nd-best allele (0) [0.00]"""
    indel_size: int = 110
    """Approximate maximum indel size considered [110]"""
    indels_2: bool = False
    """New EXPERIMENTAL indel calling model (diploid reference consensus)"""
    indels_cns: bool = False
    """New EXPERIMENTAL indel calling model with edlib"""
    seqq_offset: int = 120
    """Indel-cns tuning for indel seq-qual scores [120]"""
    no_indels_cns: bool = False
    """Disable CNS mode, to use after a -X profile"""
    poly_mqual: bool = False
    """(Edlib mode) Use minimum quality within homopolymers"""

    # bam_list: str | None = None
    # """List of input BAM filenames, one per line"""
    # fasta_ref: str | None = None
    # """Faidx indexed reference sequence file"""
    # read_groups: str | None = None
    # """Select or exclude read groups listed in the file"""


class BcfToolsPolysomy(SnappyModel, BcfToolsMixinLocus):
    verbose: bool = False
    """Verbosity level"""
    peak_size: float = 0.1
    """Minimum peak size (0-1, larger is stricter) [0.1]"""
    cn_penalty: float = 0.7
    """Penalty for increasing CN (0-1, larger is stricter) [0.7]"""
    fit_th: float = 3.3
    """Goodness of fit threshold (>0, smaller is stricter) [3.3]"""
    include_aa: bool = False
    """Include the AA peak in CN2 and CN3 evaluation"""
    min_fraction: float = 0.1
    """Minimum distinguishable fraction of aberrant cells [0.1]"""
    peak_symmetry: float = 0.5
    """Peak symmetry threshold (0-1, larger is stricter) [0.5]"""

    # output_dir: PATH


class BcfToolsRoh(SnappyModel, BcfToolsMixinOutput, BcfToolsMixinLocus, BcfToolsMixinSample):
    AF_dflt: float | None = None
    """if AF is not known, use this allele frequency [skip]"""
    AF_tag: str | None = None
    """use TAG for allele frequency"""
    buffer_size: list[int] = [0]
    """
    Buffer size and the number of overlapping sites, 0 for unlimited [0]
    If the first number is negative, it is interpreted as the maximum memory to use, in MB. The default overlap is set to roughly 1% of the buffer size.
    """
    GTs_only: float = 30.0
    """use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes. Safe value to use is 30 to account for GT errors."""
    ignore_homref: bool = False
    """skip hom-ref genotypes (0/0)"""
    include_noalt: bool = False
    """include sites with no ALT allele (ignored by default)"""
    skip_indels: bool = False
    """skip indels as their genotypes are enriched for errors"""
    rec_rate: float | None = None
    """constant recombination rate per bp"""
    hw_to_az: float = 6.7e-8
    """P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state [6.7e-8]"""
    az_to_hw: float = 5e-9
    """P(HW|AZ) transition probability from AZ to HW state [5e-9]"""
    viterbi_training: float = 1e-10
    """estimate HMM parameters, <float> is the convergence threshold, e.g. 1e-10 (experimental)"""

    # AF_file: <file>
    # """read allele frequencies from file (CHR\tPOS\tREF,ALT\tAF)"""
    # genetic_map: <file>
    # """genetic map in IMPUTE2 format, single file or mask, where string '{CHROM}' is replaced with chromosome name"""
    # estimate_AF: [TAG],<file>
    # """
    # Estimate AF from FORMAT/TAG (GT or PL) of all samples ("-") or samples listed in <file>.
    # If TAG is not given, the frequency is estimated from GT by default
    # """


class BcfToolsStats(SnappyModel, BcfToolsMixinFilter, BcfToolsMixinLocus, BcfToolsMixinSample):
    af_bins: list[float] = [0.1, 0.5, 1.0]
    """Allele frequency bins, a list (0.1,0.5,1) or a file (0.1\n0.5\n1)"""
    af_tag: str | None = None
    """Allele frequency tag to use, by default estimated from AN,AC or GT"""
    first_allele_only: bool = False
    """Include only 1st allele at multiallelic sites"""
    collapse: str
    """Treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]"""
    depth: tuple[int, int, int] = (0, 500, 1)
    """Depth distribution: min,max,bin size [0,500,1]"""
    apply_filters: list = []
    """Require at least one of the listed FILTER strings (e.g. "PASS,.")"""
    split_by_ID: bool = False
    """Collect stats for sites with ID separately (known vs novel)"""
    user_tstv: str | None = None
    """
    TAG[:min:max:n]: Collect Ts/Tv stats for any tag using the given binning [0:1:100]
    A subfield can be selected as e.g. 'PV4[0]', here the first value of the PV4 tag
    """
    verbose: bool = False
    """Produce verbose per-site and per-sample output"""

    # exons: FILE.gz
    # """Tab-delimited file with exons for indel frameshifts (chr,beg,end; 1-based, inclusive, bgzip compressed)"""
    # fasta_ref: str | None = None  # FILE
    # """Faidx indexed reference sequence file to determine INDEL context"""
