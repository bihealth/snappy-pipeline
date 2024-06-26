from enum import Enum

from pydantic import Field
from typing_extensions import Annotated

from snappy_pipeline.models import EnumField, SnappyModel, SnappyStepModel
from snappy_pipeline.models.validators import ToolsMixin


class Tool(Enum):
    BBDUK = "bbduk"
    FASTP = "fastp"


class Interleaved(Enum):
    AUTO = "auto"


class Qin(Enum):
    AUTO = "auto"
    FIELD_33 = "33"
    FIELD_64 = "64"


class Qout(Enum):
    AUTO = "auto"
    FIELD_33 = "33"
    FIELD_64 = "64"


class Statscolumns(Enum):
    INTEGER_3 = 3
    INTEGER_5 = 5


class Gcbins(Enum):
    AUTO = "auto"


class Maxhistlen(Enum):
    AUTO = "auto"


class Idbins(Enum):
    AUTO = "auto"


class Ktrim(Enum):
    F = "f"
    R = "r"
    L = "l"


class Qtrim(Enum):
    RL = "rl"
    F = "f"
    R = "r"
    L = "l"
    W = "w"


class Barcodefilter(Enum):
    T = "t"
    F = "f"
    CRASH = "crash"


class Entropytrim(Enum):
    F = "f"
    """Do not entropy-trim"""

    R = "r"
    """Trim low entropy on the right end only."""

    L = "l"
    """Trim low entropy on the left end only."""

    RL = "rl"
    """Trim low entropy on both ends."""


class Entropymask(Enum):
    F = "f"
    T = "t"
    LC = "lc"


class UmiLoc(Enum):
    INDEX1 = "index1"
    INDEX2 = "index2"
    READ1 = "read1"
    READ2 = "read2"
    PER_INDEX = "per_index"
    PER_READ = "per_read"
    NONE = ""


class Fastp(SnappyModel):
    num_threads: int = 0
    trim_front1: int = 0
    """
    trimming how many bases in front for read1, default is 0 (int [=0])
    """

    trim_tail1: int = 0
    """
    trimming how many bases in tail for read1, default is 0 (int [=0])
    """

    max_len1: int = 0
    """
    if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1.
    Default 0 means no limitation (int [=0])
    """

    trim_front2: int = 0
    """
    trimming how many bases in front for read2.
    If it's not specified, it will follow read1's settings (int [=0])
    """

    trim_tail2: int = 0
    """
    trimming how many bases in tail for read2.
    If it's not specified, it will follow read1's settings (int [=0])
    """

    max_len2: int = 0
    """
    if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2.
    Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
    """

    dedup: bool = False
    """
    enable deduplication to drop the duplicated reads/pairs
    """

    dup_calc_accuracy: Annotated[int, Field(0, ge=0, le=6)]
    """
    accuracy level to calculate duplication (1~6),
    higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G).
    Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
    """

    dont_eval_duplication: bool = True
    """
    don't evaluate duplication rate to save time and use less memory.
    """

    trim_poly_g: bool = True
    """
    force polyG tail trimming,
    by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
    """

    poly_g_min_len: int = 8
    """
    the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
    """

    trim_poly_x: bool = False
    """
    enable polyX trimming in 3' ends.
    """

    poly_x_min_len: int = 10
    """
    the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
    """

    cut_front: bool = False
    """
    move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
    """

    cut_tail: bool = False
    """
    move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
    """

    cut_right: bool = False
    """
    move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
    """

    cut_front_window_size: int = 4
    """
    the window size option of cut_front, default to cut_window_size if not specified (int [=4])
    """

    cut_front_mean_quality: int = 20
    """
    the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
    """

    cut_tail_window_size: int = 4
    """
    the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
    """

    cut_tail_mean_quality: int = 20
    """
    the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
    """

    cut_right_window_size: int = 4
    """
    the window size option of cut_right, default to cut_window_size if not specified (int [=4])
    """

    cut_right_mean_quality: int = 20
    """
    the mean quality requirement option for cut_right,
    default to cut_mean_quality if not specified (int [=20])
    """

    disable_quality_filtering: bool = False
    """
    quality filtering is enabled by default.
    If this option is specified, quality filtering is disabled
    """

    qualified_quality_phred: int = 15
    """
    the quality value that a base is qualified.
    Default 15 means phred quality >=Q15 is qualified. (int [=15])
    """

    unqualified_percent_limit: int = 40
    """
    how many percents of bases are allowed to be unqualified (0~100).
    Default 40 means 40% (int [=40])
    """

    n_base_limit: int = 5
    """
    if one read's number of N base is >n_base_limit, then this read/pair is discarded.
    Default is 5 (int [=5])
    """

    average_qual: int = 0
    """
    if one read's average quality score <avg_qual, then this read/pair is discarded.
    Default 0 means no requirement (int [=0])
    """

    disable_length_filtering: bool = False
    """
    length filtering is enabled by default.
    If this option is specified, length filtering is disabled
    """

    length_required: int = 15
    """
    reads shorter than length_required will be discarded, default is 15. (int [=15])
    """

    length_limit: int = 0
    """
    reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
    """

    low_complexity_filter: bool = False
    """
    enable low complexity filter.
    The complexity is defined as the percentage of base that is different from its next base
    (base[i] != base[i+1]).
    """

    complexity_threshold: int = 30
    """
    the threshold for low complexity filter (0~100).
    Default is 30, which means 30% complexity is required. (int [=30])
    """

    filter_by_index1: str = ""
    """
    specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line
    (string [=])
    """

    filter_by_index2: str = ""
    """
    specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line
    (string [=])
    """

    filter_by_index_threshold: int = 0
    """
    the allowed difference of index barcode for index filtering,
    default 0 means completely identical. (int [=0])
    """

    correction: bool = False
    """
    enable base correction in overlapped regions (only for PE data), default is disabled
    """

    overlap_len_require: int = 30
    """
    the minimum length to detect overlapped region of PE reads.
    This will affect overlap analysis based PE merge, adapter trimming and correction.
    30 by default. (int [=30])
    """

    overlap_diff_limit: int = 5
    """
    the maximum number of mismatched bases to detect overlapped region of PE reads.
    This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.
    (int [=5])
    """

    overlap_diff_percent_limit: int = 20
    """
    the maximum percentage of mismatched bases to detect overlapped region of PE reads.
    This will affect overlap analysis based PE merge, adapter trimming and correction.
    Default 20 means 20%. (int [=20])
    """

    umi: bool = False
    """
    enable unique molecular identifier (UMI) preprocessing
    """

    umi_loc: UmiLoc = UmiLoc.NONE
    """
    specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read/""),
    default is "" (none).
    """

    umi_len: int = 0
    """
    if the UMI is in read1/read2, its length should be provided (int [=0])
    """

    umi_prefix: Annotated[str, Field("", examples=["UMI"])]
    """
    if specified, an underline will be used to connect prefix and UMI
    (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG).
    No prefix by default (string [=])
    """

    umi_skip: int = 0
    """
    if the UMI is in read1/read2, fastp can skip several bases following UMI,
    default is 0 (int [=0])
    """

    overrepresentation_analysis: bool = False
    """
    enable overrepresented sequence analysis.
    """


class Bbduk(SnappyModel):
    """Note: The author recommends setting tpe=t & tbo=t when adapter trimming paired reads."""

    adapter_sequences: Annotated[
        list[str],
        Field(
            examples=[
                [
                    "/fast/work/groups/cubi/projects/biotools/static_data/app_support/"
                    "bbtools/39.01/resources/adapters.fa",
                    "/fast/work/groups/cubi/projects/biotools/static_data/app_support/"
                    "bbtools/39.01/resources/phix174_ill.ref.fa.gz",
                ]
            ]
        ),
    ]

    num_threads: int = 8
    interleaved: bool | Interleaved = "auto"
    """
    (int) t/f overrides interleaved autodetection.
    """

    qin: Qin = "auto"
    """
    Input quality offset: 33 (Sanger), 64, or auto.
    """

    copyundefined: bool = False
    """
    (cu) Process non-AGCT IUPAC reference bases by making all possible unambiguous copies.
    Intended for short motifs or adapter barcodes, as time/memory use is exponential.
    """

    nzo: bool = True
    """
    Only write statistics about ref sequences with nonzero hits.
    """

    qout: Qout = "auto"
    """
    Output quality offset: 33 (Sanger), 64, or auto.
    """

    statscolumns: Statscolumns = 3
    """
    (cols) Number of columns for stats output, 3 or 5. 5 includes base counts.
    """

    rename: bool = False
    """
    Rename reads to indicate which sequences they matched.
    """

    refnames: bool = False
    """
    Use names of reference files rather than scaffold IDs.
    """

    trd: bool = False
    """
    Truncate read and ref names at the first whitespace.
    """

    ordered: bool = False
    """
    Set to true to output reads in same order as input.
    """

    gcbins: int | Gcbins = "auto"
    """
    Number gchist bins.  Set to 'auto' to use read length.
    """

    maxhistlen: int | Maxhistlen = 6000
    """
    Set an upper bound for histogram lengths; higher uses more memory.
    The default is 6000 for some histograms and 80000 for others.
    """

    histbefore: bool = True
    """
    Calculate histograms from reads before processing.
    """

    idbins: int | Idbins = 100
    """
    Number idhist bins.  Set to 'auto' to use read length.
    """

    k: Annotated[int, Field(21, ge=1)]
    """
    Kmer length used for finding contaminants.
    Contaminants shorter than k will not be found.
    k must be at least 1. bbduk default: 27
    """

    rcomp: bool = True
    """
    Look for reverse-complements of kmers in addition to forward kmers.
    """

    maskmiddle: bool = True
    """
    (mm) Treat the middle base of a kmer as a wildcard,
    to increase sensitivity in the presence of errors.
    """

    minkmerhits: int = 1
    """
    (mkh) Reads need at least this many matching kmers to be considered as matching the reference.
    """

    minkmerfraction: float = 0
    """
    (mkf) A reads needs at least this fraction of its total kmers to hit a ref,
    in order to be considered a match.
    If this and minkmerhits are set, the greater is used.
    """

    mincovfraction: float = 0
    """
    (mcf) A reads needs at least this fraction of its total bases to be covered by ref kmers
    to be considered a match.
    If specified, mcf overrides mkh and mkf.
    """

    hammingdistance: int = 1
    """
    (hdist) Maximum Hamming distance for ref kmers (subs only).
    Memory use is proportional to (3*K)^hdist. bbduk default: 0
    """

    qhdist: int = 0
    """
    Hamming distance for query kmers; impacts speed, not memory.
    """

    editdistance: int = 0
    """
    (edist) Maximum edit distance from ref kmers (subs and indels).
    Memory use is proportional to (8*K)^edist.
    """

    hammingdistance2: int = 0
    """
    (hdist2) Sets hdist for short kmers, when using mink.
    """

    qhdist2: int = 0
    """
    Sets qhdist for short kmers, when using mink.
    """

    editdistance2: int = 0
    """
    (edist2) Sets edist for short kmers, when using mink.
    """

    forbidn: bool = False
    """
    (fn) Forbids matching of read kmers containing N.
    By default, these will match a reference 'A' if hdist>0 or edist>0, to increase sensitivity.
    """

    removeifeitherbad: bool = True
    """
    (rieb) Paired reads get sent to 'outmatch' if either is match
    (or either is trimmed shorter than minlen).
    Set to false to require both.
    """

    trimfailures: bool = False
    """
    Instead of discarding failed reads, trim them to 1bp.
    This makes the statistics a bit odd.
    """

    findbestmatch: bool = False
    """
    (fbm) If multiple matches, associate read with sequence sharing most kmers. Reduces speed.
    """

    skipr1: bool = False
    """
    Don't do kmer-based operations on read 1.
    """

    skipr2: bool = False
    """
    Don't do kmer-based operations on read 2.
    """

    ecco: bool = False
    """
    For overlapping paired reads only.
    Performs error- correction with BBMerge prior to kmer operations.
    """

    ktrim: Ktrim = "r"
    """
    Trim reads to remove bases matching reference kmers. Values:
      f (don't trim), [bbduk default]
      r (trim to the right),
      l (trim to the left)
    """

    kmask: str = ""
    """
    Replace bases matching ref kmers with another symbol.
    Allows any non-whitespace character, and processes short kmers on both ends if mink is set.
    'kmask: lc' will convert masked bases to lowercase.
    """

    maskfullycovered: bool = False
    """
    (mfc) Only mask bases that are fully covered by kmers.
    """

    ksplit: bool = False
    """
    For single-ended reads only.
    Reads will be split into pairs around the kmer.
    If the kmer is at the end of the read, it will be trimmed instead.
    Singletons will go to out, and pairs will go to outm.
    Do not use ksplit with other operations such as quality-trimming or filtering.
    """

    mink: int = 11
    """
    Look for shorter kmers at read tips down to this length, when k-trimming or masking.
    0 means disabled.
    Enabling this will disable maskmiddle. bbduk default: 0 (disabled)
    """

    qtrim: Qtrim = "rl"
    """
    Trim read ends to remove bases with quality below trimq.
    Performed AFTER looking for kmers.  Values:
      rl (trim both ends),
      f (neither end),  [bbduk default]
      r (right end only),
      l (left end only),
      w (sliding window).
    """

    trimq: float = 25
    """
    Regions with average quality BELOW this will be trimmed,
    if qtrim is set to something other than f.
    Can be a floating-point number like 7.3.
    Very strict quality threshold, bbduk default: 6
    """

    minlength: int = 35
    """
    (ml) Reads shorter than this after trimming will be discarded.
    Pairs will be discarded if both are shorter.
    bbduk default: 10
    """

    mlf: int = 0
    """
    (minlengthfraction) Reads shorter than this fraction of original length after trimming
    will be discarded.
    """

    minavgquality: int = 0
    """
    (maq) Reads with average quality (after trimming) below this will be discarded.
    """

    maqb: int = 0
    """
    If positive, calculate maq from this many initial bases.
    """

    minbasequality: int = 0
    """
    (mbq) Reads with any base below this quality (after trimming) will be discarded.
    """

    maxns: int = -1
    """
    If non-negative, reads with more Ns than this (after trimming) will be discarded.
    """

    mcb: int = 0
    """
    (minconsecutivebases) Discard reads without at least this many consecutive called bases.
    """

    ottm: bool = False
    """
    (outputtrimmedtomatch) Output reads trimmed to shorter than minlength to outm rather than discarding.
    """

    tp: int = 0
    """
    (trimpad) Trim this much extra around matching kmers.
    """

    tbo: bool = False
    """
    (trimbyoverlap) Trim adapters based on where paired reads overlap.
    Note: The author recommends setting tpe=t & tbo=t when adapter trimming paired reads.
    """

    strictoverlap: bool = True
    """
    Adjust sensitivity for trimbyoverlap mode.
    """

    minoverlap: int = 14
    """
    Require this many bases of overlap for detection.
    """

    mininsert: int = 40
    """
    Require insert size of at least this for overlap.
    Should be reduced to 16 for small RNA sequencing.
    """

    tpe: bool = False
    """
    (trimpairsevenly) When kmer right-trimming, trim both reads to the minimum length of either.
    Note: The author recommends setting tpe=t & tbo=t when adapter trimming paired reads.
    """

    forcetrimleft: int = 0
    """
    (ftl) If positive, trim bases to the left of this position (exclusive, 0-based).
    """

    forcetrimright: int = 0
    """
    (ftr) If positive, trim bases to the right of this position (exclusive, 0-based).

    """
    forcetrimright2: int = 0
    """
    (ftr2) If positive, trim this many bases on the right end.
    """

    forcetrimmod: int = 5
    """
    (ftm) If positive, right-trim length to be equal to zero, modulo this number. bbduk default: 0
    """

    restrictleft: int = 0
    """
    If positive, only look for kmer matches in the leftmost X bases.
    """

    restrictright: int = 0
    """
    If positive, only look for kmer matches in the rightmost X bases.
    """

    mingc: float = 0
    """
    Discard reads with GC content below this.
    """

    maxgc: float = 1
    """
    Discard reads with GC content above this.
    """

    gcpairs: bool = True
    """
    Use average GC of paired reads.    Deprecated option? Also affects gchist.
    """

    tossjunk: bool = False
    """
    Discard reads with invalid characters as bases.
    """

    swift: bool = False
    """
    Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2.
    """

    chastityfilter: bool = False
    """
    (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
    """

    barcodefilter: Barcodefilter = "f"
    """
    Remove reads with unexpected barcodes if barcodes is set,
    or barcodes containing 'N' otherwise.
    A barcode must be the last part of the read header.
    Values:
      t:     Remove reads with bad barcodes.
      f:     Ignore barcodes.
      crash: Crash upon encountering bad barcodes.
    """

    barcodes: str = ""
    """
    File of barcodes.
    """

    xmin: int = -1
    """
    If positive, discard reads with a lesser X coordinate.
    """

    ymin: int = -1
    """
    If positive, discard reads with a lesser Y coordinate.
    """

    xmax: int = -1
    """
    If positive, discard reads with a greater X coordinate.
    """

    ymax: int = -1
    """
    If positive, discard reads with a greater Y coordinate.
    """

    trimpolya: int = 0
    """
    If greater than 0, trim poly-A or poly-T tails of at least this length on either end of reads.
    """

    trimpolygleft: int = 0
    """
    If greater than 0, trim poly-G prefixes of at least this length on the left end of reads.
    Does not trim poly-C.
    """

    trimpolygright: int = 8
    """
    If greater than 0, trim poly-G tails of at least this length on the right end of reads.
    Does not trim poly-C. bbduk default: don't trim polyG (trimpolyg=0)
    """

    trimpolyg: int = 0
    """
    This sets both left and right at once.
    """

    filterpolyg: int = 8
    """
    If greater than 0, remove reads with a poly-G prefix of at least this length (on the left).
    Note: there are also equivalent poly-C flags.
    """

    entropy: float = -1
    """
    Set between 0 and 1 to filter reads with entropy below that value.
    Higher is more stringent.
    """

    entropywindow: int = 50
    """
    Calculate entropy using a sliding window of this length.
    """

    entropyk: int = 5
    """
    Calculate entropy using kmers of this length.
    """

    minbasefrequency: float = 0
    """
    Discard reads with a minimum base frequency below this.
    """

    entropytrim: Entropytrim = "f"
    """
    Values:
      f:  (false) Do not entropy-trim.
      r:  (right) Trim low entropy on the right end only.
      l:  (left) Trim low entropy on the left end only.
      rl: (both) Trim low entropy on both ends.
    NOTE: If set, entropytrim overrides entropymask.
    """

    entropymask: Entropymask = "f"
    """
    Values:
      f:  (filter) Discard low-entropy sequences.
      t:  (true) Mask low-entropy parts of sequences with N.
      lc: Change low-entropy parts of sequences to lowercase.
    """

    entropymark: bool = False
    """
    Mark each base with its entropy value.
    This is on a scale of 0-41 and is reported as quality scores,
    so the output should be fastq or fasta+qual. NOTE: If set, entropytrim overrides entropymask.
    """

    cardinality: bool = False
    """
    (loglog) Count unique kmers using the LogLog algorithm.
    """

    cardinalityout: bool = False
    """
    (loglogout) Count unique kmers in output reads.
    """

    loglogk: Annotated[int, Field(31, gt=0)]
    """
    Use this kmer length for counting.
    """

    loglogbuckets: Annotated[int, Field(2048, gt=0)]
    """
    Use this many buckets for counting.
    """


class AdapterTrimming(SnappyStepModel, ToolsMixin):
    path_link_in: str | None = None
    """Override data set configuration search paths for FASTQ files"""

    tools: Annotated[list[Tool], EnumField(Tool, min_length=1, default=["bbduk", "fastp"])]
    bbduk: Bbduk | None = None
    fastp: Fastp | None = None
