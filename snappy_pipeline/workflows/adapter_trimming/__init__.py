# -*- coding: utf-8 -*-
"""Implementation of the ``adapter_trimming`` step

The adapter_trimming step performs adapter & quality trimming of reads (DNA or RNA).
The tools are highly configurable, and provide feedback of the success of the operation.

==========
Step Input
==========

For each library defined in all sample sheets, the instances of this step will search for the input
files according to the configuration.  The found read files will be linked into
``work/input_links/{library_name}`` (status quo, not a output path, thus path not guaranteed
to be stable between minor versions).

The search paths can be overridden using the step configuration option ``path_link_in``.
``path_link_in`` is a general features that enables pre-processing steps, typically before mapping.

----------------------
Data Set Configuration
----------------------

Consider the following data set definition from the main configuration file.

.. code-block:: yaml

    data_sets:
      first_batch:
        file: 01_first_batch.tsv
        search_patterns:
          # Note that currently only "left" and "right" key known
          - {'left': '*/L???/*_R1.fastq.gz', 'right': '*/L???/*_R2.fastq.gz'}
        search_paths: ['../input/01_first_batch']

Here, the data set ``first_batch`` is defined.  The sample sheet file is named
``01_first_batch.tsv`` and looked for in the relative path to the configuration file.  The input
search will be start in the (one, but could be more than one) path ``../input/01_first_batch``
(relative to the directory containing the configuration file).  The sample sheet provides a
``folderName`` ``extraInfo`` entry for each NGS library.  This folder name is searched for (e.g.,
``P001-N1-DNA1-WES``).  Once such a folder is found, the patterns in the values of the dict
``search_patterns`` are used for locating the paths of the actual files.

Currently, the only supported keys in the ``search_patterns`` dict are ``"left"`` and ``"right""``
(the latter can be omitted when only searching for single-end reads).

Consider the following example:

::

  ../input/
    `-- 01_first_batch
        |-- P001-N1-DNA1-WES1
        |   `-- 42KF5AAXX
        |       `-- L001
        |           |-- P001-N1-DNA1-WES1_R1.fastq.gz
        |           |-- P001-N1-DNA1-WES1_R1.fastq.gz.md5
        |           |-- P001-N1-DNA1-WES1_R2.fastq.gz
        |           `-- P001-N1-DNA1-WES1_R2.fastq.gz.md5
        [...]

Here, the folder ``01_first_batch`` will be searched for a directory named ``P001-N1-DNA1-WES``.
After finding, the relative paths ``42KF5AAXX/L001/P001-N1-DNA1-WES1_R1.fastq.gz`` and
``42KF5AAXX/L001/P001-N1-DNA1-WES1_R2.fastq.gz`` will be found and used for the left/right parts of
a paired read set.

------------------------------------------------------
Overriding data set confguration with ``path_link_in``
------------------------------------------------------

When the config option ``path_link_in`` is set, it takes precedence on the search paths defined in the
data set configuration.

The searching for input files will follow the same rules as defined in the data set configuration,
except that the base path for the search provided by one single path defined in the configuration of the
step.


Mixing Single-End and Paired-End Reads
======================================

By default, it is checked that for each ``search_pattern``, the same number of matching files
has to be found, otherwise directories are ignored.  The reason is to reduce the number of
possible errors when linking in files.  You can change this behaviour by specifying
``mixed_se_pe: True`` in the data set information.  Then, it will be allowed to have the matches
for the ``right`` entry to be empty.  However, you will need to consistently have either SE or
PE data for each library; it is allowed to mix SE and PE libraries within one project but not
to have PE and SE data for one library.

Note that mixing single-end and paired-end reads is not (yet) supported when overriding the data set
configuration by setting a value to the configuration option ``path_link_in``.


===========
Step Output
===========

Adapter trimming will be performed for all NGS libraries in all sample sheets.  For each combination
of tool library, a directory ``{tool}/{lib_name}-{lib_pk}/out`` will be created.
Therein, trimmed fastq files will be created.

The input structure and file names will be maintained on output.  For example, it might look as
follows for the example from above:

::

    output/
    +-- bbduk
    |   `-- out
    |       `-- P001-N1-DNA1-WES1
    |           |-- 42KF5AAXX
    |           |   `-- L001
    |           |       |-- P001-N1-DNA1-WES1_R1.fastq.gz
    |           |       |-- P001-N1-DNA1-WES1_R1.fastq.gz.md5
    |           |       |-- P001-N1-DNA1-WES1_R2.fastq.gz
    |           |       `-- P001-N1-DNA1-WES1_R2.fastq.gz.md5
    |           `-- .done
    [...]


=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_adapter_trimming.rst

================================
Available Adapter Trimming Tools
================================

The following adpter trimming tools are currently available

- ``"bbduk"``
- ``"fastp"``

"""

from collections import OrderedDict
import os

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    ResourceUsage,
    get_ngs_library_folder_name,
)

#: Adatper trimming tools
TRIMMERS = ("bbduk", "fastp")

#: Default configuration for the hla_typing schema
DEFAULT_CONFIG = r"""
# Default configuration adapter_trimming
step_config:
  adapter_trimming:
    path_link_in: ""  # OPTIONAL Override data set configuration search paths for FASTQ files
    tools: [bbduk, fastp]  # REQUIRED, available: 'bbduk' and 'fastp'.
    bbduk:
      adapter_sequences: []  # REQUIRED
      # - /fast/work/groups/cubi/projects/biotools/static_data/app_support/bbtools/39.01/resources/adapters.fa
      # - /fast/work/groups/cubi/projects/biotools/static_data/app_support/bbtools/39.01/resources/phix174_ill.ref.fa.gz
      # Note: The author recommends setting tpe=t & tbo=t when adapter trimming paired reads.
      num_threads: 8

      # Non-default parameters from https://www.biostars.org/p/268221/
      # & https://github.com/ewels/MultiQC/issues/1146#issuecomment-607980076

      # Input parameters:
      interleaved: auto    # (int) t/f overrides interleaved autodetection.
      qin: auto            # Input quality offset: 33 (Sanger), 64, or auto.
      copyundefined: f     # (cu) Process non-AGCT IUPAC reference bases by making all
                           # possible unambiguous copies.  Intended for short motifs
                           # or adapter barcodes, as time/memory use is exponential.

      # Output parameters:
      nzo: t               # Only write statistics about ref sequences with nonzero hits.
      qout: auto           # Output quality offset: 33 (Sanger), 64, or auto.
      statscolumns: 3      # (cols) Number of columns for stats output, 3 or 5.
                           # 5 includes base counts.
      rename: f            # Rename reads to indicate which sequences they matched.
      refnames: f          # Use names of reference files rather than scaffold IDs.
      trd: f               # Truncate read and ref names at the first whitespace.
      ordered: f           # Set to true to output reads in same order as input.

      # Histogram output parameters:
      gcbins: auto         # Number gchist bins.  Set to 'auto' to use read length.
      maxhistlen: 6000     # Set an upper bound for histogram lengths; higher uses
                           # more memory.  The default is 6000 for some histograms
                           # and 80000 for others.

      # Histograms for mapped sam/bam files only:
      histbefore: t        # Calculate histograms from reads before processing.
      idbins: 100          # Number idhist bins.  Set to 'auto' to use read length.

      # Processing parameters:
      k: 21                # Kmer length used for finding contaminants.  Contaminants
                           # shorter than k will not be found.  k must be at least 1.
                           # bbduk default: 27
      rcomp: t             # Look for reverse-complements of kmers in addition to
                           # forward kmers.
      maskmiddle: t        # (mm) Treat the middle base of a kmer as a wildcard, to
                           # increase sensitivity in the presence of errors.
      minkmerhits: 1       # (mkh) Reads need at least this many matching kmers
                           # to be considered as matching the reference.
      minkmerfraction: 0.0 # (mkf) A reads needs at least this fraction of its total
                           # kmers to hit a ref, in order to be considered a match.
                           # If this and minkmerhits are set, the greater is used.
      mincovfraction: 0.0  # (mcf) A reads needs at least this fraction of its total
                           # bases to be covered by ref kmers to be considered a match.
                           # If specified, mcf overrides mkh and mkf.
      hammingdistance: 1   # (hdist) Maximum Hamming distance for ref kmers (subs only).
                           # Memory use is proportional to (3*K)^hdist.
                           # bbduk default: 0
      qhdist: 0            # Hamming distance for query kmers; impacts speed, not memory.
      editdistance: 0      # (edist) Maximum edit distance from ref kmers (subs
                           # and indels).  Memory use is proportional to (8*K)^edist.
      hammingdistance2: 0  # (hdist2) Sets hdist for short kmers, when using mink.
      qhdist2: 0           # Sets qhdist for short kmers, when using mink.
      editdistance2: 0     # (edist2) Sets edist for short kmers, when using mink.
      forbidn: f           # (fn) Forbids matching of read kmers containing N.
                           # By default, these will match a reference 'A' if
                           # hdist>0 or edist>0, to increase sensitivity.
      removeifeitherbad: t # (rieb) Paired reads get sent to 'outmatch' if either is
                           # match (or either is trimmed shorter than minlen).
                           # Set to false to require both.
      trimfailures: f      # Instead of discarding failed reads, trim them to 1bp.
                           # This makes the statistics a bit odd.
      findbestmatch: f     # (fbm) If multiple matches, associate read with sequence
                           # sharing most kmers.  Reduces speed.
      skipr1: f            # Don't do kmer-based operations on read 1.
      skipr2: f            # Don't do kmer-based operations on read 2.
      ecco: f              # For overlapping paired reads only.  Performs error-
                           # correction with BBMerge prior to kmer operations.

      # Trimming/Filtering/Masking parameters:
      # Note - if ktrim, kmask, and ksplit are unset, the default behavior is kfilter.
      # All kmer processing modes are mutually exclusive.
      # Reads only get sent to 'outm' purely based on kmer matches in kfilter mode.

      ktrim: r             # Trim reads to remove bases matching reference kmers.
                           # Values:
                           #   f (don't trim), [bbduk default]
                           #   r (trim to the right),
                           #   l (trim to the left)
      kmask:               # Replace bases matching ref kmers with another symbol.
                           # Allows any non-whitespace character, and processes short
                           # kmers on both ends if mink is set.  'kmask: lc' will
                           # convert masked bases to lowercase.
      maskfullycovered: f  # (mfc) Only mask bases that are fully covered by kmers.
      ksplit: f            # For single-ended reads only.  Reads will be split into
                           # pairs around the kmer.  If the kmer is at the end of the
                           # read, it will be trimmed instead.  Singletons will go to
                           # out, and pairs will go to outm.  Do not use ksplit with
                           # other operations such as quality-trimming or filtering.
      mink: 11             # Look for shorter kmers at read tips down to this length,
                           # when k-trimming or masking.  0 means disabled.  Enabling
                           # this will disable maskmiddle.
                           # bbduk default: 0 (disabled)
      qtrim: rl            # Trim read ends to remove bases with quality below trimq.
                           # Performed AFTER looking for kmers.  Values:
                           #   rl (trim both ends),
                           #   f (neither end),  [bbduk default]
                           #   r (right end only),
                           #   l (left end only),
                           #   w (sliding window).
      trimq: 25            # Regions with average quality BELOW this will be trimmed,
                           # if qtrim is set to something other than f.  Can be a
                           # floating-point number like 7.3.
                           # Very strict quality threshold, bbduk default: 6
      minlength: 35        # (ml) Reads shorter than this after trimming will be
                           # discarded.  Pairs will be discarded if both are shorter.
                           # bbduk default: 10
      mlf: 0               # (minlengthfraction) Reads shorter than this fraction of
                           # original length after trimming will be discarded.
      minavgquality: 0     # (maq) Reads with average quality (after trimming) below
                           # this will be discarded.
      maqb: 0              # If positive, calculate maq from this many initial bases.
      minbasequality: 0    # (mbq) Reads with any base below this quality (after
                           # trimming) will be discarded.
      maxns: -1            # If non-negative, reads with more Ns than this
                           # (after trimming) will be discarded.
      mcb: 0               # (minconsecutivebases) Discard reads without at least
                           # this many consecutive called bases.
      ottm: f              # (outputtrimmedtomatch) Output reads trimmed to shorter
                           # than minlength to outm rather than discarding.
      tp: 0                # (trimpad) Trim this much extra around matching kmers.
      tbo: f               # (trimbyoverlap) Trim adapters based on where paired
                           # reads overlap.
      strictoverlap: t     # Adjust sensitivity for trimbyoverlap mode.
      minoverlap: 14       # Require this many bases of overlap for detection.
      mininsert: 40        # Require insert size of at least this for overlap.
                           # Should be reduced to 16 for small RNA sequencing.
      tpe: f               # (trimpairsevenly) When kmer right-trimming, trim both
                           # reads to the minimum length of either.
      forcetrimleft: 0     # (ftl) If positive, trim bases to the left of this position
                           # (exclusive, 0-based).
      forcetrimright: 0    # (ftr) If positive, trim bases to the right of this position
                           # (exclusive, 0-based).
      forcetrimright2: 0   # (ftr2) If positive, trim this many bases on the right end.
      forcetrimmod: 5      # (ftm) If positive, right-trim length to be equal to zero,
                           # modulo this number.
                           # bbduk default: 0
      restrictleft: 0      # If positive, only look for kmer matches in the
                           # leftmost X bases.
      restrictright: 0     # If positive, only look for kmer matches in the
                           # rightmost X bases.
      mingc: 0             # Discard reads with GC content below this.
      maxgc: 1             # Discard reads with GC content above this.
      # gcpairs: t           # Use average GC of paired reads.    Deprecated option?
      #                      # Also affects gchist.
      tossjunk: f          # Discard reads with invalid characters as bases.
      swift: f             # Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2.

      # Header-parsing parameters - these require Illumina headers:
      chastityfilter: f    # (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
      barcodefilter: f     # Remove reads with unexpected barcodes if barcodes is set,
                           # or barcodes containing 'N' otherwise.  A barcode must be
                           # the last part of the read header.  Values:
                           #   t:     Remove reads with bad barcodes.
                           #   f:     Ignore barcodes.
                           #   crash: Crash upon encountering bad barcodes.
      barcodes: ""         # File of barcodes.
      xmin: -1             # If positive, discard reads with a lesser X coordinate.
      ymin: -1             # If positive, discard reads with a lesser Y coordinate.
      xmax: -1             # If positive, discard reads with a greater X coordinate.
      ymax: -1             # If positive, discard reads with a greater Y coordinate.

      # Polymer trimming:
      trimpolya: 0         # If greater than 0, trim poly-A or poly-T tails of
                           # at least this length on either end of reads.
      trimpolygleft: 0     # If greater than 0, trim poly-G prefixes of at least this
                           # length on the left end of reads.  Does not trim poly-C.
      trimpolygright: 0    # If greater than 0, trim poly-G tails of at least this
                           # length on the right end of reads.  Does not trim poly-C.
      trimpolyg: 8         # This sets both left and right at once.
                           # bbduk default: don't trim polyG (trimpolyg=0)
      filterpolyg: 0       # If greater than 0, remove reads with a poly-G prefix of
                           # at least this length (on the left).
      # Note: there are also equivalent poly-C flags.

      # Entropy/Complexity parameters:
      entropy: -1          # Set between 0 and 1 to filter reads with entropy below
                           # that value.  Higher is more stringent.
      entropywindow: 50    # Calculate entropy using a sliding window of this length.
      entropyk: 5          # Calculate entropy using kmers of this length.
      minbasefrequency: 0  # Discard reads with a minimum base frequency below this.
      entropytrim: f       # Values:
                           #    f:  (false) Do not entropy-trim.
                           #    r:  (right) Trim low entropy on the right end only.
                           #    l:  (left) Trim low entropy on the left end only.
                           #    rl: (both) Trim low entropy on both ends.
      entropymask: f       # Values:
                           #    f:  (filter) Discard low-entropy sequences.
                           #    t:  (true) Mask low-entropy parts of sequences with N.
                           #    lc: Change low-entropy parts of sequences to lowercase.
      entropymark: f       # Mark each base with its entropy value.  This is on a scale
                           # of 0-41 and is reported as quality scores, so the output
                           # should be fastq or fasta+qual.
      # NOTE: If set, entropytrim overrides entropymask.

      # Cardinality estimation:
      cardinality: f       # (loglog) Count unique kmers using the LogLog algorithm.
      cardinalityout: f    # (loglogout) Count unique kmers in output reads.
      loglogk: 31          # Use this kmer length for counting.
      loglogbuckets: 2048  # Use this many buckets for counting.

    fastp:
      num_threads: 4

      trim_front1: 0                      # trimming how many bases in front for read1, default is 0 (int [=0])
      trim_tail1: 0                       # trimming how many bases in tail for read1, default is 0 (int [=0])
      max_len1: 0                         # if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
      trim_front2: 0                      # trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
      trim_tail2: 0                       # trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
      max_len2: 0                         # if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
      dedup: False                        # enable deduplication to drop the duplicated reads/pairs
      dup_calc_accuracy: 0                # accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
      dont_eval_duplication: True         # don't evaluate duplication rate to save time and use less memory.
      trim_poly_g: True                   # force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
      poly_g_min_len: 8                   # the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
      trim_poly_x: False                  # enable polyX trimming in 3' ends.
      poly_x_min_len: 10                  # the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
      cut_front: False                    # move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
      cut_tail: False                     # move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
      cut_right: False                    # move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
      cut_front_window_size: 4            # the window size option of cut_front, default to cut_window_size if not specified (int [=4])
      cut_front_mean_quality: 20          # the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
      cut_tail_window_size: 4             # the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
      cut_tail_mean_quality: 20           # the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
      cut_right_window_size: 4            # the window size option of cut_right, default to cut_window_size if not specified (int [=4])
      cut_right_mean_quality: 20          # the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
      disable_quality_filtering: False    # quality filtering is enabled by default. If this option is specified, quality filtering is disabled
      qualified_quality_phred: 15         # the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
      unqualified_percent_limit: 40       # how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
      n_base_limit: 5                     # if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
      average_qual: 0                     # if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
      disable_length_filtering: False     # length filtering is enabled by default. If this option is specified, length filtering is disabled
      length_required: 15                 # reads shorter than length_required will be discarded, default is 15. (int [=15])
      length_limit: 0                     # reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
      low_complexity_filter: False        # enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
      complexity_threshold: 30            # the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
      filter_by_index1: ""                # specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      filter_by_index2: ""                # specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      filter_by_index_threshold: 0        # the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
      correction: False                   # enable base correction in overlapped regions (only for PE data), default is disabled
      overlap_len_require: 30             # the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
      overlap_diff_limit: 5               # the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
      overlap_diff_percent_limit: 20      # the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
      umi: False                          # enable unique molecular identifier (UMI) preprocessing
      umi_loc: ""                         # specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      umi_len: 0                          # if the UMI is in read1/read2, its length should be provided (int [=0])
      umi_prefix: ""                      # if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      umi_skip: 0                         # if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])
      overrepresentation_analysis: False  # enable overrepresented sequence analysis.

""".lstrip()


class AdapterTrimmingStepPart(BaseStepPart):
    """Adapter trimming common features"""

    #: Step name
    name = ""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{trimmer}.{{library_name}}"
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config["path_link_in"],
        )

    @dictify
    def get_input_files(self, action):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        yield "done", "work/input_links/{library_name}/.done"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        # Validate action
        self._validate_action(action)
        return (
            ("out_done", self.base_path_out.format(trimmer=self.name) + "/out/.done"),
            ("report_done", self.base_path_out.format(trimmer=self.name) + "/report/.done"),
            ("rejected_done", self.base_path_out.format(trimmer=self.name) + "/rejected/.done"),
        )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)
        _ = action
        prefix = "work/{trimmer}.{{library_name}}/log/{trimmer}.{{library_name}}".format(
            trimmer=self.__class__.name
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        yield "done", "work/{trimmer}.{{library_name}}/log/.done".format(
            trimmer=self.__class__.name
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
            if self.config["path_link_in"]:
                folder_name = wildcards.library_name
            reads_left = self._collect_reads(wildcards, folder_name, "")
            reads_right = self._collect_reads(wildcards, folder_name, "right-")
            return {
                "library_name": wildcards.library_name,
                "input": {
                    "reads_left": {key: reads_left[key] for key in sorted(reads_left.keys())},
                    "reads_right": {key: reads_right[key] for key in sorted(reads_right.keys())},
                },
            }

        # Validate action
        self._validate_action(action)
        return args_function

    def _collect_reads(self, wildcards, folder_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        path_info = {}
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            input_path = os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)
            assert input_path not in path_info.keys()
            paths = {
                "relative_path": path_infix,
                "filename": filename,
            }
            path_info[input_path] = paths
        return path_info


class BbdukStepPart(AdapterTrimmingStepPart):
    """bbduk adapter trimming"""

    name = "bbduk"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["bbduk"]["num_threads"],
            time="12:00:00",  # 40 hours
            memory="24000M",
        )


class FastpStepPart(AdapterTrimmingStepPart):
    """fastp adapter trimming"""

    #: Step name
    name = "fastp"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["fastp"]["num_threads"],
            time="12:00:00",  # 60 hours
            memory="24000M",
        )


class LinkOutFastqStepPart(BaseStepPart):
    """Link out the trimming results (all fastqs)"""

    name = "link_out_fastq"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/{wildcards.trimmer}.{wildcards.library_name}/{{sub_dir}}/.done"
        self.base_path_out = "output/{{trimmer}}/{{library_name}}/{sub_dir}/.done"
        self.sub_dirs = ["log", "report", "out"]

    def get_input_files(self, action):
        """Return required input files"""

        def input_function(wildcards):
            """Helper wrapper function"""
            return expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)

        # Validate action
        self._validate_action(action)
        return input_function

    def get_output_files(self, action):
        """Return output files that are generated by snappy-gatk_post_bam"""
        # Validate action
        self._validate_action(action)
        return expand(self.base_path_out, sub_dir=self.sub_dirs)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking out postprocessed (or not) files"""
        # Validate action
        self._validate_action(action)
        ins = expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)
        outs = [s.format(**wildcards) for s in expand(self.base_path_out, sub_dir=self.sub_dirs)]
        assert len(ins) == len(outs)

        cmd = "din_=$(dirname {in_}) ; dout=$(dirname {out})"
        cmd = cmd + " ; fns=$(find $din_ -type f -printf '%P\\n')"
        cmd = (
            cmd
            + " ; for fn in $fns ; do"
            + "     if [[ ! -L $din_/$fn ]] ; then"
            + "       mkdir -p $(dirname $dout/$fn) ; ln -sr $din_/$fn $dout/$fn"
            + "   ; fi"
            + " ; done"
        )
        return "\n".join((cmd.format(in_=in_, out=out) for in_, out in zip(ins, outs)))

    def _validate_action(self, action):
        assert action == "run"


class AdapterTrimmingWorkflow(BaseStep):
    """Perform adapter & quality-based trimming"""

    #: Step name
    name = "adapter_trimming"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.register_sub_step_classes(
            (BbdukStepPart, FastpStepPart, LinkInStep, LinkOutFastqStepPart)
        )
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    @listify
    def get_result_files(self):
        """Return list of fixed name result files for the adapter trimming workflow"""
        tpls = (
            "output/{trimmer}/{ngs_library_name}/out/.done",
            "output/{trimmer}/{ngs_library_name}/report/.done",
            "output/{trimmer}/{ngs_library_name}/log/.done",
        )
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                for tool in self.config["tools"]:
                    for tpl in tpls:
                        yield tpl.format(trimmer=tool, ngs_library_name=ngs_library.name)
