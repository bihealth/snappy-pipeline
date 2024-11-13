from enum import StrEnum

from snappy_pipeline.models import SnappyModel


class Annotation(StrEnum):
    AS_BASEQUALITYRANKSUMTEST = 'AS_BaseQualityRankSumTest'
    AS_FISHERSTRAND = 'AS_FisherStrand'
    AS_INBREEDINGCOEFF = 'AS_InbreedingCoeff'
    AS_MAPPINGQUALITYRANKSUMTEST = 'AS_MappingQualityRankSumTest'
    AS_QUALBYDEPTH = 'AS_QualByDepth'
    AS_RMSMAPPINGQUALITY = 'AS_RMSMappingQuality'
    AS_READPOSRANKSUMTEST = 'AS_ReadPosRankSumTest'
    AS_STRANDBIASMUTECTANNOTATION = 'AS_StrandBiasMutectAnnotation'
    AS_STRANDODDSRATIO = 'AS_StrandOddsRatio'
    ALLELEFRACTION = 'AlleleFraction'
    ALLELEPSEUDODEPTH = 'AllelePseudoDepth'
    ASSEMBLYCOMPLEXITY = 'AssemblyComplexity'
    BASEQUALITY = 'BaseQuality'
    BASEQUALITYHISTOGRAM = 'BaseQualityHistogram'
    BASEQUALITYRANKSUMTEST = 'BaseQualityRankSumTest'
    CHROMOSOMECOUNTS = 'ChromosomeCounts'
    CLIPPINGRANKSUMTEST = 'ClippingRankSumTest'
    COUNTNS = 'CountNs'
    COVERAGE = 'Coverage'
    CYCLESKIPSTATUS = 'CycleSkipStatus'
    DEPTHPERALLELEBYSAMPLE = 'DepthPerAlleleBySample'
    DEPTHPERSAMPLEHC = 'DepthPerSampleHC'
    EXCESSHET = 'ExcessHet'
    FEATURIZEDREADSETS = 'FeaturizedReadSets'
    FISHERSTRAND = 'FisherStrand'
    FRAGMENTDEPTHPERALLELEBYSAMPLE = 'FragmentDepthPerAlleleBySample'
    FRAGMENTLENGTH = 'FragmentLength'
    GCCONTENT = 'GcContent'
    GENOTYPESUMMARIES = 'GenotypeSummaries'
    HAPLOTYPEFILTERINGANNOTATION = 'HaplotypeFilteringAnnotation'
    HMERINDELLENGTH = 'HmerIndelLength'
    HMERINDELNUC = 'HmerIndelNuc'
    HMERMOTIFS = 'HmerMotifs'
    INBREEDINGCOEFF = 'InbreedingCoeff'
    INDELCLASSIFY = 'IndelClassify'
    INDELLENGTH = 'IndelLength'
    LIKELIHOODRANKSUMTEST = 'LikelihoodRankSumTest'
    MAPPINGQUALITY = 'MappingQuality'
    MAPPINGQUALITYRANKSUMTEST = 'MappingQualityRankSumTest'
    MAPPINGQUALITYZERO = 'MappingQualityZero'
    ORIENTATIONBIASREADCOUNTS = 'OrientationBiasReadCounts'
    ORIGINALALIGNMENT = 'OriginalAlignment'
    POSSIBLEDENOVO = 'PossibleDeNovo'
    QUALBYDEPTH = 'QualByDepth'
    RMSMAPPINGQUALITY = 'RMSMappingQuality'
    RAWGTCOUNT = 'RawGtCount'
    READPOSRANKSUMTEST = 'ReadPosRankSumTest'
    READPOSITION = 'ReadPosition'
    REFERENCEBASES = 'ReferenceBases'
    SAMPLELIST = 'SampleList'
    STRANDBIASBYSAMPLE = 'StrandBiasBySample'
    STRANDODDSRATIO = 'StrandOddsRatio'
    TANDEMREPEAT = 'TandemRepeat'
    UNIQUEALTREADCOUNT = 'UniqueAltReadCount'
    VARIANTTYPE = 'VariantType'



class AnnotationGroup(StrEnum):
    AS_STANDARDANNOTATION = 'AS_StandardAnnotation'
    ALLELESPECIFICANNOTATION = 'AlleleSpecificAnnotation'
    GENOTYPEANNOTATION = 'GenotypeAnnotation'
    INFOFIELDANNOTATION = 'InfoFieldAnnotation'
    JUMBOGENOTYPEANNOTATION = 'JumboGenotypeAnnotation'
    JUMBOINFOANNOTATION = 'JumboInfoAnnotation'
    REDUCIBLEANNOTATION = 'ReducibleAnnotation'
    STANDARDANNOTATION = 'StandardAnnotation'
    STANDARDFLOWBASEDANNOTATION = 'StandardFlowBasedAnnotation'
    STANDARDHCANNOTATION = 'StandardHCAnnotation'
    STANDARDMUTECTANNOTATION = 'StandardMutectAnnotation'
    VARIANTANNOTATION = 'VariantAnnotation'



class AnnotationExclude(StrEnum):
    AS_STRANDBIASMUTECTANNOTATION = 'AS_StrandBiasMutectAnnotation'
    BASEQUALITY = 'BaseQuality'
    COVERAGE = 'Coverage'
    DEPTHPERALLELEBYSAMPLE = 'DepthPerAlleleBySample'
    DEPTHPERSAMPLEHC = 'DepthPerSampleHC'
    FRAGMENTDEPTHPERALLELEBYSAMPLE = 'FragmentDepthPerAlleleBySample'
    FRAGMENTLENGTH = 'FragmentLength'
    MAPPINGQUALITY = 'MappingQuality'
    ORIENTATIONBIASREADCOUNTS = 'OrientationBiasReadCounts'
    READPOSITION = 'ReadPosition'
    STRANDBIASBYSAMPLE = 'StrandBiasBySample'
    TANDEMREPEAT = 'TandemRepeat'



class DisableReadFilter(StrEnum):
    GOODCIGARREADFILTER = 'GoodCigarReadFilter'
    MAPPEDREADFILTER = 'MappedReadFilter'
    MAPPINGQUALITYAVAILABLEREADFILTER = 'MappingQualityAvailableReadFilter'
    MAPPINGQUALITYNOTZEROREADFILTER = 'MappingQualityNotZeroReadFilter'
    MAPPINGQUALITYREADFILTER = 'MappingQualityReadFilter'
    NONCHIMERICORIGINALALIGNMENTREADFILTER = 'NonChimericOriginalAlignmentReadFilter'
    NONZEROREFERENCELENGTHALIGNMENTREADFILTER = 'NonZeroReferenceLengthAlignmentReadFilter'
    NOTDUPLICATEREADFILTER = 'NotDuplicateReadFilter'
    NOTSECONDARYALIGNMENTREADFILTER = 'NotSecondaryAlignmentReadFilter'
    PASSESVENDORQUALITYCHECKREADFILTER = 'PassesVendorQualityCheckReadFilter'
    READLENGTHREADFILTER = 'ReadLengthReadFilter'
    WELLFORMEDREADFILTER = 'WellformedReadFilter'



class IntervalMergingRule(StrEnum):
    ALL = 'ALL'
    OVERLAPPING_ONLY = 'OVERLAPPING_ONLY'



class IntervalSetRule(StrEnum):
    INTERSECTION = 'INTERSECTION'
    UNION = 'UNION'



class ReadFilter(StrEnum):
    ALIGNMENTAGREESWITHHEADERREADFILTER = 'AlignmentAgreesWithHeaderReadFilter'
    ALLOWALLREADSREADFILTER = 'AllowAllReadsReadFilter'
    AMBIGUOUSBASEREADFILTER = 'AmbiguousBaseReadFilter'
    CIGARCONTAINSNONOPERATOR = 'CigarContainsNoNOperator'
    EXCESSIVEENDCLIPPEDREADFILTER = 'ExcessiveEndClippedReadFilter'
    FIRSTOFPAIRREADFILTER = 'FirstOfPairReadFilter'
    FLOWBASEDTPATTRIBUTESYMETRICREADFILTER = 'FlowBasedTPAttributeSymetricReadFilter'
    FLOWBASEDTPATTRIBUTEVALIDREADFILTER = 'FlowBasedTPAttributeValidReadFilter'
    FRAGMENTLENGTHREADFILTER = 'FragmentLengthReadFilter'
    GOODCIGARREADFILTER = 'GoodCigarReadFilter'
    HASREADGROUPREADFILTER = 'HasReadGroupReadFilter'
    HMERQUALITYSYMETRICREADFILTER = 'HmerQualitySymetricReadFilter'
    INTERVALOVERLAPREADFILTER = 'IntervalOverlapReadFilter'
    JEXLEXPRESSIONREADTAGVALUEFILTER = 'JexlExpressionReadTagValueFilter'
    LIBRARYREADFILTER = 'LibraryReadFilter'
    MAPPEDREADFILTER = 'MappedReadFilter'
    MAPPINGQUALITYAVAILABLEREADFILTER = 'MappingQualityAvailableReadFilter'
    MAPPINGQUALITYNOTZEROREADFILTER = 'MappingQualityNotZeroReadFilter'
    MAPPINGQUALITYREADFILTER = 'MappingQualityReadFilter'
    MATCHINGBASESANDQUALSREADFILTER = 'MatchingBasesAndQualsReadFilter'
    MATEDIFFERENTSTRANDREADFILTER = 'MateDifferentStrandReadFilter'
    MATEDISTANTREADFILTER = 'MateDistantReadFilter'
    MATEONSAMECONTIGORNOMAPPEDMATEREADFILTER = 'MateOnSameContigOrNoMappedMateReadFilter'
    MATEUNMAPPEDANDUNMAPPEDREADFILTER = 'MateUnmappedAndUnmappedReadFilter'
    METRICSREADFILTER = 'MetricsReadFilter'
    NONCHIMERICORIGINALALIGNMENTREADFILTER = 'NonChimericOriginalAlignmentReadFilter'
    NONZEROFRAGMENTLENGTHREADFILTER = 'NonZeroFragmentLengthReadFilter'
    NONZEROREFERENCELENGTHALIGNMENTREADFILTER = 'NonZeroReferenceLengthAlignmentReadFilter'
    NOTDUPLICATEREADFILTER = 'NotDuplicateReadFilter'
    NOTOPTICALDUPLICATEREADFILTER = 'NotOpticalDuplicateReadFilter'
    NOTPROPERLYPAIREDREADFILTER = 'NotProperlyPairedReadFilter'
    NOTSECONDARYALIGNMENTREADFILTER = 'NotSecondaryAlignmentReadFilter'
    NOTSUPPLEMENTARYALIGNMENTREADFILTER = 'NotSupplementaryAlignmentReadFilter'
    OVERCLIPPEDREADFILTER = 'OverclippedReadFilter'
    PAIREDREADFILTER = 'PairedReadFilter'
    PASSESVENDORQUALITYCHECKREADFILTER = 'PassesVendorQualityCheckReadFilter'
    PLATFORMREADFILTER = 'PlatformReadFilter'
    PLATFORMUNITREADFILTER = 'PlatformUnitReadFilter'
    PRIMARYLINEREADFILTER = 'PrimaryLineReadFilter'
    PROPERLYPAIREDREADFILTER = 'ProperlyPairedReadFilter'
    READGROUPBLACKLISTREADFILTER = 'ReadGroupBlackListReadFilter'
    READGROUPHASFLOWORDERREADFILTER = 'ReadGroupHasFlowOrderReadFilter'
    READGROUPREADFILTER = 'ReadGroupReadFilter'
    READLENGTHEQUALSCIGARLENGTHREADFILTER = 'ReadLengthEqualsCigarLengthReadFilter'
    READLENGTHREADFILTER = 'ReadLengthReadFilter'
    READNAMEREADFILTER = 'ReadNameReadFilter'
    READSTRANDFILTER = 'ReadStrandFilter'
    READTAGVALUEFILTER = 'ReadTagValueFilter'
    SAMPLEREADFILTER = 'SampleReadFilter'
    SECONDOFPAIRREADFILTER = 'SecondOfPairReadFilter'
    SEQISSTOREDREADFILTER = 'SeqIsStoredReadFilter'
    SOFTCLIPPEDREADFILTER = 'SoftClippedReadFilter'
    VALIDALIGNMENTENDREADFILTER = 'ValidAlignmentEndReadFilter'
    VALIDALIGNMENTSTARTREADFILTER = 'ValidAlignmentStartReadFilter'
    WELLFORMEDFLOWBASEDREADFILTER = 'WellformedFlowBasedReadFilter'
    WELLFORMEDREADFILTER = 'WellformedReadFilter'



class ValidationStringency(StrEnum):
    LENIENT = 'LENIENT'
    SILENT = 'SILENT'
    STRICT = 'STRICT'



class LogLevel(StrEnum):
    DEBUG = 'DEBUG'
    ERROR = 'ERROR'
    INFO = 'INFO'
    WARNING = 'WARNING'



class WriterType(StrEnum):
    ALL_POSSIBLE_HAPLOTYPES = 'ALL_POSSIBLE_HAPLOTYPES'
    CALLED_HAPLOTYPES = 'CALLED_HAPLOTYPES'
    CALLED_HAPLOTYPES_NO_READS = 'CALLED_HAPLOTYPES_NO_READS'
    NO_HAPLOTYPES = 'NO_HAPLOTYPES'



class ReferenceConfidenceMode(StrEnum):
    BP_RESOLUTION = 'BP_RESOLUTION'
    GVCF = 'GVCF'
    NONE = 'NONE'



class FlowMode(StrEnum):
    ADVANCED = 'ADVANCED'
    NONE = 'NONE'
    STANDARD = 'STANDARD'



class Implementation(StrEnum):
    FLOWBASED = 'FlowBased'
    FLOWBASEDHMM = 'FlowBasedHMM'
    PAIRHMM = 'PairHMM'



class PairHMMImplementation(StrEnum):
    AVX_LOGLESS_CACHING = 'AVX_LOGLESS_CACHING'
    AVX_LOGLESS_CACHING_OMP = 'AVX_LOGLESS_CACHING_OMP'
    EXACT = 'EXACT'
    FASTEST_AVAILABLE = 'FASTEST_AVAILABLE'
    LOGLESS_CACHING = 'LOGLESS_CACHING'
    ORIGINAL = 'ORIGINAL'



class PCRErrorModel(StrEnum):
    AGGRESSIVE = 'AGGRESSIVE'
    CONSERVATIVE = 'CONSERVATIVE'
    HOSTILE = 'HOSTILE'
    NONE = 'NONE'



class SmithWatermanImplementation(StrEnum):
    AVX_ENABLED = 'AVX_ENABLED'
    FASTEST_AVAILABLE = 'FASTEST_AVAILABLE'
    JAVA = 'JAVA'



class Mutect2(SnappyModel):
    # Required arguments
    # input: bams tumor + normal
    # output: raw vcf
    # reference: fasta

    # Arguments actually used

    genotype_germline_sites: bool = False
    """Call all apparent germline site even though they will ultimately be filtered"""
    # germline_resource: FeatureInput | None = None  # No options for class FeatureInput
    germline_resource: str | None = None  # No options for class FeatureInput
    """Population vcf of germline sequencing containing allele fractions"""

    # Arguments that must be set by derived classes (pon & calling)

    # Panel of normals arguments


    # Calling-specific arguments

    # panel_of_normals: str | None = None  # Was class FeatureInput
    # """VCF file of sites observed in normal"""
    # genotype_pon_sites: bool = False
    # """Call sites in the PoN even though they will ultimately be filtered"""
    # mitochondria_mode: bool = False
    # """Mitochondria mode sets emission and initial LODs to 0"""
    # add_output_sam_program_record: bool = True
    # """If true, adds a PG tag to created SAM/BAM/CRAM files"""
    # assembly_region_out: str | None = None
    # """Output the assembly region to this IGV formatted file"""
    # bam_writer_type: WriterType = WriterType.CALLED_HAPLOTYPES
    # """Which haplotypes should be written to the BAM"""
    # enable_all_annotations: bool = False
    # """Use all possible annotations (not for the faint of heart)"""
    # alleles: str | None = None  # Was class FeatureInput
    # """The set of alleles to force-call regardless of evidence"""
    # bam_output: bool = False  # Was class str
    # """Write assembled haplotypes"""
    # pair_hmm_results_file: bool = False  # Was class GATKPath
    # """Write exact pairHMM inputs/outputs to for debugging purposes"""


    # Optional arguments

    add_output_vcf_command_line: bool = True
    """If true, adds a command line header line to created VCF files"""
    af_of_alleles_not_in_resource: float = -1.0
    """Population allele fraction assigned to alleles not found in germline resource.  Please see docs/mutect/mutect2.pdf fora derivation of the default value"""
    annotation: list[Annotation] = []
    """One or more specific annotations to add to variant calls"""
    annotation_group: list[AnnotationGroup] = []
    """One or more groups of annotations to apply to variant calls"""
    annotations_to_exclude: list[AnnotationExclude] = []
    """One or more specific annotations to exclude from variant calls"""
    arguments_file: str | None = None  # was class File
    """read one or more arguments files and add them to the command line"""
    assembly_region_padding: int = 100
    """Number of additional bases of context to include around each assembly region"""
    base_quality_score_threshold: int = 18
    """Base qualities below this threshold will be reduced to the minimum (6)"""
    callable_depth: int = 10
    """Minimum depth to be considered callable for Mutect stats.  Does not affect genotyping"""
    disable_bam_index_caching: bool = False
    """If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified"""
    disable_read_filter: list[DisableReadFilter] = []
    """Read filters to be disabled before analysis"""
    dont_use_dragstr_pair_hmm_scores: bool = False
    """disable DRAGstr pair-hmm score even when dragstr-params-path was provided"""
    downsampling_stride: int = 1
    """Downsample a pool of reads starting within a range of one or more bases"""
    dragstr_het_hom_ratio: int = 2
    """het to hom prior ratio use with DRAGstr on"""
    enable_dynamic_read_disqualification_for_genotyping: bool = False
    """Will enable less strict read disqualification low base quality reads"""
    exclude_intervals: list[str] = []
    """One or more genomic intervals to exclude from processing"""
    f1r2_max_depth: int = 200
    """sites with depth higher than this value will be grouped"""
    f1r2_median_mq: int = 50
    """skip sites with median mapping quality below this value"""
    f1r2_min_bq: int = 20
    """exclude bases below this quality from pileup"""
    flow_order_for_annotations: list[str] = []
    """flow order used for this annotations. [readGroup:]flowOrder"""
    founder_id: list[str] = []
    """Samples representing the population 'founders'"""
    gatk_config_file: str | None = None
    """A configuration file to use with the GATK"""
    ignore_itr_artifacts: bool = False
    """Turn off read transformer that clips artifacts associated with end repair insertions near inverted tandem repeats"""
    initial_tumor_lod: float = 2.0
    """Log 10 odds threshold to consider pileup active"""
    interval_exclusion_padding: int = 0
    """Amount of padding (in bp) to add to each interval you are excluding"""
    interval_merging_rule: IntervalMergingRule = IntervalMergingRule.ALL
    """Interval merging rule for abutting intervals"""
    interval_padding: int = 0
    """Amount of padding (in bp) to add to each interval you are including"""
    interval_set_rule: IntervalSetRule = IntervalSetRule.UNION
    """Set merging approach to use for combining interval inputs"""
    intervals: list[str] = []
    """One or more genomic intervals over which to operate"""
    lenient: bool = False
    """Lenient processing of VCF files"""
    max_assembly_region_size: int = 300
    """Maximum size of an assembly region"""
    max_population_af: float = 0.01
    """Maximum population allele frequency in tumor-only mode"""
    max_reads_per_alignment_start: int = 50
    """Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable"""
    max_variants_per_shard: int = 0
    """If non-zero, partitions VCF output into shards, each containing up to the given number of records"""
    min_assembly_region_size: int = 50
    """Minimum size of an assembly region"""
    min_base_quality_score: int = 10
    """Minimum base quality required to consider a base for calling"""
    native_pair_hmm_use_double_precision: bool = False
    """use double precision in the native pairHmm. This is slower but matches the java implementation better"""
    normal_lod: float = 2.2
    """Log 10 odds threshold for calling normal variant non-germline"""
    pcr_indel_qual: int = 40
    """Phred-scaled PCR indel qual for overlapping fragments"""
    pcr_snv_qual: int = 40
    """Phred-scaled PCR SNV qual for overlapping fragments"""
    read_filter: list[ReadFilter] = []
    """Read filters to be applied before analysis"""
    read_validation_stringency: ValidationStringency = ValidationStringency.SILENT
    """Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded"""
    sites_only_vcf_output: bool = False
    """If true, don't emit genotype fields when writing vcf file output"""
    tumor_lod_to_emit: float = 3.0
    """Log 10 odds threshold to emit variant to VCF"""
    use_jdk_deflater: bool = False
    """Whether to use the JdkDeflater (as opposed to IntelDeflater)"""
    use_jdk_inflater: bool = False
    """Whether to use the JdkInflater (as opposed to IntelInflater)"""
    verbosity: LogLevel = LogLevel.INFO
    """Control verbosity of logging"""
    QUIET: bool = False
    """Whether to suppress job-summary info on System.err"""

    # Advanced arguments

    active_probability_threshold: float = 0.002
    """Minimum probability for a locus to be considered active"""
    adaptive_pruning_initial_error_rate: float = 0.001
    """Initial base error rate estimate for adaptive pruning"""
    allele_informative_reads_overlap_margin: int = 2
    """Likelihood and read-based annotations will only take into consideration reads that overlap the variant or any base no further than this distance expressed in base pairs"""
    allow_non_unique_kmers_in_ref: bool = False
    """Allow graphs that have non-unique kmers in the reference"""
    debug_assembly: bool = False
    """Print out verbose debug information about each assembly region"""
    disable_adaptive_pruning: bool = False
    """Disable the adaptive algorithm for pruning paths in the graph"""
    disable_cap_base_qualities_to_map_quality: bool = False
    """If false this disables capping of base qualities in the HMM to the mapping quality of the read"""
    disable_symmetric_hmm_normalizing: bool = False
    """Toggle to revive legacy behavior of asymmetrically normalizing the arguments to the reference haplotype"""
    disable_tool_default_annotations: bool = False
    """Disable all tool default annotations"""
    disable_tool_default_read_filters: bool = False
    """Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)"""
    dont_increase_kmer_sizes_for_cycles: bool = False
    """Disable iterating over kmer sizes when graph cycles are detected"""
    dragstr_params_path: str | None = None  # Was class GATKPath
    """location of the DRAGstr model parameters for STR error correction used in the Pair HMM. When provided, it overrides other PCR error correcting mechanisms"""
    emit_ref_confidence: ReferenceConfidenceMode = ReferenceConfidenceMode.NONE
    """Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature)"""
    expected_mismatch_rate_for_read_disqualification: float = 0.02
    """Error rate used to set expectation for post HMM read disqualification based on mismatches"""
    flow_assembly_collapse_partial_mode: bool = False
    """Collapse long flow-based hmers only up to difference in reference"""
    flow_disallow_probs_larger_than_call: bool = False
    """Cap probabilities of error to 1 relative to base call"""
    flow_fill_empty_bins_value: float = 0.001
    """Value to fill the zeros of the matrix with"""
    flow_filter_alleles: bool = False
    """pre-filter alleles before genotyping"""
    flow_filter_alleles_qual_threshold: float = 30.0
    """Threshold for prefiltering alleles on quality"""
    flow_filter_alleles_sor_threshold: float = 3.0
    """Threshold for prefiltering alleles on SOR"""
    flow_filter_lone_alleles: bool = False
    """Remove also lone alleles during allele filtering"""
    flow_lump_probs: bool = False
    """Should all probabilities of insertion or deletion in the flow be combined together"""
    flow_matrix_mods: str | None = None
    """Modifications instructions to the read flow matrix. Format is src,dst{,src,dst}+. Example: 10,12,11,12 - these instructions will copy element 10 into 11 and 12"""
    flow_mode: FlowMode = FlowMode.NONE
    """Single argument for enabling the bulk of Flow Based features. NOTE: THIS WILL OVERWRITE PROVIDED ARGUMENT CHECK TOOL INFO TO SEE WHICH ARGUMENTS ARE SET)"""
    flow_probability_scaling_factor: int = 10
    """probability scaling factor for (phred=10) for probability quantization"""
    flow_probability_threshold: float = 0.003
    """Lowest probability ratio to be used as an option"""
    flow_quantization_bins: int = 121
    """Number of bins for probability quantization"""
    flow_remove_non_single_base_pair_indels: bool = False
    """Should the probabilities of more then 1 indel be used"""
    flow_remove_one_zero_probs: bool = False
    """Remove probabilities of basecall of zero from non-zero genome"""
    flow_report_insertion_or_deletion: bool = False
    """Report either insertion or deletion, probability, not both"""
    flow_retain_max_n_probs_base_format: bool = False
    """Keep only hmer/2 probabilities (like in base format)"""
    flow_symmetric_indel_probs: bool = False
    """Should indel probabilities be symmetric in flow"""
    flow_use_t0_tag: bool = False
    """Use t0 tag if exists in the read to create flow matrix"""
    force_active: bool = False
    """If provided, all regions will be marked as active"""
    force_call_filtered_alleles: bool = False
    """Force-call filtered alleles included in the resource specified by --alleles"""
    graph_output: bool = False  # Was class str
    """Write debug assembly graph information to this file"""
    gvcf_lod_band: list[float] = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]
    """Exclusive upper bounds for reference confidence LOD bands (must be specified in increasing order)"""
    independent_mates: bool = False
    """Allow paired reads to independently support different haplotypes.  Useful for validations with ill-designed synthetic data"""
    keep_boundary_flows: bool = False
    """prevent spreading of boundary flows"""
    kmer_size: list[int] = [10, 25]
    """Kmer size to use in the read threading assembler"""
    likelihood_calculation_engine: Implementation = Implementation.PAIRHMM
    """What likelihood calculation engine to use to calculate the relative likelihood of reads vs haplotypes"""
    linked_de_bruijn_graph: bool = False
    """If enabled, the Assembly Engine will construct a Linked De Bruijn graph to recover better haplotypes"""
    max_mnp_distance: int = 1
    """Two or more phased substitutions separated by this distance or less are merged into MNPs"""
    max_num_haplotypes_in_population: int = 128
    """Maximum number of haplotypes to consider for your population"""
    max_prob_propagation_distance: int = 50
    """Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions"""
    max_suspicious_reads_per_alignment_start: int = 0
    """Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable"""
    max_unpruned_variants: int = 100
    """Maximum number of variants in graph the adaptive pruner will allow"""
    min_dangling_branch_length: int = 4
    """Minimum length of a dangling branch to attempt recovery"""
    min_pruning: int = 2
    """Minimum support to not prune paths in the graph"""
    minimum_allele_fraction: float = 0.0
    """Lower bound of variant allele fractions to consider when calculating variant LOD"""
    num_pruning_samples: int = 1
    """Number of samples that must pass the minPruning threshold"""
    pair_hmm_gap_continuation_penalty: int = 10
    """Flat gap continuation penalty for use in the Pair HMM"""
    pair_hmm_implementation: PairHMMImplementation = PairHMMImplementation.FASTEST_AVAILABLE
    """The PairHMM implementation to use for genotype likelihood calculations"""
    pcr_indel_model: PCRErrorModel = PCRErrorModel.CONSERVATIVE
    """The PCR indel model to use"""
    pedigree: str | None = None  # Was class GATKPath
    """Pedigree file for determining the population 'founders'"""
    phred_scaled_global_read_mismapping_rate: int = 45
    """The global assumed mismapping rate for reads"""
    pileup_detection: bool = False
    """If enabled, the variant caller will create pileup-based haplotypes in addition to the assembly-based haplotype generation"""
    pruning_lod_threshold: float = 2.302585092994046
    """Ln likelihood ratio threshold for adaptive pruning algorithm"""
    pruning_seeding_lod_threshold: float = 9.210340371976184
    """Ln likelihood ratio threshold for seeding subgraph of good variation in adaptive pruning algorithm"""
    recover_all_dangling_branches: bool = False
    """Recover all dangling branches"""
    reference_model_deletion_quality: int = 30
    """The quality of deletion in the reference model"""
    smith_waterman: SmithWatermanImplementation = SmithWatermanImplementation.JAVA
    """Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice"""
    smith_waterman_dangling_end_gap_extend_penalty: int = -6
    """Smith-Waterman gap-extend penalty for dangling-end recovery"""
    smith_waterman_dangling_end_gap_open_penalty: int = -110
    """Smith-Waterman gap-open penalty for dangling-end recovery"""
    smith_waterman_dangling_end_match_value: int = 25
    """Smith-Waterman match value for dangling-end recovery"""
    smith_waterman_dangling_end_mismatch_penalty: int = -50
    """Smith-Waterman mismatch penalty for dangling-end recovery"""
    smith_waterman_haplotype_to_reference_gap_extend_penalty: int = -11
    """Smith-Waterman gap-extend penalty for haplotype-to-reference alignment"""
    smith_waterman_haplotype_to_reference_gap_open_penalty: int = -260
    """Smith-Waterman gap-open penalty for haplotype-to-reference alignment"""
    smith_waterman_haplotype_to_reference_match_value: int = 200
    """Smith-Waterman match value for haplotype-to-reference alignment"""
    smith_waterman_haplotype_to_reference_mismatch_penalty: int = -150
    """Smith-Waterman mismatch penalty for haplotype-to-reference alignment"""
    smith_waterman_read_to_haplotype_gap_extend_penalty: int = -5
    """Smith-Waterman gap-extend penalty for read-to-haplotype alignment"""
    smith_waterman_read_to_haplotype_gap_open_penalty: int = -30
    """Smith-Waterman gap-open penalty for read-to-haplotype alignment"""
    smith_waterman_read_to_haplotype_match_value: int = 10
    """Smith-Waterman match value for read-to-haplotype alignment"""
    smith_waterman_read_to_haplotype_mismatch_penalty: int = -15
    """Smith-Waterman mismatch penalty for read-to-haplotype alignment"""
    soft_clip_low_quality_ends: bool = False
    """If enabled will preserve low-quality read ends as softclips (used for DRAGEN-GATK BQD genotyper model)"""

    maximum_mapping_quality: int | None = None
    """Maximum mapping quality to keep (inclusive)"""
    minimum_mapping_quality: int = 20
    """Minimum mapping quality to keep (inclusive)"""
    max_read_length: int = 2147483647
    """Keep only reads with length at most equal to the specified value"""
    min_read_length: int = 30
    """Keep only reads with length at least equal to the specified value"""

    # Arguments omitted

    # f1r2_tar_gz: File | None = None  # No options for class File
    # """If specified, collect F1R2 counts and output files into this tar.gz file"""
    # gcs_max_retries: int = 20
    # """If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection"""
    # gcs_project_for_requester_pays: str = "."
    # """Project to bill when accessing "requester pays" buckets. If unset, these buckets cannot be accessed.  User must have storage.buckets.get permission on the bucket being accessed"""
    # mutect3_alt_downsample: int = 20
    # """Downsample alt reads to this count for Mutect3 training datasets"""
    # mutect3_dataset: File | None = None  # No options for class File
    # """Destination for Mutect3 data collection"""
    # mutect3_non_artifact_ratio: int = 20
    # """Number of non-artifact data per artifact datum in Mutect3 training"""
    # mutect3_ref_downsample: int = 10
    # """Downsample ref reads to this count when generating a Mutect3 dataset"""
    # mutect3_training_mode: bool = False
    # """Collect Mutect3 data for learning"""
    # mutect3_training_truth: FeatureInput | None = None  # No options for class FeatureInput
    # """VCF file of known variants for labeling Mutect3 training data"""
    # native_pair_hmm_threads: int = 4
    # """How many threads should a native pairHMM implementation use"""
    # read_index: list[GATKPath] = []  # No options for class GATKPath
    # """Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically"""
    # sequence_dictionary: GATKPath | None = None  # No options for class GATKPath
    # """Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file"""
    # tmp_dir: GATKPath | None = None  # No options for class GATKPath
    # """Temp directory to use"""
    # cloud_index_prefetch_buffer: int = -1
    # """Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset"""
    # cloud_prefetch_buffer: int = 40
    # """Size of the cloud-only prefetch buffer (in MB; 0 to disable)"""
    # create_output_bam_index: bool = True
    # """If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file"""
    # create_output_bam_md5: bool = False
    # """If true, create a MD5 digest for any BAM/SAM/CRAM file created"""
    # create_output_variant_index: bool = True
    # """If true, create a VCF index when writing a coordinate-sorted VCF file"""
    # create_output_variant_md5: bool = False
    # """If true, create a a MD5 digest any VCF file created"""
    # disable_sequence_dictionary_validation: bool = False
    # """If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!"""
    # help: bool = False
    # """display the help message"""
    # seconds_between_progress_updates: float = 10.0
    # """Output traversal statistics every time this many seconds elapse"""
    # version: bool = False
    # """display the version number for this tool"""
    # showHidden: bool = False
    # """display hidden arguments"""
    # normal_sample: list[str] = []
    # """BAM sample name of normal(s), if any.  May be URL-encoded as output by GetSampleName with -encode argument"""
    # tumor_sample: str | None = None
    # """BAM sample name of tumor.  May be URL-encoded as output by GetSampleName with -encode argument"""

