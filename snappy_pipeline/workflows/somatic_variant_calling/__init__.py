# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_calling`` step

The ``somatic_variant_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs somatic variant calling.  The result are variant files
with somatic variants (bgzip-ed and indexed VCF files).

Usually, the somatic variant calling step is followed by the ``somatic_variant_annotation`` step.

==========
Step Input
==========

The somatic variant calling step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk`` and each read mapper
``mapper`` that the library has been aligned with, and the variant caller ``var_caller``, the
pipeline step will create a directory ``output/{mapper}.{var_caller}.{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.mutect.P001-N1-DNA1-WES1-4
    |   `-- out
    |       |-- bwa.mutect.P001-N1-DNA1-WES1-4.vcf.gz
    |       |-- bwa.mutect.P001-N1-DNA1-WES1-4.vcf.gz.tbi
    |       |-- bwa.mutect.P001-N1-DNA1-WES1-4.vcf.gz.md5
    |       `-- bwa.mutect.P001-N1-DNA1-WES1-4.vcf.gz.tbi.md5
    [...]

Generally, these files will be unfiltered, i.e., contain low-quality variants and also variants
flagged as being non-somatic.

====================
Global Configuration
====================

- If the somatic variant caller MuTect is used, then the global settings
  ``static_data_config/dbsnp`` and ``static_data_config/cosmic`` must be given
  as MuTect uses this in its algorithm.
- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_calling.rst

=================================
Available Somatic Variant Callers
=================================

The following somatic variant callers are currently available

- ``"mutect"``
- ``"scalpel"``

=======
Reports
=======

Currently, no reports are generated.
"""

from collections import OrderedDict
from itertools import chain
import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

EXTS_MATCHED = {
    "mutect": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "tbi": ".vcf.gz.tbi",
        "tbi_md5": ".vcf.gz.tbi.md5",
        "full": ".full.vcf.gz",
        "full_md5": ".full.vcf.gz.md5",
        "full_tbi": ".full.vcf.gz.tbi",
        "full_tbi_md5": ".full.vcf.gz.tbi.md5",
        "txt": ".txt",
        "txt_md5": ".txt.md5",
        "wig": ".wig",
        "wig_md5": ".wig.md5",
    },
    "scalpel": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "tbi": ".vcf.gz.tbi",
        "tbi_md5": ".vcf.gz.tbi.md5",
        "full": ".full.vcf.gz",
        "full_md5": ".full.vcf.gz.md5",
        "full_tbi": ".full.vcf.gz.tbi",
        "full_tbi_md5": ".full.vcf.gz.tbi.md5",
        "tar": ".tar.gz",
    },
    "mutect2": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "tbi": ".vcf.gz.tbi",
        "tbi_md5": ".vcf.gz.tbi.md5",
        "full": ".full.vcf.gz",
        "full_md5": ".full.vcf.gz.md5",
        "full_tbi": ".full.vcf.gz.tbi",
        "full_tbi_md5": ".full.vcf.gz.tbi.md5",
        "stats": ".vcf.stats",
        "stats_md5": ".vcf.stats.md5",
        "f1r2": ".f1r2_tar.tar.gz",
        "f1r2_md5": ".f1r2_tar.tar.gz.md5",
    },
}

#: Available somatic variant callers assuming matched samples.
SOMATIC_VARIANT_CALLERS_MATCHED = ("mutect", "mutect2", "scalpel")

#: Available somatic variant callers that just call all samples from one donor together.
SOMATIC_VARIANT_CALLERS_JOINT = (
    "bcftools_joint",
    "platypus_joint",
    "gatk_hc_joint",
    "gatk_ug_joint",
    "varscan_joint",
)

#: Available somatic variant callers
SOMATIC_VARIANT_CALLERS = tuple(
    chain(SOMATIC_VARIANT_CALLERS_MATCHED, SOMATIC_VARIANT_CALLERS_JOINT)
)

#: Available somatic variant callers assuming matched samples.
SOMATIC_VARIANT_CALLERS_MATCHED = ("mutect", "mutect2", "scalpel", "strelka2")

#: Available somatic variant callers that just call all samples from one donor together.
SOMATIC_VARIANT_CALLERS_JOINT = (
    "bcftools_joint",
    "platypus_joint",
    "gatk_hc_joint",
    "gatk_ug_joint",
    "varscan_joint",
)

#: Available somatic variant callers
SOMATIC_VARIANT_CALLERS = tuple(
    chain(SOMATIC_VARIANT_CALLERS_MATCHED, SOMATIC_VARIANT_CALLERS_JOINT)
)

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = r"""
# Default configuration somatic_variant_calling
step_config:
  somatic_variant_calling:
    drmaa_snippet: ''  # default, you can override by step below
    tools: ['mutect', 'scalpel']
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    # Configuration for joint calling with samtools+bcftools.
    bcftools_joint:
      max_depth: 4000
      max_indel_depth: 4000
      window_length: 10000000
      num_threads: 16
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
    # Configuration for joint calling with Platypus.
    platypus_joint:
      split_complex_mnvs: true  # whether or not to split complex and MNV variants
      num_threads: 16
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
    # VCF annotation databases are given as mapping from name to
    #   {'file': '/path.vcf.gz',
    #    'info_tag': 'VCF_TAG',
    #    'description': 'VCF header description'}
    # Configuration for MuTect
    mutect:
      # Parallelization configuration
      drmaa_snippet: ''          # value to pass in as additional DRMAA arguments
      num_cores: 2               # number of cores to use locally
      window_length: 3500000     # split input into windows of this size, each triggers a job
      num_jobs: 500              # number of windows to process in parallel
      use_drmaa: true            # use drmaa for parallel processing
      restart_times: 5           # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2     # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0      # truncation to first N tokens (0 for none)
      keep_tmpdir: never         # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1         # memory multiplier
      job_mult_time: 1           # running time multiplier
      merge_mult_memory: 1       # memory multiplier for merging
      merge_mult_time: 1         # running time multiplier for merging
      ignore_chroms:             # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
    # Configuration for MuTect 2
    mutect2:
      panel_of_normals: ''      # Set path to panel of normals vcf if required
      germline_resource: REQUIRED # Germline variants resource (same as panel of normals)
      common_variants: REQUIRED # Common germline variants for contamination estimation
      # Parallelization configuration
      drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
      num_cores: 2              # number of cores to use locally
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_drmaa: true           # use DRMAA for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2    # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
    # Configuration for Scalpel
    scalpel:
      path_target_regions: REQUIRED  # REQUIRED
    # Configuration for strelka2
    strelka2:
      path_target_regions: ""   # For exomes: include a bgzipped bed file with tabix index. That also triggers the --exome flag
    gatk_hc_joint:
      # Parallelization configuration
      drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
      num_cores: 2              # number of cores to use locally
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_drmaa: true           # use DRMAA for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
      # GATK HC--specific configuration
      allow_seq_dict_incompatibility: false
      annotations:
      - BaseQualityRankSumTest
      - FisherStrand
      - GCContent
      - HaplotypeScore
      - HomopolymerRun
      - MappingQualityRankSumTest
      - MappingQualityZero
      - QualByDepth
      - ReadPosRankSumTest
      - RMSMappingQuality
      - DepthPerAlleleBySample
      - Coverage
      - ClippingRankSumTest
      - DepthPerSampleHC
    gatk_ug_joint:
      # Parallelization configuration
      drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
      num_cores: 2              # number of cores to use locally
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_drmaa: true           # use DRMAA for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
      # GATK UG--specific configuration
      downsample_to_coverage: 250
      allow_seq_dict_incompatibility: false
      annotations:
      - BaseQualityRankSumTest
      - FisherStrand
      - GCContent
      - HaplotypeScore
      - HomopolymerRun
      - MappingQualityRankSumTest
      - MappingQualityZero
      - QualByDepth
      - ReadPosRankSumTest
      - RMSMappingQuality
      - DepthPerAlleleBySample
      - Coverage
      - ClippingRankSumTest
      - DepthPerSampleHC
    varscan_joint:
      # Parallelization configuration
      drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
      num_cores: 2              # number of cores to use locally
      window_length: 5000000    # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_drmaa: true           # use drmaa for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2    # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
      # Configuration for samtools mpileup
      max_depth: 4000
      max_indel_depth: 4000
      min_bq: 13
      no_baq: True
      # Configuration for Varscan
      min_coverage: 8
      min_reads2: 2
      min_avg_qual: 15
      min_var_freq: 0.01
      min_freq_for_hom: 0.75
      p_value: 99e-02
"""


class SomaticVariantCallingStepPart(BaseStepPart):
    """Base class for somatic variant calling step parts

    Variant calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/out/"
            "{{mapper}}.{var_caller}.{{tumor_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            tumor_base_path = (
                "output/{mapper}.{tumor_library}/out/" "{mapper}.{tumor_library}"
            ).format(**wildcards)
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            }

        assert action == "run", "Unsupported actions"
        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_tumor_lib_name(self, wildcards):
        """Return name of tumor library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.tumor_sample.dna_ngs_library.name

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        assert action == "run"
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, var_caller=[self.name], ext=EXT_VALUES))
        )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        prefix = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/log/"
            "{{mapper}}.{var_caller}.{{tumor_library}}"
        ).format(var_caller=self.__class__.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class MutectBaseStepPart(SomaticVariantCallingStepPart):
    """Base class for Mutect 1 and 2 step parts"""

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "cosmic", "path"),
            "COSMIC not configured but required for %s" % (self.name,),
        )
        self.parent.ensure_w_config(
            ("static_data_config", "dbsnp", "path"),
            "Path to dbSNP not configured but required for %s" % (self.name,),
        )
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_output_files(self, action):
        output_files = {}
        for k, v in EXTS_MATCHED[self.name].items():
            output_files[k] = self.base_path_out.format(var_caller=self.name, ext=v)
        return output_files

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_%s_run" % self.name] = {
            "mem": int(3.7 * 1024 * 2),
            "time": "72:00",
            "ntasks": 2,
        }


class MutectStepPart(MutectBaseStepPart):
    """Somatic variant calling with MuTect"""

    name = "mutect"


class Mutect2StepPart(MutectBaseStepPart):
    """Somatic variant calling with MuTect 2"""

    name = "mutect2"

    def __init__(self, parent):
        super().__init__(parent)
        self.actions = ("run", "filter", "contamination", "pileup_normal", "pileup_tumor")

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_input_files(self, action):
        def input_function_run(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            tumor_base_path = (
                "output/{mapper}.{tumor_library}/out/" "{mapper}.{tumor_library}"
            ).format(**wildcards)
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            }

        def input_function_filter(wildcards):
            base_path = (
                "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}".format(
                    **wildcards
                )
            )
            return {
                "raw": base_path + ".raw.vcf.gz",
                "stats": base_path + ".raw.vcf.stats",
                "f1r2": base_path + ".raw.f1r2_tar.tar.gz",
                "table": base_path + ".contamination.tbl",
                "segments": base_path + ".segments.tbl",
            }

        def input_function_contamination(wildcards):
            base_path = (
                "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}".format(
                    **wildcards
                )
            )
            return {"normal": base_path + ".normal.pileup", "tumor": base_path + ".tumor.pileup"}

        def input_function_pileup_normal(wildcards):
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            base_path = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                normal_library=self.get_normal_lib_name(wildcards), **wildcards
            )
            return {"bam": ngs_mapping(base_path + ".bam"), "bai": ngs_mapping(base_path + ".bam")}

        def input_function_pileup_tumor(wildcards):
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            base_path = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}".format(
                **wildcards
            )
            return {"bam": ngs_mapping(base_path + ".bam"), "bai": ngs_mapping(base_path + ".bam")}

        assert action in self.actions
        if action == "run":
            return input_function_run
        if action == "filter":
            return input_function_filter
        if action == "contamination":
            return input_function_contamination
        if action == "pileup_normal":
            return input_function_pileup_normal
        if action == "pileup_tumor":
            return input_function_pileup_tumor

    def get_output_files(self, action):
        assert action in self.actions
        if action == "run":
            exts = {
                "raw": ".raw.vcf.gz",
                "raw_md5": ".raw.vcf.gz.md5",
                "raw_tbi": ".raw.vcf.gz.tbi",
                "raw_tbi_md5": ".raw.vcf.gz.tbi.md5",
                "stats": ".raw.vcf.stats",
                "stats_md5": ".raw.vcf.stats.md5",
                "f1r2": ".raw.f1r2_tar.tar.gz",
                "f1r2_md5": ".raw.f1r2_tar.tar.gz.md5",
            }
        if action == "filter":
            exts = {
                "full": ".full.vcf.gz",
                "full_md5": ".full.vcf.gz.md5",
                "full_tbi": ".full.vcf.gz.tbi",
                "full_tbi_md5": ".full.vcf.gz.tbi.md5",
                "vcf": ".vcf.gz",
                "vcf_md5": ".vcf.gz.md5",
                "tbi": ".vcf.gz.tbi",
                "tbi_md5": ".vcf.gz.tbi.md5",
            }
        if action == "contamination":
            exts = {
                "table": ".contamination.tbl",
                "table_md5": ".contamination.tbl.md5",
                "segments": ".segments.tbl",
                "segments_md5": ".segments.tbl.md5",
            }
        if action == "pileup_normal":
            exts = {"pileup": ".normal.pileup", "pileup_md5": ".normal.pileup.md5"}
        if action == "pileup_tumor":
            exts = {"pileup": ".tumor.pileup", "pileup_md5": ".tumor.pileup.md5"}
        output_files = {}
        for k, v in exts.items():
            output_files[k] = self.base_path_out.format(var_caller=self.name, ext=v)
        return output_files

    def get_log_file(self, action):
        assert action in self.actions
        postfix = ""
        if action != "run":
            postfix = "." + action
        prefix = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/log/"
            "{{mapper}}.{var_caller}.{{tumor_library}}{postfix}"
        ).format(var_caller=self.__class__.name, postfix=postfix)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        log_files = {}
        for key, ext in key_ext:
            log_files[key] = prefix + ext
            log_files[key + "_md5"] = prefix + ext + ".md5"
        return log_files

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_mutect2_run"] = {
            "h_vmem": "4g",
            "h_rt": "120:00:00",
            "pe": "smp 2",
        }
        cluster_config["somatic_variant_calling_mutect2_filter"] = {
            "h_vmem": "8g",
            "h_rt": "3:59:00",
            "pe": "smp 2",
        }
        cluster_config["somatic_variant_calling_mutect2_contamination"] = {
            "h_vmem": "8g",
            "h_rt": "3:59:00",
            "pe": "smp 2",
        }
        cluster_config["somatic_variant_calling_mutect2_pileup_normal"] = {
            "h_vmem": "8g",
            "h_rt": "3:59:00",
            "pe": "smp 2",
        }
        cluster_config["somatic_variant_calling_mutect2_pileup_tumor"] = {
            "h_vmem": "8g",
            "h_rt": "3:59:00",
            "pe": "smp 2",
        }


class ScalpelStepPart(SomaticVariantCallingStepPart):
    """Somatic variant calling with Scalpel"""

    name = "scalpel"

    def check_config(self):
        if "scalpel" not in self.config["tools"]:
            return  # scalpel not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for Scalpel",
        )

    def get_output_files(self, action):
        result = super().get_output_files(action)
        somatic_ext_names = expand("full_{name}", name=EXT_NAMES)
        somatic_ext_values = expand(".full{ext}", ext=EXT_VALUES)
        result["tar"] = self.base_path_out.format(var_caller="scalpel", ext=".tar.gz")
        result.update(
            dict(
                zip(
                    somatic_ext_names,
                    expand(self.base_path_out, var_caller=[self.name], ext=somatic_ext_values),
                )
            )
        )
        return result

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_scalpel_run"] = {
            "mem": 5 * 1024 * 16,
            "time": "48:00",
            "ntasks": 16,
        }


class Strelka2StepPart(SomaticVariantCallingStepPart):
    """Somatic variant calling with strelka2/manta"""

    name = "strelka2"

    extensions = {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "tbi": ".vcf.gz.tbi",
        "tbi_md5": ".vcf.gz.tbi.md5",
        "full": ".full.vcf.gz",
        "full_md5": ".full.vcf.gz.md5",
        "full_tbi": ".full.vcf.gz.tbi",
        "full_tbi_md5": ".full.vcf.gz.tbi.md5",
        "stats": ".tsv",
        "stats_md5": ".tsv.md5",
        "report": ".xml",
        "report_md5": ".xml.md5",
        "bed": ".bed.gz",
        "bed_md5": ".bed.gz.md5",
        "bed_tbi": ".bed.gz.tbi",
        "bed_tbi_md5": ".bed.gz.tbi.md5",
    }

    def get_output_files(self, action):
        output_files = {}
        for k, v in self.extensions.items():
            output_files[k] = self.base_path_out.format(var_caller=self.name, ext=v)
        return output_files

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_strelka2_run"] = {
            "h_vmem": "4g",
            "h_rt": "24:00:00",
            "pe": "smp 8",
        }


class JointCallingStepPart(BaseStepPart):
    """Base class for joint calling."""

    #: Name of the step, to be overridden in sub class.
    name = None

    def __init__(self, parent):
        super().__init__(parent)
        if self.__class__.name is None:
            raise RuntimeError("Step name not given, override in sub class")
        self.base_path_out = (
            "work/{{mapper}}.{name}.{{donor_name}}/out/" "{{mapper}}.{name}.{{donor_name}}{ext}"
        )
        # Build shortcut from donor name to donor.
        self.donor_by_name = OrderedDict()
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for donor in sheet.donors:
                self.donor_by_name[donor.name] = donor

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            donor = self.donor_by_name[wildcards.donor_name]
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Return paths of NGS libraries.
            input_base_path = "output/{mapper}.{library.name}/out/{mapper}.{library.name}{ext}"
            result = {"bam": [], "bai": []}
            for bio_sample in donor.bio_samples.values():
                for test_sample in bio_sample.test_samples.values():
                    for ngs_library in test_sample.ngs_libraries.values():
                        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
                            result[key].append(
                                ngs_mapping(
                                    input_base_path.format(
                                        library=ngs_library, ext=ext, **wildcards
                                    )
                                )
                            )
            return result

        assert action == "run", "Unsupported actions"
        return input_function

    def get_output_files(self, action):
        """Return resulting files generated by this step."""
        assert action == "run"
        return dict(zip(EXT_NAMES, expand(self.base_path_out, name=[self.name], ext=EXT_VALUES)))

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        prefix = (
            "work/{{mapper}}.{var_caller}.{{donor_name}}/log/"
            "{{mapper}}.{var_caller}.{{donor_name}}"
        ).format(var_caller=self.__class__.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_args(self, action):
        def arg_function(wildcards):
            donor = self.donor_by_name[wildcards.donor_name]
            result = {}
            result["sample_list"] = [
                ngs_library.name
                for bio_sample in donor.bio_samples.values()
                for test_sample in bio_sample.test_samples.values()
                for ngs_library in test_sample.ngs_libraries.values()
            ]
            if "ignore_chroms" in self.parent.config[self.name]:
                ignore_chroms = self.parent.config[self.name]["ignore_chroms"]
                result["ignore_chroms"] = ignore_chroms
            return result

        return arg_function


class BcftoolsJointStepPart(JointCallingStepPart):
    """Somatic variant calling with Samtools in "joint" mode.

    Simply pass all samples of one donor through Samtools. The resulting files are then to be
    subfiltered.
    """

    name = "bcftools_joint"

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_bcftools_joint_run"] = {
            "mem": 1024 * self.parent.config["bcftools_joint"]["num_threads"],
            "time": "48:00",
            "ntasks": self.parent.config["bcftools_joint"]["num_threads"],
        }


class VarscanJointStepPart(JointCallingStepPart):
    """Somatic variant calling with Varscan in "joint" mode."""

    name = "varscan_joint"

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_varscan_joint_run"] = {
            "mem": 1024,
            "time": "48:00",
            "ntasks": 1,
        }


class PlatypusJointStepPart(JointCallingStepPart):
    """Somatic variant calling with Platypus in "joint" mode.

    Simply pass all samples of one donor through Platypus. The resulting files are then to be
    subfiltered.
    """

    name = "platypus_joint"

    def update_cluster_config(self, cluster_config):
        cluster_config["somatic_variant_calling_platypus_joint_run"] = {
            "mem": int(3.75 * 1024 * self.parent.config["platypus_joint"]["num_threads"]),
            "time": "48:00",
            "ntasks": self.parent.config["platypus_joint"]["num_threads"],
        }


class GatkHcJointStepPart(JointCallingStepPart):
    """Somatic variant calling with GATK HC in "joint" mode.

    Simply pass all samples of one donor through GATK HC. The resulting files are then to be
    subfiltered.
    """

    name = "gatk_hc_joint"

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_calling_{}_run".format(self.__class__.name)] = {
            "mem": 2 * 1024,
            "time": "80:00",
            "ntasks": 1,
        }


class GatkUgJointStepPart(JointCallingStepPart):
    """Somatic variant calling with GATK UG in "joint" mode.

    Simply pass all samples of one donor through GATK UG. The resulting files are then to be
    subfiltered.
    """

    name = "gatk_ug_joint"

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_calling_{}_run".format(self.__class__.name)] = {
            "mem": 2 * 1024,
            "time": "80:00:00",
            "ntasks": 1,
        }


class SomaticVariantCallingWorkflow(BaseStep):
    """Perform somatic variant calling"""

    name = "somatic_variant_calling"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow,
            config,
            cluster_config,
            config_lookup_paths,
            config_paths,
            workdir,
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                MutectStepPart,
                Mutect2StepPart,
                ScalpelStepPart,
                Strelka2StepPart,
                BcftoolsJointStepPart,
                GatkHcJointStepPart,
                GatkUgJointStepPart,
                PlatypusJointStepPart,
                VarscanJointStepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.{caller}.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=set(self.config["tools"]) & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=set(self.config["tools"]) & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )
        # Panel of normals
        # joint calling
        name_pattern = "{mapper}.{caller}.{donor.name}"
        yield from self._yield_result_files_joint(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=set(self.config["tools"]) & set(SOMATIC_VARIANT_CALLERS_JOINT),
            ext=EXT_VALUES,
        )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl, tumor_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def _yield_result_files_joint(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the joint somatic variant callers such as
        "Bcftools joint".
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for donor in sheet.donors:
                yield from expand(tpl, donor=[donor], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "somatic_variant_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for somatic variant calling",
        )
