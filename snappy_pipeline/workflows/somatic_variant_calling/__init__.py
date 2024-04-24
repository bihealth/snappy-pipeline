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

Two ``vcf`` files are produced:

- ``{mapper}.{var_caller}.{lib_name}.vcf.gz`` which contains only the variants that have passed all filters, or that were protected, and
- ``{mapper}.{var_caller}.{lib_name}.full.vcf.gz`` which contains all variants, with the reason for rejection in the ``FILTER`` column.

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.mutect2.P001-N1-DNA1-WES1-4
    |   `-- out
    |       |-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz
    |       |-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz.tbi
    |       |-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz.md5
    |       `-- bwa.mutect2.P001-N1-DNA1-WES1-4.vcf.gz.tbi.md5
    [...]

Generally, the callers are set up to return many variants to avoid missing clinically important ones.
They are likely to contain low-quality variants and for some callers variants flagged as being non-somatic.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_calling.rst

=================================
Available Somatic Variant Callers
=================================

The following somatic variant callers are currently available

- ``mutect2`` is the recommended caller
- ``mutect`` & ``scalpel`` are deprecated
- the joint variant callers ``bcftools``, ``platypus``, ``gatk`` & ``varscan`` are unsupported.

==========================
Efficient usage of mutect2
==========================

The recommended caller ``mutect2`` is implemented using a parallelisation mechanism to reduce execution time.
During parallelisation, the genome is split into small chunks, and parallel jobs perform the somatic variant calling in their region only.
All results are then assembled to generate the final output files.

There is at least one chunk by contig defined in the reference genome, with the chunk size upper limit given by the ``window_length`` configuration option.
So large chromosomes can be split into several chunks, while smaller contigs are left in one chunk.
Even for large chunk size, this parallelisation can create hundreds of jobs when the reference genome contains many contigs
(unplaced or unlocalized contigs, viral sequences, decoys, HLA alleles, ...).
Somatic variant calling is generally meaningless for many of these small contigs.
It is possible to configure the somatic variant calling to avoid the contigs irrelevant for downstream analysis, for example:

.. code-block:: yaml

  mutect2:
    ignore_chroms: ['*_random', 'chrUn_*', '*_decoy', "EBV", "HPV*", "HBV", "HCV-*", "HIV-*", "HTLV-1", "CMV", "KSHV", "MCV", "SV40"] # GRCh38.d1.vd1
    window_length: 300000000   # OK for WES, too large for WGS
    keep_tmpdir: onerror       # Facilitates debugging
    ...


=======
Reports
=======

Currently, no reports are generated.
"""

import os
import sys
from collections import OrderedDict
from itertools import chain

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

EXT_MATCHED = {
    "mutect": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
        "txt": ".txt",
        "txt_md5": ".txt.md5",
        "wig": ".wig",
        "wig_md5": ".wig.md5",
    },
    "scalpel": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
        "tar": ".tar.gz",
        "tar_md5": ".tar.gz.md5",
    },
    "mutect2": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
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

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = r"""
# Default configuration somatic_variant_calling
step_config:
  somatic_variant_calling:
    tools: ['mutect', 'scalpel']  # REQUIRED, examples: 'mutect' and 'scalpel'.
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    ignore_chroms:            # patterns of chromosome names to ignore
    - NC_007605    # herpes virus
    - hs37d5       # GRCh37 decoy
    - chrEBV       # Eppstein-Barr Virus
    - '*_decoy'    # decoy contig
    - 'HLA-*'      # HLA genes
    - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
    # Configuration for joint calling with samtools+bcftools.
    bcftools_joint:
      max_depth: 4000
      max_indel_depth: 4000
      window_length: 10000000
      num_threads: 16
    # Configuration for joint calling with Platypus.
    platypus_joint:
      split_complex_mnvs: true  # whether or not to split complex and MNV variants
      num_threads: 16
    # VCF annotation databases are given as mapping from name to
    #   {'file': '/path.vcf.gz',
    #    'info_tag': 'VCF_TAG',
    #    'description': 'VCF header description'}
    # Configuration for MuTect
    mutect:
      # Parallelization configuration
      num_cores: 2               # number of cores to use locally
      window_length: 3500000     # split input into windows of this size, each triggers a job
      num_jobs: 500              # number of windows to process in parallel
      use_profile: true          # use Snakemake profile for parallel processing
      restart_times: 5           # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2     # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0      # truncation to first N tokens (0 for none)
      keep_tmpdir: never         # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1         # memory multiplier
      job_mult_time: 1           # running time multiplier
      merge_mult_memory: 1       # memory multiplier for merging
      merge_mult_time: 1         # running time multiplier for merging
    # Configuration for MuTect 2
    mutect2:
      panel_of_normals: ''      # Set path to panel of normals vcf if required
      germline_resource: ''     # Germline variants resource (same as panel of normals)
      common_variants: ''       # Common germline variants for contamination estimation
      extra_arguments: []       # List additional Mutect2 arguments
                                # Each additional argument xust be in the form:
                                # "--<argument name> <argument value>"
                                # For example, to filter reads prior to calling & to
                                # add annotations to the output vcf:
                                # - "--read-filter CigarContainsNoNOperator"
                                # - "--annotation AssemblyComplexity BaseQuality"
      # Parallelization configuration
      num_cores: 2              # number of cores to use locally
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2    # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
    # Configuration for Scalpel
    scalpel:
      path_target_regions: REQUIRED  # REQUIRED
    # Configuration for strelka2
    strelka2:
      path_target_regions: ""   # For exomes: include a bgzipped bed file with tabix index. That also triggers the --exome flag
    gatk_hc_joint:
      # Parallelization configuration
      num_cores: 2              # number of cores to use locally
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
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
      num_cores: 2              # number of cores to use locally
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
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
      num_cores: 2              # number of cores to use locally
      window_length: 5000000    # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2    # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
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
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shortcut to Snakemake module for ngs_mapping
            ngs_mapping = self.parent.modules["ngs_mapping"]

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
        # Validate action
        self._validate_action(action)
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, var_caller=[self.name], ext=EXT_VALUES))
        )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)

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
        for k, v in EXT_MATCHED[self.name].items():
            output_files[k] = self.base_path_out.format(var_caller=self.name, ext=v)
        return output_files

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="3-00:00:00",  # 3 days
            memory=f"{int(3.7 * 1024 * 2)}M",
        )


class MutectStepPart(MutectBaseStepPart):
    """Somatic variant calling with MuTect"""

    #: Step name
    name = "mutect"

    #: Class available actions
    actions = ("run",)


class Mutect2StepPart(MutectBaseStepPart):
    """Somatic variant calling with Mutect2"""

    #: Step name
    name = "mutect2"

    #: Class available actions
    actions = ["run", "filter"]  # "contamination", "pileup_normal", "pileup_tumor")

    #: Class resource usage dictionary. Key: action (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "run": ResourceUsage(
            threads=2,
            time="5-00:00:00",
            memory="3584M",
        ),
        "filter": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="15872M",
        ),
        "contamination": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="7680M",
        ),
        "pileup_normal": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="8000M",
        ),
        "pileup_tumor": ResourceUsage(
            threads=2,
            time="03:59:00",
            memory="8000M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )
        if self.config[self.name]["common_variants"]:
            self.actions.extend(["contamination", "pileup_normal", "pileup_tumor"])

    def get_input_files(self, action):
        """Return input function for Mutect2 rules.

        :param action: Action (i.e., step) in the workflow, examples: 'run', 'filter',
        'contamination'.
        :type action: str

        :return: Returns input function for Mutect2 rules based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        # Return requested function
        return getattr(self, "_get_input_files_{}".format(action))

    def _get_input_files_run(self, wildcards):
        """Get input files for rule ``run``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'run', BAM and BAI files.
        """
        # Get shortcut to Snakemake module for ngs_mapping
        ngs_mapping = self.parent.modules["ngs_mapping"]

        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        normal_base_path = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
            normal_library=self.get_normal_lib_name(wildcards), **wildcards
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

    def _get_input_files_filter(self, wildcards):
        """Get input files for rule ``filter``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'filter'.
        """
        base_path = (
            "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}".format(
                **wildcards
            )
        )
        input_files = {
            "raw": base_path + ".raw.vcf.gz",
            "stats": base_path + ".raw.vcf.stats",
            "f1r2": base_path + ".raw.f1r2_tar.tar.gz",
        }
        if "contamination" in self.actions:
            input_files["table"] = base_path + ".contamination.tbl"
            input_files["segments"] = base_path + ".segments.tbl"
        return input_files

    def _get_input_files_pileup_normal(self, wildcards):
        """Get input files for rule ``pileup_normal``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'pileup_normal', BAM and BAI files.
        """
        # Get shortcut to Snakemake module for ngs_mapping
        ngs_mapping = self.parent.modules["ngs_mapping"]

        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        base_path = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
            normal_library=self.get_normal_lib_name(wildcards), **wildcards
        )
        return {"bam": ngs_mapping(base_path + ".bam"), "bai": ngs_mapping(base_path + ".bam")}

    def _get_input_files_pileup_tumor(self, wildcards):
        """Get input files for rule ``pileup_tumor``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'pileup_tumor', BAM and BAI files.
        """
        # Get shortcut to Snakemake module for ngs_mapping
        ngs_mapping = self.parent.modules["ngs_mapping"]

        base_path = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}".format(
            **wildcards
        )
        return {"bam": ngs_mapping(base_path + ".bam"), "bai": ngs_mapping(base_path + ".bam")}

    @staticmethod
    def _get_input_files_contamination(wildcards):
        """Get input files for rule ``contamination``.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'tumor_library' (e.g., 'P001-T1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns dictionary with input files for rule 'contamination', Normal and Tumor
        pileup files.
        """
        base_path = (
            "work/{mapper}.mutect2.{tumor_library}/out/{mapper}.mutect2.{tumor_library}".format(
                **wildcards
            )
        )
        return {"normal": base_path + ".normal.pileup", "tumor": base_path + ".tumor.pileup"}

    def get_output_files(self, action):
        """Get output files for Mutect2 rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns dictionary with expected output files based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Initialise variables
        exts = {}
        output_files = {}

        # Validate action
        self._validate_action(action)

        # Set expected extensions based on action
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
                "full_vcf": ".full.vcf.gz",
                "full_vcf_md5": ".full.vcf.gz.md5",
                "full_vcf_tbi": ".full.vcf.gz.tbi",
                "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
                "vcf": ".vcf.gz",
                "vcf_md5": ".vcf.gz.md5",
                "vcf_tbi": ".vcf.gz.tbi",
                "vcf_tbi_md5": ".vcf.gz.tbi.md5",
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

        # Define output dictionary
        for k, v in exts.items():
            output_files[k] = self.base_path_out.format(var_caller=self.name, ext=v)
        return output_files

    def get_log_file(self, action):
        """Get log files for Mutect2 rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns dictionary with expected log files based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Initialise variables
        postfix = ""
        log_files = {}
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )

        # Validate action
        self._validate_action(action)

        # Set expected format based on action
        if action != "run":
            postfix = "." + action
        prefix = (
            "work/{{mapper}}.{var_caller}.{{tumor_library}}/log/"
            "{{mapper}}.{var_caller}.{{tumor_library}}{postfix}"
        ).format(var_caller=self.name, postfix=postfix)

        # Define output dictionary
        for key, ext in key_ext:
            log_files[key] = prefix + ext
            log_files[key + "_md5"] = prefix + ext + ".md5"
        return log_files

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return self.resource_usage_dict.get(action)


class ScalpelStepPart(SomaticVariantCallingStepPart):
    """Somatic variant calling with Scalpel"""

    #: Step name
    name = "scalpel"

    #: Class available actions
    actions = ("run",)

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
        result["tar_md5"] = self.base_path_out.format(var_caller="scalpel", ext=".tar.gz.md5")
        result.update(
            dict(
                zip(
                    somatic_ext_names,
                    expand(self.base_path_out, var_caller=[self.name], ext=somatic_ext_values),
                )
            )
        )
        return result

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=16,  # TODO: Make Scalpel number of thread configurable.
            time="2-00:00:00",  # 2 days
            memory=f"{5 * 1024 * 16}M",
        )


class Strelka2StepPart(SomaticVariantCallingStepPart):
    """Somatic variant calling with strelka2/manta"""

    #: Step name
    name = "strelka2"

    #: Class available actions
    actions = ("run",)

    # Output extension files dictionary. Key: output type (string); Value: extension (string)
    extensions = {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
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

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="1-00:00:00",  # 1 day
            memory="4G",
        )


class JointCallingStepPart(BaseStepPart):
    """Base class for joint calling."""

    #: Name of the step, to be overridden in sub class.
    name = None

    #: Class available actions
    actions = ("run",)

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
        # Validate action
        self._validate_action(action)

        # Get shortcut to Snakemake module for ngs_mapping
        ngs_mapping = self.parent.modules["ngs_mapping"]

        def input_function(wildcards):
            """Helper wrapper function"""
            donor = self.donor_by_name[wildcards.donor_name]

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

        return input_function

    def get_output_files(self, action):
        """Return resulting files generated by this step."""
        # Validate action
        self._validate_action(action)
        return dict(zip(EXT_NAMES, expand(self.base_path_out, name=[self.name], ext=EXT_VALUES)))

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        _ = action
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
        # Validate action
        self._validate_action(action)

        def arg_function(wildcards):
            donor = self.donor_by_name[wildcards.donor_name]
            result = {
                "sample_list": [
                    ngs_library.name
                    for bio_sample in donor.bio_samples.values()
                    for test_sample in bio_sample.test_samples.values()
                    for ngs_library in test_sample.ngs_libraries.values()
                ]
            }
            if "ignore_chroms" in self.parent.config:
                ignore_chroms = self.parent.config["ignore_chroms"]
                result["ignore_chroms"] = ignore_chroms
            return result

        return arg_function


class BcftoolsJointStepPart(JointCallingStepPart):
    """Somatic variant calling with Samtools in "joint" mode.

    Simply pass all samples of one donor through Samtools. The resulting files are then to be
    subfiltered.
    """

    #: Step name
    name = "bcftools_joint"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        mem_mb = 1024 * self.parent.config["bcftools_joint"]["num_threads"]
        return ResourceUsage(
            threads=self.parent.config["bcftools_joint"]["num_threads"],
            time="2-00:00:00",  # 2 days
            memory=f"{mem_mb}M",
        )


class VarscanJointStepPart(JointCallingStepPart):
    """Somatic variant calling with Varscan in "joint" mode."""

    #: Step name
    name = "varscan_joint"

    #: Class available actions
    actions = ("run", "call_pedigree")

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="2-00:00:00",  # 2 days
            memory="1024M",
        )


class PlatypusJointStepPart(JointCallingStepPart):
    """Somatic variant calling with Platypus in "joint" mode.

    Simply pass all samples of one donor through Platypus. The resulting files are then to be
    subfiltered.
    """

    #: Step name
    name = "platypus_joint"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        mem_mb = int(3.75 * 1024 * self.parent.config["platypus_joint"]["num_threads"])
        return ResourceUsage(
            threads=self.parent.config["platypus_joint"]["num_threads"],
            time="2-00:00:00",  # 2 days
            memory=f"{mem_mb}M",
        )


class GatkHcJointStepPart(JointCallingStepPart):
    """Somatic variant calling with GATK HC in "joint" mode.

    Simply pass all samples of one donor through GATK HC. The resulting files are then to be
    subfiltered.
    """

    #: Step name
    name = "gatk_hc_joint"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="3-08:00:00",  # 3 days and 8 hours
            memory=f"{2 * 1024}M",
        )


class GatkUgJointStepPart(JointCallingStepPart):
    """Somatic variant calling with GATK UG in "joint" mode.

    Simply pass all samples of one donor through GATK UG. The resulting files are then to be
    subfiltered.
    """

    #: Step name
    name = "gatk_ug_joint"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="3-08:00:00",  # 3 days and 8 hours
            memory=f"{2 * 1024}M",
        )


class SomaticVariantCallingWorkflow(BaseStep):
    """Perform somatic variant calling"""

    #: Workflow name
    name = "somatic_variant_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            (),
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
        self.register_module("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.{caller}.{tumor_library.name}"
        for caller in set(self.config["tools"]) & set(SOMATIC_VARIANT_CALLERS_MATCHED):
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=caller,
                ext=EXT_MATCHED[caller].values() if caller in EXT_MATCHED else EXT_VALUES,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=caller,
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
