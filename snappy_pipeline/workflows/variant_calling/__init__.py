# -*- coding: utf-8 -*-
"""Implementation of the ``variant_calling`` step

The ``variant_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs germline variant calling.  The result are variant files
with germline variants (bgzip-ed and indexed VCF files).

Usually, the variant calling step is followed by the ``variant_annotation`` step.

==========
Stability
==========

The HaplotypeCaller and UnifiedGenotyper from the Genome Analysis Toolkit (GATK) are considered
stable.

The other supported callers are still in the experimental stage and may not be stable.

==========
Step Input
==========

The variant calling step uses Snakemake sub workflows for using the result of the ``ngs_mapping``
step.

===========
Step Output
===========

For all pedigrees, variant calling will be performed on the primary DNA NGS libraries of all
members, separately for each configured read mapper and variant caller.  The name of the primary
DNA NGS library of the index will be used as an identification token in the output file.  For each
read mapper, variant caller, and pedigree, the following files will be generated:

- ``{mapper}.{var_caller}.{lib_name}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.gatk_hc.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.gatk_hc.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk_hc.P001-N1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.gatk_hc.P001-N1-DNA1-WES1.vcf.gz.md5
    |       `-- bwa.gatk_hc.P001-N1-DNA1-WES1.vcf.gz.tbi.md5
    [...]

Generally, these files will be unfiltered, i.e., contain low-quality variants.

====================
Global Configuration
====================

- If GATK HaplotypeCaller or GATK UnifiedGenotyper are activated then
  ``static_data_config/dbsnp/path`` must be properly configured
- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_calling.rst

==================================
Available Germline Variant Callers
==================================

The following germline variant callers are currently available

- ``"bcftools"``  -- samtools mpileup plus bcftools
- ``"gatk_hc"`` -- GATK HaplotypeCaller
- ``"gatk_ug"`` -- GATK UnifiedGenotyper

=======
Reports
=======

Currently, the following reports are generated (and are linked from the output directory):

- bcftools_stats (txt) is always generated by default.
  Within this file, stats are broken down separately for known and novel events.
  Report contents depend on the version of bcftools used. With version 1.3.1 the report includes
  the following details:

    - SN, Summary numbers
    - TSTV, Transitions/transversions
    - SiS, Singleton stats
    - AF, Stats by non-reference allele frequency
    - QUAL, Stats by quality
    - IDD, InDel distribution
    - ST, Substitution types
    - DP, Depth distribution
    - PSC, Per-sample counts
    - PSI, Per-Sample indels
    - HWE, Hardy-Weinberg equilibrium

- jannovar_statistics (txt) is always generated by default.
  Requires jannovar_statistics/path_ser to be set to a ".ser" file.
  Using jannovar-cli and htslib version 1.3.2 the report includes the following details:

     - putative_impacts (counts by type)
     - variant_effects (counts by type)
     - genome_regions (counts by type)
     - ts_tv_count (count TS, count TV)
     - alt_allele_count (counts for each number of alleles)
     - filter_count (counts by type, if any)
     - is_filtered_count (count passed and failed, if filtering is used)
     - contig_counts (count events per chromosome)

.. _variant_calling_parallel_execution:

==================
Parallel Execution
==================

For many of the variant callers, cluster-parallel execution has been implemented (indicated by
having a ``use_profile`` configuration setting).  Here, a temporary directory with a Snakemake
workflow is written out and then executed.  The default behaviour that the temporary files are
removed in the case of an error.  This behaviour can be changed by setting the ``keep_tmpdir``
setting to ``"onerror"`` or ``"always"``.  Further, for debugging, the number of windows to
create can be limited using ``debug_trunc_tokens`` (the default of ``0``) leads to the processing
of all windows.  Resource requirements in terms of memory or running time can be boosted using
``job_mult_memory`` and ``job_mult_time`` (similarly for the joining step and ``merge_mult_*``).

When the temporary directory is kept, a failed execution can be restarted by calling ``snakemake``
in the temporary directory with the command line written to the file ``snakemake_call.sh``.
"""

from collections import OrderedDict
from itertools import chain
import re
import sys
import warnings

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, flatten, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Available germline variant callers
VARIANT_CALLERS = (
    "bcftools",
    "gatk_hc",
    "gatk_ug",
    "gatk4_hc_gvcf",
    "gatk4_hc_joint",
)

#: Default configuration for the variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_calling
step_config:
  variant_calling:
    # common configuration
    path_ngs_mapping: ../ngs_mapping  # REQUIRED

    # report generation
    baf_file_generation:
      enabled: true
      min_dp: 10  # minimal DP of variant, must be >=1
    bcftools_stats:
      enabled: true
    jannovar_stats:
      enabled: true
      path_ser: REQUIRED  # REQUIRED

    # variant calling tools and their configuration
    tools: ['gatk4_hc_gvcf']  # REQUIRED, examples: 'gatk_hc', 'gatk_ug'
    ignore_chroms:
    - NC_007605  # herpes virus
    - hs37d5     # GRCh37 decoy
    - chrEBV     # Eppstein-Barr Virus
    - '*_decoy'  # decoy contig
    - 'HLA-*'    # HLA genes

    bcftools:
      max_depth: 250
      max_indel_depth: 250
      window_length: 10000000
      num_threads: 16
    gatk4_hc_joint:
      allow_seq_dict_incompatibility: false  # REQUIRED
      annotations: []  # REQUIRED
      annotation_groups:  # REQUIRED
      - StandardAnnotation
      - StandardHCAnnotation
    gatk4_hc_gvcf:
      allow_seq_dict_incompatibility: false  # REQUIRED
      annotations: []  # REQUIRED
      annotation_groups:  # REQUIRED
      - AlleleSpecificAnnotation
      - GenotypeAnnotation
      - InfoFieldAnnotation
      - StandardAnnotation
      - StandardHCAnnotation
      - StandardFlowBasedAnnotation
    gatk_hc:
      # Parallelization configuration
      num_threads: 2            # number of cores to use locally
      window_length: 5000000    # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 0          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10  # throttling of status jobs
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
    gatk_ug:
      # Parallelization configuration
      num_threads: 2            # number of cores to use locally
      window_length: 5000000    # split input into windows of this size, each triggers a job
      num_jobs: 500             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 0          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10  # throttling of status jobs
      debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
      keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1        # memory multiplier
      job_mult_time: 1          # running time multiplier
      merge_mult_memory: 1      # memory multiplier for merging
      merge_mult_time: 1        # running time multiplier for merging
      # GATK UG--specific configuration
      allow_seq_dict_incompatibility: false
      downsample_to_coverage: 250
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
"""


class InconsistentPedigreeWarning(UserWarning):
    """Raised on inconsistencies with pedigree."""


class GetResultFilesMixin:
    """Mixin to provide ``get_result_files()`` function.

    The function will use ``get_output_files()`` for all actions to obtain
    the files in the ``output/directory`` to expect from the given step.
    """

    @listify
    def get_result_files(self):
        """Return concrete result file paths"""

        def strip_tpl(tpl):
            return tpl.replace(r",[^\.]+", "")

        index_dna_ngs_libraries = self._get_index_dna_ngs_libraries()
        for action in self.actions:
            output_files = self.get_output_files(action).values()
            result_paths_tpls = [
                strip_tpl(p) for p in flatten(output_files) if p.startswith("output/")
            ]
            for path_tpl in result_paths_tpls:
                for index_library_name, member_library_names in index_dna_ngs_libraries.items():
                    kwargs = {
                        "mapper": self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                    }
                    if "index_ngs_library" in path_tpl:
                        kwargs["index_ngs_library"] = [index_library_name]
                        kwargs["donor_ngs_library"] = [member_library_names]
                    else:
                        kwargs["library_name"] = [index_library_name]
                    if "{var_caller}" in path_tpl:
                        kwargs["var_caller"] = self.parent.config["tools"]
                    yield from expand(
                        path_tpl,
                        **kwargs,
                    )

    @dictify
    def _get_index_dna_ngs_libraries(self):
        """Return ``dict`` that maps the index DNA library name to a list of all pedigree
        member's DNA library names.
        """
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if self._is_pedigree_good(pedigree):
                    index = pedigree.index.dna_ngs_library.name
                    donors = [
                        donor.dna_ngs_library.name
                        for donor in pedigree.donors
                        if donor.dna_ngs_library
                    ]
                    yield index, donors

    def _is_pedigree_good(self, pedigree) -> bool:
        """Check pedigrees for inconsistencies and issue warning for any.

        :return: ``True`` if there was no inconsistency reported
        """
        msg = None
        donor_names = list(sorted(d.name for d in pedigree.donors))
        if not pedigree.index:  # pragma: no cover
            msg = f"INFO: pedigree without index (name: {donor_names})"
        elif not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = f"INFO: pedigree index without DNA NGS library (names: {donor_names})"
        if msg:
            warnings.warn(InconsistentPedigreeWarning(msg))
        return not msg


class GetLogFileMixin:
    """Mixin to provide the generic ``get_log_file()`` function"""

    @dictify
    def get_log_file(self, action):
        """Return dict of log files in the "log" directory."""
        _ = action
        token = f"{{mapper}}.{self.name}.{{library_name}}"
        prefix = f"work/{token}/log/{token}.{self.name}_{action}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, f"{prefix}{ext}"
            yield f"{key}_md5", f"{prefix}{ext}.md5"


class VariantCallingStepPart(GetResultFilesMixin, GetLogFileMixin, BaseStepPart):
    """Base class for germline variant calling step parts

    Variant calling is performed on a per-pedigree level.  The (one) index individual is used
    for naming the output file.
    """

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{index_library_name}}/out/"
            "{{mapper}}.{var_caller}.{{index_library_name}}{ext}"
        )
        self.base_path_tmp = self.base_path_out.replace("/out/", "/tmp/")
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return required input files"""
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]

        if not pedigree.index or not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = "INFO: pedigree without index (names: {})"
            donor_names = list(sorted(d.name for d in pedigree.donors))
            print(msg.format(donor_names), file=sys.stderr)
            yield "bam", []
            yield "ped", None
        else:
            library_name = pedigree.index.dna_ngs_library.name
            yield "ped", f"work/write_pedigree.{library_name}/out/{library_name}.ped"

            bams = []
            for donor in pedigree.donors:
                if not donor.dna_ngs_library:
                    continue  # skip
                infix = f"{wildcards.mapper}.{donor.dna_ngs_library.name}"
                bams.append(ngs_mapping(f"output/{infix}/out/{infix}.bam"))
            yield "bam", bams

    def get_output_files(self, action):
        """Return step part output files"""
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def _get_output_files_run(self):
        token = f"{{mapper}}.{self.name}.{{library_name}}"
        work_files = {
            "vcf": f"work/{token}/out/{token}.vcf.gz",
            "vcf_md5": f"work/{token}/out/{token}.vcf.gz.md5",
            "vcf_tbi": f"work/{token}/out/{token}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{token}/out/{token}.vcf.gz.tbi.md5",
        }
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]


class BcftoolsStepPart(VariantCallingStepPart):
    """Germline variant calling with bcftools"""

    #: Step name
    name = "bcftools"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            time="2-00:00:00",  # 2 days
            memory=f"{int(3.75 * 1024 * 16)}M",
        )


class GatkCallerStepPartBase(VariantCallingStepPart):
    """Germlin variant calling with GATK caller"""

    def check_config(self):
        if self.__class__.name not in self.config["tools"]:
            return  # caller not enabled, skip  # pragma: no cover
        self.parent.ensure_w_config(
            ("static_data_config", "dbsnp", "path"),
            "dbSNP not configured but required for {}".format(self.__class__.name),
        )

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
            memory=f"{14 * 1024}M",
        )


class GatkHaplotypeCallerStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK HaplotypeCaller

    This triggers the cluster-parallel variant calling with gatk_ug, GATK3.
    """

    #: Step name
    name = "gatk_hc"


class GatkUnifiedGenotyperStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK UnifiedGenotyper

    This triggers the cluster-parallel variant calling with gatk_ug, GATK3.
    """

    #: Step name
    name = "gatk_ug"


class Gatk4CallerStepPartBase(VariantCallingStepPart):
    """Base class for germline variant calling with GATK4

    In contrast to our GATK 3 wrappers, we do not perform parallelization of the variant calling.
    Note that this step will generate both GVCF and VCF files for the pedigree.
    """

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # caller not enabled, skip  # pragma: no cover
        self.parent.ensure_w_config(
            ("static_data_config", "dbsnp", "path"),
            "dbSNP not configured but required for {}".format(self.__class__.name),
        )

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="2-00:00:00",
            memory="4G",
        )


class Gatk4HaplotypeCallerJointStepPart(Gatk4CallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller doing joint calling per pedigree."""

    name = "gatk4_hc_joint"

    actions = ("run",)


class Gatk4HaplotypeCallerGvcfStepPart(Gatk4CallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller and gVCF workflow.

    In contrast to our GATK 3 wrappers, we do not perform parallelization of the variant calling.
    This will generate GVCF files for each individual and a joint GVCF file as a VCF file for the
    whole pedigree.
    """

    name = "gatk4_hc_gvcf"

    actions = ("discover", "combine_gvcfs", "genotype")

    @dictify
    def _get_input_files_discover(self, wildcards):
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        bam_path = f"output/{infix}/out/{infix}.bam"
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        yield "bam", ngs_mapping(bam_path)

    @dictify
    def _get_input_files_combine_gvcfs(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]

        if not pedigree.index or not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = "INFO: pedigree without index (names: {})"
            donor_names = list(sorted(d.name for d in pedigree.donors))
            print(msg.format(donor_names), file=sys.stderr)
            yield "ped", None
            yield "gvcf", []
        else:
            library_name = pedigree.index.dna_ngs_library.name
            yield "ped", f"work/write_pedigree.{library_name}/out/{library_name}.ped"

            gvcfs = []
            for donor in pedigree.donors:
                if not donor.dna_ngs_library:
                    continue  # skip
                infix = f"{wildcards.mapper}.gatk4_hc_gvcf_discover.{donor.dna_ngs_library.name}"
                gvcfs.append(f"work/{infix}/out/{infix}.g.vcf.gz")
            yield "gvcf", gvcfs

    @dictify
    def _get_input_files_genotype(self, wildcards):
        infix = f"{wildcards.mapper}.gatk4_hc_gvcf_combine_gvcfs.{wildcards.library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_md5", f"work/{infix}/out/{infix}.g.vcf.gz.md5"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"
        yield "gvcf_tbi_md5", f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5"

    @dictify
    def _get_output_files_discover(self):
        infix = "{mapper}.gatk4_hc_gvcf_discover.{library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_md5", f"work/{infix}/out/{infix}.g.vcf.gz.md5"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"
        yield "gvcf_tbi_md5", f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5"
        yield "output_links", []

    @dictify
    def _get_output_files_combine_gvcfs(self):
        infix = "{mapper}.gatk4_hc_gvcf_combine_gvcfs.{library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_md5", f"work/{infix}/out/{infix}.g.vcf.gz.md5"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"
        yield "gvcf_tbi_md5", f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5"
        yield "output_links", []

    @dictify
    def _get_output_files_genotype(self):
        infix = "{mapper}.gatk4_hc_gvcf.{library_name}"
        result = {
            "gvcf": f"work/{infix}/out/{infix}.g.vcf.gz",
            "gvcf_md5": f"work/{infix}/out/{infix}.g.vcf.gz.md5",
            "gvcf_tbi": f"work/{infix}/out/{infix}.g.vcf.gz.tbi",
            "gvcf_tbi_md5": f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5",
            "vcf": f"work/{infix}/out/{infix}.vcf.gz",
            "vcf_md5": f"work/{infix}/out/{infix}.vcf.gz.md5",
            "vcf_tbi": f"work/{infix}/out/{infix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{infix}/out/{infix}.vcf.gz.tbi.md5",
        }
        yield from result.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(result.values(), self.get_log_file("genotype").values())
        ]


class BcftoolsStatsStepPart(GetResultFilesMixin, GetLogFileMixin, BaseStepPart):
    """Base class for VCF statistics computation with "bcftools stats"

    Statistics are computed overall and per-sample
    """

    # TODO: maybe we need to use "--stats" anyway and can handle pedigree VCF files then...

    #: Step name
    name = "bcftools_stats"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.{index_ngs_library}/report/bcftools_stats/"
            "{mapper}.{var_caller}.{index_ngs_library}.{donor_ngs_library}"
        )

    @dictify
    def _get_input_files_run(self):
        yield "vcf", (
            "work/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}.vcf.gz"
        )

    @dictify
    def _get_output_files_run(self):
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        work_files = {key: f"{self.base_path_out}{ext}" for key, ext in ext_names.items()}
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), [self.get_log_file("run").values()])
        ]

    def get_log_file(self, action):
        """Return dict of log files in the "log" directory."""
        self._validate_action(action)
        token = f"{{mapper}}.{self.name}.{{index_library_name}}"
        prefix = f"work/{token}/log/{token}.{{donor_ngs_library}}.{self.name}_{action}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, f"{prefix}{ext}"
            yield f"{key}_md5", f"{prefix}{ext}.md5"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="02:00:00",
            memory="1024M",
        )


class JannovarStatisticsStepPart(GetResultFilesMixin, GetLogFileMixin, BaseStepPart):
    """Base class for VCF statistics computation with "jannovar statistics"

    Statistics are computed overall and per-sample
    """

    #: Step name
    name = "jannovar_stats"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.{index_ngs_library}/report/jannovar_stats/"
            "{mapper}.{var_caller}.{index_ngs_library}"
        )

    @dictify
    def get_input_files(self, action):
        """Return path to input files"""
        # Validate action
        self._validate_action(action)
        # Return path to input VCF file
        yield "vcf", (
            "work/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}.vcf.gz"
        )

    @dictify
    def get_output_files(self, action):
        """Return output files that all germline variant calling sub steps must return (VCF +
        TBI file)
        """
        # Validate action
        self._validate_action(action)
        ext_names = {"report": ".txt", "report_md5": ".txt.md5"}
        work_files = {}
        for key, ext in ext_names.items():
            work_files[key] = self.base_path_out + ext
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), [self.get_log_file("run")])
        ]

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="04:00:00",  # 4 hours
            memory=f"{int(3.75 * 1024 * 2)}M",
        )


class BafFileGenerationStepPart(GetResultFilesMixin, GetLogFileMixin, BaseStepPart):
    """Class for computing B allele fraction files.

    One file is generated per sample in the output VCF files.
    """

    #: Step name
    name = "baf_file_generation"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.{index_ngs_library}/report/baf/"
            r"{mapper}.{var_caller}.{index_ngs_library}.{donor_ngs_library,[^\.]+}.baf"
        )

    @dictify
    def get_input_files(self, action):
        """Return path to input files"""
        # Validate action
        self._validate_action(action)
        # Return path to input VCF file
        yield "vcf", (
            "work/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}.vcf.gz"
        )

    @dictify
    def get_output_files(self, action):
        """Return output files that all germline variant calling sub steps must return (VCF +
        TBI file)
        """
        # Validate action
        self._validate_action(action)
        ext_names = {"bw": ".bw", "bw_md5": ".bw.md5"}
        work_files = {}
        for key, ext in ext_names.items():
            work_files[key] = self.base_path_out + ext
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]

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
            time="02:00:00",
            memory="1024M",
        )


class VariantCallingWorkflow(BaseStep):
    """Perform germline variant calling"""

    name = "variant_calling"
    sheet_shortcut_class = GermlineCaseSheet

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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                BcftoolsStepPart,
                GatkHaplotypeCallerStepPart,
                GatkUnifiedGenotyperStepPart,
                Gatk4HaplotypeCallerJointStepPart,
                Gatk4HaplotypeCallerGvcfStepPart,
                BcftoolsStatsStepPart,
                JannovarStatisticsStepPart,
                BafFileGenerationStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        for tool in self.config["tools"]:
            yield from self.sub_steps[tool].get_result_files()
        for name in ("baf_file_generation", "bcftools_stats", "jannovar_stats"):
            if self.w_config["step_config"]["variant_calling"][name]["enabled"]:
                yield from self.sub_steps[name].get_result_files()

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for variant calling",
        )
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant calling",
        )
        # Check that only valid tools are selected
        selected = set(self.w_config["step_config"]["variant_calling"]["tools"])
        invalid = selected - set(VARIANT_CALLERS)
        if invalid:
            raise Exception(
                "Invalid variant callers selected: {}".format(  # pragma: no cover
                    list(sorted(invalid))
                )
            )
