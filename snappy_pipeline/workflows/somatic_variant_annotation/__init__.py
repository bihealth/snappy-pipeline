# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_annotation`` step

The ``somatic_variant_annotation`` step takes as the input the results of the
``somatic_variant_calling`` step (bgzip-ed and indexed VCF files) and performs annotation of the
somatic variants.  The result are annotated versions of the somatic variant VCF files (again
bgzip-ed and indexed VCF files).

==========
Step Input
==========

The somatic variant annotation step uses Snakemake sub workflows for using the result of the
``somatic_variant_calling`` step.

The main assumption is that each VCF file contains the two matched normal and tumor samples.

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``variant_calling`` step.

===========
Step Output
===========

Users can annotate all genes & transcripts overlapping with the variant locus, or
they can select one representative gene and transcript for annotation.
In the latter case, the output vcf file will only contain one annotation per variant, while
in the former case, there might be over 100 annotations for each variant.

The ordering of features driving the representative annotation choice is under user control.
The default order is:

1. ``biotype``: protein coding genes come first, it is unclear what is the order for other types of genes
2. ``mane``: the `MANE transcript <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ is selected before other transcripts
3. ``appris``: the `APPRIS principal isoform <https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii89/6701991>`_ is selected before alternates
4. ``tsl``: `Transcript Support Level <http://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html>`_ values in increasing order
5. ``ccds``: Transcripts with `CCDS <https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi>`_ ids are selected before those without
6. ``canonical``: ENSEMBL canonical transcripts are selected before the others
7. ``rank``: VEP internal ranking is used
8. ``length``: longer transcripts are preferred to shorter ones

This order is (hopefully) suitable for cBioPortal export, as well defined transcripts from protein-coding genes are selected when possible.
However, it is recommended to check the full annotation for variants in or nearby disease-relevant genes.

All annotators generate a vcf with one annotation per transcript, and some annotators
(only ENSEMBL's Variant Effect Predictor in the current implementation) can also produce another
output containing all annotations.
The single annotation vcf is named ``<mapper>.<caller>.<annotator>.vcf.gz`` and
the full annotation output is named ``<mapper>.<caller>.<annotator>.full.vcf.gz``

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_annotation.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import os
import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_JOINT,
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Names of the annotator tools
ANNOTATION_TOOLS = ("jannovar", "vep")

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_annotation
step_config:
  somatic_variant_annotation:
    tools: ["jannovar", "vep"]
    path_somatic_variant_calling: ../somatic_variant_calling   # REQUIRED
    tools_ngs_mapping: []      # default to those configured for ngs_mapping
    tools_somatic_variant_calling: []  # default to those configured for somatic_variant_calling
    jannovar:
      path_jannovar_ser: REQUIRED                # REQUIRED
      flag_off_target: False  # REQUIRED
      dbnsfp:  # configuration for default genome release, needs change if differing
        col_contig: 1
        col_pos: 2
        columns: []
      annotation_tracks_bed: []
      annotation_tracks_tsv: []
      annotation_tracks_vcf: []
      window_length: 50000000   # split input into windows of this size, each triggers a job
      num_jobs: 100             # number of windows to process in parallel
      use_profile: true         # use Snakemake profile for parallel processing
      restart_times: 5          # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 10   # throttling of job creation
      max_status_checks_per_second: 10   # throttling of status checks
      ignore_chroms:            # patterns of chromosome names to ignore
      - NC_007605  # herpes virus
      - hs37d5     # GRCh37 decoy
      - chrEBV     # Eppstein-Barr Virus
      - 'GL*'      # problematic unplaced loci
      - '*_decoy'  # decoy contig
      - 'HLA-*'    # HLA genes
    vep:
      cache_dir: ""             # Defaults to $HOME/.vep Not a good idea on the cluster
      species: homo_sapiens
      assembly: GRCh38
      cache_version: 102        # WARNING- this must match the wrapper's vep version!
      tx_flag: "gencode_basic"  # The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged".
      pick_order: ["biotype", "mane", "appris", "tsl", "ccds", "canonical", "rank", "length"]
      num_threads: 8
      buffer_size: 1000
      output_options:
      - everything
"""


class AnnotateSomaticVcfStepPart(BaseStepPart):
    """Annotate VCF file from somatic calls

    .. note:

        The ``tumor_library`` wildcard can actually be the name of a donor!
    """

    #: Only creates vcf with one annotation per variant
    has_full = False

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancre sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor

    @dictify
    def get_input_files(self, action):
        """Return path to somatic vcf input file"""
        # Validate action
        self._validate_action(action)
        tpl = (
            "output/{mapper}.{var_caller}.{tumor_library}/out/"
            "{mapper}.{var_caller}.{tumor_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["somatic_variant_calling"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{{mapper}}.{{var_caller}}.{annotator}.{{tumor_library}}/out/"
            "{{mapper}}.{{var_caller}}.{annotator}.{{tumor_library}}"
        ).format(annotator=self.annotator)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        if self.has_full:
            key_ext["full"] = ".full.vcf.gz"
            key_ext["full_tbi"] = ".full.vcf.gz.tbi"
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        """Return mapping of log files."""
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{{mapper}}.{{var_caller}}.{annotator}.{{tumor_library}}/log/"
            "{{mapper}}.{{var_caller}}.{annotator}.{{tumor_library}}"
        ).format(annotator=self.annotator)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_params(self, action):
        """Return arguments to pass down."""
        _ = action

        def params_function(wildcards):
            if wildcards.tumor_library not in self.donors:
                return {
                    "tumor_library": wildcards.tumor_library,
                    "normal_library": self.get_normal_lib_name(wildcards),
                }
            else:
                return {}

        return params_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name


class JannovarAnnotateSomaticVcfStepPart(AnnotateSomaticVcfStepPart):
    """Annotate VCF file from somatic calling using "Jannovar annotate-vcf" """

    #: Step name
    name = "jannovar"

    #: Annotator name to construct output paths
    annotator = "jannovar"

    #: Class available actions
    actions = ("annotate_somatic_vcf",)

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
            time="4-04:00:00",  # 4 days and 4 hours
            memory=f"{8 * 1024 * 2}M",
        )

    def check_config(self):
        if self.name not in self.config["tools"]:
            return
        self.parent.ensure_w_config(
            ("step_config", "somatic_variant_annotation", "jannovar", "path_jannovar_ser"),
            "Path to serialized Jannovar database",
        )


class VepAnnotateSomaticVcfStepPart(AnnotateSomaticVcfStepPart):
    """Annotate VCF file from somatic calling using ENSEMBL's VEP"""

    #: Step name
    name = "vep"

    #: Annotator name to construct output paths
    annotator = "vep"

    #: Class available actions
    actions = ("run",)

    #: Also creates vcf with all annotations
    has_full = True

    #: Allowed keywords for pick order
    PICK_ORDER = ("biotype", "mane", "appris", "tsl", "ccds", "canonical", "rank", "length")

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["vep"]["num_threads"],
            time="24:00:00",  # 24 hours
            memory=f"{16 * 1024 * 1}M",
        )

    def check_config(self):
        if self.name not in self.config["tools"]:
            return
        if self.config["vep"]["tx_flag"] not in ("merged", "refseq", "gencode_basic"):
            raise InvalidConfiguration("tx_flag must be 'gencode_basic', or 'merged' or 'refseq'")
        if not all([x in self.PICK_ORDER for x in self.config["vep"]["pick_order"]]):
            raise InvalidConfiguration(
                "pick order keywords must be in {}".format(", ".join(self.PICK_ORDER))
            )


class SomaticVariantAnnotationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    name = "somatic_variant_annotation"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            (SomaticVariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (JannovarAnnotateSomaticVcfStepPart, VepAnnotateSomaticVcfStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow(
            "somatic_variant_calling", self.config["path_somatic_variant_calling"]
        )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_somatic_variant_calling"]:
            self.config["tools_somatic_variant_calling"] = self.w_config["step_config"][
                "somatic_variant_calling"
            ]["tools"]

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        annotators = set(self.config["tools"]) & set(ANNOTATION_TOOLS)
        callers = set(self.config["tools_somatic_variant_calling"])
        name_pattern = "{mapper}.{caller}.{annotator}.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            annotator=annotators,
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            annotator=annotators,
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )
        # Annotators with full output
        full = list(
            filter(
                lambda x: self.sub_steps[x].has_full,
                set(self.config["tools"]) & set(ANNOTATION_TOOLS),
            ),
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + ".full{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            annotator=full,
            ext=EXT_VALUES,
        )
        # joint calling
        name_pattern = "{mapper}.{caller}.{annotator}.{donor.name}"
        yield from self._yield_result_files_joint(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_JOINT),
            annotator=annotators,
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
            ("step_config", "somatic_variant_annotation", "path_somatic_variant_calling"),
            "Path to variant calling not configured but required for somatic variant annotation",
        )
