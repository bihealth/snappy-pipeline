# -*- coding: utf-8 -*-
"""Implementation of purity and ploidy checking for somatic NGS samples

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_purity_ploidy_estimate.rst

"""

from collections import OrderedDict
import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Tools for estimating purity and ploidy.
PURITY_PLOIDY_TOOLS = "ascat"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  somatic_purity_ploidy_estimate:
    tools: ['ascat']
    tool_cnv_calling: cnvetti
    # Configuration with read mapper and path to mapping output. Will use this
    # for generating a pileup using samtools for obtaining the b allele
    # fraction and computing coverage.
    tool_ngs_mapping: bwa
    path_ngs_mapping: ../ngs_mapping
    # Configuration of ASCAT method.
    ascat:
      # BED file with loci for B allele frequency.
      b_af_loci: REQUIRED  # REQUIRED
""".lstrip()


class AscatStepPart(BaseStepPart):
    """Estimation of purity and ploidy using ASCAT.

    Notes
    -----

    - Although we name the virtual probes "SNP${num}", they are not guaranteed
      to
    """

    # Name of the step.
    name = "ascat"
    # The actions for generating BAF/CNV files for the tumor and/nor normal
    # sample, and finally to run the ASCAT pipeline.
    actions = (
        "baf_tumor",
        "baf_normal",
        "cnv_tumor",
        "cnv_normal",
        "cnv_tumor_wes",
        "cnv_normal_wes",
        "run_ascat",
    )

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library_name]
        return pair.normal_sample.dna_ngs_library.name

    def get_input_files(self, action):
        """Return input files"""
        assert action in self.actions
        return getattr(self, "_get_input_files_{}".format(action))()

    def _get_input_files_baf_tumor(self):
        """Return input files for generating BAF file for the tumor."""

        def func(wildcards):
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            base_path = (
                "output/{mapper}.{tumor_library_name}/out/{mapper}.{tumor_library_name}"
            ).format(**wildcards)
            return {
                "bam": ngs_mapping(base_path + ".bam"),
                "bai": ngs_mapping(base_path + ".bam.bai"),
            }

        return func

    def _get_input_files_baf_normal(self):
        """Return input files for generating BAF file for the normal."""

        def func(wildcards):
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            base_path = (
                "output/{mapper}.{normal_library_name}/out/{mapper}.{normal_library_name}"
            ).format(**wildcards)
            return {
                "bam": ngs_mapping(base_path + ".bam"),
                "bai": ngs_mapping(base_path + ".bam.bai"),
            }

        return func

    def _get_input_files_cnv_tumor(self):
        """Return input files for generating BAF file for the tumor."""
        return self._get_input_files_baf_tumor()

    def _get_input_files_cnv_normal(self):
        """Return input files for generating CNV file for the normal."""
        return self._get_input_files_baf_normal()

    def _get_input_files_cnv_tumor_wes(self):
        """Return input files for generating CNV file from copywriter for tumor."""

        def func(wildcards):
            wgs_cnv_calling = self.parent.sub_workflows["somatic_targeted_seq_cnv_calling"]
            base_path = (
                "work/{mapper}.copywriter.{tumor_library_name}/"
                "out/{mapper}.copywriter.{tumor_library_name}"
            ).format(**wildcards)
            return {"bins": wgs_cnv_calling(base_path + "_bins.txt")}

        return func

    def _get_input_files_cnv_normal_wes(self):
        """Return input files for generating CNV file from copywriter for normal."""

        def func(wildcards):
            wgs_cnv_calling = self.parent.sub_workflows["somatic_targeted_seq_cnv_calling"]
            # look up tumor to normal
            for k, v in self.tumor_ngs_library_to_sample_pair.items():
                if v.normal_sample.dna_ngs_library.name == wildcards["normal_library_name"]:
                    tumor_library = k
                    # break
            base_path = (
                "work/{mapper}.copywriter.{tumor_library_name}/"
                "out/{mapper}.copywriter.{tumor_library_name}"
            ).format(tumor_library_name=tumor_library, **wildcards)
            return {"bins": wgs_cnv_calling(base_path + "_bins.txt")}

        return func

    def _get_input_files_run_ascat(self):
        """Return input files for actually running ASCAT."""

        @dictify
        def func(wildcards):
            result = {
                "baf_tumor": (
                    "work/{mapper}.ascat_baf_tumor.{tumor_library_name}/out/"
                    "{mapper}.ascat_baf_tumor.{tumor_library_name}.txt"
                ),
                "baf_normal": (
                    "work/{mapper}.ascat_baf_normal.{normal_library_name}/out/"
                    "{mapper}.ascat_baf_normal.{normal_library_name}.txt"
                ),
                "cnv_tumor": (
                    "work/{mapper}.ascat_cnv_tumor.{tumor_library_name}/out/"
                    "{mapper}.ascat_cnv_tumor.{tumor_library_name}.txt"
                ),
                "cnv_normal": (
                    "work/{mapper}.ascat_cnv_normal.{normal_library_name}/out/"
                    "{mapper}.ascat_cnv_normal.{normal_library_name}.txt"
                ),
            }
            normal_library_name = self.get_normal_lib_name(wildcards)
            for key, value in result.items():
                yield key, value.format(normal_library_name=normal_library_name, **wildcards)

        return func

    def get_output_files(self, action):
        """Return output files"""
        assert action in self.actions
        return getattr(self, "_get_output_files_{}".format(action))()

    def _get_output_files_baf_tumor(self):
        """Return output files for generating BAF file for the tumor."""
        return {
            "txt": (
                "work/{mapper}.ascat_baf_tumor.{tumor_library_name}/out/{mapper}."
                "ascat_baf_tumor.{tumor_library_name}.txt"
            )
        }

    def _get_output_files_baf_normal(self):
        """Return output files for generating BAF file for the normal."""
        return {
            "txt": (
                "work/{mapper}.ascat_baf_normal.{normal_library_name}/out/{mapper}."
                "ascat_baf_normal.{normal_library_name}.txt"
            )
        }

    def _get_output_files_cnv_tumor(self):
        """Return output files for generating BAF file for the tumor."""
        return {
            "txt": (
                "work/{mapper}.ascat_cnv_tumor.{tumor_library_name}/out/{mapper}."
                "ascat_cnv_tumor.{tumor_library_name}.txt"
            )
        }

    def _get_output_files_cnv_normal(self):
        """Return output files for generating CNV file for the normal."""
        return {
            "txt": (
                "work/{mapper}.ascat_cnv_normal.{normal_library_name}/out/{mapper}."
                "ascat_cnv_normal.{normal_library_name}.txt"
            )
        }

    @dictify
    def _get_output_files_run_ascat(self):
        """Return output files for actually running ASCAT."""
        yield "done", touch("work/{mapper}.ascat.{tumor_library_name}/out/.done")
        INFIXES = ("goodness_of_fit", "ploidy", "segments", "segments_raw")
        for infix in INFIXES:
            path = (
                "work/{mapper}.ascat.{tumor_library_name}/out/{tumor_library_name}_%s.txt"
            ) % infix
            yield infix, path

    def get_log_file(self, action):
        """Return path to log file"""
        assert action in self.actions
        log_dict = {
            "baf_tumor": (
                "work/{mapper}.ascat_baf_tumor.{tumor_library_name}/log/"
                "{mapper}.ascat_baf_tumor.{tumor_library_name}.log"
            ),
            "baf_normal": (
                "work/{mapper}.ascat_baf_normal.{normal_library_name}/log/"
                "{mapper}.ascat_baf_normal.{normal_library_name}.log"
            ),
            "cnv_tumor": (
                "work/{mapper}.ascat_cnv_tumor.{tumor_library_name}/log/"
                "{mapper}.ascat_cnv_tumor.{tumor_library_name}.log"
            ),
            "cnv_normal": (
                "work/{mapper}.ascat_cnv_normal.{normal_library_name}/log/"
                "{mapper}.ascat_cnv_normal.{normal_library_name}.log"
            ),
            "run_ascat": (
                "work/{mapper}.ascat.{tumor_library_name}/log/"
                "{mapper}.ascat.{tumor_library_name}.log"
            ),
        }
        return {"log": log_dict[action]}

    def update_cluster_config(self, cluster_config):
        for action in self.actions:
            cluster_config["somatic_purity_ploidy_estimate_ascat_%s" % action] = {
                "mem": 10 * 1024 * 8,
                "time": "48:00",
                "ntasks": 8,
            }


class SomaticPurityPloidyEstimateWorkflow(BaseStep):
    """Perform purity and ploidy estimation"""

    name = "somatic_purity_ploidy_estimate"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
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
        self.register_sub_step_classes((AscatStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        if self.config["tool_cnv_calling"] == "copywriter":
            self.register_sub_workflow(
                "somatic_targeted_seq_cnv_calling",
                self.config["path_somatic_targeted_seq_cnv_calling"],
            )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        name_pattern = "{mapper}.{tool}.{ngs_library.name}"
        for tool in self.config["tools"]:
            for sheet in self.shortcut_sheets:
                for donor in sheet.donors:
                    # Skip all donors that do not have a non-tumor bio sample, estimation only
                    # implemented for matched samples at the moment.
                    has_normal = any(not s.is_tumor for s in donor.bio_samples.values())
                    if not has_normal:
                        continue
                    for bio_sample in donor.bio_samples.values():
                        if not bio_sample.is_tumor:
                            continue
                        for test_sample in bio_sample.test_samples.values():
                            ngs_library = bio_sample.dna_ngs_library
                            name_pattern_value = name_pattern.format(
                                mapper=self.config["tool_ngs_mapping"],
                                tool=tool,
                                ngs_library=ngs_library,
                            )
                            yield os.path.join("output", name_pattern_value, "out", ".done")
