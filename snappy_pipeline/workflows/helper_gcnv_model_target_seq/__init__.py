# -*- coding: utf-8 -*-
"""Implementation of the ``helper_gcnv_model_target_seq`` step

The ``helper_gcnv_model_wgs`` step takes as the input the results of the ``ngs_mapping`` step
(aligned germline reads) and builds a model that can be used by GATK4 gCNV for a particular
library kit.

==========
Step Input
==========

The step uses Snakemake sub workflows for the result of the ``ngs_mapping``
(aligned reads BAM files).

===========
Step Output
===========

All donors will be used to generate the two parts of the required gCNV model, specifically:
``ploidy-model`` and ``cnv_calls-model``. Both are required to execute gCNV in CASE mode.

For example, the relevant directories might look as follows:

::

    work/
    +-- bwa.gcnv_contig_ploidy.<library_kit_name>
        `-- out
            `-- bwa.gcnv_contig_ploidy.<library_kit_name>
                |-- SAMPLE_0
                |   |-- contig_ploidy.tsv
                |   |-- global_read_depth.tsv
                |   |-- mu_psi_s_log__.tsv
                |   |-- sample_name.txt
                |   `-- std_psi_s_log__.tsv
                |-- [...]
                `-- bwa.gcnv_contig_ploidy.<library_kit_name>
                    `-- ploidy-model
                        |-- contig_ploidy_prior.tsv
                        |-- gcnvkernel_version.json
                        |-- interval_list.tsv
                        |-- mu_mean_bias_j_lowerbound__.tsv
                        |-- mu_psi_j_log__.tsv
                        |-- ploidy_config.json
                        |-- std_mean_bias_j_lowerbound__.tsv
                        `-- std_psi_j_log__.tsv
    +-- bwa.gcnv_call_cnvs.<library_kit_name>.***_of_***
        `-- out
            `-- bwa.gcnv_call_cnvs.<library_kit_name>.***_of_***
                |-- cnv_calls-calls
                |   |-- SAMPLE_0
                |       `-- [...]
                |   |-- [...]
                |-- cnv_calls-model
                |  |-- denoising_config.json
                |  |-- gcnvkernel_version.json
                |  |-- interval_list.tsv
                |  |-- log_q_tau_tk.tsv
                |  |-- mu_W_tu.tsv
                |  |-- mu_ard_u_log__.tsv
                |  |-- mu_log_mean_bias_t.tsv
                |  |-- mu_psi_t_log__.tsv
                |  |-- std_W_tu.tsv
                |  |-- std_ard_u_log__.tsv
                |  |-- std_log_mean_bias_t.tsv
                |  `-- std_psi_t_log__.tsv
                `-- cnv_calls-tracking
                    `-- [...]

====================
Global Configuration
====================

- At the moment, no global configuration is used.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_helper_gcnv_model_target_seq.rst

"""
import os
import re

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand, glob_wildcards

from snappy_pipeline.utils import DictQuery, dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.targeted_seq_cnv_calling import GcnvStepPart

#: Default configuration for the helper_gcnv_model_wgs schema
DEFAULT_CONFIG = r"""
# Default configuration helper_gcnv_model_target_seq
step_config:
  helper_gcnv_model_target_seq:
    path_ngs_mapping: ../ngs_mapping  # REQUIRED

    gcnv:
      path_uniquely_mapable_bed: null  # REQUIRED - path to BED file with uniquely mappable regions.
      path_target_interval_list_mapping: []  # REQUIRED - define one or more set of target intervals.
      # The following will match both the stock IDT library kit and the ones
      # with spike-ins seen fromr Yale genomics.  The path above would be
      # mapped to the name "default".
      # - name: IDT_xGen_V1_0
      #   pattern: "xGen Exome Research Panel V1\\.0*"
      #   path: "path/to/targets.bed"
"""


class BuildGcnvModelStepPart(GcnvStepPart):
    """Build model for GATK4 gCNV"""

    #: Step name
    name = "gcnv"

    #: Class available actions
    actions = (
        "preprocess_intervals",
        "annotate_gc",
        "filter_intervals",
        "scatter_intervals",
        "coverage",
        "contig_ploidy",
        "scatter_intervals",
        "call_cnvs_cohort_mode",
        "post_germline_calls",
        "post_germline_calls_cohort_mode",
        "merge_cohort_vcfs",
    )

    def validate_request(self):
        """Override validation method.

        Build model always in COHORT MODE.
        """
        return "cohort_mode"

    @dictify
    def _get_input_files_post_germline_calls_cohort_mode(self, wildcards, checkpoints):
        checkpoint = checkpoints.build_gcnv_model_scatter_intervals
        library_kit = self.ngs_library_to_kit.get(wildcards.library_name)
        scatter_out = checkpoint.get(library_kit=library_kit, **wildcards).output[0]
        shards = list(
            map(
                os.path.basename,
                glob_wildcards(os.path.join(scatter_out, "temp_{shard}/{file}")).shard,
            )
        )
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield "calls", [
            "work/{name_pattern}.{shard}/out/{name_pattern}.{shard}/.done".format(
                name_pattern=name_pattern, shard=shard
            )
            for shard in shards
        ]
        ext = "ploidy"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield ext, "work/{name_pattern}/out/{name_pattern}/.done".format(name_pattern=name_pattern)

    @dictify
    def _build_ngs_library_to_kit(self):
        gcnv_config = DictQuery(self.w_config).get("step_config/helper_gcnv_model_target_seq/gcnv")
        if not gcnv_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.parent.all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"

        # Build mapping
        regexes = {
            item["pattern"]: item["name"]
            for item in gcnv_config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result


class HelperBuildTargetSeqGcnvModelWorkflow(BaseStep):
    """Perform gCNV model building for WES samples by library kit"""

    #: Workflow name
    name = "helper_gcnv_model_target_seq"

    #: Default biomed sheet class
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
        self.register_sub_step_classes((BuildGcnvModelStepPart,))
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Build mapping from NGS DNA library to library kit
        self.ngs_library_to_kit = self.sub_steps["gcnv"].ngs_library_to_kit

    @listify
    def get_result_files(self):
        """Return list of result files for the gCNV build model workflow."""
        ext_values = (".ratio.tsv", ".interval.vcf.gz", ".vcf.gz")
        donors = {donor.name for donor in self.all_donors() if donor.dna_ngs_library}
        name_pattern = "{mapper}.gcnv_post_germline_calls.{index.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("work", name_pattern, "out", name_pattern + "{ext}"),
            donors,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            ext=ext_values,
        )

    def _yield_result_files(self, tpl, donors, **kwargs):
        """Build output paths from path template and extension list.

        Will only yield the result files for pedigrees where the index is in ``donors``.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if pedigree.index.name in donors:
                    yield from expand(tpl, index=[pedigree.index], **kwargs)

    @listify
    def all_donors(self, include_background=True):
        """Get all donors.

        :param include_background: Boolean flag to defined if background should be included or not.
        Default: True, i.e., background will be included.

        :return: Returns list of all donors in sample sheet.
        """
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    def pick_kits_and_donors(self):
        """Return ``(library_kits, donors)`` with the donors with a matching kit and the kits with a
        matching donor.
        """
        kit_counts = {name: 0 for name in self.ngs_library_to_kit.values()}
        for name in self.ngs_library_to_kit.values():
            kit_counts[name] += 1
        donors = [
            donor
            for donor in self.all_donors()
            if donor.dna_ngs_library and donor.dna_ngs_library.name in self.ngs_library_to_kit
        ]
        return list(sorted(set(self.ngs_library_to_kit.values()))), donors, kit_counts

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        self.ensure_w_config(
            ("step_config", "helper_gcnv_model_target_seq", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for gCNV model building.",
        )
