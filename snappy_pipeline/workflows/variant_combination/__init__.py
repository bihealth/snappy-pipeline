# -*- coding: utf-8 -*-
"""Implementation of the germline ``variant_combination`` step

This step allows to combine results from ``variant_filtration`` and ``wgs_sv_filtration`` in
different ways to obtain candidates for heterozygous compound causative variants.

Currently, this only properly works for the case of trios.  For any other, empty variant
lists will be generated.

==========
Step Input
==========

The variant calling step uses Snakemake sub workflows for using the result of the
``variant_filtration`` and ``wgs_sv_filtration`` steps.

===========
Step Output
===========

.. note:: TODO

====================
Global Configuration
====================

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_combination.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import os
import os.path
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnknownFiltrationSourceException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow
from snappy_pipeline.workflows.wgs_cnv_filtration import WgsCnvFiltrationWorkflow
from snappy_pipeline.workflows.wgs_sv_filtration import WgsSvFiltrationWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Default configuration for the variant_combination step
DEFAULT_CONFIG = r"""
step_config:
  variant_combination:
    path_variant_filtration: ../variant_filtration
    path_wgs_sv_filtration: ../wgs_sv_filtration
    path_wgs_cnv_filtration: ../wgs_cnv_filtration
    tools_ngs_mapping:
    - bwa
    tools_variant_calling:
    - gatk_hc
    tools_wgs_sv_calling:
    - delly2
    tools_wgs_cnv_calling:
    - erds_sv2
    combinations:
#    - name: all
#      operation: vars_intersect
#      left: wgs_sv_filtration:no_filter.all.whole_genome
#      right: variant_filtration:no_filter.de_novo.freq_all.whole_genome.score_all.passthrough
    # SV and a coding small variant share a common TAD.
    - name: small_coding_sv_common_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed
      left: wgs_sv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.coding.\
              passthrough"
    # SV and a conserved small variant share a common TAD.
    - name: small_conserved_sv_common_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed
      left: wgs_sv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.conserved.\
              passthrough"
    # SV and a coding small variant share a common *limb* TAD.
    - name: small_coding_sv_common_limb_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed
      left: wgs_sv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.coding.\
              passthrough"
    # SV and a conserved small variant share a common *limb* TAD.
    - name: small_conserved_sv_common_limb_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed
      left: wgs_sv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.conserved.\
              passthrough"
    # CNV and a coding small variant share a common TAD.
    - name: small_coding_cnv_common_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed
      left: wgs_cnv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.coding.\
              passthrough"
    # CNV and a conserved small variant share a common TAD.
    - name: small_conserved_cnv_common_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed
      left: wgs_cnv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.conserved.\
              passthrough"
    # CNV and a coding small variant share a common *limb* TAD.
    - name: small_coding_cnv_common_limb_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed
      left: wgs_cnv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.coding.\
              passthrough"
    # CNV and a conserved small variant share a common *limb* TAD.
    - name: small_conserved_cnv_common_limb_tad
      operation: vars_share_interval
      args:
        intervals_bed: /fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed
      left: wgs_cnv_filtration:conservative.dominant.whole_genome
      right: "variant_filtration:conservative.dominant.recessive_freq.whole_genome.conserved.\
              passthrough"
"""

#: Ways of combining variants.
COMBINATORS = ("vars_intersect", "vars_share_interval")

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")


class CombineVariantsStepPartBase(BaseStepPart):
    """Base class for the different combinations."""

    #: Name of the step (e.g., for rule names)
    name = None
    #: Token to use in file name
    name_pattern = None

    def __init__(self, parent):
        super().__init__(parent)
        assert self.name_pattern is not None, "Set into class..."
        name_pattern = self.name_pattern
        # Note that the call to ``replace`` here removes the wildcard constraint that is needed
        # in the patern.  The constraint *could* also be moved with the ``wildcard_constraints``
        # keyword recently introduced at some later point at the cost of longer Snakefiles
        # and another set of additional functions.
        self.base_path_out = os.path.join(
            "work", name_pattern, "out", name_pattern.replace(r",[^\.]+", "") + "{ext}"
        )
        self.path_log = os.path.join(
            "work", name_pattern, "out", name_pattern.replace(r",[^\.]+", "") + ".log"
        )

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_combination_{}_run".format(self.name)] = {
            "mem": int(3.7 * 1024 * 2),
            "time": "01:00",
            "ntasks": 2,
        }

    def get_input_files(self, action):
        assert action == "run"

        @dictify
        def input_function(wildcards):
            yield "ped", "work/write_pedigree.{index_library}/out/{index_library}.ped".format(
                **wildcards
            )
            info = self.parent.combinations[self.name][wildcards.combination]
            for side in ("left", "right"):
                source, name_pattern = info[side].split(":", 1)
                for name, ext in zip(EXT_NAMES, EXT_VALUES):
                    yield "{}_{}".format(side, name), self._get_path(
                        wildcards.mapper,
                        wildcards["{}_caller".format(side)],
                        wildcards.index_library,
                        source,
                        name_pattern,
                    ) + ext

        return input_function

    def get_args(self, action):
        """Get "args" from step configuration for the given combination."""
        assert action == "run"

        def args_function(wildcards):
            info = self.parent.combinations[self.name][wildcards.combination]
            return info.get("args", {})

        return args_function

    def _get_path(self, mapper, caller, index_library, source, name_pattern):
        """Generate path to input file"""
        # Initialise variables
        valid_filtration_sources = ["variant_filtration", "wgs_sv_filtration", "wgs_cnv_filtration"]

        # Validate filtration source
        if source not in valid_filtration_sources:
            valid_filtration_sources_str = ", ".join(valid_filtration_sources)
            error_message = (
                "User requested unknown source '{source}'. Valid sources: {valid}.".format(
                    source=source, valid=valid_filtration_sources_str
                )
            )
            raise UnknownFiltrationSourceException(error_message)

        if source == "variant_filtration":
            return self._get_path_variant_filtration(mapper, caller, index_library, name_pattern)
        elif source == "wgs_sv_filtration":
            return self._get_path_wgs_sv_filtration(mapper, caller, index_library, name_pattern)
        elif source == "wgs_cnv_filtration":
            return self._get_path_wgs_cnv_filtration(mapper, caller, index_library, name_pattern)

        # Keep consistent return statements
        return None

    def _get_path_variant_filtration(self, mapper, caller, index_library, name_pattern):
        workflow = self.parent.sub_workflows["variant_filtration"]
        tpl = "{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library}.{name_pattern}"
        chunk = tpl.format(
            mapper=mapper, caller=caller, index_library=index_library, name_pattern=name_pattern
        )
        return workflow(os.path.join("output", chunk, "out", chunk))

    def _get_path_wgs_sv_filtration(self, mapper, caller, index_library, name_pattern):
        workflow = self.parent.sub_workflows["wgs_sv_filtration"]
        tpl = "{mapper}.{caller}.annotated.filtered.{index_library}.{name_pattern}"
        chunk = tpl.format(
            mapper=mapper, caller=caller, index_library=index_library, name_pattern=name_pattern
        )
        return workflow(os.path.join("output", chunk, "out", chunk))

    def _get_path_wgs_cnv_filtration(self, mapper, caller, index_library, name_pattern):
        workflow = self.parent.sub_workflows["wgs_cnv_filtration"]
        tpl = "{mapper}.{caller}.annotated.filtered.{index_library}.{name_pattern}"
        chunk = tpl.format(
            mapper=mapper, caller=caller, index_library=index_library, name_pattern=name_pattern
        )
        return workflow(os.path.join("output", chunk, "out", chunk))

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{ext}", ext)

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log


class VarsIntersectStepPart(CombineVariantsStepPartBase):
    name = "vars_intersect"
    name_pattern = (
        "{mapper}.combined_variants.vars_intersect.{index_library}.{combination}."
        "{left_caller}.{right_caller}"
    )


class VarsShareIntervalStepPart(CombineVariantsStepPartBase):
    name = "vars_share_interval"
    name_pattern = (
        "{mapper}.combined_variants.vars_share_interval.{index_library}."
        "{combination}.{left_caller}.{right_caller}"
    )


class VariantCombinationWorkflow(BaseStep):
    """Combination of germline structural and/or sequence variants."""

    name = "variant_combination"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
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
            (VariantFiltrationWorkflow, WgsSvFiltrationWorkflow, WgsCnvFiltrationWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                VarsIntersectStepPart,
                VarsShareIntervalStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("variant_filtration", self.config["path_variant_filtration"])
        self.register_sub_workflow("wgs_sv_filtration", self.config["path_wgs_sv_filtration"])
        self.register_sub_workflow("wgs_cnv_filtration", self.config["path_wgs_cnv_filtration"])
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_variant_calling"]:
            self.config["tools_variant_calling"] = self.w_config["step_config"]["variant_calling"][
                "tools"
            ]
        if not self.config["tools_wgs_sv_calling"]:
            self.config["tools_wgs_sv_calling"] = self.w_config["step_config"]["wgs_sv_calling"][
                "tools"
            ]["dna"]
        self.combinations = {key: {} for key in COMBINATORS}
        for entry in self.config["combinations"]:
            self.combinations[entry["operation"]][entry["name"]] = entry

    @listify
    def get_result_files(self):
        """Return list of result files for the variant filtration workflow."""
        # Generate output paths without extracting individuals.
        name_pattern_tpl = (
            "{mapper}.combined_variants.{combinator}.{index_library.name}."
            "{combination}.{left_caller}.{right_caller}"
        )
        for combinator in COMBINATORS:
            name_pattern = name_pattern_tpl.replace("{combinator}", combinator)
            for combination in self.combinations[combinator].values():
                left = combination["left"].split(":", 1)
                assert left[0] in ("wgs_sv_filtration", "variant_filtration", "wgs_cnv_filtration")
                if left[0] == "wgs_sv_filtration":
                    left_callers = self.config["tools_wgs_sv_calling"]
                elif left[0] == "wgs_cnv_filtration":
                    left_callers = self.config["tools_wgs_cnv_calling"]
                else:
                    left_callers = self.config["tools_variant_calling"]

                right = combination["right"].split(":", 1)
                assert right[0] in ("wgs_sv_filtration", "variant_filtration", "wgs_cnv_filtration")
                if right[0] == "wgs_sv_filtration":
                    right_callers = self.config["tools_wgs_sv_calling"]
                elif right[0] == "wgs_cnv_filtration":
                    right_callers = self.config["tools_wgs_cnv_calling"]
                else:
                    right_callers = self.config["tools_variant_calling"]

                yield from self._yield_result_files(
                    os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                    mapper=self.config["tools_ngs_mapping"],
                    combination=[combination["name"]],
                    left_caller=left_callers,
                    right_caller=right_callers,
                    ext=EXT_VALUES,
                )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(msg.format(pedigree), file=sys.stderr)
                    continue
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_combination", "path_variant_filtration"),
            ("Path to variant_annotation not configured but required for variant_combination"),
        )
        self.ensure_w_config(
            ("step_config", "variant_combination", "path_wgs_sv_filtration"),
            ("Path to wgs_sv_filtration not configured but required for variant_combination"),
        )
        self.ensure_w_config(
            ("step_config", "variant_combination", "path_wgs_cnv_filtration"),
            ("Path to wgs_cnv_filtration not configured but required for variant_combination"),
        )
