# -*- coding: utf-8 -*-
"""Implementation of the ``targeted_seq_cnv_annotation`` step

The ``targeted_seq_cnv_annotation`` step takes as the input the results of the
``targeted_seq_cnv_calling`` step (called germline CNVs) and
``variant_calling`` (called small germline variants) and performs annotation
and filtration of the structural variants.

Such filters include:

- quality filter for removing calls with low support,
- inheritance compatibility filter that checks for compatibility of inheritance (only works
  for trios currently),
- various filters related to lower the false discovery rate in rare/de novo variant calling,
  e.g., counting number of affected individuals in cohort outside the index' family.

.. note::

    Status: not implemented yet

==========
Step Input
==========

The variant annotation step uses the output of the following CUBI pipeline steps:

- ``targeted_seq_cnv_calling``
- ``variant_annotation``

===========
Step Output
===========

.. note: TODO

====================
Global Configuration
====================

.. note: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_targeted_seq_cnv_annotation.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import os
import re
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import DictQuery, dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.targeted_seq_cnv_calling import TargetedSeqCnvCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Default configuration of the targeted_seq_cnv_filtration step
DEFAULT_CONFIG = r"""
# Default configuration targeted_seq_cnv_annotation
step_config:
  targeted_seq_cnv_annotation:
    path_ngs_mapping: ../ngs_mapping
    path_targeted_seq_cnv_calling: ../targeted_seq_cnv_calling
    tools_ngs_mapping:
    - bwa
    tools_targeted_seq_cnv_calling:
    - xhmm
    bed_files: []
"""


class VcfCnvFilterStepPart(BaseStepPart):
    """Annotate VCF using targeted_seq_cnv_filter.py script."""

    name = "vcf_cnv_filter"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/.done"
        self.log_path = (
            "work/{mapper}.{caller}.annotated.{index_ngs_library}/"
            "log/snakemake.targeted_seq_cnv_filter.log"
        )
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        xhmm_config = DictQuery(self.w_config).get("step_config/targeted_seq_cnv_calling/xhmm")
        if not xhmm_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.parent._all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"

        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in xhmm_config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.parent._all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    def get_input_files(self, action):
        """Return input function returning input file dict."""

        @dictify
        def input_function(wildcards):
            tpl = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
            yield "ped", tpl.format(**wildcards)
            # yield "ped", os.path.realpath(
            #    "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
            # ).format(**wildcards)
            tpl = (
                "output/{mapper}.{caller}.{index_ngs_library}/out/"
                "{mapper}.{caller}.{index_ngs_library}"
            )
            KEY_EXT = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
            # SVs
            targeted_seq_cnv_calling = self.parent.sub_workflows["targeted_seq_cnv_calling"]
            if wildcards.caller == "xhmm":
                library_kit = self.ngs_library_to_kit[wildcards.index_ngs_library]
                tpl = tpl.replace(".{index_ngs_library}", "_genotype.%s" % library_kit)
            elif wildcards.caller == "gcnv":
                library_kit = self.ngs_library_to_kit[wildcards.index_ngs_library]
                tpl = tpl.replace(".{index_ngs_library}", "_merge_cohort_vcfs.%s" % library_kit)
            elif wildcards.caller == "cnvetti_hom":
                library_kit = self.ngs_library_to_kit[wildcards.index_ngs_library]
                tpl = tpl.replace(".{index_ngs_library}", "_merge_cohort_vcfs.%s" % library_kit)
            for key, ext in KEY_EXT.items():
                yield key, targeted_seq_cnv_calling(tpl + ext).format(**wildcards)
            return

        assert action == "run"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        assert action == "run"
        prefix = (
            "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/"
            "{mapper}.{caller}.annotated.{index_ngs_library}"
        )
        KEY_EXT = {
            "vcf": ".vcf.gz",
            "tbi": ".vcf.gz.tbi",
            "vcf_md5": ".vcf.gz.md5",
            "tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in KEY_EXT.items():
            yield key, prefix + ext

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        prefix = (
            "work/{mapper}.{caller}.annotated.{index_ngs_library}/log/"
            "{mapper}.{caller}.annotated.{index_ngs_library}"
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["targeted_seq_cnv_annotation_targeted_seq_cnv_filter"] = {
            "mem": 5 * 1024 * 2,
            "time": "100:00",
            "tasks": 2,
        }


class TargetedSeqCnvAnnotationWorkflow(BaseStep):
    """Perform germline targeted sequencing CNV annotation"""

    name = "targeted_seq_cnv_annotation"
    sheet_shortcut_class = GermlineCaseSheet

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
            (TargetedSeqCnvCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VcfCnvFilterStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow(
            "targeted_seq_cnv_calling", self.config["path_targeted_seq_cnv_calling"]
        )
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Copy over "tools" setting from targeted_seq_cnv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_targeted_seq_cnv_calling"]:
            self.config["tools_targeted_seq_cnv_calling"] = self.w_config["step_config"][
                "variant_calling"
            ]["tools"]
        # Build mapping from NGS DNA library to library kit.
        self.ngs_library_to_kit = self.sub_steps["vcf_cnv_filter"].ngs_library_to_kit

    @listify
    def _all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            filter(is_not_background, sheets)
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    def _pick_kits_and_donors(self):
        """Return ``(library_kits, donors)`` with the donors with a matching kit and the kits with a
        matching donor.
        """
        kit_counts = {name: 0 for name in self.ngs_library_to_kit.values()}
        for name in self.ngs_library_to_kit.values():
            kit_counts[name] += 1
        donors = [
            donor
            for donor in self._all_donors()
            if donor.dna_ngs_library and donor.dna_ngs_library.name in self.ngs_library_to_kit
        ]
        return list(sorted(set(self.ngs_library_to_kit.values()))), donors, kit_counts

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        # Get list of library kits and donors to use.
        library_kits, donors, kit_counts = self._pick_kits_and_donors()
        # Actually yield the result files.
        token = "{mapper}.{caller}.annotated.{index_library.name}"
        if "xhmm" in self.config["tools_targeted_seq_cnv_calling"]:
            min_kit_usages = 10
            chosen_kits = [kit for kit in library_kits if kit_counts.get(kit, 0) > min_kit_usages]
            chosen_donors = [
                donor.name
                for donor in donors
                if self.ngs_library_to_kit.get(donor.dna_ngs_library.name) in chosen_kits
            ]
            yield from self._yield_result_files(
                os.path.join("output", token, "log", token + "{ext}"),
                chosen_donors,
                mapper=self.config["tools_ngs_mapping"],
                caller=["xhmm"],
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )
        if "gcnv" in self.config["tools_targeted_seq_cnv_calling"]:
            min_kit_usages = 10
            chosen_kits = [kit for kit in library_kits if kit_counts.get(kit, 0) > min_kit_usages]
            chosen_donors = [
                donor.name
                for donor in donors
                if self.ngs_library_to_kit.get(donor.dna_ngs_library.name) in chosen_kits
            ]
            yield from self._yield_result_files(
                os.path.join("output", token, "log", token + "{ext}"),
                chosen_donors,
                mapper=self.config["tools_ngs_mapping"],
                caller=["gcnv"],
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )

    def _yield_result_files(self, tpl, donors, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if pedigree.index.name not in donors:
                    msg = "INFO: pedigree index {}, not included in {}"
                    print(msg.format(pedigree.index.name, donors), file=sys.stderr)
                elif not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(msg.format(donors), file=sys.stderr)
                elif not pedigree.index.dna_ngs_library:  # pragma: no cover
                    msg = "INFO: pedigree index without DNA NGS library (names: {})"
                    print(msg.format(donors), file=sys.stderr)
                elif pedigree.index.dna_ngs_library.name not in self.ngs_library_to_kit:
                    msg = "INFO: pedigree index DNA NGS library kit not in mapping, skipping: {}"
                    print(msg.format(pedigree.index.dna_ngs_library.name), file=sys.stderr)
                else:
                    yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "targeted_seq_cnv_annotation", "tools_ngs_mapping"),
            (
                "NGS mapping tools not configured but required for targeted sequencing "
                "CNV annotation"
            ),
        )
        self.ensure_w_config(
            ("step_config", "targeted_seq_cnv_annotation", "tools_targeted_seq_cnv_calling"),
            (
                "targeted sequencing CNV calling tools not configured but required for targeted "
                "sequencing CNV annotation"
            ),
        )
