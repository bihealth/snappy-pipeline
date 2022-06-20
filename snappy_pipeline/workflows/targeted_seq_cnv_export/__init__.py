# -*- coding: utf-8 -*-
"""Implementation of the ``targeted_seq_cnv_export`` step.

The ``targeted_seq_cnv_export`` step takes as the input the results of the
``targeted_seq_cnv_annotation`` step and uses ``varfish-annotator-cli annotate`` commmand to
create files fit for import into VarFish Server.

==========
Stability
==========

TODO

==========
Step Input
==========

The WGS SV export step uses Snakemake sub workflows for using the result of the
``targeted_seq_cnv_export`` step.

===========
Step Output
===========

TODO

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

.. include:: DEFAULT_CONFIG_targeted_seq_cnv.rst

==================
Parallel Execution
==================

Parallel execution is not performed currently.
"""

from collections import OrderedDict
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
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.targeted_seq_cnv_annotation import TargetedSeqCnvAnnotationWorkflow
from snappy_pipeline.workflows.targeted_seq_cnv_calling import TargetedSeqCnvCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")

#: Infixes to use for file name generation
INFIXES = ("gts", "feature-effects", "db-infos")

#: Extensions of files to create as main payload (VCF)
VCF_EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Default configuration for the targeted_seq_cnv_export step
DEFAULT_CONFIG = r"""
# Default configuration targeted_seq_cnv_export.
step_config:
  targeted_seq_cnv_export:
    path_ngs_mapping: ../ngs_mapping
    path_targeted_seq_cnv_annotation: ../targeted_seq_cnv_annotation
    tools_ngs_mapping: null
    tools_targeted_seq_cnv_calling: null
    path_refseq_ser: REQUIRED    # REQUIRED: path to RefSeq .ser file
    path_ensembl_ser: REQUIRED   # REQUIRED: path to ENSEMBL .ser file
    path_db: REQUIRED            # REQUIRED: spath to annotator DB file to use
"""


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    #: step name
    name = "varfish_annotator"

    #: Class available actions
    actions = ("annotate",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/.done"
        )
        # Build shortcut from index library name to donor
        self.index_ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_donor.update(sheet.index_ngs_library_to_donor)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        xhmm_config = DictQuery(self.w_config).get("step_config/targeted_seq_cnv_calling/xhmm")
        if not xhmm_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.parent.all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"

        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in xhmm_config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    @dictify
    def get_input_files(self, action):
        """Return path to pedigree input file"""
        # Validate action
        self._validate_action(action)
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        tpl = (
            "output/{mapper}.{var_caller}.annotated.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.annotated.{index_ngs_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        targeted_seq_cnv_annotation = self.parent.sub_workflows["targeted_seq_cnv_annotation"]
        for key, ext in key_ext.items():
            yield key, targeted_seq_cnv_annotation(tpl + ext)

    @dictify
    def _build_ngs_library_to_kit(self):
        xhmm_config = DictQuery(self.w_config).get("step_config/targeted_seq_cnv_calling/xhmm")
        if not xhmm_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.parent.all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"
        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in xhmm_config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
        )
        for infix in INFIXES:
            key = infix.replace("-", "_")
            yield key, prefix + ".%s.tsv.gz" % infix
            yield key + "_md5", prefix + ".%s.tsv.gz.md5" % infix

    @dictify
    def _get_log_file(self, action):
        assert action == "annotate"
        prefix = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
            "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
        )

        key_ext = (
            ("wrapper", ".wrapper.py"),
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

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
            time="4-04:00:00",  # 4 days and 4 hours
            memory=f"{7 * 1024 * 2}M",
        )

    def get_params(self, action):
        assert action == "annotate"

        def get_params_func(wildcards):
            result = {"is_wgs": False, "step_name": "targeted_seq_cnv_export"}
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            for donor in pedigree.donors:
                if (
                    donor.dna_ngs_library
                    and donor.dna_ngs_library.extra_infos.get("libraryType") == "WGS"
                ):
                    result["is_wgs"] = True
                    return result
            return result

        return get_params_func


class TargetedSeqCnvExportWorkflow(BaseStep):
    """Perform germline Targeted Sequencing CNV export"""

    name = "targeted_seq_cnv_export"
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
            (TargetedSeqCnvAnnotationWorkflow, TargetedSeqCnvCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VarfishAnnotatorAnnotateStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow(
            "targeted_seq_cnv_annotation", self.config["path_targeted_seq_cnv_annotation"]
        )
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Copy over "tools" setting from targeted_seq_cnv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_targeted_seq_cnv_calling"]:
            # Remove plain ERDS as it does not do multi-sample genotypeing
            tools = self.w_config["step_config"]["targeted_seq_cnv_calling"]["tools"]
            self.config["tools_targeted_seq_cnv_calling"] = [t for t in tools if t != "erds"]
        # Build mapping from NGS DNA library to library kit.
        self.ngs_library_to_kit = self.sub_steps["varfish_annotator"].ngs_library_to_kit

    @listify
    def all_donors(self, include_background=True):
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
            for donor in self.all_donors()
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
        name_pattern = "{mapper}.{caller}.varfish_annotated.{index_library.name}"
        if "xhmm" in self.config["tools_targeted_seq_cnv_calling"]:
            min_kit_usages = 10
            chosen_kits = {kit for kit in library_kits if kit_counts.get(kit, 0) > min_kit_usages}
            chosen_donors = [
                donor.name
                for donor in donors
                if self.ngs_library_to_kit.get(donor.dna_ngs_library.name) in chosen_kits
            ]
            yield from self._yield_result_files(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
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
            chosen_kits = {kit for kit in library_kits if kit_counts.get(kit, 0) > min_kit_usages}
            chosen_donors = [
                donor.name
                for donor in donors
                if self.ngs_library_to_kit.get(donor.dna_ngs_library.name) in chosen_kits
            ]
            yield from self._yield_result_files(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
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
            ("step_config", "targeted_seq_cnv_export", "path_targeted_seq_cnv_annotation"),
            "Path to WGS SV annotation not configured but required for WGS SV export",
        )
