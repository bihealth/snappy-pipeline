# -*- coding: utf-8 -*-
"""Implementation of the ``variant_export`` step.

The ``variant_export`` step takes as the input the results of the ``variant_calling`` step and
uses ``varfish-annotator-cli annotate`` commmand to create files fit for import into VarFish
Server.

==========
Stability
==========

TODO

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``variant_calling`` step.

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

The default configuration is as follows. Note that the ``path_jannovar_ser`` parameter must
be set to point to the desired transcript annotations db as generated by ``jannovar download``.

.. include:: DEFAULT_CONFIG_variant_export.rst

==================
Parallel Execution
==================

Parallel execution is not performed currently.
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
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")

#: Infixes to use for file name generation
INFIXES = ("gts", "db-infos", "bam-qc")

# TODO: the number of restart times is high because tabix in HTSJDK/Jannovar is flaky...

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_export.
step_config:
  variant_export:
    path_ngs_mapping: ../ngs_mapping
    path_variant_calling: ../variant_calling
    tools_ngs_mapping: null
    tools_variant_calling: null
    path_exon_bed: REQUIRED      # REQUIRED: exon BED file to use when handling WGS data
    path_refseq_ser: REQUIRED    # REQUIRED: path to RefSeq .ser file
    path_ensembl_ser: REQUIRED   # REQUIRED: path to ENSEMBL .ser file
    path_db: REQUIRED            # REQUIRED: spath to annotator DB file to use
"""


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    #: Step name
    name = "varfish_annotator"

    #: Class available actions
    actions = ("annotate", "bam_qc")

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/.done"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return path to pedigree input file"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_%s" % action)

    @dictify
    def _get_input_files_annotate(self, wildcards):
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        tpl = (
            "output/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["variant_calling"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl + ext).format(
                mapper=wildcards.mapper,
                var_caller=wildcards.var_caller,
                index_ngs_library=wildcards.index_ngs_library,
            )

    def _get_input_files_bam_qc(self, wildcards):
        vals = {"mapper": wildcards.mapper, "var_caller": wildcards.var_caller}
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Get names of primary libraries of the selected pedigree.  The pedigree is selected
        # by the primary DNA NGS library of the index.
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "cov_qc": []}
        for donor in pedigree.donors:
            if not donor.dna_ngs_library:
                continue
            tpl = (
                "output/{mapper}.{index_ngs_library}/report/bam_qc/"
                "{mapper}.{index_ngs_library}.bam.%s.txt"
            ).format(**vals, index_ngs_library=donor.dna_ngs_library.name)
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(ngs_mapping(tpl % key))
            if not donor.dna_ngs_library.name in self.parent.ngs_library_to_kit:
                continue
            path = (
                "output/{mapper}.{index_ngs_library}/report/cov_qc/"
                "{mapper}.{index_ngs_library}.txt"
            ).format(**vals, index_ngs_library=donor.dna_ngs_library.name)
            result["cov_qc"].append(ngs_mapping(path))
        return result

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)

        infixes = None
        if action == "annotate":
            infixes = ("gts", "db-infos")
        elif action == "bam_qc":
            infixes = ("bam-qc",)
        prefix = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
        )
        for infix in infixes:
            key = infix.replace("-", "_")
            yield key, prefix + ".%s.tsv.gz" % infix
            yield key + "_md5", prefix + ".%s.tsv.gz.md5" % infix

    @dictify
    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
            "{mapper}.{var_caller}.varfish_annotated.%s.{index_ngs_library}"
        ) % action
        key_ext = (
            ("wrapper", ".wrapper.py"),
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )

        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

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
            time="4-04:00:00",  # 4 day and 4 hours
            memory=f"{7 * 1024 * 2}M",
        )

    def get_params(self, action):
        assert action == "annotate", f"Option only valid for action 'annotate' (used: '{action}')."

        def get_params_func(wildcards):
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            for donor in pedigree.donors:
                if (
                    donor.dna_ngs_library
                    and donor.dna_ngs_library.extra_infos.get("libraryType") == "WGS"
                ):
                    return {"is_wgs": True}
            return {"is_wgs": False}

        return get_params_func


class VariantExportWorkflow(BaseStep):
    """Perform germline variant export"""

    #: Workflow name
    name = "variant_export"

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
            (VariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VarfishAnnotatorAnnotateStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("variant_calling", self.config["path_variant_calling"])
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_variant_calling"]:
            self.config["tools_variant_calling"] = self.w_config["step_config"]["variant_calling"][
                "tools"
            ]
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        cov_config = DictQuery(self.w_config).get("step_config/ngs_mapping/target_coverage_report")
        if not cov_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for sheet in self.shortcut_sheets:
                for donor in sheet.donors:
                    for bio_sample in donor.bio_samples.values():
                        for test_sample in bio_sample.test_samples.values():
                            for library in test_sample.ngs_libraries.values():
                                yield library.name, "default"

        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in cov_config["path_target_interval_list_mapping"]
            if item["name"] != "__default__"
        }
        result = {}
        for sheet in self.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    for test_sample in bio_sample.test_samples.values():
                        for library in test_sample.ngs_libraries.values():
                            if library.extra_infos.get("libraryKit"):
                                library_kit = library.extra_infos.get("libraryKit")
                                for pattern, name in regexes.items():
                                    if re.match(pattern, library_kit):
                                        yield library.name, name
        return result

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_pattern = "{mapper}.{caller}.varfish_annotated.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + ".{infix}{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
            infix=INFIXES,
            ext=EXTS,
        )
        for action in ("annotate", "bam_qc"):
            name_pattern_action = name_pattern.replace(
                ".varfish_annotated.", ".varfish_annotated.%s." % action
            )
            prefix = os.path.join("output", name_pattern, "log", name_pattern_action + "{ext}")
            yield from self._yield_result_files(
                prefix,
                mapper=self.config["tools_ngs_mapping"],
                caller=self.config["tools_variant_calling"],
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:  # pragma: no cover
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                elif not pedigree.index.dna_ngs_library:  # pragma: no cover
                    msg = "INFO: pedigree index without DNA NGS library (names: {})"
                    print(
                        msg.format(  # pragma: no cover
                            list(sorted(d.name for d in pedigree.donors))
                        ),
                        file=sys.stderr,
                    )
                    continue  # pragma: no cover
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_export", "path_variant_calling"),
            "Path to variant calling not configured but required for variant annotation",
        )
