# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_sv_export`` step.

The ``wgs_sv_export`` step takes as the input the results of the ``wgs_sv_annotation`` step and
uses ``varfish-annotator-cli annotate`` commmand to create files fit for import into VarFish
Server.

==========
Stability
==========

TODO

==========
Step Input
==========

The WGS SV export step uses Snakemake sub workflows for using the result of the
``wgs_sv_export`` step.

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

.. include:: DEFAULT_CONFIG_wgs_sv.rst

==================
Parallel Execution
==================

Parallel execution is not performed currently.
"""

import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.wgs_sv_annotation import WgsSvAnnotationWorkflow
from snappy_pipeline.workflows.wgs_sv_calling import WgsSvCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")

#: Infixes to use for file name generation
INFIXES = ("gts", "feature-effects", "db-infos")

#: Default configuration for the wgs_sv_export step
DEFAULT_CONFIG = r"""
# Default configuration wgs_sv_export.
step_config:
  wgs_sv_export:
    path_wgs_sv_annotation: ../wgs_sv_annotation
    tools_ngs_mapping: null
    tools_wgs_sv_calling: null
    path_refseq_ser: REQUIRED    # REQUIRED: path to RefSeq .ser file
    path_ensembl_ser: REQUIRED   # REQUIRED: path to ENSEMBL .ser file
    path_db: REQUIRED            # REQUIRED: spath to annotator DB file to use
"""


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    name = "varfish_annotator"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/.done"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    @dictify
    def get_input_files(self, action):
        """Return path to pedigree input file"""
        assert action == "annotate"
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        tpl = (
            "output/{mapper}.{var_caller}.annotated.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.annotated.{index_ngs_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        wgs_sv_annotation = self.parent.sub_workflows["wgs_sv_annotation"]
        for key, ext in key_ext.items():
            yield key, wgs_sv_annotation(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        assert action == "annotate"
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

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["wgs_sv_export_varfish_annotator_annotate_svs"] = {
            "mem": 7 * 1024 * 2,
            "time": "100:00",
            "ntasks": 2,
        }

    def get_params(self, action):
        assert action == "annotate"

        def get_params_func(wildcards):
            result = {"is_wgs": True, "step_name": "wgs_sv_export"}
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


class WgsSvExportWorkflow(BaseStep):
    """Perform germline WGS SV export"""

    name = "wgs_sv_export"
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
            (WgsSvAnnotationWorkflow, WgsSvCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VarfishAnnotatorAnnotateStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("wgs_sv_annotation", self.config["path_wgs_sv_annotation"])
        # Copy over "tools" setting from wgs_sv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_wgs_sv_calling"]:
            self.config["tools_wgs_sv_calling"] = self.w_config["step_config"]["wgs_sv_calling"][
                "tools"
            ]["dna"]

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_pattern = "{mapper}.{caller}.varfish_annotated.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + ".{infix}{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_wgs_sv_calling"],
            infix=INFIXES,
            ext=EXTS,
        )
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_wgs_sv_calling"],
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
            ("step_config", "wgs_sv_export", "path_wgs_sv_annotation"),
            ("Path to WGS SV annotation not configured but required for WGS SV export"),
        )
