# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_cnv_export`` step.

The ``wgs_cnv_export`` step takes as the input the results of the ``wgs_cnv_annotation`` step and
uses ``varfish-annotator-cli annotate`` commmand to create files fit for import into VarFish
Server.

==========
Stability
==========

This step is considered is considered stable for short Illumina reads.

==========
Step Input
==========

The WGS CNV export step uses Snakemake sub workflows for using the result of the
``wgs_cnv_calling`` step.

===========
Step Output
===========

For all pedigrees, annotation will be performed on the merged VCF file contained all samples.
The name of the primary DNA NGS library of the index will be used as an identification token
in the output file.  For each read mapper, sv caller, and pedigree, the following is an example of
the files that will be generated:

::
    output/
    +-- varfish_annotated.P001-N1-DNA1-WGS1
    |   `-- out
    |   |   |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz
    |   |   |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz.md5
    |   |   |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz
    |   |   |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz.md5
    |   |   |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz
    |   |   +-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz.md5
    |   |
    |   +-- log
    |       |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.conda_info.txt
    |       |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.conda_info.txt.md5
    |       |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.conda_list.txt
    |       |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.conda_list.txt.md5
    |       |-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.log
    |       +-- bwa.gcnv.varfish_annotated.P001-N1-DNA1-WGS1.log.md5
    |
    [...]

====================
Global Configuration
====================

Not applicable.

=====================
Default Configuration
=====================

.. include:: DEFAULT_CONFIG_wgs_cnv.rst

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
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.wgs_cnv_annotation import WgsCnvAnnotationWorkflow
from snappy_pipeline.workflows.wgs_cnv_calling import WgsCnvCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")

#: Infixes to use for file name generation
INFIXES = ("gts", "feature-effects", "db-infos")

#: Default configuration for the wgs_cnv_export step
DEFAULT_CONFIG = r"""
# Default configuration wgs_cnv_export.
step_config:
  wgs_cnv_export:
    path_wgs_cnv_annotation: ../wgs_cnv_annotation
    tools_ngs_mapping: null
    tools_wgs_cnv_calling: null
    release: GRCh37              # REQUIRED: default 'GRCh37'
    path_refseq_ser: REQUIRED    # REQUIRED: path to RefSeq .ser file
    path_ensembl_ser: REQUIRED   # REQUIRED: path to ENSEMBL .ser file
    path_db: REQUIRED            # REQUIRED: spath to annotator DB file to use
"""


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    #: Step name
    name = "varfish_annotator"

    #: Class available actions
    actions = ("annotate",)

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
        # Validate action
        self._validate_action(action)
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        tpl = (
            "output/{mapper}.{var_caller}.annotated.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.annotated.{index_ngs_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        wgs_cnv_annotation = self.parent.sub_workflows["wgs_cnv_annotation"]
        for key, ext in key_ext.items():
            yield key, wgs_cnv_annotation(tpl + ext)

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
        # Validate action
        self._validate_action(action)
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
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="4-04:00:00",  # 4 days and 4 hours
            memory=f"{7 * 1024 * 2}M",
        )

    def get_params(self, action):
        # Validate action
        self._validate_action(action)

        def get_params_func(wildcards):
            result = {"is_wgs": True, "step_name": "wgs_cnv_export"}
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


class WgsCnvExportWorkflow(BaseStep):
    """Perform germline WGS SV export"""

    #: Workflow name
    name = "wgs_cnv_export"

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
            (WgsCnvAnnotationWorkflow, WgsCnvCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VarfishAnnotatorAnnotateStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("wgs_cnv_annotation", self.config["path_wgs_cnv_annotation"])
        # Copy over "tools" setting from wgs_cnv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_wgs_cnv_calling"]:
            # Remove plain ERDS as it does not do multi-sample genotypeing
            tools = self.w_config["step_config"]["wgs_cnv_calling"]["tools"]
            self.config["tools_wgs_cnv_calling"] = [t for t in tools if t != "erds"]

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_patternn = "{mapper}.{caller}.varfish_annotated.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_patternn, "out", name_patternn + ".{infix}{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_wgs_cnv_calling"],
            infix=INFIXES,
            ext=EXTS,
        )
        yield from self._yield_result_files(
            os.path.join("output", name_patternn, "log", name_patternn + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_wgs_cnv_calling"],
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
            ("step_config", "wgs_cnv_export", "path_wgs_cnv_annotation"),
            "Path to WGS SV annotation not configured but required for WGS SV export",
        )
