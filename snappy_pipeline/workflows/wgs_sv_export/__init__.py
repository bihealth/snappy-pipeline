# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_sv_export`` step.

The ``wgs_sv_export`` step takes as the input the results of the ``wgs_sv_annotation`` step and
uses ``varfish-annotator-cli annotate-sv`` commmand to create files fit for import into VarFish
Server.

==========
Stability
==========

This step is considered is considered stable for short Illumina reads.

==========
Step Input
==========

The WGS SV export step uses Snakemake sub workflows for using the result of the
``wgs_sv_calling`` step.

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
    |   |   |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz
    |   |   |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz.md5
    |   |   |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz
    |   |   |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz.md5
    |   |   |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz
    |   |   +-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz.md5
    |   |
    |   +-- log
    |       |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.conda_info.txt
    |       |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.conda_info.txt.md5
    |       |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.conda_list.txt
    |       |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.conda_list.txt.md5
    |       |-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.log
    |       +-- bwa.delly2.varfish_annotated.P001-N1-DNA1-WGS1.log.md5
    |
    [...]


====================
Global Configuration
====================

Not applicable.

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
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.wgs_sv_annotation import WgsSvAnnotationWorkflow
from snappy_pipeline.workflows.wgs_sv_calling import WgsSvCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

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
    path_wgs_sv_calling: ../wgs_sv_calling
    tools_ngs_mapping: null
    tools_wgs_sv_calling: null
    release: GRCh37              # REQUIRED: default 'GRCh37'
    path_refseq_ser: REQUIRED    # REQUIRED: path to RefSeq .ser file
    path_ensembl_ser: REQUIRED   # REQUIRED: path to ENSEMBL .ser file
    path_db: REQUIRED            # REQUIRED: spath to annotator DB file to use
    varfish_server_compatibility: false # OPTIONAL: build output compatible with varfish-server v1.2 (Anthenea) and early versions of the v2 (Bollonaster)
"""


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate-sv"."""

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

    def get_input_files(self, action):
        """Return path to pedigree input file"""
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            tpl = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
            yield "ped", tpl.format(**wildcards)
            if wildcards.var_caller == "popdel":
                tpl = (
                    "output/{mapper}.{var_caller}.annotated.{index_ngs_library}/out/"
                    "{mapper}.{var_caller}.annotated.{index_ngs_library}"
                )
                subworkflow = self.parent.sub_workflows["wgs_sv_annotation"]
            else:
                tpl = (
                    "output/{mapper}.{var_caller}.{index_ngs_library}/out/"
                    "{mapper}.{var_caller}.{index_ngs_library}"
                )
                subworkflow = self.parent.sub_workflows["wgs_sv_calling"]
            key_ext = {
                "vcf": ".vcf.gz",
                "vcf_md5": ".vcf.gz.md5",
                "tbi": ".vcf.gz.tbi",
                "tbi_md5": ".vcf.gz.tbi.md5",
            }
            for key, ext in key_ext.items():
                yield key, subworkflow(tpl.format(**wildcards) + ext)

        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files for the export"""
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
            result = {
                "is_wgs": True,
                "step_name": "wgs_sv_export",
                "varfish_server_compatibility": self.config["varfish_server_compatibility"],
            }
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

    #: Workflow name
    name = "wgs_sv_export"

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
            (WgsSvAnnotationWorkflow, WgsSvCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VarfishAnnotatorAnnotateStepPart, LinkOutStepPart)
        )
        # Copy over "tools" setting from wgs_sv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_wgs_sv_calling"]:
            self.config["tools_wgs_sv_calling"] = self.w_config["step_config"]["wgs_sv_calling"][
                "tools"
            ]["dna"]
        # Register sub workflows
        if "popdel" in self.config["tools_wgs_sv_calling"]:
            self.register_sub_workflow("wgs_sv_annotation", self.config["path_wgs_sv_annotation"])
        if (
            "delly2" in self.config["tools_wgs_sv_calling"]
            or "sniffles2" in self.config["tools_wgs_sv_calling"]
        ):
            self.register_sub_workflow("wgs_sv_calling", self.config["path_wgs_sv_calling"])

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
            "Path to WGS SV annotation not configured but required for WGS SV export",
        )
