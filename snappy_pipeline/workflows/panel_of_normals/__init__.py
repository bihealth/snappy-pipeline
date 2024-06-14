# -*- coding: utf-8 -*-
"""Implementation of the ``panel_of_normals`` step

The ``panel_of_normals`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and creates background information for somatic variant calling.
and/or somatic copy number calling. This background information is summarized as a
panel of normals.

Usually, the ``panel_of_normals`` step is required by somatic variant calling or
somatic copy number calling tools.

==========
Step Input
==========

The somatic variant calling step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step.

By default, all normal DNA samples in the ``ngs_mapping`` step are using to create the panel of normals.
However, the user can select of subset of those using the ``path_normals_list`` configuration option
(which can be different for the different tools).
In this case, the libraries listed in the file will be used, **even if they are not flagged as corresponding to normal samples**.

===========
Step Output
===========

For each panel of normals tool, the step outputs one set of files describing the panel.
For example, the ``mutect2`` panel of normal generates ``{mapper}.mutect2.pon.vcf.gz``
and associated files (md5 sums indices).

The normals that have been used, as well as the individual files (for example
vcf files for each normal) are kept in the ``work`` directory. This enables the
augmentation of the panel by new files when they become available.

.. warning::

    Panel of normals are powerful tools to reduce systematic bias in the analysis of sequencing data.
    However, they should be built using data generated as similarily as possible.
    In particular, a panel of normals should only contain data collected with the **same** exome enrichment kit.
    It is also essential to use such a panel on tumor samples collected in the same way.

================================
Notes on the ``cnvkit`` workflow
================================

``cnvkit`` is a set of tools originally designed to call somatic copy number alterations from exome data.
Its design is modular, which enables its use for whole genome and amplicon data.

Provided that sufficient normal samples are available, the ``cnvkit`` `documentation <https://cnvkit.readthedocs.io/en/stable/>`_
recommends the creation of a panel of normal (called ``reference``) for exome and whole genome data.

.. note::

    ``cnvkit`` provides a tool to encapsulate common practice workflows (``batch``), depending on the type of data, and on the availability of optional inputs.
    The actual workflow to generate this reference is slightly different between exome and whole genome data.
    The current implementation recapitulates the common practice, while still dispaching computations on multiple cluster nodes.

-----------
Access file
-----------

``cnvkit`` can use a bed file describing the accessible regions for coverage computations.
The ``cnvkit`` distribution provides it for the ``GRCh37`` human genome release, but incompletely only for ``GRCh38``.
Therefore, a tentative ``access`` tool has been added, to generate this bed file when the user knows which locii should be excluded from coverage.
Its output (``output/cnvkit.access/out/cnvkit.access.bed``) is optional, but its presence impacts of the way the target and antitarget regions are computed in whole genome mode.

.. note::

    In a nutshell, for exome data, the accessibility file is only used to create antitarget regions.
    For genome data, it is used by the ``autobin`` tool to compute the average target size used during target regions creation.
    If it is present, the target size is computed in amplicon mode, and when it is absent,
    an accessibility file is created with default settings, which value is used by ``autobin`` is whole genome mode.

To generate the access file from a bed file containing regions to exclude from further coverage computations,
the user must proceed in two steps:

First, she needs to run the ``access`` tool to create the desired access file

.. code-block:: yaml

    panel_of_normals:
        tools: [access]
        access:
            exclude: <absolute path to excluded locii bed file>

This will create ``output/cnvkit.access/out/cnvkit.access.bed`` from the genomic sequence & excluded regions.

------------------------
Panel of normal creation
------------------------

If the user wants to create her own access file, then the panel of normal can only be created after the ``access`` tool has been run.
If she decides that the access file provided in the ``cnvkit`` distribution is suitable (no excluded region),
then she can skip the ``access`` tool step and directly creates her panel of normals.

In both cases, the configuration might read:

.. code-block:: yaml

    panel_of_normals:
        tools: [cnvkit]                                               # , access]
        access: <absolute path to access file>                        # Even when created by the ``access`` tool.
        path_target_regions: <absolute path to baits>                 # Keep empty for WGS data
        path_normals_list: <absolute path to list of normal samples>  # Keep empty to use all available normals

Note that there is no provision (yet) to automatically create separate panel of normals for males & females.
If the number of samples collected in the same fashion is large enough, it is nevertheless the way to achieve best results.

-------
Reports
-------

Report tables can be found in the ``output/{mapper}.cnvkit/report`` directory.
Two tables are produced, grouping results for all normal samples together:

- ``metrics.txt``: coverage metrics over target and antitarget regions.
- ``sex.txt``: prediction of the donor's gender based on the coverage of chromosome X & Y target and antitarget regions.

The cnvkit authors recommend to check these reports to ensure that all data is suitable for panel of normal creation.

================
Notes ``purecn``
================

In the current implementation, the ``purecn`` panel of normals is required when calling somatic copy numbers in the ``somatic_targeted_seq_cnv_calling`` step.
In turn, the ``purecn`` panel of normals requires the availability of a ``mutect2`` panel of normals.
This is because ``mutect2`` is used as somatic variant caller, rather than the older ``mutect`` which is the ``PureCN`` default.

The ``PureCN`` docker container is used, rather than conda environments, because of the complexity of PureCN R packages requirements
(including github-only changes to older packages).

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_panel_of_normals.rst

=====================================
Panel of normals generation for tools
=====================================

- Panel of normal for ``mutect2`` somatic variant caller
- Panel of normal for ``cvnkit`` somatic Copy Number Alterations caller

"""

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import PanelOfNormals as PanelOfNormalsConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Names of the tools that might use panel of normals
TOOLS = ("mutect2", "cnvkit", "access", "purecn")

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = PanelOfNormalsConfigModel.default_config_yaml_string()


class PanelOfNormalsStepPart(BaseStepPart):
    """Base class for panel of normals step parts

    Two steps: the preparation is done on each normal samples separately, and the panel creation
    merges all the individual results in the the panel.
    """

    #: Step name (default, must be overwritten)
    name = None

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.normal_libraries = list(self._get_normal_libraries())
        if self.name and (cfg := self.config.get(self.name)):
            if path := cfg.get("path_normals_list"):
                self.normal_libraries = []
                with open(path, "rt") as f:
                    for line in f:
                        self.normal_libraries.append(line.strip())

    def _get_normal_libraries(self):
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                for _, bio_sample in donor.bio_samples.items():
                    if bio_sample.is_tumor:
                        continue
                    for _, test_sample in bio_sample.test_samples.items():
                        extraction_type = test_sample.extra_infos.get("extractionType", "DNA")
                        if extraction_type.lower() == "dna":
                            for _, ngs_library in test_sample.ngs_libraries.items():
                                yield ngs_library.name

    @staticmethod
    @dictify
    def _get_log_file(tpl):
        """Return all log files files"""
        ext_dict = {
            "conda_list": "conda_list.txt",
            "conda_list_md5": "conda_list.txt.md5",
            "conda_info": "conda_info.txt",
            "conda_info_md5": "conda_info.txt.md5",
            "log": "log",
            "log_md5": "log.md5",
        }
        for key, ext in ext_dict.items():
            yield key, tpl + "." + ext


class PureCnStepPart(PanelOfNormalsStepPart):
    """Creating a panel of normals with GC-corrected coverage"""

    #: Step name
    name = "purecn"

    #: Actions
    actions = ("install", "prepare", "coverage", "create_panel")

    #: Resources
    resource_usage = {
        "install": ResourceUsage(
            threads=1,
            time="01:00:00",
            memory="24G",
        ),
        "prepare": ResourceUsage(
            threads=1,
            time="04:00:00",  # 4 hours
            memory="24G",
        ),
        "coverage": ResourceUsage(
            threads=1,
            time="04:00:00",  # 4 hours
            memory="24G",
        ),
        "create_panel": ResourceUsage(
            threads=1,
            time="12:00:00",  # 12 hours
            memory="32G",
        ),
    }

    def get_input_files(self, action):
        self._validate_action(action)
        self.ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        if action == "prepare":
            return {"container": "work/containers/out/purecn.simg"}
        if action == "coverage":
            return self._get_input_files_coverage
        if action == "create_panel":
            return self._get_input_files_create

    @dictify
    def _get_input_files_coverage(self, wildcards):
        yield "container", "work/containers/out/purecn.simg"
        yield (
            "intervals",
            "work/purecn/out/{}_{}.list".format(
                self.config.purecn.enrichment_kit_name,
                self.config.purecn.genome_name,
            ),
        )
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
        yield "bam", self.ngs_mapping(tpl.format(**wildcards))

    @dictify
    def _get_input_files_create(self, wildcards):
        yield "container", "work/containers/out/purecn.simg"
        tpl = "work/{mapper}.purecn/out/{mapper}.purecn.{library_name}_coverage_loess.txt.gz"
        yield (
            "normals",
            [
                tpl.format(mapper=wildcards.mapper, library_name=lib)
                for lib in self.normal_libraries
            ],
        )

    def get_output_files(self, action):
        self._validate_action(action)

        if action == "install":
            return {"container": "work/containers/out/purecn.simg"}
        if action == "prepare":
            base_out = "{}_{}".format(
                self.config.purecn.enrichment_kit_name,
                self.config.purecn.genome_name,
            )
            return {
                "intervals": "work/purecn/out/" + base_out + ".list",
                "optimized": "work/purecn/out/" + base_out + ".bed.gz",
                "tbi": "work/purecn/out/" + base_out + ".bed.gz.tbi",
                "intervals_md5": "work/purecn/out/" + base_out + ".list.md5",
                "optimized_md5": "work/purecn/out/" + base_out + ".bed.gz.md5",
                "tbi_md5": "work/purecn/out/" + base_out + ".bed.gz.tbi.md5",
            }
        if action == "coverage":
            return {
                "coverage": "work/{mapper}.purecn/out/{mapper}.purecn.{library_name,.+-DNA[0-9]+-WES[0-9]+}_coverage_loess.txt.gz"
            }
        if action == "create_panel":
            return {
                "db": "work/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.rds",
                "db_md5": "work/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.rds.md5",
                "mapbias": "work/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.rds",
                "mapbias_md5": "work/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.rds.md5",
                "lowcov": "work/{mapper}.purecn/out/{mapper}.purecn.low_coverage_targets.bed",
                "hq": "work/{mapper}.purecn/out/{mapper}.purecn.hq_sites.bed",
                "plot": "work/{mapper}.purecn/out/{mapper}.purecn.interval_weights.png",
            }

    def get_log_file(self, action):
        tpls = {
            "install": "work/containers/log/purecn",
            "prepare": "work/purecn/log/{}_{}".format(
                self.config.purecn.enrichment_kit_name,
                self.config.purecn.genome_name,
            ),
            "coverage": "work/{mapper}.purecn/log/{mapper}.purecn.{library_name,.+-DNA[0-9]+-WES[0-9]+}",
            "create_panel": "work/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals",
        }
        assert action in self.actions
        return self._get_log_file(tpls[action])


class Mutect2StepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    #: Step name
    name = "mutect2"

    #: Class available actions
    actions = ("prepare_panel", "create_panel")

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage = {
        "prepare_panel": ResourceUsage(
            threads=2,
            time="3-00:00:00",  # 3 days
            memory="8G",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            time="48:00:00",  # 48 hours
            memory="30G",
        ),
    }

    def check_config(self):
        if self.name not in self.config.tools:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_input_files(self, action):
        """Return input files for mutect2 variant calling"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "prepare_panel": self._get_input_files_prepare_panel,
            "create_panel": self._get_input_files_create_panel,
        }
        return mapping[action]

    def _get_input_files_prepare_panel(self, wildcards):
        """Helper wrapper function for single sample panel preparation"""
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        return {"normal_bam": bam, "normal_bai": bam + ".bai"}

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for merging individual results & panel creation"""
        paths = []
        tpl = "work/{mapper}.{tool}/out/{mapper}.{tool}.{normal_library}.prepare.vcf.gz"
        for normal in self.normal_libraries:
            paths.append(tpl.format(normal_library=normal, tool=self.name, **wildcards))
        return {"normals": paths}

    @dictify
    def get_output_files(self, action):
        """Return panel of normal files"""
        self._validate_action(action)
        ext_dict = {
            "vcf": "vcf.gz",
            "vcf_md5": "vcf.gz.md5",
            "vcf_tbi": "vcf.gz.tbi",
            "vcf_tbi_md5": "vcf.gz.tbi.md5",
        }
        tpls = {
            "prepare_panel": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare",
            "create_panel": "work/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals",
        }
        for key, ext in ext_dict.items():
            yield key, tpls[action] + "." + ext
        if action == "create_panel":
            yield "db", "work/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB.tar.gz"
            yield "db_md5", "work/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB.tar.gz.md5"

    @classmethod
    def get_log_file(cls, action):
        """Return panel of normal files"""
        tpls = {
            "prepare_panel": "work/{mapper}.mutect2/log/{mapper}.mutect2.{normal_library}.prepare",
            "create_panel": "work/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals",
        }
        assert action in cls.actions
        return cls._get_log_file(tpls[action])


class CnvkitStepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "target",
        "antitarget",
        "coverage",
        "create_panel",
        "report",
    )

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage = {
        "target": ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8G",
        ),
        "antitarget": ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8G",
        ),
        "coverage": ResourceUsage(
            threads=8,
            time="02:00:00",  # 2 hours
            memory="16G",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="16G",
        ),
        "report": ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="16G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.is_wgs = self.config.cnvkit.path_target_regions == ""

    def check_config(self):
        if self.name not in self.config.tools:
            return  # cnvkit not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_args(self, action):
        self._validate_action(action)
        if self.is_wgs:
            method = "wgs"
        else:
            method = "hybrid"
        return {"method": method, "flat": (len(self.normal_libraries) == 0)}

    def get_input_files(self, action):
        """Return input files for cnvkit panel of normals creation"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "target": self._get_input_files_target,
            "antitarget": self._get_input_files_antitarget,
            "coverage": self._get_input_files_coverage,
            "create_panel": self._get_input_files_create_panel,
            "report": self._get_input_files_report,
            "access": None,
        }
        return mapping[action]

    def _get_input_files_target(self, wildcards):
        """Helper wrapper function to estimate target average size in wgs mode"""
        if not self.is_wgs:
            return {}
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        bams = [
            ngs_mapping(tpl.format(mapper=wildcards["mapper"], normal_library=x))
            for x in self.normal_libraries
        ]
        bais = [x + ".bai" for x in bams]
        input_files = {"bams": bams, "bais": bais}
        return input_files

    def _get_input_files_antitarget(self, wildcards):
        """Helper wrapper function for computing antitarget locations"""
        if self.is_wgs:
            return {}
        return {
            "target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed".format(**wildcards),
        }

    def _get_input_files_coverage(self, wildcards):
        """Helper wrapper function for computing coverage"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        return {
            "target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed".format(**wildcards),
            "antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed".format(
                **wildcards
            ),
            "bam": bam,
            "bai": bam + ".bai",
        }

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for computing panel of normals"""
        tpl = "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.targetcoverage.cnn"
        targets = [
            tpl.format(mapper=wildcards["mapper"], normal_library=x) for x in self.normal_libraries
        ]
        tpl = "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.antitargetcoverage.cnn"
        antitargets = [
            tpl.format(mapper=wildcards["mapper"], normal_library=x) for x in self.normal_libraries
        ]
        tpl = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.{normal_library}.coverage.{ext}"
        logs = [
            tpl.format(mapper=wildcards["mapper"], normal_library=x, ext=ext)
            for x in self.normal_libraries
            for ext in ("log", "conda_list.txt", "conda_info.txt")
        ]
        return {
            "target": (
                targets
                if targets
                else "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed".format(**wildcards)
            ),
            "antitarget": (
                antitargets
                if antitargets
                else "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed".format(**wildcards)
            ),
            "logs": logs if targets or antitargets else [],
        }

    def _get_input_files_report(self, wildcards):
        """Helper wrapper function for the panel of normals report"""
        tpl = "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.targetcoverage.cnn"
        targets = [
            tpl.format(mapper=wildcards["mapper"], normal_library=x) for x in self.normal_libraries
        ]
        tpl = "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.antitargetcoverage.cnn"
        antitargets = [
            tpl.format(mapper=wildcards["mapper"], normal_library=x) for x in self.normal_libraries
        ]
        return {
            "target": targets,
            "antitarget": antitargets,
        }

    def get_output_files(self, action):
        """Return panel of normal files"""
        if action == "target":
            return self._get_output_files_target()
        elif action == "antitarget":
            return self._get_output_files_antitarget()
        elif action == "coverage":
            return self._get_output_files_coverage()
        elif action == "create_panel":
            return self._get_output_files_create_panel()
        elif action == "report":
            return self._get_output_files_report()
        elif action == "access":
            return self._get_output_files_access()
        else:
            self._validate_action(action)

    def _get_output_files_target(self):
        return {
            "target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed",
            "target_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed.md5",
        }

    def _get_output_files_antitarget(self):
        return {
            "antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed",
            "antitarget_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed.md5",
        }

    def _get_output_files_coverage(self):
        return {
            "target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.targetcoverage.cnn",
            "target_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.targetcoverage.cnn.md5",
            "antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.antitargetcoverage.cnn",
            "antitarget_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{normal_library}.antitargetcoverage.cnn.md5",
        }

    def _get_output_files_create_panel(self):
        return {
            "panel": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.cnn",
            "panel_md5": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.cnn.md5",
            "log": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.merged.tar.gz",
            "log_md5": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.merged.tar.gz.md5",
        }

    def _get_output_files_report(self):
        return {
            "sex": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.sex.tsv",
            "sex_md5": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.sex.tsv.md5",
            "metrics": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.metrics.tsv",
            "metrics_md5": "work/{mapper}.cnvkit/report/{mapper}.cnvkit.metrics.tsv.md5",
        }

    def _get_output_files_access(self):
        return {
            "access": "work/cnvkit.access/out/cnvkit.access.bed",
            "access_md5": "work/cnvkit.access/out/cnvkit.access.bed.md5",
        }

    @classmethod
    def get_log_file(cls, action):
        """Return panel of normal files"""
        tpls = {
            "target": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.target",
            "antitarget": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.antitarget",
            "coverage": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.{normal_library}.coverage",
            "create_panel": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.panel_of_normals",
            "report": "work/{mapper}.cnvkit/log/{mapper}.cnvkit.report",
            "access": "work/cnvkit.access/log/cnvkit.access",
        }
        assert action in cls.actions
        return cls._get_log_file(tpls[action])


class AccessStepPart(PanelOfNormalsStepPart):
    """Utility to create access file for cnvkit"""

    name = "access"
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8G",
        )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        return None

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        tpl = "work/cnvkit.access/out/cnvkit.access.bed"
        return {"access": tpl, "access_md5": tpl + ".md5"}

    @classmethod
    def get_log_file(cls, action):
        """Return log files"""
        assert action in cls.actions
        return cls._get_log_file("work/cnvkit.access/log/cnvkit.access")


class PanelOfNormalsWorkflow(BaseStep):
    """Creates a panel of normals"""

    # Workflow name
    name = "panel_of_normals"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

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
            config_model_class=PanelOfNormalsConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                Mutect2StepPart,
                CnvkitStepPart,
                AccessStepPart,
                PureCnStepPart,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        result_files = []

        log_ext_list = [
            "log",
            "log.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
        ]

        if "mutect2" in set(self.config.tools) & set(TOOLS):
            tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals.{ext}"
            ext_list = ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB.{ext}"
            ext_list = ("tar.gz", "tar.gz.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if "cnvkit" in set(self.config.tools) & set(TOOLS):
            tpls = [
                ("output/{mapper}.cnvkit/out/{mapper}.cnvkit.target.{ext}", ("bed", "bed.md5")),
                ("output/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.{ext}", ("bed", "bed.md5")),
                (
                    "output/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.{ext}",
                    ("cnn", "cnn.md5"),
                ),
                (
                    "output/{mapper}.cnvkit/report/{mapper}.cnvkit.sex.{ext}",
                    ("tsv", "tsv.md5"),
                ),
                (
                    "output/{mapper}.cnvkit/report/{mapper}.cnvkit.metrics.{ext}",
                    ("tsv", "tsv.md5"),
                ),
            ]
            for tpl, ext_list in tpls:
                result_files.extend(self._expand_result_files(tpl, ext_list))
            tpls = [
                "output/{mapper}.cnvkit/log/{mapper}.cnvkit.target.{ext}",
                "output/{mapper}.cnvkit/log/{mapper}.cnvkit.antitarget.{ext}",
                "output/{mapper}.cnvkit/log/{mapper}.cnvkit.panel_of_normals.{ext}",
                "output/{mapper}.cnvkit/log/{mapper}.cnvkit.report.{ext}",
            ]
            for tpl in tpls:
                result_files.extend(self._expand_result_files(tpl, log_ext_list))
            tpl = "output/{mapper}.cnvkit/log/{mapper}.cnvkit.merged.tar.gz{ext}"
            result_files.extend(self._expand_result_files(tpl, ("", ".md5")))

        if "access" in set(self.config.tools) & set(TOOLS):
            tpl = "output/cnvkit.access/out/cnvkit.access.bed"
            result_files.extend([tpl + md5 for md5 in ("", ".md5")])
            tpl = "output/cnvkit.access/log/cnvkit.access.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if "purecn" in set(self.config.tools) & set(TOOLS):
            tpl = "output/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))
            tpl = "output/purecn/out/{}_{}.{{ext}}".format(
                self.config.purecn.enrichment_kit_name,
                self.config.purecn.genome_name,
            )
            ext_list = ("list", "list.md5", "bed.gz", "bed.gz.md5", "bed.gz.tbi", "bed.gz.tbi.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/purecn/log/{}_{}.{{ext}}".format(
                self.config.purecn.enrichment_kit_name,
                self.config.purecn.genome_name,
            )
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        return result_files

    def _expand_result_files(self, tpl, ext_list):
        for mapper in self.w_config.step_config["ngs_mapping"].tools.dna:
            for ext in ext_list:
                yield tpl.format(mapper=mapper, ext=ext)
