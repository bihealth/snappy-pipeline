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
    These regions are essentially the accessible regions minus the target regions (with edge effect correction).

Access files can be generated from the genome reference ``fasta`` file, and optionally ``bed`` file(s) containing regions to exclude from further computations.
In this case, the user must proceed in two steps:

First, she needs to run the ``access`` tool to create the desired access file

.. code-block:: yaml

    panel_of_normals:
        tools: [access]
        access:
            exclude: <absolute path to excluded locii bed file>

This will create ``output/cnvkit.access/out/cnvkit.access.bed`` from the genomic sequence & excluded regions.

When there are no exclusion regions, the access file is automatically created using only the reference genome, and removing masked regions.

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

----------------------
Notes on the algorithm
----------------------

The choice of steps depends on the library type: whole exome sequencing is different from whole genome sequencing (panel not implemented yet).

The reference is assembled on coverage computed for all normal samples.
The coverage is always computed on target regions, and separately on antitarget regions only for exome data, not for whole genome or panel data.

For exome and panel data, target regions are obtained from the baits bed file, adding gene information & edge effects correction in the case of exome data.
For WGS data, the target regions are the full accessible regions in the genome. The user can define those accessible region (using ``access``).
But when she has left this option empty, the accessible regions are automatically defined based on the reference genome.

To create the target regions from the baits (or from the accessible regions), the target average bin size must be set.
There is a reasonable default value for exome data, but an additional ``autobin`` step is required for the whole genome data.
In ``batch`` mode, this value is computed from the coverage over the full genome
.. note::

    The ``cnvkit batch`` command also allows the creation of a flat reference, when there are no normal samples.
    This is not implemented in the ``panel_of_normals`` step, for obvious reasons.
    Using a flat reference for CNV computations is nevertheless possible, it is implemented in the ``somatic_cnv_calling`` step.

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

import re

from enum import StrEnum

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import Wildcards, InputFiles

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from snappy_pipeline.models.common import SexOrigin, SexValue

from .model import PanelOfNormals as PanelOfNormalsConfigModel
from .model import PureCn as PureCnConfig
from .model import CnvKit as CnvkitConfig

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Names of the tools that might use panel of normals
TOOLS = ("mutect2", "cnvkit", "purecn")

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = PanelOfNormalsConfigModel.default_config_yaml_string()


#: Known library types
class LibraryType(StrEnum):
    WES = "WES"
    WGS = "WGS"
    Panel = "Panel"


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
        known_libraries = self._get_normal_libraries()
        self.normal_libraries = known_libraries
        if self.name and (cfg := self.config.get(self.name)):
            if path := cfg.get("path_normals_list"):
                self.normal_libraries = {}
                with open(path, "rt") as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        library_name = line.strip()
                        assert (
                            line in known_libraries.keys()
                        ), f"Unknown requested library {library_name}"
                        self.normal_libraries[library_name] = known_libraries[library_name]
        self.libraryType, self.libraryKit = self._validate_normal_libraries()

        self.ignored = []
        if len(self.config.get("ignore_chroms", [])) > 0:
            self.ignored += self.config.ignore_chroms

    def _get_normal_libraries(self) -> dict[str, dict[str, str]]:
        normal_libraries = {}
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    if bio_sample.is_tumor:
                        continue
                    for test_sample in bio_sample.test_samples.values():
                        extraction_type = test_sample.extra_infos.get("extractionType", "DNA")
                        if extraction_type.lower() == "dna":
                            for library in test_sample.ngs_libraries.values():
                                normal_libraries[library.name] = self._get_extra_info(library)
        return normal_libraries

    def _validate_normal_libraries(self) -> tuple[str, str]:
        libraries = self.normal_libraries
        libraryType = None
        libraryKit = None
        for library in libraries:
            assert (
                libraryType is None or libraryType == libraries[library]["libraryType"]
            ), "Panel of normal cannot be built from multiple library types"
            libraryType = libraries[library]["libraryType"]
            if libraryType == LibraryType.WES:
                assert (
                    libraryKit is None or libraryKit == libraries[library]["libraryKit"]
                ), "Panel of normal cannot be built from multiple library kits"
                libraryKit = libraries[library]["libraryKit"]
        return (libraryType, libraryKit)

    @staticmethod
    def _get_extra_info(library) -> dict[str, str]:
        extra_info = {}
        assert "libraryType" in library.extra_infos, f"Undefined type of library {library.name}"
        assert (
            library.extra_infos.get("libraryType") in LibraryType
        ), f"Unknown library type {library.extra_infos.get('libraryType')}"
        extra_info["libraryType"] = library.extra_infos.get("libraryType")
        if extra_info["libraryType"] == LibraryType.WES:
            assert (
                "libraryKit" in library.extra_infos
            ), f"Undefined exome kit for library {library.name}"
            extra_info["libraryKit"] = library.extra_infos.get("libraryKit", "__default__")
        extra_info["sex"] = library.parent.parent.parent.extra_infos.get("sex", None)
        return extra_info

    @staticmethod
    @dictify
    def _get_log_file(tpl, has_sh=False):
        """Return all log files files"""
        ext_dict = {"conda_list": "conda_list.txt", "conda_info": "conda_info.txt", "log": "log"}
        if has_sh:
            ext_dict["sh"] = "sh"
        for key, ext in ext_dict.items():
            yield key, tpl + "." + ext
            yield key + "_md5", tpl + "." + ext + ".md5"


class PureCnStepPart(PanelOfNormalsStepPart):
    """Creating a panel of normals with GC-corrected coverage"""

    #: Step name
    name = "purecn"

    #: Actions
    actions = ("install", "prepare", "coverage", "create_panel")

    #: Resources
    resource_usage = {
        "install": ResourceUsage(threads=1, time="01:00:00", memory="24G"),
        "prepare": ResourceUsage(threads=1, time="04:00:00", memory="24G"),
        "coverage": ResourceUsage(threads=1, time="04:00:00", memory="24G"),
        "create_panel": ResourceUsage(threads=1, time="12:00:00", memory="32G"),
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
                self.config.purecn.path_target_interval_list_mapping[0].name,
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
                for lib in self.normal_libraries.keys()
            ],
        )

    def get_output_files(self, action):
        self._validate_action(action)
        if self.name not in self.config.tools:
            return {}

        if action == "install":
            return {"container": "work/containers/out/purecn.simg"}
        if action == "prepare":
            base_out = "{}_{}".format(
                self.config.purecn.path_target_interval_list_mapping[0].name,
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
        if self.name not in self.config.tools:
            return {}

        tpls = {
            "install": "work/containers/log/purecn",
            "prepare": "work/purecn/log/{}_{}".format(
                self.config.purecn.path_target_interval_list_mapping[0].name,
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
        "prepare_panel": ResourceUsage(threads=2, time="3-00:00:00", memory="8G"),
        "create_panel": ResourceUsage(threads=2, time="48:00:00", memory="30G"),
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
        for normal in self.normal_libraries.keys():
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
    """Build reference covergage for cnvkit"""

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "access",
        "autobin",
        "target",
        "antitarget",
        "coverage",
        "create_panel",
        "sex",
    )

    # Overwrite defaults
    default_resource_usage = ResourceUsage(threads=1, time="03:59:59", memory="7680M")  # 4h
    resource_usage = {"coverage": ResourceUsage(threads=8, time="11:59:59", memory="7680M")}

    def __init__(self, parent):
        super().__init__(parent)

        if self.name in self.config.tools:
            self.is_wgs = self.libraryType == LibraryType.WGS
            self.is_wes = self.libraryType == LibraryType.WES

            self.cfg: CnvkitConfig = self.config.get(self.name)

            self.ignored += self.cfg.ignore_chroms
            self.ignored = set(self.ignored)

            self._set_cnvkit_pipeline_logic()

            self.path_baits = self._get_path_baits()

            self.base_out = "work/{mapper}.cnvkit/out/{mapper}.cnvkit."
            self.base_out_lib = (
                "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}."
            )

    def _set_cnvkit_pipeline_logic(self):
        """
        Creates instance variables to choose path in cnvkit pipeline

        Access: regions accessible for CNV calling (unmasked)
            path_access or when missing build from genome reference + optional list of excluded region

        Target: regions of good coverage
            From baits (WES) or accessible regions (WGS) + estimate of target size from config or autobin step

        Antitarget: regions of low coverage
            antitarget = access - target, only WES, otherwise empty

        Reference:
            Flat: based on targets & antitargets only
            Cohort: from panel_of_normals step
            File: from another cohort or public data (reference + target + antitarget [WES only])
            Paired: reference built from the target & antitarget coverage of one normal sample only (paired with the tumor)
        """
        self.compute_avg_target_size = self.is_wgs and self.cfg.target.avg_size is None
        self.create_access = not self.cfg.path_access
        self.plain_access = (
            not self.cfg.path_access
            and len(self.cfg.access.exclude) == 0
            and self.cfg.access.min_gap_size is None
        )

    def _get_cohort_sex(self) -> SexValue | None:
        match self.cfg.sample_sex.source:
            case SexOrigin.CONFIG:
                return self.cfg.sample_sex.default
            case SexOrigin.AUTOMATIC:
                return None
            case SexOrigin.SAMPLESHEET:
                sex = None
                for library, extra_info in self.normal_libraries.items():
                    if extra_info.get("sex", None) is None:
                        assert sex is None, f"Sex of library {library} not defined in samplesheet"
                    else:
                        if sex is None:
                            sex = SexValue(extra_info.get("sex"))
                        else:
                            assert sex == SexValue(
                                extra_info.get("sex")
                            ), "Multiple sex in the cohort, use 'auto' in sex source"
                return sex

    def _get_path_baits(self) -> str | None:
        if not self.is_wes:
            return None
        default = None
        for item in self.cfg.path_target_interval_list_mapping:
            if item.name == self.libraryKit:
                return item.path
            elif item.name == "__default__":
                default = item.path
        if default is None:
            raise ValueError(f"Missing library kit definition for {self.libraryKit}")
        return default

    def get_input_files(self, action):
        """Return input files for cnvkit panel of normals creation"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action.replace("/", "_")))

    def get_args(self, action):
        """Return parameters input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_args_{}".format(action.replace("/", "_")))

    @dictify
    def get_output_files(self, action):
        """Return panel of normal output files"""
        self._validate_action(action)
        output_files = {}
        match action:
            case "access":
                output_files = {"access": self.base_out + "access.bed"}
            case "autobin":
                output_files = {"result": self.base_out + "autobin.txt"}
            case "target":
                output_files = {"target": self.base_out + "target.bed"}
            case "antitarget":
                output_files = {"antitarget": self.base_out + "antitarget.bed"}
            case "coverage":
                output_files = {"coverage": self.base_out_lib + "{region,(target|antitarget)}.cnn"}
            case "create_panel":
                output_files = {"reference": self.base_out + "panel_of_normals.cnn"}
            case "sex":
                output_files = {"sex": self.base_out + "sex.tsv"}

        for k, v in output_files.items():
            yield k, v
            yield k + "_md5", v + ".md5"

    @dictify
    def get_log_file(self, action):
        """Return panel of normal log files"""
        # Validate action
        self._validate_action(action)

        base_log = "work/{mapper}.cnvkit/log/{mapper}.cnvkit."
        base_log_lib = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}."
        if action in ("access", "autobin", "target", "antitarget", "create_panel", "sex"):
            tpl = base_log + action
        elif action in ("coverage",):
            tpl = base_log_lib + "{region,(target|antitarget)}"
        else:
            raise ValueError(f"Logs of action '{action}' not implemented yet")

        for key, ext in (
            ("conda_list", ".conda_list.txt"),
            ("conda_info", ".conda_info.txt"),
            ("log", ".log"),
            ("sh", ".sh"),
        ):
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    @listify
    def get_result_files(self) -> list[str]:
        if self.name not in self.config.tools:
            return []

        result_files = []

        result_files += list(self.get_output_files("create_panel").values())
        result_files += list(self.get_log_file("create_panel").values())

        result_files += list(self.get_output_files("target").values())
        result_files += list(self.get_log_file("target").values())

        if self.libraryType == LibraryType.WES:
            result_files += list(self.get_output_files("antitarget").values())
            result_files += list(self.get_log_file("antitarget").values())

        result_files += list(self.get_output_files("sex").values())
        result_files += list(self.get_log_file("sex").values())

        return filter(lambda x: not x.endswith(".md5"), result_files)

    def _get_input_files_access(self, wildcards: Wildcards) -> dict[str, str]:
        assert self.create_access, "Access shouldn't be created, already available"
        return {}

    def _get_args_access(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        assert self.create_access, "Access shouldn't be created, already available"
        return dict(input) | {
            "reference": self.w_config.static_data_config.reference.path,
            "min-gap-size": self.cfg.access.min_gap_size,
            "exclude": self.cfg.access.exclude,
            "ignore_chroms": list(self.ignored),
        }

    def _get_input_files_autobin(self, wildcards: Wildcards) -> dict[str, str]:
        assert (
            self.compute_avg_target_size
        ), "Trying to estimate average target size for non-WGS samples"
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        input_files = {
            "bams": [
                ngs_mapping(tpl.format(mapper=wildcards["mapper"], normal_library=x))
                for x in self.normal_libraries.keys()
            ]
        }
        if self.create_access:
            if self.plain_access:
                input_files["access"] = self.base_out.format(**wildcards) + "access.bed"
            else:
                input_files["target"] = self.base_out.format(**wildcards) + "access.bed"
        return input_files

    def _get_args_autobin(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        assert (
            self.compute_avg_target_size
        ), "Trying to estimate average target size for non-WGS samples"
        args = dict(input) | {"bp-per-bin": 50000}
        if self.plain_access:
            args["method"] = "wgs"
        else:
            args["method"] = "amplicon"
            if "target" not in args:
                args["target"] = self.cfg.path_access
        return args

    def _get_input_files_target(self, wildcards: Wildcards) -> dict[str, str]:
        """Helper wrapper function to estimate target average size in wgs mode"""
        input_files = {}
        if self.is_wgs:
            if self.create_access:
                input_files["interval"] = self.base_out.format(**wildcards) + "access.bed"
            if self.compute_avg_target_size:
                input_files["avg-size"] = self.base_out.format(**wildcards) + "autobin.txt"
        return input_files

    def _get_args_target(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        if self.libraryType == LibraryType.WES:
            args = {
                "avg-size": self.cfg.target.avg_size,
                "split": self.cfg.target.split,
                "interval": self.path_baits,
            }
        else:
            assert self.is_wgs, "Panel not implemented yet"
            args = dict(input) | {"split": self.cfg.target.split}
            if args.get("avg-size", None) is not None:
                args["avg-size"] = self._read_autobin_output(args["avg-size"])
            elif self.cfg.target.avg_size is not None:
                args["avg-size"] = self.cfg.target.avg_size
            else:
                args["avg-size"] = 5000
        if self.w_config.static_data_config.get("features", None):
            args["annotate"] = self.w_config.static_data_config.features.path
            args["short-names"] = self.cfg.target.short_names
        return args

    def _get_input_files_antitarget(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {"target": self.base_out.format(**wildcards) + "target.bed"}
        if self.create_access:
            input_files["access"] = self.base_out.format(**wildcards) + "access.bed"
        return input_files

    def _get_args_antitarget(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "avg-size": self.cfg.antitarget.avg_size,
            "min-size": self.cfg.antitarget.min_size,
        }
        if "access" not in args:
            args["access"] = self.cfg.path_access
        return args

    def _get_input_files_coverage(self, wildcards: Wildcards) -> dict[str, str]:
        """Helper wrapper function for computing coverage"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        return {
            "intervals": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.{region}.bed".format(
                **wildcards
            ),
            "bam": bam,
            "bai": bam + ".bai",
        }

    def _get_args_coverage(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        return dict(input) | {
            "reference": self.w_config.static_data_config.reference.path,
            "min-mapq": self.cfg.coverage.min_mapq,
            "count": self.cfg.coverage.count,
        }

    def _get_input_files_create_panel(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = self.base_out_lib + "target.cnn"
        targets = [
            tpl.format(mapper=wildcards["mapper"], library_name=x)
            for x in self.normal_libraries.keys()
        ]
        if self.libraryType == LibraryType.WES:
            tpl = self.base_out_lib + "antitarget.cnn"
            antitargets = [
                tpl.format(mapper=wildcards["mapper"], library_name=x)
                for x in self.normal_libraries.keys()
            ]
        else:
            antitargets = []
        return {"normals": targets + antitargets}

    def _get_args_create_panel(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(input) | {
            "reference": self.w_config.static_data_config.reference.path,
            "cluster": self.cfg.cluster,
            "no-gc": not self.cfg.gc,
            "no-rmask": not self.cfg.rmask,
            "no-edge": not self.cfg.get("edge", self.is_wes),
            "diploid-parx-genome": self.cfg.diploid_parx_genome,
        }
        if self.cfg.cluster:
            args["min-cluster-size"] = self.cfg.min_cluster_size
        sample_sex = self._get_cohort_sex()
        if sample_sex is not None:
            args["sample-sex"] = str(sample_sex)
        args["male-reference"] = self.cfg.male_reference
        return args

    def _get_input_files_sex(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = self.base_out_lib + "target.cnn"
        coverages = [
            tpl.format(mapper=wildcards["mapper"], library_name=x)
            for x in self.normal_libraries.keys()
        ]
        if self.is_wes:
            tpl = self.base_out_lib + "antitarget.cnn"
            coverages += [
                tpl.format(mapper=wildcards["mapper"], library_name=x)
                for x in self.normal_libraries.keys()
            ]
        return {"coverages": coverages}

    def _get_args_sex(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        return dict(input) | {"diploid-parx-genome": self.cfg.diploid_parx_genome}

    def _read_autobin_output(self, filename: str) -> int:
        nb = r"([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)"
        pattern = re.compile("^Target:[ \t]+" + nb + "[ \t]+" + nb + "$")
        with open(filename) as f:
            for line in f:
                m = pattern.match(line)
                if m:
                    return int(float(m.groups()[4]))
        return -1


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

        log_ext_list = ["log", "conda_list.txt", "conda_info.txt"]

        if "mutect2" in set(self.config.tools) & set(TOOLS):
            tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals.{ext}"
            ext_list = ("vcf.gz", "vcf.gz.tbi")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB.{ext}"
            ext_list = ("tar.gz",)
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if "cnvkit" in set(self.config.tools) & set(TOOLS):
            cnvkit_files = self.sub_steps["cnvkit"].get_result_files()
            for work in cnvkit_files:
                output = work.replace("work/", "output/", 1)
                result_files.extend(self._expand_result_files(output, ("",)))

        if "purecn" in set(self.config.tools) & set(TOOLS):
            tpl = "output/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.{ext}"
            ext_list = ("rds",)
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.{ext}"
            ext_list = ("rds",)
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))
            tpl = "output/purecn/out/{}_{}.{{ext}}".format(
                # TODO: select enrichment kit
                self.config.purecn.path_target_interval_list_mapping[0].name,
                self.config.purecn.genome_name,
            )
            ext_list = ("list", "bed.gz", "bed.gz.tbi")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/purecn/log/{}_{}.{{ext}}".format(
                self.config.purecn.path_target_interval_list_mapping[0].name,
                self.config.purecn.genome_name,
            )
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        return result_files

    def _expand_result_files(self, tpl, ext_list):
        for mapper in self.w_config.step_config["ngs_mapping"].tools.dna:
            for ext in ext_list:
                yield tpl.format(mapper=mapper, ext=ext)
                yield tpl.format(mapper=mapper, ext=ext) + ".md5"
