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
For example, the ``mutect2`` panel of normal generates ``mutect2.pon.vcf.gz``
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
        path_access: <absolute path to access file>                   # Even when created by the ``access`` tool.
        path_target: <absolute path to baits>                         # Keep empty for WGS data
        path_normals_list: <absolute path to list of normal samples>  # Keep empty to use all available normals

Note that there is no provision (yet) to automatically create separate panel of normals for males & females.
If the number of samples collected in the same fashion is large enough, it is nevertheless the way to achieve best results.

-------
Reports
-------

Report tables can be found in the ``output/cnvkit/report`` directory.
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

from snappy_pipeline.models.cnvkit import Gender as CnvKitGender
from snappy_pipeline.models.cnvkit import PanelOfNormals as CnvKitModel
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import PanelOfNormals as PanelOfNormalsConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

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
            runtime="1h",
            mem="24GB",
        ),
        "prepare": ResourceUsage(
            threads=1,
            runtime="4h",  # 4 hours
            mem="24GB",
        ),
        "coverage": ResourceUsage(
            threads=1,
            runtime="4h",  # 4 hours
            mem="24GB",
        ),
        "create_panel": ResourceUsage(
            threads=1,
            runtime="12h",  # 12 hours
            mem="32GB",
        ),
    }

    def get_input_files(self, action):
        if self.name != self.config.tool:
            return {}
        self._validate_action(action)
        self.ngs_mapping = self.parent.modules["ngs_mapping"]
        if action == "prepare":
            return {
                "container": "work/containers/out/purecn.simg",
                "reference": self.w_config.static_data_config.reference.path,
            }
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
        tpl = "output/{library_name}/out/{library_name}.bam"
        yield "bam", self.ngs_mapping(tpl.format(**wildcards))

    @dictify
    def _get_input_files_create(self, wildcards):
        yield "container", "work/containers/out/purecn.simg"
        tpl = "work/purecn/out/{library_name}_coverage_loess.txt.gz"
        yield "normals", [tpl.format(library_name=lib) for lib in self.normal_libraries]
        # The Mutect2 genomicsDB is the output of the upstream panel_of_normals (mutect2) task;
        # resolve it through the registered module so Snakemake tracks it as a real dependency.
        pon_module = self.parent.modules["panel_of_normals"]
        yield "genomicsdb", pon_module("work/mutect2/out/mutect2.genomicsDB.tar.gz")

    def get_output_files(self, action):
        if self.name != self.config.tool:
            return {}
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
            return {"coverage": "work/purecn/out/{library_name}_coverage_loess.txt.gz"}
        if action == "create_panel":
            return {
                "db": "work/purecn/out/purecn.panel_of_normals.rds",
                "db_md5": "work/purecn/out/purecn.panel_of_normals.rds.md5",
                "mapbias": "work/purecn/out/purecn.mapping_bias.rds",
                "mapbias_md5": "work/purecn/out/purecn.mapping_bias.rds.md5",
                "lowcov": "work/purecn/out/purecn.low_coverage_targets.bed",
                "hq": "work/purecn/out/purecn.hq_sites.bed",
                "plot": "work/purecn/out/purecn.interval_weights.png",
            }

    def get_args(self, action):
        self._validate_action(action)
        if action == "coverage":
            return getattr(self, f"_get_args_{action}")
        else:
            return {"config": self.config.get(self.name).model_dump(by_alias=True)}

    def _get_args_coverage(self, wildcards):
        mapper = str(self.parent.get_task_config("ngs_mapping").tool)
        return {
            "config": self.config.get(self.name).model_dump(by_alias=True),
            "mapper": mapper,
            "library_name": wildcards.library_name,
        }

    def get_log_file(self, action):
        if self.name != self.config.tool:
            return {}
        tpls = {
            "install": "work/containers/log/purecn",
            "prepare": "work/purecn/log/{}_{}".format(
                self.config.purecn.enrichment_kit_name,
                self.config.purecn.genome_name,
            ),
            "coverage": "work/purecn/log/{library_name}",
            "create_panel": "work/purecn/log/purecn.panel_of_normals",
        }
        assert action in self.actions
        return self._get_log_file(tpls[action])


class Mutect2StepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    #: Step name
    name = "mutect2"

    #: Class available actions
    actions = ("scatter", "prepare_panel", "gather", "create_panel")

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage = {
        "scatter": ResourceUsage(
            threads=1,
            runtime="2m",
            mem="1000MB",
        ),
        "prepare_panel": ResourceUsage(
            threads=2,
            runtime="3d",  # 3 days
            mem="8GB",
        ),
        "gather": ResourceUsage(
            threads=1,
            runtime="4h",
            mem="8192MB",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            runtime="48h",  # 48 hours
            mem="30GB",
        ),
    }

    def get_input_files(self, action):
        """Return input files for mutect2 variant calling"""
        # Validate action
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def _get_input_files_scatter(self, wildcards):
        return {"fai": self.w_config.static_data_config.reference.path + ".fai"}

    def _get_input_files_prepare_panel(self, wildcards):
        """Helper wrapper function for single sample panel preparation"""
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.modules["ngs_mapping"]
        tpl = "output/{normal_library}/out/{normal_library}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        scatteritem_base_path = "work/{normal_library}/par/scatter/{scatteritem}.region.bed"
        return {
            "normal_bam": bam,
            "normal_bai": bam + ".bai",
            "region": scatteritem_base_path.format(**wildcards),
            "reference": self.w_config.static_data_config.reference.path,
        }

    def _get_input_files_gather(self, wildcards):
        gather = self.parent.workflow.globals.get("gather")
        gather = getattr(gather, self.name)
        tpl = "work/{normal_library}/par/run/{{scatteritem}}.vcf.gz".format(**wildcards)
        return {"vcf": gather(tpl)}

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for merging individual results & panel creation"""
        paths = []
        tpl = "work/{normal_library}/out/{normal_library}.prepare.vcf.gz"
        for normal in self.normal_libraries:
            paths.append(tpl.format(normal_library=normal, **wildcards))
        return {
            "normals": paths,
            "reference": self.w_config.static_data_config.reference.path,
            "germline_resource": self.config.mutect2.germline_resource,
        }

    def get_output_files(self, action):
        """Return panel of normal files"""
        self._validate_action(action)

        if action == "scatter":
            scatter = self.parent.workflow.globals.get("scatter")
            scatter = getattr(scatter, self.name)
            tpl = "work/{{normal_library}}/par/scatter/{scatteritem}.region.bed"
            return {"regions": scatter(tpl)}

        ext_dict = {
            "vcf": "vcf.gz",
            "vcf_md5": "vcf.gz.md5",
            "vcf_tbi": "vcf.gz.tbi",
            "vcf_tbi_md5": "vcf.gz.tbi.md5",
        }

        tpls = {
            "prepare_panel": "work/{normal_library}/par/run/{scatteritem}",
            "gather": "work/{normal_library}/out/{normal_library}.prepare",
            "create_panel": "work/mutect2/out/mutect2.panel_of_normals",
        }
        output_files = {}
        for key, ext in ext_dict.items():
            output_files[key] = tpls[action] + "." + ext
        if action == "create_panel":
            output_files["db"] = "work/mutect2/out/mutect2.genomicsDB.tar.gz"
            output_files["db_md5"] = "work/mutect2/out/mutect2.genomicsDB.tar.gz.md5"
        return output_files

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_scatter(self, wildcards):
        return {
            "ignore_chroms": self.config.ignore_chroms,
            "padding": self.config.mutect2.padding,
        }

    def _get_args_prepare_panel(self, wildcards):
        return {
            "max_mnp_distance": 0,
            "java_options": self.config.mutect2.java_options,
            "extra_arguments": self.config.mutect2.extra_arguments,
        }

    def _get_args_gather(self, wildcards):
        return {}

    def _get_args_create_panel(self, wildcards):
        return self.config.mutect2.genomicsdb.model_dump(by_alias=True)

    def get_log_file(self, action):
        """Get log files for Mutect2 rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns dictionary with expected log files based on inputted action.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)

        # Set expected format based on action
        tpl = "{normal_library}"
        match action:
            case "gather":
                postfix = ""
            case "prepare_panel":
                postfix = ".{scatteritem}"
            case "create_panel":
                tpl = self.name
                postfix = ".panel_of_normals"
            case _:
                postfix = "." + action

        prefix = f"work/{tpl}/log/{tpl}{postfix}"

        return self._get_log_file(prefix)


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
            runtime="2h",  # 2 hours
            mem="8GB",
        ),
        "antitarget": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="8GB",
        ),
        "coverage": ResourceUsage(
            threads=8,
            runtime="2h",  # 2 hours
            mem="16GB",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="16GB",
        ),
        "report": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="16GB",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        if self.name == self.config.tool:
            self.is_wgs = self.config.cnvkit.path_target == ""

    def check_config(self):
        if self.name != self.config.tool:
            return None  # cnvkit not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_args(self, action):
        if self.name != self.config.tool:
            return None  # cnvkit not enabled, skip
        self._validate_action(action)
        cfg: CnvKitModel = self.config.get(self.name)
        if action == "create_panel":
            action = "reference"
        if args := getattr(cfg, action, {}):
            args = args.model_dump(by_alias=True)
        if action == "reference":
            if cfg.gender != CnvKitGender.guess:
                args["gender"] = cfg.gender
            if cfg.male_reference:
                args["male_reference"] = True
            args["flat"] = len(self.normal_libraries) == 0
            if not cfg.path_target:
                args["edge_correction"] = False
        if action == "target":
            args["bp_per_bin"] = cfg.bp_per_bin
        return args

    def get_input_files(self, action):
        """Return input files for cnvkit panel of normals creation"""
        if self.name != self.config.tool:
            return None  # cnvkit not enabled, skip
        # Validate action
        self._validate_action(action)
        mapping = {
            "target": self._get_input_files_target,
            "antitarget": self._get_input_files_antitarget,
            "coverage": self._get_input_files_coverage,
            "create_panel": self._get_input_files_create_panel,
            "report": self._get_input_files_report,
            "access": self._get_input_files_access,
        }
        return mapping[action]

    def _get_input_files_access(self, wildcards):
        return {"reference": self.w_config.static_data_config.reference.path}

    def _get_input_files_target(self, wildcards):
        """Helper wrapper function to estimate target average size in wgs mode"""
        if not self.is_wgs:
            input_files = {"target": self.config.cnvkit.path_target}
            if self.config.cnvkit.path_annotation:
                input_files["annotate"] = self.config.cnvkit.path_annotation
            return input_files
        ngs_mapping = self.parent.modules["ngs_mapping"]
        tpl = "output/{normal_library}/out/{normal_library}.bam"
        bams = [ngs_mapping(tpl.format(normal_library=x)) for x in self.normal_libraries]
        bais = [x + ".bai" for x in bams]
        input_files = {
            "bams": bams,
            "bais": bais,
            "reference": self.w_config.static_data_config.reference.path,
        }
        if self.config.cnvkit.path_access:
            input_files["access"] = self.config.cnvkit.path_access
        if self.config.cnvkit.path_annotation:
            input_files["annotate"] = self.config.cnvkit.path_annotation
        return input_files

    def _get_input_files_antitarget(self, wildcards):
        """Helper wrapper function for computing antitarget locations"""
        if self.is_wgs:
            return {}
        result = {
            "target": "work/cnvkit/out/cnvkit.target.bed".format(**wildcards),
        }
        if path_access := self.config.cnvkit.path_access:
            result["access"] = path_access
        return result

    def _get_input_files_coverage(self, wildcards):
        """Helper wrapper function for computing coverage"""
        ngs_mapping = self.parent.modules["ngs_mapping"]
        tpl = "output/{normal_library}/out/{normal_library}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        return {
            "target": "work/cnvkit/out/cnvkit.target.bed".format(**wildcards),
            "antitarget": "work/cnvkit/out/cnvkit.antitarget.bed".format(**wildcards),
            "bam": bam,
            "bai": bam + ".bai",
            "reference": self.w_config.static_data_config.reference.path,
        }

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for computing panel of normals"""
        tpl = "work/cnvkit/out/{normal_library}.targetcoverage.cnn"
        targets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        tpl = "work/cnvkit/out/{normal_library}.antitargetcoverage.cnn"
        antitargets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        tpl = "work/cnvkit/log/{normal_library}.coverage.{ext}"
        logs = [
            tpl.format(normal_library=x, ext=ext)
            for x in self.normal_libraries
            for ext in ("log", "conda_list.txt", "conda_info.txt")
        ]
        return {
            "target": (
                targets if targets else "work/cnvkit/out/cnvkit.target.bed".format(**wildcards)
            ),
            "antitarget": (
                antitargets
                if antitargets
                else "work/cnvkit/out/cnvkit.antitarget.bed".format(**wildcards)
            ),
            "logs": logs if targets or antitargets else [],
            "reference": self.w_config.static_data_config.reference.path,
        }

    def _get_input_files_report(self, wildcards):
        """Helper wrapper function for the panel of normals report"""
        tpl = "work/cnvkit/out/{normal_library}.targetcoverage.cnn"
        targets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        tpl = "work/cnvkit/out/{normal_library}.antitargetcoverage.cnn"
        antitargets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        return {
            "target": targets,
            "antitarget": antitargets,
        }

    def get_output_files(self, action):
        """Return panel of normal files"""
        if self.name != self.config.tool:
            return {}  # cnvkit not enabled, skip
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
            "target": "work/cnvkit/out/cnvkit.target.bed",
            "target_md5": "work/cnvkit/out/cnvkit.target.bed.md5",
        }

    def _get_output_files_antitarget(self):
        return {
            "antitarget": "work/cnvkit/out/cnvkit.antitarget.bed",
            "antitarget_md5": "work/cnvkit/out/cnvkit.antitarget.bed.md5",
        }

    def _get_output_files_coverage(self):
        return {
            "target": "work/cnvkit/out/{normal_library}.targetcoverage.cnn",
            "target_md5": "work/cnvkit/out/{normal_library}.targetcoverage.cnn.md5",
            "antitarget": "work/cnvkit/out/{normal_library}.antitargetcoverage.cnn",
            "antitarget_md5": "work/cnvkit/out/{normal_library}.antitargetcoverage.cnn.md5",
        }

    def _get_output_files_create_panel(self):
        return {
            "panel": "work/cnvkit/out/cnvkit.panel_of_normals.cnn",
            "panel_md5": "work/cnvkit/out/cnvkit.panel_of_normals.cnn.md5",
            "log": "work/cnvkit/log/cnvkit.merged.tar.gz",
            "log_md5": "work/cnvkit/log/cnvkit.merged.tar.gz.md5",
        }

    def _get_output_files_report(self):
        return {
            "sex": "work/cnvkit/report/cnvkit.sex.tsv",
            "sex_md5": "work/cnvkit/report/cnvkit.sex.tsv.md5",
            "metrics": "work/cnvkit/report/cnvkit.metrics.tsv",
            "metrics_md5": "work/cnvkit/report/cnvkit.metrics.tsv.md5",
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
            "target": "work/cnvkit/log/cnvkit.target",
            "antitarget": "work/cnvkit/log/cnvkit.antitarget",
            "coverage": "work/cnvkit/log/{normal_library}.coverage",
            "create_panel": "work/cnvkit/log/cnvkit.panel_of_normals",
            "report": "work/cnvkit/log/cnvkit.report",
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
            runtime="2h",  # 2 hours
            mem="8GB",
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
    consumes = {DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})): True}
    produces = [DataSignature(DataType.MODELS, frozenset({"pon"}))]

    config_model_class = PanelOfNormalsConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str | None = None,
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=(NgsMappingWorkflow,),
            task_name=task_name,
            **kwargs,
        )
        # Initialize sub-workflows
        self.register_module("ngs_mapping")
        # When tool=purecn, the genomicsDB produced by an upstream mutect2 PON task is a
        # tracked Snakemake input; register it so paths are resolved relative to that task.
        if self.config.tool == "purecn":
            self.register_module("panel_of_normals")
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

        if self.config.tool == "mutect2":
            tpl = "output/mutect2/out/mutect2.panel_of_normals.{ext}"
            ext_list = ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/mutect2/out/mutect2.genomicsDB.{ext}"
            ext_list = ("tar.gz", "tar.gz.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/mutect2/log/mutect2.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if self.config.tool == "cnvkit":
            tpls = [
                ("output/cnvkit/out/cnvkit.target.{ext}", ("bed", "bed.md5")),
                ("output/cnvkit/out/cnvkit.antitarget.{ext}", ("bed", "bed.md5")),
                (
                    "output/cnvkit/out/cnvkit.panel_of_normals.{ext}",
                    ("cnn", "cnn.md5"),
                ),
                (
                    "output/cnvkit/report/cnvkit.sex.{ext}",
                    ("tsv", "tsv.md5"),
                ),
                (
                    "output/cnvkit/report/cnvkit.metrics.{ext}",
                    ("tsv", "tsv.md5"),
                ),
            ]
            for tpl, ext_list in tpls:
                result_files.extend(self._expand_result_files(tpl, ext_list))
            tpls = [
                "output/cnvkit/log/cnvkit.target.{ext}",
                "output/cnvkit/log/cnvkit.antitarget.{ext}",
                "output/cnvkit/log/cnvkit.panel_of_normals.{ext}",
                "output/cnvkit/log/cnvkit.report.{ext}",
            ]
            for tpl in tpls:
                result_files.extend(self._expand_result_files(tpl, log_ext_list))
            tpl = "output/cnvkit/log/cnvkit.merged.tar.gz{ext}"
            result_files.extend(self._expand_result_files(tpl, ("", ".md5")))

        if self.config.tool == "access":
            tpl = "output/cnvkit.access/out/cnvkit.access.bed"
            result_files.extend([tpl + md5 for md5 in ("", ".md5")])
            tpl = "output/cnvkit.access/log/cnvkit.access.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if self.config.tool == "purecn":
            tpl = "output/purecn/out/purecn.panel_of_normals.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/purecn/out/purecn.mapping_bias.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/purecn/log/purecn.panel_of_normals.{ext}"
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
        for ext in ext_list:
            yield tpl.format(ext=ext)
