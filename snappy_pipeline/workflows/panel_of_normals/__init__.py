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

================================
Notes on the ``cnvkit`` workflow
================================

``cnvkit`` is a set of tools originally designed to call somatic copy number alterations from exome data.
Its design is modular, which enables its use for whole genome and amplicon data.

``cnvkit`` provides a tool to encapsulate common practice workflows (``batch``), depending on the type of data, and on the availability of optional inputs.
The current implementation recapitulates the common practice, while still dispaching computations on multiple cluster nodes.

For exome and whole genome data, the ``cnvkit`` `documentation <https://cnvkit.readthedocs.io/en/stable/>`_
recommends the creation of a panel of normal (called ``reference``).
The actual workflow to generate this reference is slightly different between exome and whole genome data,
and also changes depending whether an accessibility file is provided by the user or not.

Therefore, the ``cnvkit`` tool to generate such accessibility file is implemented as a separate tool.
If a user wants to create this accessibility file with ``cnvkit`` tools, then she must first run the ``access`` tool.
Only after it has been created can she use it to generate the panel of normals.
For that, she will need to modify the configuration file, adding ``cnvkit`` in the list of tools, and setting the ``access`` parameter to the output of the ``access`` tool.

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

``access`` is used to create a genome accessibility file that can be used for ``cnvkit`` panel of normals creation.
Its output (``output/cnvkit.access/out/cnvkit.access.bed``) is optional, but its presence impacts of the way the target and antitarget regions are computed in whole genome mode.

In a nutshell, for exome data, the accessibility file is only used to create antitarget regions.
For genome data, it is used by the ``autobin`` tool to compute the average target size used during target regions creation.
If it is present, the target size is computed in amplicon mode, and when it is absent,
an accessibility file is created with default settings, which value is used by ``autobin`` is whole genome mode.

This follows the internal ``batch`` code of ``cnvkit``.

=======
Reports
=======

Report tables can be found in the ``output/{mapper}.cnvkit/report`` directory.
Two tables are produced, grouping results for all normal samples together:

- ``metrics.txt``: coverage metrics over target and antitarget regions.
- ``sex.txt``: prediction of the donor's gender based on the coverage of chromosome X & Y target and antitarget regions.

The cnvkit authors recommend to check these reports to ensure that all data is suitable for panel of normal creation.
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

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Names of the tools that might use panel of normals
TOOLS = ("mutect2", "cnvkit", "access", "purecn")

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = r"""
# Default configuration somatic_variant_calling
step_config:
  panel_of_normals:
    tools: ['mutect2']  # REQUIRED - available: 'mutect2'
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    ignore_chroms: []   # patterns of chromosome names to ignore
                        # hs37d5: [NC_007605, hs37d5, chrEBV, '*_decoy', 'HLA-*', 'GL000220.*']
                        # GRCh38.d1.vd1: [chrEBV, 'HPV*', CMV, HBV, 'HCV-*', 'HIV-*', KSHV, 'HTLV-1', MCV, '*_decoy', 'chrUn_GL00220*', SV40]
    # Configuration for mutect2
    mutect2:
      path_normals_list: null    # Optional file listing libraries to include in panel
      germline_resource: REQUIRED
      # Java options
      java_options: ' -Xmx16g '
      # Parallelization configuration
      num_cores: 2               # number of cores to use locally
      window_length: 100000000   # split input into windows of this size, each triggers a job
      num_jobs: 500              # number of windows to process in parallel
      use_profile: true          # use Snakemake profile for parallel processing
      restart_times: 5           # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2     # throttling of job creation
      max_status_checks_per_second: 10 # throttling of status checks
      debug_trunc_tokens: 0      # truncation to first N tokens (0 for none)
      keep_tmpdir: never         # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1         # memory multiplier
      job_mult_time: 1           # running time multiplier
      merge_mult_memory: 1       # memory multiplier for merging
      merge_mult_time: 1         # running time multiplier for merging
    cnvkit:
      path_normals_list: ""       # Optional file listing libraries to include in panel
      path_target_regions: ""     # Bed files of targetted regions (Missing when creating a panel of normals for WGS data)
      access: ""                  # Access bed file (output/cnvkit.access/out/cnvkit.access.bed when create_cnvkit_acces was run)
      annotate: ""                # [target] Optional targets annotations
      target_avg_size: 0          # [target] Average size of split target bins (0: use default value)
      bp_per_bin: 50000           # [autobin] Expected base per bin
      split: True                 # [target] Split large intervals into smaller ones
      antitarget_avg_size: 0      # [antitarget] Average size of antitarget bins (0: use default value)
      min_size: 0                 # [antitarget] Min size of antitarget bins (0: use default value)
      min_mapq: 0                 # [coverage] Mininum mapping quality score to count a read for coverage depth
      count: False                # [coverage] Alternative couting algorithm
      min_cluster_size: 0         # [reference] Minimum cluster size to keep in reference profiles. 0 for no clustering
      gender: ""                  # [reference] Specify the chromosomal sex of all given samples as male or female. Guess when missing
      male_reference: False       # [reference & sex] Create male reference
      gc_correction: True         # [reference] Use GC correction
      edge_correction: True       # [reference] Use edge correction
      rmask_correction: True      # [reference] Use rmask correction
      drop_low_coverage: False    # [metrics] Drop very-low-coverage bins before calculations
    access:                       # Creates access file for cnvkit, based on genomic sequence & excluded regions (optionally)
      exclude: []                 # [access] Bed file of regions to exclude (mappability, blacklisted, ...)
      min_gap_size: 0             # [access] Minimum gap size between accessible sequence regions (0: use default value)
    purecn:
      path_normals_list: ""       # Optional file listing libraries to include in panel
      path_bait_regions: REQUIRED # Bed files of enrichment kit sequences (MergedProbes for Agilent SureSelect), recommended by PureCN author
      path_genomicsDB: REQUIRED   # Mutect2 genomicsDB created during panel_of_normals
      genome_name: "unknown"      # Must be one from hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6, canFam3
      enrichment_kit_name: "unknown" # For filename only...
      mappability: "" # GRCh38: /fast/work/groups/cubi/projects/biotools/static_data/app_support/PureCN/hg38/mappability.bw
      reptiming: ""   # Nothing for GRCh38
      seed: 1234567
"""


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
        if self.name and self.config[self.name].get("path_normals_list"):
            self.normal_libraries = []
            with open(self.config[self.name]["path_normals_list"], "rt") as f:
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
            time="04:00:00",  # 4 hours
            memory="24G",
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
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        yield "intervals", "work/purecn/out/{}_{}.list".format(
            self.config["purecn"]["enrichment_kit_name"],
            self.config["purecn"]["genome_name"],
        )
        yield "bam", self.ngs_mapping(tpl.format(**wildcards))

    @dictify
    def _get_input_files_create(self, wildcards):
        yield "container", "work/containers/out/purecn.simg"
        tpl = "work/{mapper}.purecn/out/{mapper}.{normal_library}_coverage_loess.txt.gz"
        yield "normals", [
            tpl.format(mapper=wildcards.mapper, normal_library=lib) for lib in self.normal_libraries
        ]

    def get_output_files(self, action):
        self._validate_action(action)

        if action == "install":
            return {"container": "work/containers/out/purecn.simg"}
        if action == "prepare":
            base_out = "{}_{}".format(
                self.config["purecn"]["enrichment_kit_name"],
                self.config["purecn"]["genome_name"],
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
                "coverage": "work/{mapper}.purecn/out/{mapper}.{normal_library,.+-DNA[0-9]+-WES[0-9]+}_coverage_loess.txt.gz"
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
                self.config["purecn"]["enrichment_kit_name"],
                self.config["purecn"]["genome_name"],
            ),
            "coverage": "work/{mapper}.purecn/log/{mapper}.{normal_library,.+-DNA[0-9]+-WES[0-9]+}",
            "create_panel": "work/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals",
        }
        assert action in self.actions
        return self._get_log_file(tpls[action])

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # PureCN not enabled, skip
        self.parent.ensure_w_config(
            ("step_config", "panel_of_normals", self.name, "path_bait_regions"),
            "Path to exome panel bait regions not defined for tool {}".format(self.name),
        )


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
        if self.name not in self.config["tools"]:
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
        self.is_wgs = self.config["cnvkit"]["path_target_regions"] == ""

    def check_config(self):
        if self.name not in self.config["tools"]:
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
        return {
            "target": targets
            if targets
            else "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed".format(**wildcards),
            "antitarget": antitargets
            if antitargets
            else "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed".format(**wildcards),
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

    def get_resource_usage(self, action):
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
            (NgsMappingWorkflow,),
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
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

        if "mutect2" in set(self.config["tools"]) & set(TOOLS):
            tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals.{ext}"
            ext_list = ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB.{ext}"
            ext_list = ("tar.gz", "tar.gz.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if "cnvkit" in set(self.config["tools"]) & set(TOOLS):
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

        if "access" in set(self.config["tools"]) & set(TOOLS):
            tpl = "output/cnvkit.access/out/cnvkit.access.bed"
            result_files.extend([tpl + md5 for md5 in ("", ".md5")])
            tpl = "output/cnvkit.access/log/cnvkit.access.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        if "purecn" in set(self.config["tools"]) & set(TOOLS):
            tpl = "output/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals.{ext}"
            result_files.extend(self._expand_result_files(tpl, log_ext_list))
            tpl = "output/purecn/out/{}_{}.{{ext}}".format(
                self.config["purecn"]["enrichment_kit_name"],
                self.config["purecn"]["genome_name"],
            )
            ext_list = ("list", "list.md5", "bed.gz", "bed.gz.md5", "bed.gz.tbi", "bed.gz.tbi.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/purecn/log/{}_{}.{{ext}}".format(
                self.config["purecn"]["enrichment_kit_name"],
                self.config["purecn"]["genome_name"],
            )
            result_files.extend(self._expand_result_files(tpl, log_ext_list))

        return result_files

    def _expand_result_files(self, tpl, ext_list):
        for mapper in self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"]:
            for ext in ext_list:
                yield tpl.format(mapper=mapper, ext=ext)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "panel_of_normals", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for somatic variant calling",
        )
