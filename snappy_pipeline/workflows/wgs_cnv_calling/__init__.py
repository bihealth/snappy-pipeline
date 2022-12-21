# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_cnv_calling`` step

The ``wgs_cnv_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned germline reads) and performs CNV calling on them.

The result are called CNVs in VCF format for each donor's primary NGS library.  The only supported
CNV caller is Illumina Canvas that does not support integrating CNV results from multiple
samples.

In the future, CNV callers supporting pedigrees or merging the CNV caller output should be
supported.

In contrats to WGS SV calling, CNV calling only considers copy number variants and no inversions.
However, large-range CNVs are more reliably detected with read depth signal considered by CNV
calling.

==========
Stability
==========

This step is considered experimental, use it at your own discretion.

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` (for the aligned reads BAM files).

===========
Step Output
===========

For all donors, CNV calling will be performed on the primary DNA NGS library, separately
for each configured read mapper and CNV calling tool.  The name of this primary DNA NGS library
will be used for the identification token in the output file.  For each donor, read mapper,
and CNV calling tool, the following files will be generated.

- ``{mapper}.{cnv_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{cnv_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{cnv_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{cnv_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.cnvkit.P001-N1-DNA1-WGS1-4
    |   `-- out
    |       |-- bwa.cnvkit.P001-N1-DNA1-WGS1-4.vcf.gz
    |       |-- bwa.cnvkit.P001-N1-DNA1-WGS1-4.vcf.gz.tbi
    |       |-- bwa.cnvkit.P001-N1-DNA1-WGS1-4.vcf.gz.md5
    |       `-- bwa.cnvkit.P001-N1-DNA1-WGS1-4.vcf.gz.tbi.md5
    [...]

====================
Global Configuration
====================

- At the moment, no global configuration is used

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_wgs_cnv_calling.rst

=====================
Available CNV Callers
=====================

The following CNV callers are currently available

- ``"canvas"``
- ``"Delly2"``
- ``"gCNV"``

=======
Reports
=======

Currently, no reports are generated.
"""

from collections import OrderedDict
import os

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
from snappy_pipeline.workflows.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload (VCF)
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Name/extension pairs for VCF files.
VCF_KEY_EXTS = dict(zip(EXT_NAMES, EXT_VALUES))

#: Name/extension pairs for BCF files.
BCF_KEY_EXTS = {"bcf": ".bcf", "bcf_md5": ".bcf.md5", "csi": ".bcf.csi", "csi_md5": ".bcf.csi.md5"}

#: Available WGS CNV callers
WGS_CNV_CALLERS = ("cnvetti", "delly2", "gcnv")

#: Default configuration for the wgs_cnv_calling schema
DEFAULT_CONFIG = r"""
# Default configuration wgs_cnv_calling
step_config:
  wgs_cnv_calling:
    path_ngs_mapping: ../ngs_mapping          # REQUIRED
    tools: [delly2, gcnv]                     # REQUIRED

    cnvetti:
      window_length: null  # Optional
      count_kind: null     # Optional
      segmentation: null   # Optional
      normalization: null  # Optional
      preset: deep_wgs
      presets:
        shallow_wgs:
          window_length: 20000
          count_kind: Fragments
          segmentation: HaarSeg
          normalization: MedianGcBinned
        deep_wgs:
          window_length: 200
          count_kind: Coverage
          segmentation: HaarSeg
          normalization: MedianGcBinned

    delly2:
      mappability: null  # REQUIRED - Path to mappability map
      minsize: 1000  # Merge min size - set to Delly default
      maxsize: 100000  # Merge max size - set to Delly default

    gcnv:
      # Path to gCNV model - will execute analysis in CASE MODE.
      #
      # Example of precomputed model:
      # - library: "Agilent SureSelect Human All Exon V6"  # Library name
      #   contig_ploidy: /path/to/ploidy-model         # Output from `DetermineGermlineContigPloidy`
      #   model_pattern: /path/to/model_*              # Output from `GermlineCNVCaller`
      precomputed_model_paths: []  # REQUIRED

      # Path to BED file with uniquely mappable regions.
      path_uniquely_mapable_bed: null  # REQUIRED
"""


class CnvettiStepPart(BaseStepPart):
    """Shallow WGS CNV calling with CNVetti."""

    #: Step name
    name = "cnvetti"

    #: Supported actions
    actions = (
        "coverage",
        "segment",
        "merge_segments",
        "genotype",
        "merge_genotypes",
        "reorder_vcf",
    )

    #: Directory infixes
    dir_infixes = {
        "coverage": "{mapper}.cnvetti_coverage.{library_name}",
        "segment": "{mapper}.cnvetti_segment.{library_name}",
        "merge_segments": "{mapper}.cnvetti_merge_segments",
        "genotype": "{mapper}.cnvetti_genotype.{library_name}",
        "merge_genotypes": "{mapper}.cnvetti_merge_genotypes",
        "reorder_vcf": "{mapper}.cnvetti_reorder_vcf.{library_name}",
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{index_ngs_library}}/out/"
            "{{mapper}}.{var_caller}.{{index_ngs_library}}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return input function for CNVetti action."""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_coverage(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Yield input BAM and BAI file
        bam_tpl = "output/{mapper}.{ngs_library}/out/{mapper}.{ngs_library}{ext}"
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield key, ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_input_files_segment(self, wildcards):
        name_pattern = "{mapper}.cnvetti_coverage.{ngs_library}".format(**wildcards)
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    def _get_input_files_merge_segments(self, wildcards):
        result = {}
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            ngs_libraries = [
                donor.dna_ngs_library
                for pedigree in sheet.cohort.pedigrees
                for donor in pedigree.donors
                if donor.dna_ngs_library
            ]
            for ngs_library in ngs_libraries:
                name_pattern = "{mapper}.cnvetti_segment.{ngs_library}".format(
                    ngs_library=ngs_library.name, **wildcards
                )
                for key, ext in BCF_KEY_EXTS.items():
                    result.setdefault(key, []).append(
                        os.path.join("work", name_pattern, "out", name_pattern + ".segments" + ext)
                    )
        return result

    @dictify
    def _get_input_files_genotype(self, wildcards):
        # The sites list BCF file.
        name_pattern = "{mapper}.cnvetti_merge_segments".format(**wildcards)
        for key, ext in BCF_KEY_EXTS.items():
            yield "sites_" + key, os.path.join("work", name_pattern, "out", name_pattern + ext)
        # The original coverage BCF file.
        name_pattern = "{mapper}.cnvetti_coverage.{ngs_library}".format(**wildcards)
        for key, ext in BCF_KEY_EXTS.items():
            yield "coverage_" + key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    def _get_input_files_merge_genotypes(self, wildcards):
        result = {}
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            ngs_libraries = [
                donor.dna_ngs_library
                for pedigree in sheet.cohort.pedigrees
                for donor in pedigree.donors
                if donor.dna_ngs_library
            ]
            for ngs_library in ngs_libraries:
                name_pattern = "{mapper}.cnvetti_genotype.{ngs_library}".format(
                    ngs_library=ngs_library.name, **wildcards
                )
                for key, ext in BCF_KEY_EXTS.items():
                    result.setdefault(key, []).append(
                        os.path.join("work", name_pattern, "out", name_pattern + ext)
                    )
        return result

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        name_pattern = "{mapper}.cnvetti_merge_genotypes".format(**wildcards)
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    def get_output_files(self, action):
        """Return output files that CNVetti creates for the given action."""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_coverage(self):
        name_pattern = "{mapper}.cnvetti_coverage.{ngs_library}"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_segment(self):
        name_pattern = "{mapper}.cnvetti_segment.{ngs_library}"
        modifiers = (("segments_", ".segments"), ("windows_", ".windows"))
        for kmod, emod in modifiers:
            for key, ext in BCF_KEY_EXTS.items():
                yield kmod + key, os.path.join(
                    "work", name_pattern, "out", name_pattern + emod + ext
                )

    @dictify
    def _get_output_files_merge_segments(self):
        name_pattern = "{mapper}.cnvetti_merge_segments"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_genotype(self):
        name_pattern = "{mapper}.cnvetti_genotype.{ngs_library}"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_merge_genotypes(self):
        name_pattern = "{mapper}.cnvetti_merge_genotypes"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_reorder_vcf(self):
        name_pattern = "{mapper}.cnvetti.{ngs_library}"
        for key, ext in VCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)
        if action in ("merge_segments", "merge_genotypes"):
            prefix = "work/{{mapper}}.cnvetti_{action}/log/{{mapper}}.cnvetti_{action}".format(
                action=action
            )
        else:
            name_pattern = (
                "work/{{mapper}}.cnvetti_{action}.{{ngs_library}}/log/"
                "{{mapper}}.cnvetti_{action}.{{ngs_library}}"
            )
            prefix = name_pattern.format(action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.ngs_library]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        # Return resource
        return ResourceUsage(
            threads=1,
            time="04:00:00",  # 4 hours
            memory=f"{12 * 1024}M",
        )


class Delly2StepPart(BaseStepPart):
    """WGS CNV calling with Delly."""

    #: Step name
    name = "delly2"

    #: Supported actions
    actions = ("call", "merge_calls", "genotype", "merge_genotypes", "filter", "bcf_to_vcf")

    #: Directory infixes
    dir_infixes = {
        "call": r"{mapper,[^\.]+}.delly2.call.{library_name,[^\.]+}",
        "merge_calls": r"{mapper,[^\.]+}.delly2.merge_calls.{index_ngs_library,[^\.]+}",
        "genotype": "{mapper}.delly2.genotype.{library_name}",
        "merge_genotypes": r"{mapper}.delly2.merge_genotypes.{index_ngs_library,[^\.]+}",
        "filter": r"{mapper,[^\.]+}.delly2.filter.{index_ngs_library,[^\.]+}",
        "bcf_to_vcf": r"{mapper,[^\.]+}.delly2.bcf_to_vcf.{index_ngs_library,[^\.]+}",
    }

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "cheap_action": ResourceUsage(
            threads=2,
            time="4-00:00:00",  # 4 days
            memory=f"{7 * 1024 * 2}M",
        ),
        "default": ResourceUsage(
            threads=2,
            time="7-00:00:00",  # 7 days
            memory=f"{20 * 1024 * 2}M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return input function for Delly action."""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_call(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Yield input BAM and BAI file
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield key, ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    def _get_input_files_index_dependent_rules(self, wildcards, step):
        assert step in ("call", "genotype", "filter")
        # Create path template to per-sample call/genotype BCF
        infix = self.dir_infixes[step]
        infix = infix.replace(r",[^\.]+", "")
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        # Yield paths to pedigree's per-sample call BCF files
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                yield tpl.format(library_name=donor.dna_ngs_library.name, **wildcards)

    @listify
    def _get_input_files_merge_calls(self, wildcards):
        """Return input files for "merge_calls" action"""
        yield from self._get_input_files_index_dependent_rules(wildcards, "call")

    @dictify
    def _get_input_files_genotype(self, wildcards):
        """Return input files for "genotype" action"""
        pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
        # Per-pedigree site BCF file
        infix = self.dir_infixes["merge_calls"]
        infix = infix.replace(r",[^\.]+", "")
        infix = infix.format(
            mapper=wildcards.mapper, index_ngs_library=pedigree.index.dna_ngs_library.name
        )
        yield "bcf", os.path.join("work", infix, "out", infix + ".bcf").format(**wildcards)
        # BAM files
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @listify
    def _get_input_files_merge_genotypes(self, wildcards):
        """Return input files for "merge_genotypes" action"""
        yield from self._get_input_files_index_dependent_rules(wildcards, "genotype")

    @dictify
    def _get_input_files_filter(self, wildcards):
        name_pattern = "{mapper}.delly2.merge_genotypes.{index_ngs_library}".format(**wildcards)
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_input_files_bcf_to_vcf(self, wildcards):
        name_pattern = "{mapper}.delly2.filter.{library_name}".format(**wildcards)
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    def get_output_files(self, action):
        """Return output files that Delly creates for the given action."""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_call(self):
        name_pattern = "{mapper}.delly2.call.{library_name}"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_merge_calls(self):
        for name, ext in BCF_KEY_EXTS.items():
            infix = self.dir_infixes["merge_calls"]
            infix2 = infix.replace(r",[^\.]+", "")
            yield name, "work/" + infix + "/out/" + infix2 + ext

    @dictify
    def _get_output_files_genotype(self):
        name_pattern = "{mapper}.delly2.genotype.{library_name}"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_merge_genotypes(self):
        name_pattern = "{mapper}.delly2.merge_genotypes.{index_ngs_library}"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_filter(self):
        name_pattern = "{mapper}.delly2.filter.{index_ngs_library}"
        for key, ext in BCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_bcf_to_vcf(self):
        name_pattern = "{mapper}.delly2.{library_name}"
        for key, ext in VCF_KEY_EXTS.items():
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)
        if action in ("merge_calls", "merge_genotypes", "filter"):
            prefix = (
                f"work/{{mapper}}.delly2.{action}.{{index_ngs_library}}/log/"
                f"{{mapper}}.delly2.{action}.{{index_ngs_library}}"
            )
        else:
            name_pattern = (
                "work/{{mapper}}.delly2.{action}.{{library_name}}/log/"
                "{{mapper}}.delly2.{action}.{{library_name}}"
            )
            prefix = name_pattern.format(action=action)
        key_ext = (
            ("log", ".log"),
            ("log_md5", ".log.md5"),
            ("conda_info", ".conda_info.txt"),
            ("conda_info_md5", ".conda_info.txt.md5"),
            ("conda_list", ".conda_list.txt"),
            ("conda_list_md5", ".conda_list.txt.md5"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        if action in ("merge_genotypes", "merge_calls", "bcf_to_vcf"):  # cheap actions
            return self.resource_usage_dict.get("cheap_action")
        else:
            return self.resource_usage_dict.get("default")


class RunGcnvWgsStepPart(RunGcnvStepPart):
    """WGS CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        # No mapping given as WGS, we will use the "default" one for all.
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library:
                yield donor.dna_ngs_library.name, "default"


class WgsCnvCallingWorkflow(BaseStep):
    """Perform germline WGS CNV calling"""

    #: Workflow name
    name = "wgs_cnv_calling"

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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                CnvettiStepPart,
                Delly2StepPart,
                RunGcnvWgsStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def all_donors(self, include_background=True):
        """Get all donors.

        :param include_background: Boolean flag to defined if background should be included or not.
        Default: True, i.e., background will be included.

        :return: Returns list of all donors in sample sheet.
        """
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the germline WGS CNV calling workflow"""
        name_pattern = "{mapper}.{caller}.{index.dna_ngs_library.name}"
        tools = [t for t in self.config["tools"] if t != "gcnv"]
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            index_only=True,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=tools,
            ext=EXT_VALUES,
        )
        if "gcnv" in self.config["tools"]:
            yield from self.sub_steps["gcnv"].get_result_files()

    def _yield_result_files(self, tpl, index_only, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in self.shortcut_sheets:
            for pedigree in sheet.cohort.pedigrees:
                if index_only:
                    yield from expand(tpl, index=[pedigree.index], **kwargs)
                else:
                    donors = [d for d in pedigree.donors if d.dna_ngs_library]
                    yield from expand(tpl, index=donors, **kwargs)

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        self.ensure_w_config(
            ("step_config", "wgs_cnv_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for germline variant calling",
        )
