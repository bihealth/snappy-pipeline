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
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` (for the aligned reads BAM files) and the ``variant_calling`` (for looking at the
B allele frequency).

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

=======
Reports
=======

Currently, no reports are generated.
"""

from collections import OrderedDict
import os

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload (VCF)
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Name/extension pairs for VCF files.
VCF_KEY_EXTS = dict(zip(EXT_NAMES, EXT_VALUES))

#: Name/extension pairs for BCF files.
BCF_KEY_EXTS = {"bcf": ".bcf", "bcf_md5": ".bcf.md5", "csi": ".bcf.csi", "csi.md5": ".bcf.csi.md5"}

#: Available WGS CNV callers
WGS_CNV_CALLERS = ("erds", "erds_sv2", "cnvetti")

#: Default configuration for the wgs_cnv_calling schema
DEFAULT_CONFIG = r"""
# Default configuration wgs_cnv_calling
step_config:
  wgs_cnv_calling:
    path_ngs_mapping: ../ngs_mapping          # REQUIRED
    path_variant_calling: ../variant_calling  # REQUIRED
    variant_calling_tool: REQUIRED            # REQUIRED
    tools:
    - erds_sv2
    cnvetti:
      window_length: null
      count_kind: null
      segmentation: null
      normalization: null
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
    sv2:
      path_hg19: /fast/projects/cubit/18.12/static_data/reference/hg19/ucsc/hg19.fa  # REQUIRED
      path_hg38: /fast/projects/cubit/18.12/static_data/reference/hg38/ucsc/hg38.fa  # REQUIRED
      path_mm10: /fast/projects/cubit/18.12/static_data/reference/mm10/ucsc/mm10.fa  # REQUIRED
      path_sv2_resource: /fast/work/users/holtgrem_c/cubit_incoming/SV2/v1.4.3.4   # REQUIRED
      path_sv2_models: /fast/work/users/holtgrem_c/cubit_incoming/SV2/v1.4.3.4/training_sets   # REQUIRED
"""


class CnvettiStepPart(BaseStepPart):
    """Shallow WGS CNV calling with CNVetti."""

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
        assert action in self.actions
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
        assert action in self.actions
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
        if action in ("merge_segments", "merge_genotypes"):
            prefix = "work/{{mapper}}.cnvetti_{action}/log/" "{{mapper}}.cnvetti_{action}".format(
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

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration for XHMM CNV calling"""
        for action in self.actions:
            cluster_config["wgs_cnv_calling_cnvetti_{}".format(action)] = {
                "mem": int(12 * 1024),
                "time": "04:00",
                "ntasks": 1,
            }


class ErdsStepPart(BaseStepPart):
    """WGS CNV calling with ERDS

    WGS CNV calling is performed on the primary DNA NGS library for each donor.
    """

    name = "erds"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.erds.{{library_name}}/out/{{mapper}}.erds.{{library_name}}{ext}"
        )
        # Build shortcut from index library name to donor
        self.index_ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_donor.update(sheet.index_ngs_library_to_donor)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return input function for ERDS rule"""

        @dictify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to sub workflows
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            var_calling = self.parent.sub_workflows["variant_calling"]
            # Get names of BAM and BAI file for the given donor
            bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
            for ext in (".bam", ".bam.bai"):
                yield ext.split(".")[-1], ngs_mapping(bam_tpl.format(ext=ext, **wildcards))
            # Get names of the VCF files for the given library's donor/pedigree
            pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
            vcf_base = (
                "output/{mapper}.{var_caller}.{index_library_name}/out/"
                "{mapper}.{var_caller}.{index_library_name}"
            ).format(
                var_caller=self.config["variant_calling_tool"],
                index_library_name=pedigree.index.dna_ngs_library.name,
                **wildcards
            )
            yield "vcf", var_calling(vcf_base + ".vcf.gz")
            yield "tbi", var_calling(vcf_base + ".vcf.gz.tbi")

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files that ERDS returns return (VCF + TBI file)"""
        assert action == "run"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, self.base_path_out.format(ext=ext)

    def get_log_file(self, action):
        """Return path to log file"""
        return "work/{mapper}.erds.{library_name}/log/snakemake.wgs_cnv_calling.log"

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration for ERDS WGS CNV calling"""
        cluster_config["wgs_cnv_calling_erds_run"] = {
            "mem": 32 * 1024,
            "time": "48:00",
            "ntasks": 1,
        }


class ErdsSv2StepPart(BaseStepPart):
    """WGS SV identification using ERDS+SV2

    ERDS+SV2 supports the calling based on whole cohorts.  The rough steps are as follows:

    - Perform variant calling on each sample individually with ERDS ("erds_sv2_call")
    - Merge called variants to get a cohort-wide site list ("erds_sv2_merge_calls")
    - Perform genotyping of the variants in the cohort-wide site list in each sample
      with SV2 ("erds_sv2_genotype")
    - Merge cohort-wide site list ("erds_sv2_merge_genotypes"); using bcftools
    - Reorder VCF and put pedigree in front; later on, non-pedigree variants should be removed.
    """

    name = "erds_sv2"

    #: Actions in ERDS+SV2 workflow
    actions = (
        "call",
        "merge_calls",
        "genotype",
        "info_to_format",
        "merge_genotypes",
        "reorder_vcf",
    )

    #: Directory infixes
    dir_infixes = {
        "call": "{mapper}.erds_sv2.call.{library_name}",
        "merge_calls": "{mapper}.erds_sv2.merge_calls",
        "genotype": "{mapper}.erds_sv2.genotype.{library_name}",
        "info_to_format": "{mapper}.erds_sv2.info_to_format.{library_name}",
        "merge_genotypes": "{mapper}.erds_sv2.merge_genotypes",
        "reorder_vcf": r"{mapper}.erds_sv2.{index_ngs_library,(?!call|merge_|genotype|info_to_format).*}",
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

    def get_library_extra_infos(self, wildcards):
        """Returns library extra infos for the given library name"""
        return self.library_name_to_library[wildcards.library_name].ngs_library.extra_infos

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        assert action in self.actions
        mapping = {
            "call": self._get_input_files_call,
            "merge_calls": self._get_input_files_merge_calls,
            "genotype": self._get_input_files_genotype,
            "info_to_format": self._get_input_files_info_to_format,
            "merge_genotypes": self._get_input_files_merge_genotypes,
            "reorder_vcf": self._get_input_files_reorder_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_call(self, wildcards):
        """Return input files for "call" action"""
        return ErdsStepPart(self.parent).get_input_files("run")(wildcards)

    @listify
    def _get_input_files_merge_calls(self, wildcards):
        """Return input files for "merge_calls" action"""
        tpl = os.path.join(
            "work", self.dir_infixes["call"], "out", self.dir_infixes["call"] + ".vcf.gz"
        )
        for donor in self._donors_with_dna_ngs_library():
            yield tpl.format(library_name=donor.dna_ngs_library.name, **wildcards)

    @dictify
    def _get_input_files_genotype(self, wildcards):
        """Return input files for "genotype" action"""
        # Get shorcut to sub workflows
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        var_calling = self.parent.sub_workflows["variant_calling"]
        # Get names of BAM file for the given donor
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        yield "bam", ngs_mapping(bam_tpl.format(ext=".bam", **wildcards))
        # CNV sites VCF file
        infix = self.dir_infixes["merge_calls"]
        yield "vcf_cnv", os.path.join("work", infix, "out", infix + ".vcf.gz").format(**wildcards)
        # Small variants VCF file
        pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
        vcf_base = (
            "output/{mapper}.{var_caller}.{index_library_name}/out/"
            "{mapper}.{var_caller}.{index_library_name}"
        ).format(
            var_caller=self.config["variant_calling_tool"],
            index_library_name=pedigree.index.dna_ngs_library.name,
            **wildcards
        )
        yield "vcf_small", var_calling(vcf_base + ".vcf.gz")
        # PED file
        tpl = "work/write_pedigree.{index_library_name}/out/{index_library_name}.ped"
        yield "ped", tpl.format(index_library_name=pedigree.index.dna_ngs_library.name, **wildcards)

    def _get_input_files_info_to_format(self, wildcards):
        """Return input files for "info_to_format" action"""
        infix = self.dir_infixes["genotype"]
        tpl = os.path.join("work", infix, "out", infix + ".vcf.gz")
        yield tpl.format(**wildcards)

    @listify
    def _get_input_files_merge_genotypes(self, wildcards):
        """Return input files for "merge_genotypes" action"""
        for donor in self._donors_with_dna_ngs_library():
            infix = self.dir_infixes["info_to_format"]
            tpl = os.path.join("work", infix, "out", infix + ".vcf.gz")
            yield tpl.format(**{**wildcards, "library_name": donor.dna_ngs_library.name})

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        """Return input files for "reorder_vcf" action"""
        infix = self.dir_infixes["merge_genotypes"]
        tpl = os.path.join("work", infix, "out", infix + ".vcf.gz")
        yield "vcf", tpl.format(**wildcards)

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action in self.actions
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            infix1 = self.dir_infixes[action]
            infix2 = infix1.replace(r",(?!call|merge_|genotype|info_to_format).*", "")
            yield name, "work/" + infix1 + "/out/" + infix2 + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action in self.actions
        infix1 = self.dir_infixes[action]
        return "work/" + infix1 + "/log/snakemake.log"

    def update_cluster_config(self, cluster_config):
        for action in self.actions:
            # TODO: refine for ERDS and SV2
            if action in ("info_to_format", "merge_genotypes", "merge_calls", "reorder_vcf"):
                # cheap actions
                cluster_config["wgs_cnv_calling_erds_sv2_{}".format(action)] = {
                    "mem": int(3.75 * 1024 * 2),
                    "time": "24:00",
                    "ntasks": 2,
                }
            elif action == "call":
                # ERDS is pretty memory hungry
                cluster_config["wgs_cnv_calling_erds_sv2_{}".format(action)] = {
                    "mem": 40 * 1024,
                    "time": "150:00",
                    "ntasks": 1,
                }
            else:
                cluster_config["wgs_cnv_calling_erds_sv2_{}".format(action)] = {
                    "mem": int(7.5 * 1024 * 4),
                    "time": "150:00",
                    "ntasks": 4,
                }


class WgsCnvCallingWorkflow(BaseStep):
    """Perform germline WGS CNV calling"""

    name = "wgs_cnv_calling"
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
            (NgsMappingWorkflow, VariantCallingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, CnvettiStepPart, ErdsStepPart, ErdsSv2StepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        self.register_sub_workflow("variant_calling", self.config["path_variant_calling"])

    @listify
    def get_result_files(self):
        """Return list of result files for the germline WGS CNV calling workflow"""
        name_pattern = "{mapper}.{caller}.{index.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            index_only=True,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=[t for t in self.config["tools"] if t != "erds"],
            ext=EXT_VALUES,
        )
        if "erds" in self.config["tools"]:
            yield from self._yield_result_files(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                index_only=False,
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=["erds"],
                ext=EXT_VALUES,
            )

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
        self.ensure_w_config(
            ("step_config", "wgs_cnv_calling", "path_variant_calling"),
            (
                "Path to germline (small) variant calling not configured but required for germline "
                "WGS CNV calling"
            ),
        )
        self.ensure_w_config(
            ("step_config", "wgs_cnv_calling", "variant_calling_tool"),
            "Name of the germline (small) variant calling tool",
        )
        if (
            self.w_config["step_config"]["wgs_cnv_calling"]["variant_calling_tool"]
            not in self.w_config["step_config"]["variant_calling"]["tools"]
        ):
            tpl = "Variant caller {} is not selected in variant_calling step"
            raise MissingConfiguration(
                tpl.format(self.w_config["step_config"]["wgs_cnv_calling"]["variant_calling_tool"])
            )
