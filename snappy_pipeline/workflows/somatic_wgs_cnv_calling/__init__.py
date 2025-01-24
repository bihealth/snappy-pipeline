# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_wgs_cnv_calling`` step

The ``somatic_wgs_cnv_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned NGS reads) and performs somatic CNV calling on them.  The result are called CNVs in VCF
format.

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` and ``somatic_variant_calling`` steps.  Somatic (small) variant calling is required
for b-allele based filtration.  For the somatic variant calling, one somatic (small) variant
caller must be configured of which to use the results.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk`` and each read mapper
``mapper`` that the library has been aligned with, and the variant caller ``var_caller``, the
pipeline step will create a directory ``output/{mapper}.{var_caller}.{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.canvas.P001-T1-DNA1-WGS1-4
    |   `-- out
    |       |-- bwa.canvas.P001-T1-DNA1-WGS1-4.vcf.gz
    |       |-- bwa.canvas.P001-T1-DNA1-WGS1-4.vcf.gz.tbi
    |       |-- bwa.canvas.P001-T1-DNA1-WGS1-4.vcf.gz.md5
    |       `-- bwa.canvas.P001-T1-DNA1-WGS1-4.vcf.gz.tbi.md5
    [...]

Generally, these files will be unfiltered, i.e., contain low-quality variants and also variants
flagged as being non-somatic.

====================
Global Configuration
====================

None so far

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_wgs_cnv_calling.rst

=============================
Available Somatic CNV Callers
=============================

The following somatic CNV callers are currently available

- ``"canvas"``

=======
Reports
=======

Currently, no reports are generated.
"""

import os
import sys
from collections import OrderedDict
from itertools import chain
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import SomaticWgsCnvCalling as SomaticWgsCnvCallingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
BCF_EXT_VALUES = (".bcf", ".bcf.csi", ".bcf.md5", ".bcf.csi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available somatic WGS CNV callers
SOMATIC_WGS_CNV_CALLERS = ("canvas", "cnvetti", "control_freec")

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = SomaticWgsCnvCallingConfigModel.default_config_yaml_string()


class SomaticWgsCnvCallingStepPart(BaseStepPart):
    """Base class for somatic WGS CNV calling steps

    WGS CNV calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    # TODO: unify with somatic (small) variant calling base class?

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{cancer_library}}/out/"
            "{{mapper}}.{var_caller}.{{cancer_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched tumor sample
        self.cancer_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.cancer_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build shortcut from library name to donor.
        self.donor_to_cancer_ngs_libraries = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for pair in sheet.all_sample_pairs_by_tumor_dna_ngs_library.values():
                self.donor_to_cancer_ngs_libraries.setdefault(pair.donor.name, [])
                if pair.tumor_sample.name not in self.donor_to_cancer_ngs_libraries:
                    self.donor_to_cancer_ngs_libraries[pair.donor.name].append(
                        pair.tumor_sample.dna_ngs_library.name
                    )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflows
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            cancer_base_path = (
                "output/{mapper}.{cancer_library}/out/" "{mapper}.{cancer_library}"
            ).format(**wildcards)
            yield "normal_bam", ngs_mapping(normal_base_path + ".bam")
            yield "normal_bai", ngs_mapping(normal_base_path + ".bam.bai")
            yield "tumor_bam", ngs_mapping(cancer_base_path + ".bam")
            yield "tumor_bai", ngs_mapping(cancer_base_path + ".bam.bai")

        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.cancer_ngs_library_to_sample_pair[wildcards.cancer_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        # Validate action
        self._validate_action(action)
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, var_caller=[self.name], ext=EXT_VALUES))
        )

    @dictify
    def _get_log_file(self, action):
        """Return path to log file"""
        # Validate action
        self._validate_action(action)

        name_pattern = "{{mapper}}.{var_caller}.{{cancer_library}}".format(
            var_caller=self.__class__.name
        )
        prefix = "work/{name_pattern}/log/{name_pattern}".format(name_pattern=name_pattern)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext


class CanvasSomaticWgsStepPart(SomaticWgsCnvCallingStepPart):
    """Somatic WGS SV calling with Canvas"""

    #: Step name
    name = "canvas"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{int(3.75 * 1024 * 16)}M",
        )

    def get_params(self, action):
        self._validate_action(action)

        def args_fn(_wildcards):
            return dict(self.config.canvas.path_reference)

        return args_fn


class CnvettiSomaticWgsStepPart(SomaticWgsCnvCallingStepPart):
    """Somatic WGS CNV calling with CNVetti"""

    #: Step name
    name = "cnvetti"

    #: Class available actions
    actions = ("coverage", "tumor_normal_ratio", "segment")

    #: Extension file dictionary. Key: file type (string); Value: file extension (string)
    bcf_dict = {
        "bcf": ".bcf",
        "bcf_csi": ".bcf.csi",
        "bcf_md5": ".bcf.md5",
        "bcf_csi_md5": ".bcf.csi.md5",
    }

    #: Parameters to pass to the wrapper
    params_for_action = {
        "coverage": ("window_length", "count_kind", "normalization"),
        "segment": ("segmentation",),
    }

    def get_input_files(self, action):
        """Return input function for CNVetti rule"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_coverage(self, wildcards):
        """Return input files that "cnvetti coverage" needs"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Yield input BAM and BAI file
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for ext in (".bam", ".bam.bai"):
            yield ext.split(".")[-1], ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_input_files_tumor_normal_ratio(self, wildcards):
        """Return input files that the merge step ("bcftools merge") needs"""
        # TODO: Potential bug as 'library_name' is required in the wildcards but also obtained using
        #  `get_normal_lib_name()`.
        #  Error: "TypeError: str.format() got multiple values for keyword argument 'library_name'"
        libraries = {"tumor": wildcards.library_name, "normal": self.get_normal_lib_name(wildcards)}
        for kind, library_name in libraries.items():
            key = "{}_bcf".format(kind)
            name_pattern = "{mapper}.cnvetti_coverage.{library_name}".format(
                library_name=library_name, **wildcards
            )
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=".bcf"
                ),
            )

    @dictify
    def _get_input_files_segment(self, wildcards):
        """Return input files that "cnvetti segment" needs"""
        for key, ext in self.bcf_dict.items():
            name_pattern = "{mapper}.cnvetti_tumor_normal_ratio.{library_name}".format(**wildcards)
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    def get_output_files(self, action):
        """Return output files that CNVetti creates for the given action"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_coverage(self):
        for key, ext in self.bcf_dict.items():
            name_pattern = "{mapper}.cnvetti_coverage.{library_name}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    @dictify
    def _get_output_files_tumor_normal_ratio(self):
        for key, ext in self.bcf_dict.items():
            name_pattern = "{mapper}.cnvetti_tumor_normal_ratio.{library_name}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    @dictify
    def _get_output_files_segment(self):
        for key, ext in self.bcf_dict.items():
            name_pattern = "{mapper}.cnvetti_segment.{library_name}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    def get_args(self, action: str) -> dict[str, Any]:
        """Return args (params) that CNVetti creates for the given action"""
        # Validate action
        self._validate_action(action)

        params = {}
        if action in self.params_for_action:
            cfg = getattr(self.config, self.name)
            assert cfg.preset in cfg.presets, f"Undefined preset '{cfg.preset}'"
            for k in self.params_for_action[action]:
                v = getattr(cfg, k, None)
                if v is None:
                    assert (
                        k in cfg.presets[cfg.preset]
                    ), f"Missing parameter '{k}' from preset '{cfg.preset}'"
                    v = cfg.presets[cfg.preset].get(k)
                params[k] = v

        if action == "coverage":
            params["reference"] = self.parent.w_config.static_data_config.reference.path

        return getattr(self, "_get_args_{}".format(action))

    @dictify
    def get_log_file(self, action):
        """Return path to log file"""
        name_pattern = "{{mapper}}.cnvetti_{action}.{{library_name}}".format(action=action)
        prefix = "work/{name_pattern}/log/{name_pattern}".format(name_pattern=name_pattern)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="04:00:00",  # 4 hours
            memory=f"{int(3.75 * 1024 * 4)}M",
        )


class CnvkitSomaticWgsStepPart(SomaticWgsCnvCallingStepPart):
    """Somatic WGS CNV calling with cnvkit.py"""

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "coverage",
        "fix",
        "segment",
        "call",
        "postprocess",
        "export",
        "plot",
        "report",
    )

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "plot": ResourceUsage(
            threads=1,
            time="1-00:00:00",  # 1 day
            memory=f"{30 * 1024}M",
        ),
        "coverage": ResourceUsage(
            threads=8,
            time="1-00:00:00",  # 1 day
            memory=f"{16 * 1024}M",
        ),
        "default": ResourceUsage(
            threads=1,
            time="1-00:00:00",  # 1 day
            memory=f"{int(7.5 * 1024)}M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        method_mapping = {
            "coverage": self._get_input_files_coverage,
            "call": self._get_input_files_call,
            "fix": self._get_input_files_fix,
            "segment": self._get_input_files_segment,
            "postprocess": self._get_input_files_postprocess,
            "export": self._get_input_files_export,
            "plot": self._get_input_files_plot,
            "report": self._get_input_files_report,
        }
        return method_mapping[action]

    def _get_input_files_coverage(self, wildcards):
        # BAM/BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam.bai"),
        }
        return input_files

    @staticmethod
    def _get_input_files_fix(wildcards):
        tpl_base = "{mapper}.cnvkit.{library_name}"
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".{target}coverage.cnn"
        input_files = {
            "target": tpl.format(target="target", **wildcards),
            "antitarget": tpl.format(target="antitarget", **wildcards),
        }
        return input_files

    @staticmethod
    def _get_input_files_segment(wildcards):
        cnr_pattern = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cnr"
        input_files = {"cnr": cnr_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_call(wildcards):
        segment_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.segment.cns"
        )
        input_files = {"segment": segment_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_postprocess(wildcards):
        segment_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.call.cns"
        )
        input_files = {"call": segment_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_export(wildcards):
        cns_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.call.cns"
        )
        input_files = {"cns": cns_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_plot(wildcards):
        tpl = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{ext}"
        input_files = {
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def _get_input_files_report(self, wildcards):
        tpl = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{ext}"
        input_files = {
            "target": tpl.format(ext="targetcoverage.cnn", **wildcards),
            "antitarget": tpl.format(ext="antitargetcoverage.cnn", **wildcards),
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def get_output_files(self, action):
        """Return output files for the given action"""
        # Validate action
        self._validate_action(action)
        method_mapping = {
            "coverage": self._get_output_files_coverage,
            "fix": self._get_output_files_fix,
            "call": self._get_output_files_call,
            "postprocess": self._get_output_files_postprocess,
            "segment": self._get_output_files_segment,
            "export": self._get_output_files_export,
            "plot": self._get_output_files_plot,
            "report": self._get_output_files_report,
        }
        return method_mapping[action]()

    @staticmethod
    def _get_output_files_coverage():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
            output_files[target + "_md5"] = output_files[target] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_fix():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cnr")
        return {"ratios": tpl, "ratios_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_segment():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".segment.cns")
        return {"segments": tpl, "segments_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_call():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".call.cns")
        return {"calls": tpl, "calls_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_postprocess():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cns")
        return {"final": tpl, "final_md5": tpl + ".md5"}

    @dictify
    def _get_output_files_plot(self):
        plots = (("diagram", "pdf"), ("heatmap", "pdf"), ("scatter", "png"))
        chrom_plots = (("heatmap", "pdf"), ("scatter", "png"))
        chroms = list(chain(range(1, 23), ["X", "Y"]))
        output_files = {}
        # Yield file name pairs for global plots
        tpl = (
            "work/{{mapper}}.cnvkit.{{library_name}}/report/"
            "{{mapper}}.cnvkit.{{library_name}}.{plot}.{ext}"
        )
        for plot, ext in plots:
            output_files[plot] = tpl.format(plot=plot, ext=ext)
            output_files[plot + "_md5"] = output_files[plot] + ".md5"
        # Yield file name pairs for the chromosome-wise plots
        tpl_chrom = (
            "work/{{mapper}}.cnvkit.{{library_name}}/report/"
            "{{mapper}}.cnvkit.{{library_name}}.{plot}.chr{chrom}.{ext}"
        )
        for plot, ext in chrom_plots:
            for chrom in chroms:
                key = "{plot}_chr{chrom}".format(plot=plot, chrom=chrom)
                output_files[key] = tpl_chrom.format(plot=plot, ext=ext, chrom=chrom)
                output_files[key + "_md5"] = output_files[key] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_export():
        exports = (("bed", "bed"), ("seg", "seg"), ("vcf", "vcf.gz"), ("tbi", "vcf.gz.tbi"))
        output_files = {}
        tpl = (
            "work/{{mapper}}.cnvkit.{{library_name}}/out/"
            "{{mapper}}.cnvkit.{{library_name}}.{ext}"
        )
        for export, ext in exports:
            output_files[export] = tpl.format(export=export, ext=ext)
            output_files[export + "_md5"] = output_files[export] + ".md5"
        return output_files

    @dictify
    def _get_output_files_report(self):
        reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        output_files = {}
        tpl = (
            "work/{{mapper}}.cnvkit.{{library_name}}/report/"
            "{{mapper}}.cnvkit.{{library_name}}.{report}.txt"
        )
        for report in reports:
            output_files[report] = tpl.format(report=report)
            output_files[report + "_md5"] = output_files[report] + ".md5"
        return output_files

    def get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{{mapper}}.cnvkit.{{library_name}}/log/"
            "{{mapper}}.cnvkit.{action}.{{library_name}}"
        ).format(action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        log_files = {}
        for key, ext in key_ext:
            log_files[key] = prefix + ext
            log_files[key + "_md5"] = prefix + ext + ".md5"
        return log_files

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        if action == "plot":
            return self.resource_usage_dict.get("plot")
        elif action == "coverage":
            return self.resource_usage_dict.get("coverage")
        else:
            return self.resource_usage_dict.get("default")


class ControlFreecSomaticWgsStepPart(SomaticWgsCnvCallingStepPart):
    """Somatic WGS CNV calling with Control-FreeC"""

    #: Step name
    name = "control_freec"

    #: Class available actions
    actions = ("run", "transform", "plot")

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "plot": ResourceUsage(
            threads=1,
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{2 * 30 * 1024}M",
        ),
        "transform": ResourceUsage(
            threads=1,
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{2 * 8 * 1024}M",
        ),
        "run": ResourceUsage(
            threads=8,
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{int(2 * 3.75 * 1024 * 8)}M",
        ),
    }

    def get_output_files(self, action):
        # Initialise variable
        result = {}

        # Validate action
        self._validate_action(action)

        if action == "run":
            result["ratio"] = self.base_path_out.format(var_caller=self.name, ext=".ratio.txt")
            result["ratio_md5"] = self.base_path_out.format(
                var_caller=self.name, ext=".ratio.txt.md5"
            )
        elif action == "transform":
            transform_ext_names = ("log2", "call", "segments", "cns", "cnr")
            transform_ext_values = (
                "_gene_log2.txt",
                "_gene_call.txt",
                "_segments.txt",
                ".cns",
                ".cnr",
            )
            result = dict(
                zip(
                    transform_ext_names,
                    expand(self.base_path_out, var_caller=[self.name], ext=transform_ext_values),
                )
            )
        elif action == "plot":
            plot_ext_names = ("heatmap", "scatter", "diagram")
            plot_ext_values = (".heatmap.png", ".scatter.png", ".diagram.pdf")
            result = dict(
                zip(
                    plot_ext_names,
                    expand(self.base_path_out, var_caller=[self.name], ext=plot_ext_values),
                )
            )

        return result

    def get_args(self, action: str) -> dict[str, Any]:
        # Validate action
        self._validate_action(action)
        cfg = self.config.control_freec
        if action == "run":
            return {
                "path_chrlenfile": cfg.path_chrlenfile,
                "path_mappability": cfg.path_mappability,
                "path_mappability_enabled": cfg.path_mappability_enabled,
                "window_size": cfg.window_size,
            }
        elif action == "transform":
            return {
                "org_obj": cfg.convert.org_obj,
                "tx_obj": cfg.convert.tx_obj,
                "bs_obj": cfg.convert.bs_obj,
            }
        elif action == "plot":
            return {}

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return self.resource_usage_dict.get(action)


class SomaticWgsCnvCallingWorkflow(BaseStep):
    """Perform somatic WGS CNV calling"""

    #: Workflow name
    name = "somatic_wgs_cnv_calling"

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
            config_model_class=SomaticWgsCnvCallingConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                CanvasSomaticWgsStepPart,
                CnvettiSomaticWgsStepPart,
                CnvkitSomaticWgsStepPart,
                ControlFreecSomaticWgsStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        self.register_sub_workflow(
            "somatic_variant_calling", self.config.path_somatic_variant_calling
        )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.somatic_variant_calling_tool:
            self.config.somatic_variant_calling_tool = self.w_config.step_config[
                "somatic_variant_calling"
            ].tools[0]

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.{caller}.{cancer_library.name}"
        tpl = os.path.join("output", name_pattern, "out", name_pattern + "{ext}")
        vcf_tools = [
            t for t in self.config.tools if t not in ("cnvetti", "control_freec", "cnvkit")
        ]
        bcf_tools = [t for t in self.config.tools if t in ("cnvetti",)]
        yield from self._yield_result_files(
            tpl,
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            caller=vcf_tools,
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files(
            tpl,
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            caller=bcf_tools,
            ext=BCF_EXT_VALUES,
        )
        if "control_freec" in self.config.tools:
            yield from self._yield_result_files(
                tpl,
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller="control_freec",
                ext=[
                    ".ratio.txt",
                    ".ratio.txt.md5",
                    ".gene_log2.txt",
                    ".gene_call.txt",
                    ".segments.txt",
                    ".scatter.png",
                    ".heatmap.png",
                    ".diagram.pdf",
                ],
            )
        # Plots for cnvetti
        if "cnvkit" in self.config.tools:
            exts = (".cnr", ".cns", ".bed", ".seg", ".vcf.gz", ".vcf.gz.tbi")
            yield from self._yield_result_files(
                tpl,
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller="cnvkit",
                ext=exts,
            )
            yield from self._yield_result_files(
                tpl,
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller="cnvkit",
                ext=[ext + ".md5" for ext in exts],
            )
            reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
            yield from self._yield_report_files(
                (
                    "output/{mapper}.{caller}.{cancer_library.name}/report/"
                    "{mapper}.{caller}.{cancer_library.name}.{ext}"
                ),
                [(report, "txt", False) for report in reports],
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller="cnvkit",
            )
            if self.config.cnvkit.plot:
                plots = (
                    ("diagram", "pdf", False),
                    ("heatmap", "pdf", True),
                    ("scatter", "png", True),
                )
                yield from self._yield_report_files(
                    (
                        "output/{mapper}.{caller}.{cancer_library.name}/report/"
                        "{mapper}.{caller}.{cancer_library.name}.{ext}"
                    ),
                    plots,
                    mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                    caller="cnvkit",
                )
        if "cnvetti" in bcf_tools:
            for sheet in filter(is_not_background, self.shortcut_sheets):
                for donor in sheet.donors:
                    if donor.all_pairs:
                        name_pattern = "{mapper}.cnvetti_plot.{donor}"
                        for ext in (".png", ".png.md5"):
                            yield from expand(
                                os.path.join(
                                    "output", name_pattern, "out", name_pattern + "_genome" + ext
                                ),
                                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                                donor=[donor.name],
                            )
                            yield from expand(
                                os.path.join(
                                    "output",
                                    name_pattern,
                                    "out",
                                    name_pattern + "_chr{chrom}" + ext,
                                ),
                                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                                donor=[donor.name],
                                chrom=map(str, chain(range(1, 23), ("X", "Y"))),
                            )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )  # pragma: no cover
                    print(
                        msg.format(sample_pair.tumor_sample.name), file=sys.stderr
                    )  # pragma: no cover
                    continue  # pragma: no cover
                yield from expand(
                    tpl, cancer_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def _yield_report_files(self, tpl, exts, **kwargs):
        """Build report paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )  # pragma: no cover
                    print(
                        msg.format(sample_pair.tumor_sample.name), file=sys.stderr
                    )  # pragma: no cover
                    continue  # pragma: no cover
                for plot, ext, chrom in exts:
                    yield from expand(
                        tpl,
                        cancer_library=[sample_pair.tumor_sample.dna_ngs_library],
                        ext=plot + "." + ext,
                        **kwargs,
                    )
                    yield from expand(
                        tpl,
                        cancer_library=[sample_pair.tumor_sample.dna_ngs_library],
                        ext=plot + "." + ext + ".md5",
                        **kwargs,
                    )
                    if chrom:
                        for c in map(str, chain(range(1, 23), ("X", "Y"))):
                            yield from expand(
                                tpl,
                                cancer_library=[sample_pair.tumor_sample.dna_ngs_library],
                                ext=plot + ".chr" + c + "." + ext,
                                **kwargs,
                            )
                            yield from expand(
                                tpl,
                                cancer_library=[sample_pair.tumor_sample.dna_ngs_library],
                                ext=plot + ".chr" + c + "." + ext + ".md5",
                                **kwargs,
                            )
