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

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk``
the pipeline step will create a directory ``output/{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{lib_name}-{lib_pk}.vcf.gz``
- ``{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- P001-T1-DNA1-WGS1-4
    |   `-- out
    |       |-- P001-T1-DNA1-WGS1-4.vcf.gz
    |       |-- P001-T1-DNA1-WGS1-4.vcf.gz.tbi
    |       |-- P001-T1-DNA1-WGS1-4.vcf.gz.md5
    |       `-- P001-T1-DNA1-WGS1-4.vcf.gz.tbi.md5
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

.. include:: DEFAULT_CONFIG_somatic_wgs_cnv_calling.rst

=============================
Available Somatic CNV Callers
=============================

The following somatic CNV callers are currently available

- ``"cnvetti"``
- ``"cnvkit"``
- ``"control_freec"``

=======
Reports
=======

Currently, no reports are generated.
"""

import os
from itertools import chain
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.models.cnvkit import Gender as CnvkitGender
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import SomaticWgsCnvCalling as SomaticWgsCnvCallingConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
BCF_EXT_VALUES = (".bcf", ".bcf.csi", ".bcf.md5", ".bcf.csi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available somatic WGS CNV callers
SOMATIC_WGS_CNV_CALLERS = ("cnvetti", "control_freec")

#: Default configuration for the somatic_variant_calling schema


class SomaticWgsCnvCallingStepPart(BaseStepPart):
    """Base class for somatic WGS CNV calling steps

    WGS CNV calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    # TODO: unify with somatic (small) variant calling base class?

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{{tumor_library}}/out/{{tumor_library}}{ext}"

    @dictify
    def _get_input_files_run(self, wildcards):
        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        normal = self.parent.get_upstream_paths(
            "alignments", library_name=self.get_normal_lib_name(wildcards)
        )
        tumor = self.parent.get_upstream_paths("alignments", library_name=wildcards.tumor_library)
        yield "normal_bam", normal.bam
        yield "normal_bai", normal.bai
        yield "tumor_bam", tumor.bam
        yield "tumor_bai", tumor.bai

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        lib_name = getattr(wildcards, "tumor_library", None) or getattr(
            wildcards, "library_name", None
        )
        if not lib_name:
            return None
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == lib_name]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        # Validate action
        self._validate_action(action)
        return dict(zip(EXT_NAMES, expand(self.base_path_out, ext=EXT_VALUES)))

    @dictify
    def _get_log_file(self, action):
        """Return path to log file"""
        # Validate action
        self._validate_action(action)

        name_pattern = "{{tumor_library}}".format()
        prefix = "work/{name_pattern}/log/{name_pattern}".format(name_pattern=name_pattern)
        if len(self.actions) > 1:  # one log per action, not one shared by all
            prefix += f".{action}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext


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

    @dictify
    def _get_input_files_coverage(self, wildcards):
        """Return input files that "cnvetti coverage" needs"""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        # Yield input BAM and BAI file
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    @dictify
    def _get_input_files_tumor_normal_ratio(self, wildcards):
        tumor_library = getattr(wildcards, "library_name", None) or getattr(
            wildcards, "tumor_library", None
        )
        normal_library = self.get_normal_lib_name(wildcards)
        libraries = {"tumor": tumor_library, "normal": normal_library}
        for kind, library_name in libraries.items():
            if not library_name:
                continue
            key = "{}_bcf".format(kind)
            yield (
                key,
                f"work/{library_name}/out/{library_name}.coverage.bcf",
            )

    @dictify
    def _get_input_files_segment(self, wildcards):
        """Return input files that "cnvetti segment" needs"""
        for key, ext in self.bcf_dict.items():
            name_pattern = "{tumor_library}".format(**wildcards)
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}.tumor_normal_ratio{ext}".format(
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
            name_pattern = "{library_name}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}.coverage{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    @dictify
    def _get_output_files_tumor_normal_ratio(self):
        for key, ext in self.bcf_dict.items():
            name_pattern = "{library_name}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}.tumor_normal_ratio{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    @dictify
    def _get_output_files_segment(self):
        for key, ext in self.bcf_dict.items():
            name_pattern = "{tumor_library}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    def _get_params_coverage(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._get_params("coverage")

    def _get_params_tumor_normal_ratio(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._get_params("tumor_normal_ratio")

    def _get_params_segment(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._get_params("segment")

    def _get_params(self, action: str) -> dict[str, Any]:
        """Return args (params) that CNVetti creates for the given action"""
        cfg = getattr(self.config, self.name)

        preset_cfg = cfg.presets.get(cfg.preset)
        assert preset_cfg is not None, f"Undefined preset '{cfg.preset}'"
        preset_values = (
            preset_cfg.model_dump(by_alias=True)
            if hasattr(preset_cfg, "model_dump")
            else preset_cfg
        )

        params = {}
        if action in self.params_for_action:
            for k in self.params_for_action[action]:
                v = getattr(cfg, k, None)
                if v is None:
                    assert isinstance(preset_values, dict) and k in preset_values, (
                        f"Missing parameter '{k}' from preset '{cfg.preset}'"
                    )
                    v = preset_values.get(k)
                params[k] = v

        return params

    @dictify
    def get_log_file(self, action):
        """Return path to log file"""
        wildcard_name = (
            "library_name" if action in {"coverage", "tumor_normal_ratio"} else "tumor_library"
        )
        name_pattern = f"{{{wildcard_name}}}"
        prefix = f"work/{name_pattern}/log/{name_pattern}.{action}"
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
            runtime="4h",  # 4 hours
            mem=f"{int(3.75 * 1024 * 4)}MB",
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
            runtime="1d",  # 1 day
            mem=f"{30 * 1024}MB",
        ),
        "coverage": ResourceUsage(
            threads=8,
            runtime="1d",  # 1 day
            mem=f"{16 * 1024}MB",
        ),
        "default": ResourceUsage(
            threads=1,
            runtime="1d",  # 1 day
            mem=f"{int(7.5 * 1024)}MB",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)

    def _cnvkit_params(self, action: str) -> dict[str, Any]:
        cfg = self.config.cnvkit
        args = {}
        if action != "report":
            args = getattr(cfg, action).model_dump(by_alias=True)
        else:
            args = {
                "breaks": cfg.breaks.model_dump(by_alias=True),
                "genemetrics": cfg.genemetrics.model_dump(by_alias=True),
                "segmetrics": cfg.segmetrics.model_dump(by_alias=True),
            }
        if action in ("segment", "call", "report"):
            args["drop_low_coverage"] = cfg.drop_low_coverage
        if action in ("call", "diagram"):
            if cfg.gender != CnvkitGender.guess:
                args["gender"] = cfg.gender
            if cfg.male_reference:
                args["male_reference"] = cfg.male_reference
        return args

    def _get_params_coverage(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._cnvkit_params("coverage")

    def _get_params_fix(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._cnvkit_params("fix")

    def _get_params_segment(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._cnvkit_params("segment")

    def _get_params_call(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._cnvkit_params("call")

    def _get_params_plot(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._cnvkit_params("diagram")

    def _get_params_report(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._cnvkit_params("report")

    def _get_input_files_coverage(self, wildcards):
        # BAM/BAI file
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        return {"bam": alignments.bam, "bai": alignments.bai}

    def _get_input_files_fix(self, wildcards):
        tpl_base = "{library_name}"
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".{target}coverage.cnn"
        input_files = {
            "target": tpl.format(target="target", **wildcards),
            "antitarget": tpl.format(target="antitarget", **wildcards),
        }
        return input_files

    def _get_input_files_segment(self, wildcards):
        cnr_pattern = "work/{library_name}/out/{library_name}.cnr"
        input_files = {"cnr": cnr_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_call(self, wildcards):
        segment_pattern = "work/{library_name}/out/{library_name}.segment.cns"
        input_files = {"segment": segment_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_postprocess(self, wildcards):
        segment_pattern = "work/{library_name}/out/{library_name}.call.cns"
        input_files = {"call": segment_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_export(self, wildcards):
        cns_pattern = "work/{library_name}/out/{library_name}.call.cns"
        input_files = {"cns": cns_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_plot(self, wildcards):
        tpl = "work/{library_name}/out/{library_name}.{ext}"
        input_files = {
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def _get_input_files_report(self, wildcards):
        tpl = "work/{library_name}/out/{library_name}.{ext}"
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
        name_pattern = "{library_name}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
            output_files[target + "_md5"] = output_files[target] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_fix():
        name_pattern = "{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cnr")
        return {"ratios": tpl, "ratios_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_segment():
        name_pattern = "{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".segment.cns")
        return {"segments": tpl, "segments_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_call():
        name_pattern = "{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".call.cns")
        return {"calls": tpl, "calls_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_postprocess():
        name_pattern = "{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cns")
        return {"final": tpl, "final_md5": tpl + ".md5"}

    @dictify
    def _get_output_files_plot(self):
        plots = (("diagram", "pdf"), ("heatmap", "pdf"), ("scatter", "png"))
        chrom_plots = (("heatmap", "pdf"), ("scatter", "png"))
        chroms = list(chain(range(1, 23), ["X", "Y"]))
        output_files = {}
        # Yield file name pairs for global plots
        name_pattern = "{library_name}"
        for plot, ext in plots:
            output_files[plot] = f"work/{name_pattern}/report/{name_pattern}.{plot}.{ext}"
            output_files[plot + "_md5"] = output_files[plot] + ".md5"
        # Yield file name pairs for the chromosome-wise plots
        for plot, ext in chrom_plots:
            for chrom in chroms:
                key = "{plot}_chr{chrom}".format(plot=plot, chrom=chrom)
                output_files[key] = (
                    f"work/{name_pattern}/report/{name_pattern}.{plot}.chr{chrom}.{ext}"
                )
                output_files[key + "_md5"] = output_files[key] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_export():
        exports = (
            ("bed", ".bed"),
            ("seg", ".dnacopy.seg"),
            ("vcf", ".vcf.gz"),
            ("tbi", ".vcf.gz.tbi"),
        )
        output_files = {}
        name_pattern = "{library_name}"
        for export, suffix in exports:
            output_files[export] = f"work/{name_pattern}/out/{name_pattern}{suffix}"
            output_files[export + "_md5"] = output_files[export] + ".md5"
        return output_files

    @dictify
    def _get_output_files_report(self):
        reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        output_files = {}
        name_pattern = "{library_name}"
        for report in reports:
            output_files[report] = f"work/{name_pattern}/report/{name_pattern}.{report}.txt"
            output_files[report + "_md5"] = output_files[report] + ".md5"
        return output_files

    def get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)
        prefix = f"work/{{library_name}}/log/{action}.{{library_name}}"
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
            runtime="40h",  # 1 day and 16 hours
            mem=f"{2 * 30 * 1024}MB",
        ),
        "transform": ResourceUsage(
            threads=1,
            runtime="40h",  # 1 day and 16 hours
            mem=f"{2 * 8 * 1024}MB",
        ),
        "run": ResourceUsage(
            threads=8,
            runtime="40h",  # 1 day and 16 hours
            mem=f"{int(2 * 3.75 * 1024 * 8)}MB",
        ),
    }

    def get_output_files(self, action):
        # Initialise variable
        result = {}

        # Validate action
        self._validate_action(action)

        if action == "run":
            result["ratio"] = self.base_path_out.format(ext=".ratio.txt")
            result["ratio_md5"] = self.base_path_out.format(ext=".ratio.txt.md5")
        elif action == "transform":
            transform_ext_names = ("log2", "call", "segments", "cns", "cnr")
            transform_ext_values = (
                ".gene_log2.txt",
                ".gene_call.txt",
                ".segments.txt",
                ".cns",
                ".cnr",
            )
            result = dict(
                zip(
                    transform_ext_names,
                    expand(self.base_path_out, ext=transform_ext_values),
                )
            )
        elif action == "plot":
            plot_ext_names = ("heatmap", "scatter", "diagram")
            plot_ext_values = (".heatmap.png", ".scatter.png", ".diagram.pdf")
            result = dict(
                zip(
                    plot_ext_names,
                    expand(self.base_path_out, ext=plot_ext_values),
                )
            )

        return result

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "chrlenfile", self.config.control_freec.path_chrlenfile
        if self.config.control_freec.path_mappability_enabled:
            yield "mappability", self.config.control_freec.path_mappability

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.control_freec
        return {
            "path_mappability_enabled": cfg.path_mappability_enabled,
            "window_size": cfg.window_size,
        }

    def _get_params_transform(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.control_freec
        return {
            "org_obj": cfg.convert.org_obj,
            "tx_obj": cfg.convert.tx_obj,
            "bs_obj": cfg.convert.bs_obj,
            "tumor_library": wildcards.tumor_library,
        }

    def _get_params_plot(self, wildcards: Wildcards) -> dict[str, Any]:
        return self.config.control_freec.diagram.model_dump(by_alias=True)

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
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"}))]

    config_model_class = SomaticWgsCnvCallingConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local somatic WGS CNV output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"dnacopy_seg": f"output/{lib}/out/{lib}.dnacopy.seg"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        selected_tool = self.config.tool
        match selected_tool:
            case Tool.cnvetti:
                selected_sub_step = CnvettiSomaticWgsStepPart
            case Tool.cnvkit:
                selected_sub_step = CnvkitSomaticWgsStepPart
            case Tool.control_freec:
                selected_sub_step = ControlFreecSomaticWgsStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {selected_tool}")
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                selected_sub_step,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow"""
        tool = self.config.tool
        for entity in self.output_entities:
            if tool == "cnvetti":
                yield from expand(
                    os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                    tumor_library=[entity],
                    ext=BCF_EXT_VALUES,
                )
            elif tool == "control_freec":
                yield from expand(
                    os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                    tumor_library=[entity],
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
            elif tool == "cnvkit":
                exts = (".cnr", ".cns", ".bed", ".dnacopy.seg", ".vcf.gz", ".vcf.gz.tbi")
                yield from expand(
                    os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                    tumor_library=[entity],
                    ext=exts,
                )
                yield from expand(
                    os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                    tumor_library=[entity],
                    ext=[ext + ".md5" for ext in exts],
                )
                reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
                for report in reports:
                    yield from expand(
                        os.path.join(
                            "output", "{tumor_library}", "report", "{tumor_library}.{ext}"
                        ),
                        tumor_library=[entity],
                        ext=[f"{report}.txt", f"{report}.txt.md5"],
                    )
                for plot, ext, chrom in (
                    ("diagram", "pdf", False),
                    ("heatmap", "pdf", True),
                    ("scatter", "png", True),
                ):
                    yield from expand(
                        os.path.join(
                            "output", "{tumor_library}", "report", "{tumor_library}.{ext}"
                        ),
                        tumor_library=[entity],
                        ext=[f"{plot}.{ext}", f"{plot}.{ext}.md5"],
                    )
                    if chrom:
                        for c in map(str, chain(range(1, 23), ("X", "Y"))):
                            yield from expand(
                                os.path.join(
                                    "output", "{tumor_library}", "report", "{tumor_library}.{ext}"
                                ),
                                tumor_library=[entity],
                                ext=[
                                    f"{plot}.chr{c}.{ext}",
                                    f"{plot}.chr{c}.{ext}.md5",
                                ],
                            )
            else:
                yield from expand(
                    os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                    tumor_library=[entity],
                    ext=EXT_VALUES,
                )
