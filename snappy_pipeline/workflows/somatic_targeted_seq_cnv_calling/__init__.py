# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_target_seq_cnv_calling`` step

This step allows for the detection of CNV events for cancer samples from targeted sequenced (e.g.,
exomes or large panels).  The wrapped tools start from the aligned reads (thus off ``ngs_mapping``)
and generate CNV calls for somatic variants.

The wrapped tools implement different strategies.  Some work "reference free" and just use the
somatic BAM files for their input, some work in "matched cancer normal mode" and need the cancer
and normal BAM files, others again need both normal and cancer BAM files, and additionally a
set of non-cancer BAM files for their background.

==========
Step Input
==========

Gene somatic CNV calling for targeted sequencing starts off the aligned reads, i.e.,
``ngs_mapping``.

===========
Step Output
===========

There is no widely used standard to report copy number alterations.
In absence of a better solution, all CNV tools implemented in somatic pipeline output the segmentation table loosely following the `DNAcopy format <https://bioconductor.org/packages/devel/bioc/manuals/DNAcopy/man/DNAcopy.pdf>`_.`
The copy number call may or may not be present, and the chromosome number is replaced by its name.
The segmentation output is in file ``output/<mapper>.<cnv caller>.<lib name>/out/<mapper>.<cnv caller>.<lib name>_dnacopy.seg``.

::

    output/
    +-- bwa.cnvkit.P001-N1-DNA1-WES1
    |   |-- out
    |   |   |-- bwa.cnvkitP001-N1-DNA1-WES1_dnacopy.seg
            [...]

Note that tool ``cnvetti`` doesn't follow the snappy convention above:
the tool name is followed by an underscore & the action, where the action is one of ``coverage``, ``segment`` and ``postprocess``.
For example, the output directory would contain a directory named ``bwa.cnvetti_coverage.P002-T1-DNA1-WES1``.

.. note:: Tool-Specific Output

    Each tool produces its own set of outputs, generally not in standard format.
    Some of these files are linked from ``work`` to ``output``, but not necessarily all of them.
    Some tools (for example ``cnvkit``) also produces a report, with tables and figures.


=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_targeted_seq_cnv_calling.rst

=====================================
Available Somatic Targeted CNV Caller
=====================================

- ``cnvkit``
- ``sequenza``
- ``purecn``. Note that ``purecn`` requires a panel of normals and a second set of variants called by ``mutect2``, that includes germline ones.

"""

import os
import os.path
import re
from itertools import chain
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import UnsupportedActionException
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
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import Cnvkit as CnvkitModel
from .model import SequenzaExtraArgs, SequenzaExtractExtraArgs, SequenzaFitExtraArgs
from .model import SomaticTargetedSeqCnvCalling as SomaticTargetedSeqCnvCallingConfigModel, Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the somatic_targeted_seq_cnv_calling step

#: JSON key for "isCancer"
KEY_IS_CANCER = "isCancer"

#: Value for "libraryType" is whole exome sequencing
VALUE_WES = "WES"

#: Value for "libraryType" is panel sequencing
VALUE_PANEL = "Panel-seq"

#: Values for targeted sequencing
VALUES_TARGETED_SEQ = (VALUE_WES, VALUE_PANEL)

#: Standard key/extension values for BCF files
BCF_KEY_EXTS = (
    ("bcf", ".bcf"),
    ("bcf_md5", ".bcf.md5"),
    ("bcf_csi", ".bcf.csi"),
    ("bcf_csi_md5", ".bcf.csi.md5"),
)


class SomaticTargetedSeqCnvCallingStepPart(BaseStepPart):
    """Shared code for all caller classes in somatic_targeted_seq_cnv_calling"""

    def __init__(self, parent):
        super().__init__(parent)

    def _resolve_library_name(self, library_name: str) -> str:
        df = self.parent.build_library_dataframe()
        if library_name in df["library_name"].values:
            return library_name
        if "." in library_name:
            unprefixed_name = library_name.split(".")[-1]
            if unprefixed_name in df["library_name"].values:
                return unprefixed_name
        return library_name

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        df = self.parent.build_library_dataframe()
        library_name = self._resolve_library_name(wildcards.tumor_library)
        tumor_df = df[df["library_name"] == library_name]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    @staticmethod
    @dictify
    def _get_log_file_from_prefix(prefix):
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


def expand_id(*args):
    """Returns a dict that can be passed into expand to get identity for the given values

    ::

        >> expand_id(('foo', 'bar'))
        {'foo': ['{foo}'], 'bar': ['{bar}']}

    """
    return {key: ["{{{}}}".format(key)] for key in args}


def format_id(*args):
    """Returns a dict that can be passed into format to get identity for the given values

    ::

        >> expand_id(('foo', 'bar'))
        {'foo': '{foo}', 'bar': '{bar}'}

    """
    return {key: "{{{}}}".format(key) for key in args}


class SequenzaStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using sequenza"""

    #: Step name
    name = "sequenza"

    #: Class available actions
    actions = (
        "install",
        "gcreference",
        "coverage",
        "run",
    )

    resource_usage = {
        "coverage": ResourceUsage(
            threads=1,
            runtime="24h",
            mem="24GB",
        ),
        "run": ResourceUsage(
            threads=4,
            runtime="24h",
            mem="64GB",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)

        method_mapping = {
            "coverage": self._get_input_files_coverage(),
            "run": self._get_input_files_run(),
        }
        return method_mapping[action]

    def _get_input_files_coverage(self):
        @dictify
        def input_function(wildcards):
            ngs_mapping = self.parent.upstream("ngs_mapping")
            tumor_library = self._resolve_library_name(wildcards.tumor_library)
            normal_base_path = "output/{normal_library}/out/{normal_library}".format(
                normal_library=self.get_normal_lib_name(wildcards), **wildcards
            )
            tumor_base_path = f"output/{tumor_library}/out/{tumor_library}"
            yield (
                "gc",
                "work/static_data/out/sequenza.{length}.wig.gz".format(
                    length=self.config.sequenza.length,
                ),
            )
            yield "normal_bam", ngs_mapping(normal_base_path + ".bam")
            yield "normal_bai", ngs_mapping(normal_base_path + ".bam.bai")
            yield "tumor_bam", ngs_mapping(tumor_base_path + ".bam")
            yield "tumor_bai", ngs_mapping(tumor_base_path + ".bam.bai")

        return input_function

    def _get_input_files_run(self):
        @dictify
        def input_function(wildcards):
            yield "packages", "work/R_packages/out/sequenza.done"
            name_pattern = "{tumor_library}"
            yield "seqz", f"work/{name_pattern}/out/{name_pattern}.seqz.gz"

        return input_function

    def get_output_files(self, action):
        if action == "install":
            return {"done": "work/R_packages/out/sequenza.done"}
        elif action == "gcreference":
            return {
                "gc": "work/static_data/out/sequenza.{length}.wig.gz".format(
                    length=self.config.sequenza.length,
                )
            }
        elif action == "coverage":
            name_pattern = "{tumor_library}"
            return {
                "seqz": f"work/{name_pattern}/out/{name_pattern}.seqz.gz",
                "seqz_md5": f"work/{name_pattern}/out/{name_pattern}.seqz.gz.md5",
            }
        elif action == "run":
            name_pattern = "{tumor_library}"
            return {
                "seg": f"work/{name_pattern}/out/{name_pattern}_dnacopy.seg",
                "seg_md5": f"work/{name_pattern}/out/{name_pattern}_dnacopy.seg.md5",
                "done": f"work/{name_pattern}/report/.done",
            }
        else:
            raise UnsupportedActionException(
                "Action '{action}' is not supported. Valid options: {valid}".format(
                    action=action, valid=", ".join(self.actions)
                )
            )

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    @staticmethod
    def _coerce_model(model_cls, value):
        if isinstance(value, model_cls):
            return value
        return model_cls.model_validate(value or {})

    def _get_args_coverage(self, wildcards: Wildcards) -> dict[str, Any]:
        extra_args = self._coerce_model(SequenzaExtraArgs, self.config.sequenza.extra_args)
        return {
            "reference": self.parent.w_config.static_data_config.reference.path,
            "length": self.config.sequenza.length,
            "ignore_chroms": self.config.sequenza.ignore_chroms,
            "extra_arguments": extra_args.model_dump(by_alias=True),
        }

    def _get_args_gcreference(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "reference": self.parent.w_config.static_data_config.reference.path,
            "length": self.config.sequenza.length,
        }

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        extra_args_extract = self._coerce_model(
            SequenzaExtractExtraArgs, self.config.sequenza.extra_args_extract
        )
        extra_args_fit = self._coerce_model(
            SequenzaFitExtraArgs, self.config.sequenza.extra_args_fit
        )
        return {
            "reference": self.parent.w_config.static_data_config.reference.path,
            "assembly": self.config.sequenza.assembly,
            "ignore_chroms": self.config.sequenza.ignore_chroms,
            "extra_args_extract": extra_args_extract.model_dump(by_alias=True),
            "extra_args_fit": extra_args_fit.model_dump(by_alias=True),
            "library_name": wildcards.tumor_library,
        }

    def get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)
        if action == "install":
            prefix = "work/R_packages/log/sequenza"
        elif action == "gcreference":
            prefix = "work/static_data/log/sequenza.{length}".format(
                length=self.config.sequenza.length,
            )
        else:
            name_pattern = "{tumor_library}"
            prefix = os.path.join("work", name_pattern, "log", name_pattern + "." + action)
        return self._get_log_file_from_prefix(prefix)


class PureCNStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using PureCN"""

    #: Step name
    name = "purecn"

    #: Class available actions
    actions = ("coverage", "run")

    resource_usage = {
        "coverage": ResourceUsage(
            threads=1,
            runtime="4h",
            mem="24GB",
        ),
        "run": ResourceUsage(
            threads=4,
            runtime="24h",
            mem="96GB",
        ),
    }

    def get_input_files(self, action):
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        action_mapping = {
            "coverage": self._get_input_files_coverage,
            "run": self._get_input_files_run,
        }
        return action_mapping[action]

    @dictify
    def _get_input_files_run(self, wildcards):
        name_pattern = "{tumor_library}".format(**wildcards)
        yield (
            "tumor",
            os.path.join(
                "work",
                name_pattern,
                "out",
                name_pattern + "_coverage_loess.txt.gz",
            ).format(**wildcards),
        )
        pon = self.parent.upstream("panel_of_normals")
        somatic_vcf = self.parent.get_upstream_paths(
            "somatic_variants", library_name=wildcards.tumor_library
        )
        yield "vcf", getattr(somatic_vcf, "full_vcf", None) or getattr(somatic_vcf, "vcf", None)
        purecn_cfg = self.config.purecn
        yield "normaldb", pon("output/purecn/out/purecn.panel_of_normals.rds")
        yield "mapping_bias", pon("output/purecn/out/purecn.mapping_bias.rds")
        yield (
            "intervals",
            pon(
                f"output/purecn/out/{purecn_cfg.enrichment_kit_name}_{purecn_cfg.genome_name}.list"
            ),
        )

    @dictify
    def _get_input_files_coverage(self, wildcards):
        ngs_mapping = self.parent.upstream("ngs_mapping")
        pon = self.parent.upstream("panel_of_normals")
        name_pattern = "{tumor_library}".format(**wildcards)
        base_path = os.path.join("output", name_pattern, "out", name_pattern)
        yield "bam", ngs_mapping(base_path + ".bam")
        yield "bai", ngs_mapping(base_path + ".bam.bai")
        purecn_cfg = self.config.purecn
        yield (
            "intervals",
            pon(
                f"output/purecn/out/{purecn_cfg.enrichment_kit_name}_{purecn_cfg.genome_name}.list"
            ),
        )

    def get_output_files(self, action):
        """Return output paths, dependent on rule"""
        # Validate action
        self._validate_action(action)
        name_pattern = "{tumor_library}"
        prefix = os.path.join("work", name_pattern, "out", name_pattern)
        action_mapping = {
            "coverage": {"coverage": prefix + "_coverage_loess.txt.gz"},
            "run": {
                "segments": prefix + "_dnacopy.seg",
                "ploidy": prefix + ".csv",
                "pvalues": prefix + "_amplification_pvalues.csv",
                "vcf": prefix + ".vcf.gz",
                "vcf_tbi": prefix + ".vcf.gz.tbi",
                "loh": prefix + "_loh.csv",
                "segments_md5": prefix + "_dnacopy.seg.md5",
                "ploidy_md5": prefix + ".csv.md5",
                "pvalues_md5": prefix + "_amplification_pvalues.csv.md5",
                "vcf_md5": prefix + ".vcf.gz.md5",
                "vcf_tbi_md5": prefix + ".vcf.gz.tbi.md5",
                "loh_md5": prefix + "_loh.csv.md5",
            },
        }
        return action_mapping[action]

    def get_args(self, action):
        self._validate_action(action)
        return self._get_args_all

    def _get_args_all(self, wildcards):
        mapper = str(self.parent.get_task_config("ngs_mapping").tool)
        config_dump = self.config.get(self.name).model_dump(by_alias=True)
        # Inject PON file paths resolved from the panel_of_normals dependency so that
        # the wrapper can access them via config["path_*"] as before.
        pon = self.parent.upstream("panel_of_normals")
        purecn_cfg = self.config.purecn
        config_dump["path_panel_of_normals"] = pon("output/purecn/out/purecn.panel_of_normals.rds")
        config_dump["path_mapping_bias"] = pon("output/purecn/out/purecn.mapping_bias.rds")
        config_dump["path_intervals"] = pon(
            f"output/purecn/out/{purecn_cfg.enrichment_kit_name}_{purecn_cfg.genome_name}.list"
        )
        return {
            "config": config_dump,
            "mapper": mapper,
            "library_name": wildcards.tumor_library,
        }

    def get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)

        name_pattern = "{tumor_library}"
        prefix = os.path.join("work", name_pattern, "log", name_pattern + "." + action)
        return self._get_log_file_from_prefix(prefix)


class CnvKitStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using cnvkit"""

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

    default_resource_usage = ResourceUsage(threads=1, runtime="4h", mem="7680MB")

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage = {
        "plot": ResourceUsage(
            threads=1,
            runtime="8h",
            mem=f"{30 * 1024}MB",
        ),
        "coverage": ResourceUsage(
            threads=8,
            runtime="8h",
            mem=f"{16 * 1024}MB",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: CnvkitModel = self.config.get(self.name)

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
        ngs_mapping = self.parent.upstream("ngs_mapping")
        base_path = "output/{tumor_library}/out/{tumor_library}".format(**wildcards)
        return {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam.bai"),
            "reference": self.w_config.static_data_config.reference.path,
            "target": self.config.cnvkit.path_target,
            "antitarget": self.config.cnvkit.path_antitarget,
        }

    def _get_input_files_fix(self, wildcards):
        tpl_base = "{tumor_library}"
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".{target}coverage.cnn"
        return {
            "target": tpl.format(target="target", **wildcards),
            "antitarget": tpl.format(target="antitarget", **wildcards),
            "ref": self.parent.upstream("panel_of_normals")(
                "output/cnvkit/out/cnvkit.panel_of_normals.cnn"
            ),
        }

    def _get_input_files_segment(self, wildcards):
        cnr_pattern = "work/{tumor_library}/out/{tumor_library}.cnr"
        input_files = {"cnr": cnr_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_call(self, wildcards):
        segment_pattern = "work/{tumor_library}/out/{tumor_library}.segment.cns"
        input_files = {"segment": segment_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_postprocess(self, wildcards):
        segment_pattern = "work/{tumor_library}/out/{tumor_library}.segment.cns"
        call_pattern = "work/{tumor_library}/out/{tumor_library}.call.cns"
        input_files = {
            "segment": segment_pattern.format(**wildcards),
            "call": call_pattern.format(**wildcards),
        }
        return input_files

    def _get_input_files_export(self, wildcards):
        cns_pattern = "work/{tumor_library}/out/{tumor_library}.call.cns"
        input_files = {"cns": cns_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_plot(self, wildcards):
        tpl = "work/{tumor_library}/out/{tumor_library}.{ext}"
        input_files = {
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def _get_input_files_report(self, wildcards):
        tpl = "work/{tumor_library}/out/{tumor_library}.{ext}"
        input_files = {
            "target": tpl.format(ext="targetcoverage.cnn", **wildcards),
            "antitarget": tpl.format(ext="antitargetcoverage.cnn", **wildcards),
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def get_args(self, action):
        self._validate_action(action)
        if action == "plot":
            action = "diagram"
        if args := getattr(self.cfg, action, {}):
            args = args.model_dump(by_alias=True)
        if action == "report":
            args["breaks"] = self.cfg.breaks.model_dump(by_alias=True)
            args["genemetrics"] = self.cfg.genemetrics.model_dump(by_alias=True)
            args["segmetrics"] = self.cfg.segmetrics.model_dump(by_alias=True)
        if action in ("segment", "call", "report"):
            args["drop_low_coverage"] = self.cfg.drop_low_coverage
        if action in ("call", "diagram") and self.cfg.gender != CnvkitGender.guess:
            action["gender"] = self.cfg.gender
        if action in ("call", "diagram") and self.cfg.male_reference:
            action["male_reference"] = self.cfg.male_reference
        return args

    def get_output_files(self, action):
        """Return output files for the given action"""
        if action == "coverage":
            return self._get_output_files_coverage()
        elif action == "fix":
            return self._get_output_files_fix()
        elif action == "segment":
            return self._get_output_files_segment()
        elif action == "call":
            return self._get_output_files_call()
        elif action == "postprocess":
            return self._get_output_files_postprocess()
        elif action == "export":
            return self._get_output_files_export()
        elif action == "plot":
            return self._get_output_files_plot()
        elif action == "report":
            return self._get_output_files_report()
        else:
            self._validate_action(action)

    @staticmethod
    def _get_output_files_coverage():
        name_pattern = "{tumor_library}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
            output_files[target + "_md5"] = output_files[target] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_fix():
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cnr")
        return {"ratios": tpl, "ratios_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_segment():
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".segment.cns")
        return {"segments": tpl, "segments_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_call():
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".call.cns")
        return {"calls": tpl, "calls_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_postprocess():
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + "_dnacopy.seg")
        return {
            "final": tpl,
            "final_md5": tpl + ".md5",
        }

    @dictify
    def _get_output_files_plot(self):
        plots = (("diagram", "pdf"), ("scatter", "png"))
        chrom_plots = (("scatter", "png"),)
        chroms = list(chain(range(1, 23), ["X", "Y"]))
        output_files = {}
        name_pattern = "{tumor_library}"
        # Yield file name pairs for global plots
        for plot, ext in plots:
            tpl = os.path.join("work", name_pattern, "report", name_pattern + f".{plot}.{ext}")
            output_files[plot] = tpl
            output_files[plot + "_md5"] = tpl + ".md5"
        # Yield file name pairs for the chromosome-wise plots
        for plot, ext in chrom_plots:
            for chrom in chroms:
                key = f"{plot}_chr{chrom}"
                tpl = os.path.join(
                    "work", name_pattern, "report", name_pattern + f".{plot}.chr{chrom}.{ext}"
                )
                output_files[key] = tpl
                output_files[key + "_md5"] = tpl + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_export():
        exports = (
            ("bed", "bed.gz"),
            ("bed_tbi", "bed.gz.tbi"),
            ("seg", "seg"),
            ("vcf", "vcf.gz"),
            ("vcf_tbi", "vcf.gz.tbi"),
        )
        output_files = {}
        name_pattern = "{tumor_library}"
        for export, ext in exports:
            tpl = os.path.join("work", name_pattern, "out", name_pattern + f".{ext}")
            output_files[export] = tpl
            output_files[export + "_md5"] = tpl + ".md5"
        return output_files

    @dictify
    def _get_output_files_report(self):
        reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        output_files = {}
        name_pattern = "{tumor_library}"
        for report in reports:
            tpl = os.path.join("work", name_pattern, "report", name_pattern + f".{report}.txt")
            output_files[report] = tpl
            output_files[report + "_md5"] = tpl + ".md5"
        return output_files

    def get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)
        prefix = f"work/{{tumor_library}}/log/{action}.{{tumor_library}}"
        return self._get_log_file_from_prefix(prefix)


class SomaticTargetedSeqCnvCallingWorkflow(BaseStep):
    """Perform somatic targeted sequencing CNV calling"""

    #: Workflow name
    name = "somatic_targeted_seq_cnv_calling"
    consumes = {DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})): True}
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"}))]

    config_model_class = SomaticTargetedSeqCnvCallingConfigModel

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
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local somatic targeted CNV output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{tumor_library}")
        return {"done": f"output/{lib}/out/.done"}

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
        selected_tool = self.config.tool
        match selected_tool:
            case Tool.cnvkit:
                selected_sub_step = CnvKitStepPart
            case Tool.sequenza:
                selected_sub_step = SequenzaStepPart
            case Tool.purecn:
                selected_sub_step = PureCNStepPart
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
        """Return list of result files for the somatic targeted sequencing CNV calling step"""
        tool = self.config.tool
        sub_step = self.sub_steps[tool]
        tool_actions = {
            Tool.cnvkit: ("fix", "postprocess", "report", "plot", "export"),
            Tool.sequenza: ("coverage", "run"),
            Tool.purecn: ("run",),
        }
        for action in tool_actions[tool]:
            output_files = sub_step.get_output_files(action)
            paths = (
                list(output_files.values()) if isinstance(output_files, dict) else [output_files]
            )
            for p in paths:
                if isinstance(p, str) and p.startswith("work/"):
                    out_p = re.sub(r"^work/", "output/", p)
                    yield from expand(out_p, tumor_library=self.output_entities)

            log_files = sub_step.get_log_file(action)
            log_paths = (
                list(log_files.values())
                if isinstance(log_files, dict)
                else [log_files]
                if isinstance(log_files, str)
                else []
            )
            for lp in log_paths:
                if isinstance(lp, str) and lp.startswith("work/"):
                    out_lp = re.sub(r"^work/", "output/", lp)
                    yield from expand(out_lp, tumor_library=self.output_entities)
