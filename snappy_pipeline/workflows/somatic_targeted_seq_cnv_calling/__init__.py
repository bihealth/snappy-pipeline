# -*- coding: utf-8 -*-
"""Somatic CNV calling for targeted sequencing data

The step detects copy number alterations in cancer samples from targeted sequencing (exomes or large panels).
The ``tool`` field selects the caller, which starts from the aligned reads:

- ``cnvkit`` compares the tumor coverage with a CNVkit panel of normals (``depends_on.panel_of_normals``, from a ``panel_of_normals`` task with ``tool: cnvkit``) and reads its target and antitarget regions.
- ``sequenza`` compares the tumor with its matched normal and estimates the purity and ploidy.
- ``purecn`` needs a PureCN panel of normals (``depends_on.panel_of_normals``, from a ``panel_of_normals`` task with ``tool: purecn``, which also gives the container) and somatic variants (``depends_on.variants``).
  The variants must come from Mutect2 and include the germline variants.
  See :ref:`step_panel_of_normals`.

The task gives one result per tumor DNA library.
``sequenza`` leaves out tumor libraries without a matched normal, with a warning.
All callers write the segments in a table that loosely follows the `DNAcopy format <https://bioconductor.org/packages/devel/bioc/manuals/DNAcopy/man/DNAcopy.pdf>`_.
The copy number call may be missing, and the chromosome number is replaced by its name.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
Every caller has the result file ``{tumor_library}/out/{tumor_library}.dnacopy.seg``.
The other files depend on the tool and are generally not in a standard format.

- ``cnvkit``: ``{tumor_library}/out/{tumor_library}.cnr``, ``.seg``, ``.bed.gz`` (with ``.tbi``) and ``.vcf.gz`` (with ``.tbi``).
  The directory ``{tumor_library}/report/`` holds the tables ``{tumor_library}.breaks.txt``, ``.genemetrics.txt``, ``.segmetrics.txt``, ``.sex.txt`` and ``.metrics.txt``, and the plots ``{tumor_library}.diagram.pdf`` and ``{tumor_library}.scatter.png`` (also per chromosome, ``{tumor_library}.scatter.chr1.png`` and so on).
- ``sequenza``: ``{tumor_library}/out/{tumor_library}.seqz.gz``, and the directory ``{tumor_library}/report/`` (without ``.md5`` file) with the plots and tables of ``sequenza.results``.
- ``purecn``: ``{tumor_library}/out/{tumor_library}.csv``, ``.amplification_pvalues.csv``, ``.loh.csv`` and ``.vcf.gz`` (with ``.tbi``).

Downstream tasks read the key ``dnacopy_seg``, and for ``sequenza`` also the key ``seqz``.

Reports
-------

With ``tool: cnvkit``, the plots ``{tumor_library}.diagram.pdf`` and ``{tumor_library}.scatter.png`` are part of the Snakemake report.
"""

import os
from itertools import chain
from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.models.cnvkit import Gender as CnvkitGender
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import Cnvkit as CnvkitModel
from .model import SequenzaExtraArgs, SequenzaExtractExtraArgs, SequenzaFitExtraArgs, Tool
from .model import SomaticTargetedSeqCnvCalling as SomaticTargetedSeqCnvCallingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


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

    entity_wildcard = "tumor_library"

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

    #: ``run`` is an R script without command line, ``coverage`` runs ``sequenza-utils bam2seqz``
    extra_arguments_action = "coverage"

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

    result_keys = {
        "coverage": ("seqz", "seqz_md5"),
        "run": ("seg", "seg_md5", "report"),
    }

    @dictify
    def _get_input_files_coverage(self, wildcards):
        tumor_library = self._resolve_library_name(wildcards.tumor_library)
        normal = self.parent.get_upstream_paths(
            "alignments", library_name=self.get_normal_lib_name(wildcards)
        )
        tumor = self.parent.get_upstream_paths("alignments", library_name=tumor_library)
        yield (
            "gc",
            "work/static_data/out/gc.{length}.wig.gz".format(
                length=self.config.sequenza.length,
            ),
        )
        yield "normal_bam", normal.bam
        yield "normal_bai", normal.bai
        yield "tumor_bam", tumor.bam
        yield "tumor_bai", tumor.bai
        yield from self.parent.reference_inputs().items()

    def _get_input_files_gcreference(self, wildcards):
        return self.parent.reference_inputs()

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "packages", "work/R_packages/out/packages.done"
        name_pattern = "{tumor_library}"
        yield "seqz", f"work/{name_pattern}/out/{name_pattern}.seqz.gz"
        yield from self.parent.reference_inputs().items()

    def _get_output_files_install(self):
        return {"done": "work/R_packages/out/packages.done"}

    def _get_output_files_gcreference(self):
        return {"gc": f"work/static_data/out/gc.{self.config.sequenza.length}.wig.gz"}

    def _get_output_files_coverage(self):
        name_pattern = "{tumor_library}"
        return {
            "seqz": f"work/{name_pattern}/out/{name_pattern}.seqz.gz",
            "seqz_md5": f"work/{name_pattern}/out/{name_pattern}.seqz.gz.md5",
        }

    def _get_output_files_run(self):
        name_pattern = "{tumor_library}"
        return {
            "seg": f"work/{name_pattern}/out/{name_pattern}.dnacopy.seg",
            "seg_md5": f"work/{name_pattern}/out/{name_pattern}.dnacopy.seg.md5",
            # sequenza.results writes plots and tables that the step does not name
            "report": directory(f"work/{name_pattern}/report"),
        }

    @staticmethod
    def _coerce_model(model_cls, value):
        if isinstance(value, model_cls):
            return value
        return model_cls.model_validate(value or {})

    def _get_params_coverage(self, wildcards: Wildcards) -> dict[str, Any]:
        extra_args = self._coerce_model(SequenzaExtraArgs, self.config.sequenza.extra_args)
        return {
            "length": self.config.sequenza.length,
            "ignore_contigs": self.config.sequenza.ignore_contigs,
            "bam2seqz_options": extra_args.model_dump(by_alias=True),
        }

    def _get_params_gcreference(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "length": self.config.sequenza.length,
        }

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        extra_args_extract = self._coerce_model(
            SequenzaExtractExtraArgs, self.config.sequenza.extra_args_extract
        )
        extra_args_fit = self._coerce_model(
            SequenzaFitExtraArgs, self.config.sequenza.extra_args_fit
        )
        return {
            "assembly": self.config.sequenza.assembly,
            "ignore_contigs": self.config.sequenza.ignore_contigs,
            "extra_args_extract": extra_args_extract.model_dump(by_alias=True),
            "extra_args_fit": extra_args_fit.model_dump(by_alias=True),
            "library_name": wildcards.tumor_library,
        }

    def _get_log_file_install(self):
        return {"log": "work/R_packages/log/packages.log"}

    def _get_log_file_gcreference(self):
        return {"log": f"work/static_data/log/gc.{self.config.sequenza.length}.log"}

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}


class PureCNStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using PureCN"""

    #: Step name
    name = "purecn"

    #: Class available actions
    actions = ("coverage", "run")

    extra_arguments_action = "run"

    intermediate_outputs = {"coverage": ("coverage",)}

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

    @property
    def result_keys(self):
        return {"run": tuple(self._get_output_files_run())}

    @dictify
    def _get_input_files_run(self, wildcards):
        name_pattern = "{tumor_library}".format(**wildcards)
        yield (
            "tumor",
            os.path.join(
                "work",
                name_pattern,
                "out",
                name_pattern + ".coverage_loess.txt.gz",
            ).format(**wildcards),
        )
        pon = self.parent.get_upstream_paths("panel_of_normals")
        somatic_vcf = self.parent.get_upstream_paths(
            "variants", library_name=wildcards.tumor_library
        )
        yield "vcf", somatic_vcf.vcf
        yield "normaldb", pon.panel_of_normals
        yield "mapping_bias", pon.mapping_bias
        yield "intervals", pon.intervals
        yield "container", pon.container

    @dictify
    def _get_input_files_coverage(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.tumor_library
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        pon = self.parent.get_upstream_paths("panel_of_normals")
        yield "intervals", pon.intervals
        yield "container", pon.container

    _prefix = "work/{tumor_library}/out/{tumor_library}"

    def _get_output_files_coverage(self):
        return {"coverage": self._prefix + ".coverage_loess.txt.gz"}

    def _get_output_files_run(self):
        prefix = self._prefix
        return {
            "segments": prefix + ".dnacopy.seg",
            "ploidy": prefix + ".csv",
            "pvalues": prefix + ".amplification_pvalues.csv",
            "vcf": prefix + ".vcf.gz",
            "vcf_tbi": prefix + ".vcf.gz.tbi",
            "loh": prefix + ".loh.csv",
            "segments_md5": prefix + ".dnacopy.seg.md5",
            "ploidy_md5": prefix + ".csv.md5",
            "pvalues_md5": prefix + ".amplification_pvalues.csv.md5",
            "vcf_md5": prefix + ".vcf.gz.md5",
            "vcf_tbi_md5": prefix + ".vcf.gz.tbi.md5",
            "loh_md5": prefix + ".loh.csv.md5",
        }

    def _get_params_coverage(self, wildcards):
        # The files reach the wrappers as inputs
        config_dump = self.config.get(self.name).model_dump(
            by_alias=True, exclude={"path_container"}
        )
        return {"config": config_dump}

    def _get_params_run(self, wildcards):
        return {**self._get_params_coverage(wildcards), "library_name": wildcards.tumor_library}

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}


class CnvKitStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using cnvkit"""

    #: Step name
    name = "cnvkit"

    report_outputs = {
        "plot": {
            key: ReportOutput("cnvkit", caption="cnvkit_plots.rst")
            for key in (
                "diagram",
                "scatter",
                *(f"scatter_chr{chrom}" for chrom in (*range(1, 23), "X", "Y")),
            )
        },
        "report": {
            key: ReportOutput("cnvkit", caption="cnvkit_reports.rst")
            for key in ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        },
    }

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

    extra_arguments_action = "segment"

    intermediate_outputs = {
        "coverage": ("target", "antitarget"),
        "segment": ("segments",),
        "call": ("calls",),
    }

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

    @property
    def result_keys(self):
        actions = ("fix", "postprocess", "export", "plot", "report")
        return {action: tuple(self.get_output_files(action)) for action in actions}

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: CnvkitModel = self.config.get(self.name)

    def _get_input_files_coverage(self, wildcards):
        # BAM/BAI file
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.tumor_library
        )
        pon = self.parent.get_upstream_paths("panel_of_normals")
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
            "target": pon.target,
            "antitarget": pon.antitarget,
        }

    def _get_input_files_fix(self, wildcards):
        tpl_base = "{tumor_library}"
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".{target}coverage.cnn"
        return {
            "target": tpl.format(target="target", **wildcards),
            "antitarget": tpl.format(target="antitarget", **wildcards),
            "ref": self.parent.get_upstream_paths("panel_of_normals").panel_of_normals,
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

    def get_params(self, action):
        self._validate_action(action)
        if action == "postprocess":
            return self._get_params_postprocess
        return lambda wildcards: self._cnvkit_params(action)

    def _get_params_postprocess(self, wildcards):
        return {"library_name": wildcards.tumor_library, "ploidy": self.cfg.call.ploidy}

    def _cnvkit_params(self, action):
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
            args["gender"] = self.cfg.gender
        if action in ("call", "diagram") and self.cfg.male_reference:
            args["male_reference"] = self.cfg.male_reference
        return args

    @staticmethod
    def _get_output_files_coverage():
        name_pattern = "{tumor_library}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
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
        return {"segments": tpl}

    @staticmethod
    def _get_output_files_call():
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".call.cns")
        return {"calls": tpl}

    @staticmethod
    def _get_output_files_postprocess():
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".dnacopy.seg")
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

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{action}.{{tumor_library}}.log"}


class SomaticTargetedSeqCnvCallingWorkflow(BaseStep):
    """Perform somatic targeted sequencing CNV calling"""

    #: Workflow name
    name = "somatic_targeted_seq_cnv_calling"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"}))]

    config_model_class = SomaticTargetedSeqCnvCallingConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """The somatic CNVs; the panel of normals must have been built with the task's tool."""
        pon = upstream.get("panel_of_normals")
        if config.tool in (Tool.cnvkit, Tool.purecn) and str(config.tool) not in pon.tags:
            raise ValueError(
                f"tool={config.tool} needs a panel_of_normals task with tool={config.tool}, "
                f"but depends_on.panel_of_normals is {config.depends_on.panel_of_normals!r}, "
                f"which produces {pon}"
            )
        return tuple(cls.produces)

    library_requirement = "kind == 'cancer'"

    @classmethod
    def task_library_requirement(cls, config):
        """sequenza compares each tumor with its matched normal, cnvkit and purecn with the panel
        of normals."""
        if config.tool == Tool.sequenza:
            return MATCHED_NORMAL_REQUIREMENT
        return cls.library_requirement

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local somatic targeted CNV output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        paths = {"dnacopy_seg": f"work/{lib}/out/{lib}.dnacopy.seg"}
        if config.tool == Tool.sequenza:
            paths["seqz"] = f"work/{lib}/out/{lib}.seqz.gz"
        return paths

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
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
        self.register_sub_step_classes((selected_sub_step,))
