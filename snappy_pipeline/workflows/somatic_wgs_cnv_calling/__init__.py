# -*- coding: utf-8 -*-
"""Somatic CNV calling for WGS data

The step detects copy number alterations in cancer whole genomes from the aligned reads of each tumor and its matched normal.
The ``tool`` field selects the caller:

- ``cnvetti`` writes the called CNVs as a BCF file.
- ``cnvkit`` runs the CNVkit workflow for whole genomes and exports the segments in several formats.
- ``control_freec`` runs Control-FREEC and transforms its ratios into log2 values, calls and segments.
  It needs ``depends_on.database``: an ``external_data`` task with ``{type: database, tags: [control_freec]}``, the key ``chr_len`` (the chromosome lengths) and, optionally, the key ``mappability``.

The task gives one result per tumor DNA library.
The files are unfiltered.
They contain low-quality calls.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
The files depend on the tool and are generally not in a standard format.

- ``cnvetti``: ``{tumor_library}/out/{tumor_library}.bcf`` with ``.csi``.
- ``cnvkit``: ``{tumor_library}/out/{tumor_library}.cnr``, ``.cns``, ``.bed``, ``.dnacopy.seg`` and ``.vcf.gz`` (with ``.tbi``).
  The directory ``{tumor_library}/report/`` holds the tables ``{tumor_library}.breaks.txt``, ``.genemetrics.txt``, ``.segmetrics.txt``, ``.sex.txt`` and ``.metrics.txt``, and the plots ``{tumor_library}.diagram.pdf``, ``{tumor_library}.scatter.png`` and ``{tumor_library}.heatmap.pdf`` (also per chromosome, ``{tumor_library}.heatmap.chr1.pdf`` and so on).
- ``control_freec``: ``{tumor_library}/out/{tumor_library}.ratio.txt``, ``.gene_log2.txt``, ``.gene_call.txt``, ``.segments.txt``, ``.heatmap.png``, ``.scatter.png`` and ``.diagram.pdf``.

Downstream tasks read the key ``dnacopy_seg``, which only ``cnvkit`` provides, so ``cnvetti`` and ``control_freec`` tasks do not satisfy their requirement.
``cnvetti`` also provides the keys ``bcf`` and ``bcf_csi``, and ``control_freec`` the keys ``log2``, ``call`` and ``segments``.
"""

import os
from itertools import chain
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.models.cnvkit import Gender as CnvkitGender
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import SomaticWgsCnvCalling as SomaticWgsCnvCallingConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class SomaticWgsCnvCallingStepPart(BaseStepPart):
    """Base class for somatic WGS CNV calling steps

    WGS CNV calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    # TODO: unify with somatic (small) variant calling base class?

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


class CnvettiSomaticWgsStepPart(SomaticWgsCnvCallingStepPart):
    """Somatic WGS CNV calling with CNVetti"""

    #: Step name
    name = "cnvetti"

    #: Class available actions
    actions = ("coverage", "tumor_normal_ratio", "segment")

    extra_arguments_action = "coverage"

    entity_wildcard = "tumor_library"

    published_keys = {"segment": ("bcf", "bcf_csi", "bcf_md5", "bcf_csi_md5")}

    default_resource_usage = ResourceUsage(threads=4, runtime="4h", mem=f"{int(3.75 * 1024 * 4)}MB")

    intermediate_outputs = {
        "coverage": ("bcf", "bcf_csi"),
        "tumor_normal_ratio": ("bcf", "bcf_csi"),
    }

    #: Extension file dictionary. Key: file type (string); Value: file extension (string)
    bcf_dict = {"bcf": ".bcf", "bcf_csi": ".bcf.csi"}

    #: The BCF files of the last action are published, so they get md5 files
    bcf_md5_dict = {"bcf_md5": ".bcf.md5", "bcf_csi_md5": ".bcf.csi.md5"}

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
        yield from self.parent.reference_inputs().items()

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
            yield (
                key + "_csi",
                f"work/{library_name}/out/{library_name}.coverage.bcf.csi",
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
        for key, ext in (self.bcf_dict | self.bcf_md5_dict).items():
            name_pattern = "{tumor_library}"
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    def _get_params_coverage(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._get_params("coverage")

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

    def _get_log_file(self, action):
        """Return path to log file"""
        wildcard_name = (
            "library_name" if action in {"coverage", "tumor_normal_ratio"} else "tumor_library"
        )
        name_pattern = f"{{{wildcard_name}}}"
        prefix = f"work/{name_pattern}/log/{name_pattern}.{action}"
        return {"log": prefix + ".log"}


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

    extra_arguments_action = "segment"

    intermediate_outputs = {
        "coverage": ("target", "antitarget"),
        "segment": ("segments",),
        "call": ("calls",),
    }

    resource_usage = {
        "plot": ResourceUsage(threads=1, runtime="1d", mem=f"{30 * 1024}MB"),
        "coverage": ResourceUsage(threads=8, runtime="1d", mem=f"{16 * 1024}MB"),
    }

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem=f"{int(7.5 * 1024)}MB")

    @property
    def published_keys(self):
        return {
            "fix": ("ratios", "ratios_md5"),
            "postprocess": ("final", "final_md5"),
            "export": tuple(self._get_output_files_export()),
            "plot": tuple(self._get_output_files_plot()),
            "report": tuple(self._get_output_files_report()),
        }

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

    def _get_params_postprocess(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"library_name": wildcards.library_name}

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

    @staticmethod
    def _get_output_files_coverage():
        name_pattern = "{library_name}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
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
        return {"segments": tpl}

    @staticmethod
    def _get_output_files_call():
        name_pattern = "{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".call.cns")
        return {"calls": tpl}

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

    def _get_log_file(self, action):
        """Return path to log file for the given action"""
        return {"log": f"work/{{library_name}}/log/{action}.{{library_name}}.log"}


class ControlFreecSomaticWgsStepPart(SomaticWgsCnvCallingStepPart):
    """Somatic WGS CNV calling with Control-FreeC"""

    #: Step name
    name = "control_freec"

    #: Class available actions
    actions = ("run", "transform", "plot")

    extra_arguments_action = "run"

    intermediate_outputs = {"transform": ("cns", "cnr")}

    entity_wildcard = "tumor_library"

    published_keys = {
        "run": ("ratio", "ratio_md5"),
        "transform": ("log2", "call", "segments", "log2_md5", "call_md5", "segments_md5"),
        "plot": ("heatmap", "scatter", "diagram", "heatmap_md5", "scatter_md5", "diagram_md5"),
    }

    resource_usage = {
        "plot": ResourceUsage(threads=1, runtime="40h", mem=f"{2 * 30 * 1024}MB"),
        "transform": ResourceUsage(threads=1, runtime="40h", mem=f"{2 * 8 * 1024}MB"),
        "run": ResourceUsage(threads=8, runtime="40h", mem=f"{int(2 * 3.75 * 1024 * 8)}MB"),
    }

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}

    def _get_output_files_run(self):
        base = "work/{tumor_library}/out/{tumor_library}"
        return {"ratio": base + ".ratio.txt", "ratio_md5": base + ".ratio.txt.md5"}

    def _get_output_files_transform(self):
        base = "work/{tumor_library}/out/{tumor_library}"
        result = {
            "log2": base + ".gene_log2.txt",
            "call": base + ".gene_call.txt",
            "segments": base + ".segments.txt",
            "cns": base + ".cns",
            "cnr": base + ".cnr",
        }
        for name in ("log2", "call", "segments"):  # published
            result[name + "_md5"] = result[name] + ".md5"
        return result

    def _get_output_files_plot(self):
        base = "work/{tumor_library}/out/{tumor_library}"
        result = {
            "heatmap": base + ".heatmap.png",
            "scatter": base + ".scatter.png",
            "diagram": base + ".diagram.pdf",
        }
        result.update({name + "_md5": path + ".md5" for name, path in result.items()})
        return result

    def _get_input_files_transform(self, wildcards):
        return {"ratio": "work/{tumor_library}/out/{tumor_library}.ratio.txt"}

    def _get_input_files_plot(self, wildcards):
        transform = self.get_output_files("transform").items()
        return {k: v.format(**wildcards) for k, v in transform if not k.endswith("_md5")}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        database = self.parent.get_upstream_paths("database")
        yield "chrlenfile", database.chr_len
        if database.mappability:
            yield "mappability", database.mappability

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.control_freec
        return {
            "path_mappability_enabled": bool(
                self.parent.get_upstream_paths("database").mappability
            ),
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


class SomaticWgsCnvCallingWorkflow(BaseStep):
    """Perform somatic WGS CNV calling"""

    #: Workflow name
    name = "somatic_wgs_cnv_calling"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"}))]

    config_model_class = SomaticWgsCnvCallingConfigModel

    library_requirement = "kind == 'cancer'"

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
        prefix = f"work/{lib}/out/{lib}"
        match config.tool:
            case Tool.cnvkit:
                return {"dnacopy_seg": prefix + ".dnacopy.seg"}
            case Tool.cnvetti:
                return {"bcf": prefix + ".bcf", "bcf_csi": prefix + ".bcf.csi"}
            case Tool.control_freec:
                return {
                    "log2": prefix + ".gene_log2.txt",
                    "call": prefix + ".gene_call.txt",
                    "segments": prefix + ".segments.txt",
                }
        return {}

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
        self.register_sub_step_classes((selected_sub_step,))
