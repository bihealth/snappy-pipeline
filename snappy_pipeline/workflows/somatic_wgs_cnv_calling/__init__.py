# -*- coding: utf-8 -*-
"""Somatic CNV calling for WGS data

The step detects copy number alterations in cancer whole genomes from the aligned reads of each tumor and its matched normal.
The ``tool`` field selects the caller:

- ``cnvetti`` writes the called CNVs as a BCF file.
- ``cnvkit`` runs the CNVkit workflow for whole genomes and exports the segments in several formats.
  It needs ``depends_on.panel_of_normals``: a ``panel_of_normals`` task with ``tool: cnvkit`` in whole genome mode (without ``depends_on.kits``), which also gives the target and antitarget regions.
  Without normals, the panel task builds a flat reference.
- ``control_freec`` runs Control-FREEC and transforms its ratios into log2 values, calls and segments.
  It needs ``depends_on.database``: an ``external_data`` task with ``{type: database, tags: [control_freec]}``, the key ``chr_len`` (the chromosome lengths) and, optionally, the key ``mappability``.

The task gives one result per tumor DNA library.
``control_freec`` leaves out tumor libraries without a matched normal, with a warning.
It takes the sex of each donor from the sample sheet (``XY`` for male, ``XX`` for female donors) and ``control_freec.unknown_sex`` for the others.
``cnvetti`` segments the coverage of such a tumor alone.
The files are unfiltered.
They contain low-quality calls.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
The files depend on the tool and are generally not in a standard format.

- ``cnvetti``: the segments ``{tumor_library}/out/{tumor_library}.bcf`` and the segmented windows ``{tumor_library}/out/{tumor_library}.windows.bcf``, each with ``.csi``.
- ``cnvkit``: ``{tumor_library}/out/{tumor_library}.cnr``, ``.dnacopy.seg``, ``.seg``, ``.bed.gz`` (with ``.tbi``) and ``.vcf.gz`` (with ``.tbi``).
  The directory ``{tumor_library}/report/`` holds the tables ``{tumor_library}.breaks.txt``, ``.genemetrics.txt``, ``.segmetrics.txt``, ``.sex.txt`` and ``.metrics.txt``, and the plots ``{tumor_library}.diagram.pdf`` and ``{tumor_library}.scatter.png`` (also per chromosome, ``{tumor_library}.scatter.chr1.png`` and so on).
- ``control_freec``: ``{tumor_library}/out/{tumor_library}.ratio.txt``, ``.gene_log2.txt``, ``.gene_call.txt``, ``.segments.txt``, ``.scatter.png`` and ``.diagram.pdf``.

Downstream tasks read the key ``dnacopy_seg``, which only ``cnvkit`` provides, so ``cnvetti`` and ``control_freec`` tasks do not satisfy their requirement.
``cnvetti`` also provides the keys ``bcf`` and ``bcf_csi`` (the segments), and ``control_freec`` the keys ``log2``, ``call`` and ``segments``.
"""

from functools import cached_property
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.cnvkit import SomaticCnvkitStepPart

from .model import SomaticWgsCnvCalling as SomaticWgsCnvCallingConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class CnvettiSomaticWgsStepPart(BaseStepPart):
    """Somatic WGS CNV calling with CNVetti"""

    #: Step name
    name = "cnvetti"

    #: Class available actions
    actions = ("coverage", "tumor_normal_ratio", "segment")
    wrappers = {
        "segment": "cnvetti/wgs/segment",
        "tumor_normal_ratio": "cnvetti/wgs/tumor_normal_ratio",
        "coverage": "cnvetti/wgs/coverage",
    }

    extra_arguments_action = "coverage"

    entity_wildcard = "tumor_library"

    default_resource_usage = ResourceUsage(threads=4, runtime="4h", mem=f"{int(3.75 * 1024 * 4)}MB")

    intermediate_outputs = {
        "coverage": ("bcf", "bcf_csi"),
        "tumor_normal_ratio": ("bcf", "bcf_csi"),
    }

    #: Extension file dictionary. Key: file type (string); Value: file extension (string)
    bcf_dict = {"bcf": ".bcf", "bcf_csi": ".bcf.csi"}

    #: The files that downstream tasks read, the segments: the key of ``get_output_paths`` and the
    #: suffix after ``work/{tumor_library}/out/{tumor_library}``
    output_suffixes = {"bcf": ".bcf", "bcf_csi": ".bcf.csi"}

    #: The files of ``cnvetti cmd segment``: the segments and the segmented windows
    segment_dict = {
        "segments_bcf": output_suffixes["bcf"],
        "segments_bcf_csi": output_suffixes["bcf_csi"],
        "windows_bcf": ".windows.bcf",
        "windows_bcf_csi": ".windows.bcf.csi",
    }

    result_keys = {"segment": tuple(k + md5 for k in segment_dict for md5 in ("", "_md5"))}

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
        tumor_library = wildcards.library_name
        normal_library = self.parent.matched_normal(tumor_library)
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
        prefix = "work/{tumor_library}/out/{tumor_library}"
        for key, ext in self.segment_dict.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

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


class CnvkitSomaticWgsStepPart(SomaticCnvkitStepPart):
    """Somatic WGS CNV calling with cnvkit.py"""

    resource_usage = {
        "plot": ResourceUsage(threads=1, runtime="1d", mem=f"{30 * 1024}MB"),
        "coverage": ResourceUsage(threads=8, runtime="1d", mem=f"{16 * 1024}MB"),
    }

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem=f"{int(7.5 * 1024)}MB")


class ControlFreecSomaticWgsStepPart(BaseStepPart):
    """Somatic WGS CNV calling with Control-FreeC"""

    report_outputs = {
        "plot": {
            key: ReportOutput("Control-FREEC", caption="control_freec_plots.rst")
            for key in ("diagram", "scatter")
        }
    }

    #: Step name
    name = "control_freec"

    #: Class available actions
    actions = ("run", "transform", "plot")
    wrappers = {
        "plot": "cnvkit/plot",
        "transform": "control_freec/transform",
        "run": "control_freec",
    }

    extra_arguments_action = "run"

    intermediate_outputs = {"transform": ("cns", "cnr")}

    entity_wildcard = "tumor_library"

    #: The files that downstream tasks read: the key of ``get_output_paths`` and the suffix after
    #: ``work/{tumor_library}/out/{tumor_library}``
    output_suffixes = {
        "log2": ".gene_log2.txt",
        "call": ".gene_call.txt",
        "segments": ".segments.txt",
    }

    result_keys = {
        "run": ("ratio", "ratio_md5"),
        "transform": ("log2", "call", "segments", "log2_md5", "call_md5", "segments_md5"),
        "plot": ("scatter", "diagram", "scatter_md5", "diagram_md5"),
    }

    resource_usage = {
        "plot": ResourceUsage(threads=1, runtime="40h", mem=f"{2 * 30 * 1024}MB"),
        "transform": ResourceUsage(threads=1, runtime="40h", mem=f"{2 * 8 * 1024}MB"),
        "run": ResourceUsage(threads=8, runtime="40h", mem=f"{int(2 * 3.75 * 1024 * 8)}MB"),
    }

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}

    @cached_property
    def freec_sex(self) -> dict[str, str]:
        """Control-FREEC's ``sex`` (XY, XX) of each library whose donor's sex is known"""
        df = self.parent.build_library_dataframe()
        sexes = {"male": "XY", "female": "XX"}
        return {
            library: sexes[str(sex).lower()]
            for library, sex in zip(df["library_name"], df["sex"])
            if str(sex).lower() in sexes
        }

    def _get_output_files_run(self):
        base = "work/{tumor_library}/out/{tumor_library}"
        return {"ratio": base + ".ratio.txt", "ratio_md5": base + ".ratio.txt.md5"}

    def _get_output_files_transform(self):
        base = "work/{tumor_library}/out/{tumor_library}"
        result = {
            **{key: base + suffix for key, suffix in self.output_suffixes.items()},
            "cns": base + ".cns",
            "cnr": base + ".cnr",
        }
        for name in ("log2", "call", "segments"):  # result files
            result[name + "_md5"] = result[name] + ".md5"
        return result

    def _get_output_files_plot(self):
        base = "work/{tumor_library}/out/{tumor_library}"
        result = {
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
        yield from self.parent.tumor_normal_alignments(wildcards.tumor_library).items()
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
            "ploidy": cfg.ploidy,
            "break_point_threshold": cfg.break_point_threshold,
            "coefficient_of_variation": cfg.coefficient_of_variation,
            "contamination_adjustment": cfg.contamination_adjustment,
            "contamination": cfg.contamination,
            "minimal_subclone_presence": cfg.minimal_subclone_presence,
            "unique_match": cfg.unique_match,
            "force_gc_content_normalization": cfg.force_gc_content_normalization,
            "min_cna_length": cfg.min_cna_length,
            "sex": self.freec_sex.get(wildcards.tumor_library, cfg.unknown_sex),
        }

    def _get_params_transform(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.control_freec
        return {
            "org_obj": cfg.convert.org_obj,
            "tx_obj": cfg.convert.tx_obj,
            "bs_obj": cfg.convert.bs_obj,
            "amplification": cfg.convert.amplification,
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

    tool_parts = {
        Tool.cnvetti: CnvettiSomaticWgsStepPart,
        Tool.cnvkit: CnvkitSomaticWgsStepPart,
        Tool.control_freec: ControlFreecSomaticWgsStepPart,
    }

    library_requirement = "kind == 'cancer'"

    @classmethod
    def task_library_requirement(cls, config):
        """control_freec compares each tumor with its matched normal; cnvkit uses the panel of
        normals, and cnvetti takes the coverage of a tumor without a normal alone."""
        if config.tool == Tool.control_freec:
            return MATCHED_NORMAL_REQUIREMENT
        return cls.library_requirement

    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local somatic WGS CNV output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        suffixes = cls.tool_parts[config.tool].output_suffixes
        return {key: f"work/{lib}/out/{lib}{suffix}" for key, suffix in suffixes.items()}
