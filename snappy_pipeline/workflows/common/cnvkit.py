"""Somatic CNV calling with CNVkit against a CNVkit panel of normals, shared by
``somatic_targeted_seq_cnv_calling`` and ``somatic_wgs_cnv_calling``"""

from typing import Any

from snappy_pipeline.models.cnvkit import Gender
from snappy_pipeline.workflows.abstract import BaseStepPart, ReportOutput

#: The chromosomes of the per-chromosome scatter plots
CHROMOSOMES = (*range(1, 23), "X", "Y")


def sex_params(cfg) -> dict[str, Any]:
    """CNVkit's ``gender`` and ``male_reference`` arguments, where the config sets them"""
    args = {}
    if cfg.gender != Gender.guess:
        args["gender"] = cfg.gender
    if cfg.male_reference:
        args["male_reference"] = True
    return args


class SomaticCnvkitStepPart(BaseStepPart):
    """Somatic CNV calling with CNVkit: the tumor coverage of the panel's target and antitarget
    regions, its ratios against the panel, segments, calls, exports, plots and reports.

    A step subclasses it with its ``entity_wildcard`` and resource usage.
    """

    name = "cnvkit"

    actions = ("coverage", "fix", "segment", "call", "postprocess", "export", "plot", "report")
    wrappers = {action: f"cnvkit/{action}" for action in actions}

    extra_arguments_action = "segment"

    intermediate_outputs = {
        "coverage": ("target", "antitarget"),
        "segment": ("segments",),
        "call": ("calls",),
    }

    report_outputs = {
        "plot": {
            key: ReportOutput("cnvkit", caption="cnvkit_plots.rst")
            for key in ("diagram", "scatter", *(f"scatter_chr{chrom}" for chrom in CHROMOSOMES))
        },
        "report": {
            key: ReportOutput("cnvkit", caption="cnvkit_reports.rst")
            for key in ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        },
    }

    #: The files that downstream tasks read: the key of ``get_output_paths`` and the suffix after
    #: ``work/{entity}/out/{entity}``
    output_suffixes = {"dnacopy_seg": ".dnacopy.seg"}

    @property
    def result_keys(self):
        actions = ("fix", "postprocess", "export", "plot", "report")
        return {action: tuple(self.get_output_files(action)) for action in actions}

    def _path(self, directory: str, suffix: str) -> str:
        """The path of a file of the entity below ``work/{entity}/<directory>/``"""
        entity = f"{{{self.entity_wildcard}}}"
        return f"work/{entity}/{directory}/{entity}{suffix}"

    def _out(self, wildcards, suffix: str) -> str:
        """The path of an earlier action's file of the entity in ``wildcards``"""
        return self._path("out", suffix).format(**wildcards)

    def _get_input_files_coverage(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=getattr(wildcards, self.entity_wildcard)
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
        return {
            "target": self._out(wildcards, ".targetcoverage.cnn"),
            "antitarget": self._out(wildcards, ".antitargetcoverage.cnn"),
            "ref": self.parent.get_upstream_paths("panel_of_normals").panel_of_normals,
        }

    def _get_input_files_segment(self, wildcards):
        return {"cnr": self._out(wildcards, ".cnr")}

    def _get_input_files_call(self, wildcards):
        return {"segment": self._out(wildcards, ".segment.cns")}

    def _get_input_files_postprocess(self, wildcards):
        return {
            "segment": self._out(wildcards, ".segment.cns"),
            "call": self._out(wildcards, ".call.cns"),
        }

    def _get_input_files_export(self, wildcards):
        return {"cns": self._out(wildcards, ".call.cns")}

    def _get_input_files_plot(self, wildcards):
        return {"cnr": self._out(wildcards, ".cnr"), "cns": self._out(wildcards, ".call.cns")}

    def _get_input_files_report(self, wildcards):
        return {
            "target": self._out(wildcards, ".targetcoverage.cnn"),
            "antitarget": self._out(wildcards, ".antitargetcoverage.cnn"),
            "cnr": self._out(wildcards, ".cnr"),
            "cns": self._out(wildcards, ".call.cns"),
        }

    def _with_md5(self, directory: str, suffixes: dict[str, str]) -> dict[str, str]:
        """The files of ``suffixes`` below ``work/{entity}/<directory>/``, each followed by its md5 file"""
        paths = {}
        for key, suffix in suffixes.items():
            paths[key] = self._path(directory, suffix)
            paths[key + "_md5"] = paths[key] + ".md5"
        return paths

    def _get_output_files_coverage(self):
        return {
            "target": self._path("out", ".targetcoverage.cnn"),
            "antitarget": self._path("out", ".antitargetcoverage.cnn"),
        }

    def _get_output_files_fix(self):
        return self._with_md5("out", {"ratios": ".cnr"})

    def _get_output_files_segment(self):
        return {"segments": self._path("out", ".segment.cns")}

    def _get_output_files_call(self):
        return {"calls": self._path("out", ".call.cns")}

    def _get_output_files_postprocess(self):
        return self._with_md5("out", {"final": self.output_suffixes["dnacopy_seg"]})

    def _get_output_files_export(self):
        exports = {
            "bed": ".bed.gz",
            "bed_tbi": ".bed.gz.tbi",
            "seg": ".seg",
            "vcf": ".vcf.gz",
            "vcf_tbi": ".vcf.gz.tbi",
        }
        return self._with_md5("out", exports)

    def _get_output_files_plot(self):
        plots = {"diagram": ".diagram.pdf", "scatter": ".scatter.png"}
        plots |= {f"scatter_chr{chrom}": f".scatter.chr{chrom}.png" for chrom in CHROMOSOMES}
        return self._with_md5("report", plots)

    def _get_output_files_report(self):
        reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        return self._with_md5("report", {report: f".{report}.txt" for report in reports})

    def _get_log_file(self, action):
        entity = f"{{{self.entity_wildcard}}}"
        return {"log": f"work/{entity}/log/{action}.{entity}.log"}

    def _get_params_coverage(self, wildcards):
        return self.config.cnvkit.coverage.model_dump(by_alias=True)

    def _get_params_fix(self, wildcards):
        return self.config.cnvkit.fix.model_dump(by_alias=True)

    def _get_params_segment(self, wildcards):
        cfg = self.config.cnvkit
        return {**cfg.segment.model_dump(by_alias=True), "drop_low_coverage": cfg.drop_low_coverage}

    def _get_params_call(self, wildcards):
        cfg = self.config.cnvkit
        return {
            **cfg.call.model_dump(by_alias=True),
            "drop_low_coverage": cfg.drop_low_coverage,
            **sex_params(cfg),
        }

    def _get_params_postprocess(self, wildcards):
        return {
            "library_name": getattr(wildcards, self.entity_wildcard),
            "ploidy": self.config.cnvkit.call.ploidy,
        }

    def _get_params_plot(self, wildcards):
        cfg = self.config.cnvkit
        return {**cfg.diagram.model_dump(by_alias=True), **sex_params(cfg)}

    def _get_params_report(self, wildcards):
        cfg = self.config.cnvkit
        return {
            "breaks": cfg.breaks.model_dump(by_alias=True),
            "genemetrics": cfg.genemetrics.model_dump(by_alias=True),
            "segmetrics": cfg.segmetrics.model_dump(by_alias=True),
            "drop_low_coverage": cfg.drop_low_coverage,
        }
