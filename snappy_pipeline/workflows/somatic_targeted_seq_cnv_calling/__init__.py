# -*- coding: utf-8 -*-
"""Somatic CNV calling for targeted sequencing data

The step detects copy number alterations in cancer samples from targeted sequencing (exomes or large panels).
The ``tool`` field selects the caller, which starts from the aligned reads:

- ``cnvkit`` compares the tumor coverage with a CNVkit panel of normals (``depends_on.panel_of_normals``, from a ``panel_of_normals`` task with ``tool: cnvkit``) and reads its target and antitarget regions.
- ``sequenza`` compares the tumor with its matched normal and estimates the purity and ploidy.
- ``purecn`` needs a PureCN panel of normals (``depends_on.panel_of_normals``, from a ``panel_of_normals`` task with ``tool: purecn``) and somatic variants (``depends_on.variants``).
  The variants must come from Mutect2 and include the germline variants.
  The PureCN jobs run in the container image of the panel of normals task (its ``purecn.container``), so the panel and the calls use the same PureCN.
  They also take the panel's ``purecn.genome_name`` and ``purecn.seed``.
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

from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    RuleOptions,
    MATCHED_NORMAL_RELATIONSHIP,
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.common.cnvkit import SomaticCnvkitStepPart

from .model import SomaticTargetedSeqCnvCalling as SomaticTargetedSeqCnvCallingConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: The path of a tumor library's files without their suffix
PREFIX = "work/{tumor_library}/out/{tumor_library}"


class SequenzaStepPart(BaseStepPart):
    """Perform somatic targeted CNV calling using sequenza"""

    #: Step name
    name = "sequenza"

    #: Class available actions
    actions = (
        "gcreference",
        "coverage",
        "run",
    )
    wrappers = {
        "run": "sequenza/run",
        "coverage": "sequenza/coverage",
        "gcreference": "sequenza/gcreference",
    }

    entity_wildcard = "tumor_library"

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

    #: The plots of ``sequenza.results`` in the report directory
    report_outputs = {
        "run": {
            "report": ReportOutput("sequenza", caption="sequenza.rst", patterns=("{plot}.pdf",))
        }
    }

    #: The files that downstream tasks read: the key of ``get_output_paths`` and the suffix after
    #: ``PREFIX``
    output_suffixes = {"dnacopy_seg": ".dnacopy.seg", "seqz": ".seqz.gz"}

    @dictify
    def _get_input_files_coverage(self, wildcards):
        yield (
            "gc",
            "work/static_data/out/gc.{length}.wig.gz".format(
                length=self.config.sequenza.length,
            ),
        )
        yield from self.parent.tumor_normal_alignments(wildcards.tumor_library).items()
        yield from self.parent.reference_inputs().items()

    def _get_input_files_gcreference(self, wildcards):
        return self.parent.reference_inputs()

    def _get_input_files_run(self, wildcards):
        """The pileup, and the genome index for the names of the contigs"""
        return {
            "seqz": PREFIX + self.output_suffixes["seqz"],
            "reference_fai": self.parent.get_upstream_paths("reference").fai,
        }

    def _get_output_files_gcreference(self):
        return {"gc": f"work/static_data/out/gc.{self.config.sequenza.length}.wig.gz"}

    def _get_output_files_coverage(self):
        seqz = PREFIX + self.output_suffixes["seqz"]
        return {"seqz": seqz, "seqz_md5": seqz + ".md5"}

    def _get_output_files_run(self):
        seg = PREFIX + self.output_suffixes["dnacopy_seg"]
        return {
            "seg": seg,
            "seg_md5": seg + ".md5",
            # sequenza.results writes plots and tables that the step does not name
            "report": directory("work/{tumor_library}/report"),
        }

    def _get_params_coverage(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "length": self.config.sequenza.length,
            "ignore_contigs": self.config.sequenza.ignore_contigs,
            "bam2seqz_options": self.config.sequenza.extra_args.model_dump(by_alias=True),
        }

    def _get_params_gcreference(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "length": self.config.sequenza.length,
        }

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.sequenza
        return {
            "assembly": cfg.assembly,
            "ignore_contigs": cfg.ignore_contigs,
            "extra_args_extract": cfg.extra_args_extract.model_dump(by_alias=True),
            "extra_args_fit": cfg.extra_args_fit.model_dump(by_alias=True),
            "library_name": wildcards.tumor_library,
        }

    def _get_log_file_gcreference(self):
        return {"log": f"work/static_data/log/gc.{self.config.sequenza.length}.log"}

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}


#: The plots of ``PureCN.R``: the solutions, their purity and ploidy, the segmentation and the
#: allele fractions per chromosome
PURECN_PLOTS = {
    "plot": ".pdf",
    "plot_local_optima": ".local_optima.pdf",
    "plot_segmentation": ".segmentation.pdf",
    "plot_chromosomes": ".chromosomes.pdf",
}


class PureCNStepPart(BaseStepPart):
    """Perform somatic targeted CNV calling using PureCN"""

    #: Step name
    name = "purecn"

    #: Class available actions
    actions = ("coverage", "run")
    wrappers = {"run": "purecn/run", "coverage": "purecn/coverage"}
    rule_options = {"coverage": RuleOptions(explicit=True), "run": RuleOptions(explicit=True)}

    entity_wildcard = "tumor_library"

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

    #: The plots of ``PureCN.R``, for the report
    report_outputs = {
        "run": {key: ReportOutput("PureCN", caption="purecn.rst") for key in PURECN_PLOTS}
    }

    #: The files that downstream tasks read: the key of ``get_output_paths`` and the suffix after
    #: ``PREFIX``
    output_suffixes = {"dnacopy_seg": ".dnacopy.seg"}

    @property
    def result_keys(self):
        return {"run": tuple(k for k in self._get_output_files_run() if k not in PURECN_PLOTS)}

    @property
    def panel(self):
        """The PureCN section of the panel of normals task, whose image, genome and seed the
        calls use, too"""
        return self.parent.get_task_config("panel_of_normals").purecn

    @property
    def image(self) -> str:
        """The PureCN image of the panel of normals task, so that the panel and the calls use the
        same PureCN"""
        return self.panel.container

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "tumor", PREFIX.format(**wildcards) + ".coverage_loess.txt.gz"
        pon = self.parent.get_upstream_paths("panel_of_normals")
        somatic_vcf = self.parent.get_upstream_paths(
            "variants", library_name=wildcards.tumor_library
        )
        yield "vcf", somatic_vcf.vcf
        yield "normaldb", pon.panel_of_normals
        yield "mapping_bias", pon.mapping_bias
        yield "intervals", pon.intervals

    @dictify
    def _get_input_files_coverage(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.tumor_library
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        pon = self.parent.get_upstream_paths("panel_of_normals")
        yield "intervals", pon.intervals

    def _get_output_files_coverage(self):
        return {"coverage": PREFIX + ".coverage_loess.txt.gz"}

    def _get_output_files_run(self):
        tables = {
            "segments": self.output_suffixes["dnacopy_seg"],
            "ploidy": ".csv",
            "pvalues": ".amplification_pvalues.csv",
            "loh": ".loh.csv",
        }
        paths = {key: PREFIX + suffix for key, suffix in tables.items()}
        return {
            **paths,
            **{key + "_md5": path + ".md5" for key, path in paths.items()},
            **vcf_outputs(PREFIX),
            **{key: PREFIX + suffix for key, suffix in PURECN_PLOTS.items()},
        }

    def _get_params_coverage(self, wildcards):
        return {
            "config": {
                **self.config.purecn.model_dump(by_alias=True),
                "genome_name": self.panel.genome_name,
                "seed": self.panel.seed,
            }
        }

    def _get_params_run(self, wildcards):
        return {**self._get_params_coverage(wildcards), "library_name": wildcards.tumor_library}

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}


class CnvKitStepPart(SomaticCnvkitStepPart):
    """Perform somatic targeted CNV calling using cnvkit"""

    entity_wildcard = "tumor_library"

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


class SomaticTargetedSeqCnvCallingWorkflow(BaseStep):
    """Perform somatic targeted sequencing CNV calling"""

    #: Workflow name
    name = "somatic_targeted_seq_cnv_calling"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"}))]

    config_model_class = SomaticTargetedSeqCnvCallingConfigModel

    tool_parts = {
        Tool.cnvkit: CnvKitStepPart,
        Tool.sequenza: SequenzaStepPart,
        Tool.purecn: PureCNStepPart,
    }

    @classmethod
    def task_produces(cls, config, upstream):
        """The somatic CNVs; the panel of normals must have been built with the task's tool."""
        pon = upstream.get("panel_of_normals")
        if config.tool in (Tool.cnvkit, Tool.purecn) and str(config.tool) not in pon.tags:
            raise ConfigError(
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

    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local somatic targeted CNV output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        suffixes = cls.tool_parts[config.tool].output_suffixes
        return {key: f"work/{lib}/out/{lib}{suffix}" for key, suffix in suffixes.items()}
