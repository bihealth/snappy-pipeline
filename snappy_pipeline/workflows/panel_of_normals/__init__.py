# -*- coding: utf-8 -*-
"""Panel of normals creation

The step builds the background models that somatic variant and copy number callers use to remove systematic artefacts.
The ``tool`` field selects the model:

- ``mutect2`` creates the Mutect2 panel of normals (VCF) and the GenomicsDB it is built from.
  It needs ``depends_on.germline_resource``: an ``external_data`` task with ``{type: variants, tags: [population]}`` and the key ``vcf``, for example the af-only gnomAD VCF of the GATK resource bundle.
- ``cnvkit`` creates the CNVkit reference from the target, antitarget and coverage files of the normals.
- ``purecn`` creates the PureCN normal database, the mapping bias and the optimized intervals.
- ``access`` creates the CNVkit access file (accessible regions) from the reference genome and, optionally, regions to exclude.

The task yields one result per cohort.
All DNA libraries of the normal samples in the sample sheets build the panel, those of background data sets included, unless ``library_selection`` narrows them.
For ``mutect2`` and ``purecn``, ``path_normals_list`` selects other libraries by name, one per line.
The listed libraries are used even if the sample sheet does not flag them as normal.

.. warning::

    A panel of normals only works on data generated like the normals.
    A panel for exomes must come from the same exome enrichment kit as the tumor samples it is used on.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
The individual files per normal stay in ``work/``, so the panel can be extended when more normals are available.

- ``mutect2``: ``panel_of_normals/out/panel_of_normals.vcf.gz`` (with ``.tbi``) and ``panel_of_normals/out/panel_of_normals.genomicsDB.tar.gz``.
- ``cnvkit``: ``panel_of_normals/out/panel_of_normals.cnn``, ``panel_of_normals.target.bed`` and ``panel_of_normals.antitarget.bed``, and the reports ``panel_of_normals/report/panel_of_normals.metrics.tsv`` and ``panel_of_normals.sex.tsv``.
- ``purecn``: ``panel_of_normals/out/panel_of_normals.rds``, ``panel_of_normals.mapping_bias.rds``, ``panel_of_normals.intervals.list`` and ``panel_of_normals.intervals.bed.gz`` (with ``.tbi``).
- ``access``: ``panel_of_normals/out/panel_of_normals.accessible.bed``.

Downstream tasks read these keys through their ``depends_on.panel_of_normals``:

- ``mutect2``: ``vcf`` and ``vcf_tbi`` (the panel, read by ``variant_calling``), and ``genomicsdb`` (read by the ``purecn`` task).
- ``cnvkit``: ``panel_of_normals``, ``target`` and ``antitarget``.
- ``purecn``: ``panel_of_normals``, ``mapping_bias``, ``intervals`` and ``container``.

The ``access`` task gives the key ``access``.
A ``cnvkit`` task names it in ``depends_on.access``, as it can an ``external_data`` task that produces ``{type: raw, tags: [access]}`` with the key ``access``.

The CNVkit authors recommend checking the metrics and sex reports before a panel is used.
The sex report predicts the donor sex from the coverage of the chromosomes X and Y.

Notes on CNVkit
---------------

CNVkit was designed for exomes, and its modular design also covers whole genomes and amplicons.
Its documentation recommends a reference built from enough normals, for exomes and for whole genomes.
The two data types differ in the workflow that creates the reference (CNVkit ``batch``).
The step follows the common practice and distributes the computations over several jobs.

Without ``depends_on.kits``, ``cnvkit`` runs in whole genome mode.
There, the reference has no edge correction, which models the coverage drop at bait edges, whatever ``cnvkit.reference.edge_correction`` says (as ``cnvkit.py batch --method wgs``).
With it, the baits of the common kit of the normals are the targets.
For exomes, the access file only restricts the antitarget regions.
For whole genomes, ``autobin`` uses it to compute the average target size.
Without normals, there is nothing to run ``autobin`` on, and ``cnvkit.target.default_avg_size`` is the average target size.
Without ``depends_on.access``, whole genome mode creates an access file with default settings.
The CNVkit distribution provides an access file for ``GRCh37``, but it is incomplete for ``GRCh38``.
Run a task with ``tool: access`` to create one from regions to exclude (``access.exclude``).

There is no provision yet to create separate panels for males and females.
With enough samples collected alike, separate panels give the best results.

Notes on PureCN
---------------

A ``purecn`` panel is required by ``somatic_targeted_seq_cnv_calling`` with ``tool: purecn``.
It needs a ``mutect2`` panel of normals task in ``depends_on.panel_of_normals``, because the step uses Mutect2 as variant caller (not Mutect, the PureCN default).
PureCN runs in its container (``purecn.container``), because its R package requirements are too complex for a conda environment.
The default image has no version tag, so set one to fix the PureCN version.
"""

from snappy_pipeline.models.cnvkit import Gender as CnvKitGender
from snappy_pipeline.models.cnvkit import PanelOfNormals as CnvKitModel
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ReportOutput, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import PanelOfNormals as PanelOfNormalsConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class PanelOfNormalsStepPart(BaseStepPart):
    """Base class for panel of normals step parts

    Two steps: the preparation is done on each normal samples separately, and the panel creation
    merges all the individual results in the the panel.
    """

    #: Step name (default, must be overwritten)
    name = None

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.normal_libraries = list(self._get_normal_libraries())
        if self.name and (cfg := self.config.get(self.name)):
            if path := cfg.get("path_normals_list"):
                self.normal_libraries = []
                with open(path, "rt") as f:
                    for line in f:
                        self.normal_libraries.append(line.strip())

    def _common_kit(self) -> str:
        """The one kit of the normals; the panel is built with its baits."""
        return self.parent.common_kit("kits", self.normal_libraries)

    def _get_normal_libraries(self):
        """Yield the names of the task's normal DNA libraries, including background data sets."""
        yield from self.parent.libraries["library_name"]

    #: Log file prefix by action
    log_templates: dict[str, str] = {}

    def _get_log_file(self, action):
        return {"log": self.log_templates[action] + ".log"}


class PureCnStepPart(PanelOfNormalsStepPart):
    """Creating a panel of normals with GC-corrected coverage"""

    report_outputs = {
        "create_panel": {"plot": ReportOutput("PureCN", caption="purecn_weights.rst")}
    }

    #: Step name
    name = "purecn"

    #: Actions
    actions = ("install", "prepare", "coverage", "create_panel")

    extra_arguments_action = "create_panel"

    #: Per-library coverage is only read by ``create_panel``
    intermediate_outputs = {"coverage": ("coverage",)}

    def __init__(self, parent):
        super().__init__(parent)
        self.kit = self._common_kit()

    log_templates = {
        "install": "work/containers/log/container",
        "prepare": "work/panel_of_normals/log/panel_of_normals.prepare",
        "coverage": "work/{library_name}/log/{library_name}.coverage",
        "create_panel": "work/panel_of_normals/log/panel_of_normals",
    }

    #: Resources
    resource_usage = {
        "install": ResourceUsage(
            threads=1,
            runtime="1h",
            mem="24GB",
        ),
        "prepare": ResourceUsage(
            threads=1,
            runtime="4h",  # 4 hours
            mem="24GB",
        ),
        "coverage": ResourceUsage(
            threads=1,
            runtime="4h",  # 4 hours
            mem="24GB",
        ),
        "create_panel": ResourceUsage(
            threads=1,
            runtime="12h",  # 12 hours
            mem="32GB",
        ),
    }

    @dictify
    def _get_input_files_prepare(self, wildcards):
        yield "container", "work/containers/out/container.simg"
        yield from self.parent.reference_inputs().items()
        yield "bait_regions", self.parent.get_upstream_paths("kits", library_kit=self.kit).baits
        for key, path in (
            ("mappability", self.config.purecn.mappability),
            ("reptiming", self.config.purecn.reptiming),
        ):
            if path:
                yield key, path

    @dictify
    def _get_input_files_coverage(self, wildcards):
        yield "container", "work/containers/out/container.simg"
        yield (
            "intervals",
            "work/panel_of_normals/out/panel_of_normals.intervals.list",
        )
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    @dictify
    def _get_input_files_create_panel(self, wildcards):
        yield "container", "work/containers/out/container.simg"
        tpl = "work/{library_name}/out/{library_name}.coverage_loess.txt.gz"
        yield "normals", [tpl.format(library_name=lib) for lib in self.normal_libraries]
        # The Mutect2 genomicsDB is the output of the upstream panel_of_normals (mutect2) task.
        yield "genomicsdb", self.parent.get_upstream_paths("panel_of_normals").genomicsdb

    def _get_output_files_install(self):
        return {"container": "work/containers/out/container.simg"}

    def _get_output_files_prepare(self):
        base_out = "work/panel_of_normals/out/panel_of_normals.intervals"
        return {
            "intervals": base_out + ".list",
            "optimized": base_out + ".bed.gz",
            "tbi": base_out + ".bed.gz.tbi",
            "intervals_md5": base_out + ".list.md5",
            "optimized_md5": base_out + ".bed.gz.md5",
            "tbi_md5": base_out + ".bed.gz.tbi.md5",
        }

    def _get_output_files_coverage(self):
        return {"coverage": "work/{library_name}/out/{library_name}.coverage_loess.txt.gz"}

    def _get_output_files_create_panel(self):
        return {
            "db": "work/panel_of_normals/out/panel_of_normals.rds",
            "db_md5": "work/panel_of_normals/out/panel_of_normals.rds.md5",
            "mapbias": "work/panel_of_normals/out/panel_of_normals.mapping_bias.rds",
            "mapbias_md5": "work/panel_of_normals/out/panel_of_normals.mapping_bias.rds.md5",
            "lowcov": "work/panel_of_normals/out/panel_of_normals.low_coverage_targets.bed",
            "hq": "work/panel_of_normals/out/panel_of_normals.hq_sites.bed",
            "plot": "work/panel_of_normals/out/panel_of_normals.interval_weights.png",
        }

    #: Config paths reach the wrappers as inputs, so they are left out of the params
    _config_files = {"path_normals_list", "mappability", "reptiming"}

    def _get_params_prepare(self, wildcards):
        return {
            "config": self.config.get(self.name).model_dump(
                by_alias=True, exclude=self._config_files
            )
        }

    _get_params_create_panel = _get_params_prepare
    _get_params_coverage = _get_params_prepare


class Mutect2StepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    #: Step name
    name = "mutect2"

    #: Class available actions
    actions = ("scatter", "prepare_panel", "gather", "create_panel")

    #: ``create_panel`` takes the arguments of the sub-section ``genomicsdb``
    extra_arguments_action = "prepare_panel"

    #: The genomicsDB of ``create_panel`` is read by the PureCN panel task, so it stays
    temp_outputs = {"scatter": ("regions",), "prepare_panel": ("vcf", "vcf_tbi")}

    intermediate_outputs = {"gather": ("vcf", "vcf_tbi")}

    log_templates = {
        "scatter": "work/{normal_library}/log/{normal_library}.scatter",
        "prepare_panel": "work/{normal_library}/log/{normal_library}.{scatteritem}",
        "gather": "work/{normal_library}/log/{normal_library}",
        "create_panel": "work/panel_of_normals/log/panel_of_normals",
    }

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage = {
        "scatter": ResourceUsage(
            threads=1,
            runtime="2m",
            mem="1000MB",
        ),
        "prepare_panel": ResourceUsage(
            threads=2,
            runtime="3d",  # 3 days
            mem="8GB",
        ),
        "gather": ResourceUsage(
            threads=1,
            runtime="4h",
            mem="8192MB",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            runtime="48h",  # 48 hours
            mem="30GB",
        ),
    }

    def _get_input_files_scatter(self, wildcards):
        return {"fai": self.parent.get_upstream_paths("reference").fai}

    def _get_input_files_prepare_panel(self, wildcards):
        """Helper wrapper function for single sample panel preparation"""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.normal_library
        )
        scatteritem_base_path = "work/{normal_library}/par/scatter/{scatteritem}.region.bed"
        return {
            "normal_bam": alignments.bam,
            "normal_bai": alignments.bai,
            "region": scatteritem_base_path.format(**wildcards),
            **self.parent.reference_inputs(),
        }

    def _get_input_files_gather(self, wildcards):
        gather = self.parent.workflow.globals.get("gather")
        gather = getattr(gather, self.name)
        tpl = "work/{normal_library}/par/run/{{scatteritem}}.vcf.gz".format(**wildcards)
        return {"vcf": gather(tpl), "vcf_tbi": gather(tpl + ".tbi")}

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for merging individual results & panel creation"""
        paths = []
        tpl = "work/{normal_library}/out/{normal_library}.prepare.vcf.gz"
        for normal in self.normal_libraries:
            paths.append(tpl.format(normal_library=normal, **wildcards))
        return {
            "normals": paths,
            "normals_tbi": [path + ".tbi" for path in paths],
            **self.parent.reference_inputs(),
            **self.parent.vcf_inputs("germline_resource"),
        }

    def _get_output_files_scatter(self):
        scatter = self.parent.workflow.globals.get("scatter")
        scatter = getattr(scatter, self.name)
        tpl = "work/{{normal_library}}/par/scatter/{scatteritem}.region.bed"
        return {"regions": scatter(tpl)}

    @staticmethod
    def _vcf_outputs(prefix):
        return {"vcf": prefix + ".vcf.gz", "vcf_tbi": prefix + ".vcf.gz.tbi"}

    def _get_output_files_prepare_panel(self):
        return self._vcf_outputs("work/{normal_library}/par/run/{scatteritem}")

    def _get_output_files_gather(self):
        return self._vcf_outputs("work/{normal_library}/out/{normal_library}.prepare")

    def _get_output_files_create_panel(self):
        prefix = "work/panel_of_normals/out/panel_of_normals"
        return {
            **self._vcf_outputs(prefix),
            "vcf_md5": prefix + ".vcf.gz.md5",
            "vcf_tbi_md5": prefix + ".vcf.gz.tbi.md5",
            "db": prefix + ".genomicsDB.tar.gz",
            "db_md5": prefix + ".genomicsDB.tar.gz.md5",
        }

    def _get_params_scatter(self, wildcards):
        return {
            "ignore_contigs": self.config.ignore_contigs,
            "padding": self.config.mutect2.padding,
        }

    def _get_params_prepare_panel(self, wildcards):
        return {
            "max_mnp_distance": 0,
            "java_options": self.config.mutect2.java_options,
            "extra_arguments": self.config.mutect2.extra_arguments,
        }

    def _get_params_create_panel(self, wildcards):
        return self.config.mutect2.genomicsdb.model_dump(by_alias=True)


class CnvkitStepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    report_outputs = {
        "report": {
            key: ReportOutput("cnvkit", caption="cnvkit_reports.rst") for key in ("sex", "metrics")
        }
    }

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "target",
        "antitarget",
        "coverage",
        "create_panel",
        "report",
    )

    #: ``create_panel`` runs ``cnvkit.py reference``
    extra_arguments_action = "create_panel"

    #: Per-library coverage is only read by ``create_panel``
    intermediate_outputs = {"coverage": ("target", "antitarget")}

    log_templates = {
        "target": "work/panel_of_normals/log/panel_of_normals.target",
        "antitarget": "work/panel_of_normals/log/panel_of_normals.antitarget",
        "coverage": "work/{normal_library}/log/{normal_library}.coverage",
        "create_panel": "work/panel_of_normals/log/panel_of_normals",
        "report": "work/panel_of_normals/log/panel_of_normals.report",
    }

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage = {
        "target": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="8GB",
        ),
        "antitarget": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="8GB",
        ),
        "coverage": ResourceUsage(
            threads=8,
            runtime="2h",  # 2 hours
            mem="16GB",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="16GB",
        ),
        "report": ResourceUsage(
            threads=2,
            runtime="2h",  # 2 hours
            mem="16GB",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.is_wgs = not self.config.depends_on.kits
        self.kit = None if self.is_wgs else self._common_kit()

    def _get_params_target(self, wildcards):
        return self._cnvkit_params("target")

    def _get_params_antitarget(self, wildcards):
        return self._cnvkit_params("antitarget")

    def _get_params_coverage(self, wildcards):
        return self._cnvkit_params("coverage")

    def _get_params_create_panel(self, wildcards):
        return self._cnvkit_params("create_panel")

    def _get_params_report(self, wildcards):
        return self._cnvkit_params("report")

    def _cnvkit_params(self, action):
        cfg: CnvKitModel = self.config.get(self.name)
        if action == "create_panel":
            action = "reference"
        if args := getattr(cfg, action, {}):
            args = args.model_dump(by_alias=True)
        if action == "reference":
            if cfg.gender != CnvKitGender.guess:
                args["gender"] = cfg.gender
            if cfg.male_reference:
                args["male_reference"] = True
            args["flat"] = len(self.normal_libraries) == 0
            if self.is_wgs:
                args["edge_correction"] = False
        if action == "target":
            args["bp_per_bin"] = cfg.bp_per_bin
        return args

    def _get_input_files_target(self, wildcards):
        """Helper wrapper function to estimate target average size in wgs mode"""
        if not self.is_wgs:
            input_files = {
                "target": self.parent.get_upstream_paths("kits", library_kit=self.kit).baits
            }
            if self.config.cnvkit.path_annotation:
                input_files["annotate"] = self.config.cnvkit.path_annotation
            return input_files
        alignments = [
            self.parent.get_upstream_paths("alignments", library_name=x)
            for x in self.normal_libraries
        ]
        input_files = {
            "bams": [x.bam for x in alignments],
            "bais": [x.bai for x in alignments],
            **self.parent.reference_inputs(),
        }
        if self.config.depends_on.access:
            input_files["access"] = self.parent.get_upstream_paths("access").access
        if self.config.cnvkit.path_annotation:
            input_files["annotate"] = self.config.cnvkit.path_annotation
        return input_files

    def _get_input_files_antitarget(self, wildcards):
        """Helper wrapper function for computing antitarget locations"""
        if self.is_wgs:
            return {}
        result = {
            "target": "work/panel_of_normals/out/panel_of_normals.target.bed".format(**wildcards),
        }
        if self.config.depends_on.access:
            result["access"] = self.parent.get_upstream_paths("access").access
        return result

    def _get_input_files_coverage(self, wildcards):
        """Helper wrapper function for computing coverage"""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.normal_library
        )
        return {
            "target": "work/panel_of_normals/out/panel_of_normals.target.bed".format(**wildcards),
            "antitarget": "work/panel_of_normals/out/panel_of_normals.antitarget.bed".format(
                **wildcards
            ),
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
        }

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for computing panel of normals"""
        tpl = "work/{normal_library}/out/{normal_library}.targetcoverage.cnn"
        targets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        tpl = "work/{normal_library}/out/{normal_library}.antitargetcoverage.cnn"
        antitargets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        return {
            "target": (
                targets
                if targets
                else "work/panel_of_normals/out/panel_of_normals.target.bed".format(**wildcards)
            ),
            "antitarget": (
                antitargets
                if antitargets
                else "work/panel_of_normals/out/panel_of_normals.antitarget.bed".format(**wildcards)
            ),
            **self.parent.reference_inputs(),
        }

    def _get_input_files_report(self, wildcards):
        """Helper wrapper function for the panel of normals report"""
        tpl = "work/{normal_library}/out/{normal_library}.targetcoverage.cnn"
        targets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        tpl = "work/{normal_library}/out/{normal_library}.antitargetcoverage.cnn"
        antitargets = [tpl.format(normal_library=x) for x in self.normal_libraries]
        return {
            "target": targets,
            "antitarget": antitargets,
        }

    def _get_output_files_target(self):
        return {
            "target": "work/panel_of_normals/out/panel_of_normals.target.bed",
            "target_md5": "work/panel_of_normals/out/panel_of_normals.target.bed.md5",
        }

    def _get_output_files_antitarget(self):
        return {
            "antitarget": "work/panel_of_normals/out/panel_of_normals.antitarget.bed",
            "antitarget_md5": "work/panel_of_normals/out/panel_of_normals.antitarget.bed.md5",
        }

    def _get_output_files_coverage(self):
        return {
            "target": "work/{normal_library}/out/{normal_library}.targetcoverage.cnn",
            "antitarget": "work/{normal_library}/out/{normal_library}.antitargetcoverage.cnn",
        }

    def _get_output_files_create_panel(self):
        return {
            "panel": "work/panel_of_normals/out/panel_of_normals.cnn",
            "panel_md5": "work/panel_of_normals/out/panel_of_normals.cnn.md5",
        }

    def _get_output_files_report(self):
        return {
            "sex": "work/panel_of_normals/report/panel_of_normals.sex.tsv",
            "sex_md5": "work/panel_of_normals/report/panel_of_normals.sex.tsv.md5",
            "metrics": "work/panel_of_normals/report/panel_of_normals.metrics.tsv",
            "metrics_md5": "work/panel_of_normals/report/panel_of_normals.metrics.tsv.md5",
        }


class AccessStepPart(PanelOfNormalsStepPart):
    """Utility to create access file for cnvkit"""

    name = "access"
    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=2, runtime="2h", mem="8GB")

    log_templates = {"run": "work/panel_of_normals/log/panel_of_normals.accessible"}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from self.parent.reference_inputs().items()
        if self.config.access.exclude:
            yield "exclude", self.config.access.exclude

    def _get_params_run(self, wildcards):
        return {"min_gap_size": self.config.access.min_gap_size}

    def _get_output_files_run(self):
        tpl = "work/panel_of_normals/out/panel_of_normals.accessible.bed"
        return {"access": tpl, "access_md5": tpl + ".md5"}


class PanelOfNormalsWorkflow(BaseStep):
    """Creates a panel of normals"""

    # Workflow name
    name = "panel_of_normals"

    #: The panel is built from normal DNA libraries, of the background data sets too; all of
    #: them unless ``library_selection`` narrows them
    library_requirement = "role == 'normal' and extraction_type == 'dna'"
    library_default = ""
    background_results = True
    produces = [DataSignature(DataType.MODELS, frozenset({"pon"}))]

    @classmethod
    def task_produces(cls, config, upstream):
        """The panel of normals of the task's tool; the access regions for ``tool: access``."""
        if config.tool == Tool.access:
            return (DataSignature(DataType.RAW, frozenset({"access"})),)
        return (DataSignature(DataType.MODELS, frozenset({"pon", str(config.tool)})),)

    config_model_class = PanelOfNormalsConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local panel-of-normals output paths for downstream consumers."""
        match config.tool:
            case Tool.mutect2:
                prefix = "work/panel_of_normals/out/panel_of_normals"
                return {
                    "vcf": prefix + ".vcf.gz",
                    "vcf_tbi": prefix + ".vcf.gz.tbi",
                    "genomicsdb": prefix + ".genomicsDB.tar.gz",
                }
            case Tool.cnvkit:
                prefix = "work/panel_of_normals/out/panel_of_normals"
                return {
                    "panel_of_normals": prefix + ".cnn",
                    "target": prefix + ".target.bed",
                    "antitarget": prefix + ".antitarget.bed",
                }
            case Tool.purecn:
                return {
                    "panel_of_normals": "work/panel_of_normals/out/panel_of_normals.rds",
                    "mapping_bias": "work/panel_of_normals/out/panel_of_normals.mapping_bias.rds",
                    "intervals": "work/panel_of_normals/out/panel_of_normals.intervals.list",
                    "container": "work/containers/out/container.simg",
                }
            case Tool.access:
                return {"access": "work/panel_of_normals/out/panel_of_normals.accessible.bed"}
        return {}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        match self.config.tool:
            case "mutect2":
                selected = Mutect2StepPart
            case "cnvkit":
                selected = CnvkitStepPart
            case "access":
                selected = AccessStepPart
            case "purecn":
                selected = PureCnStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes((selected,))

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        result_files = []

        if self.config.tool == "mutect2":
            tpl = "output/panel_of_normals/out/panel_of_normals.{ext}"
            ext_list = (
                "vcf.gz",
                "vcf.gz.md5",
                "vcf.gz.tbi",
                "vcf.gz.tbi.md5",
                "genomicsDB.tar.gz",
                "genomicsDB.tar.gz.md5",
            )
            result_files.extend(self._expand_result_files(tpl, ext_list))

        if self.config.tool == "cnvkit":
            tpls = [
                (
                    "output/panel_of_normals/out/panel_of_normals.target.{ext}",
                    ("bed", "bed.md5"),
                ),
                (
                    "output/panel_of_normals/out/panel_of_normals.antitarget.{ext}",
                    ("bed", "bed.md5"),
                ),
                ("output/panel_of_normals/out/panel_of_normals.{ext}", ("cnn", "cnn.md5")),
                (
                    "output/panel_of_normals/report/panel_of_normals.sex.{ext}",
                    ("tsv", "tsv.md5"),
                ),
                (
                    "output/panel_of_normals/report/panel_of_normals.metrics.{ext}",
                    ("tsv", "tsv.md5"),
                ),
            ]
            for tpl, ext_list in tpls:
                result_files.extend(self._expand_result_files(tpl, ext_list))

        if self.config.tool == "access":
            tpl = "output/panel_of_normals/out/panel_of_normals.accessible.bed"
            result_files.extend([tpl + md5 for md5 in ("", ".md5")])

        if self.config.tool == "purecn":
            tpl = "output/panel_of_normals/out/panel_of_normals.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/panel_of_normals/out/panel_of_normals.mapping_bias.{ext}"
            ext_list = ("rds", "rds.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))
            tpl = "output/panel_of_normals/out/panel_of_normals.intervals.{ext}"
            ext_list = ("list", "list.md5", "bed.gz", "bed.gz.md5", "bed.gz.tbi", "bed.gz.tbi.md5")
            result_files.extend(self._expand_result_files(tpl, ext_list))

        return result_files

    def _expand_result_files(self, tpl, ext_list):
        for ext in ext_list:
            yield tpl.format(ext=ext)
