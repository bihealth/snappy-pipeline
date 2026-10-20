# -*- coding: utf-8 -*-
"""Somatic neoepitope prediction with pVACtools

The step predicts neoepitopes of the tumor DNA libraries with pVACtools, run from its container: pVACseq from somatic small variants, pVACfuse from gene fusions, pVACsplice from splicing events.
The single ``tool`` is ``pvactools``; ``pvactools.pvacseq``, ``pvacfuse`` and ``pvacsplice`` are enabled separately (pVACseq by default) and share the preparation of the task.
The IEDB tools in the container are *not* free software; their use is restricted to academic, non-profit research unless a licence is acquired.

Inputs
------

- ``depends_on.hla_types`` (required) lists one or more HLA typing tasks.
  The task merges the HLA types of the tumor DNA library, its matched normal library and its matched RNA library from all of them, as far as each task typed these libraries.
- ``depends_on.variants`` (pVACseq, pVACsplice): somatic variants annotated by VEP (``variant_annotation`` with ``tool: vep``), with the ``Frameshift`` and ``Wildtype`` plugins of pVACtools in ``vep.plugins``.
  For all transcripts instead of one consequence per variant, it is an annotation task with ``vep.pick: false``.
- ``depends_on.fusions`` (pVACfuse): arriba fusions of the matched RNA library.
- ``depends_on.alignments_rna`` and ``depends_on.features`` (pVACsplice): STAR alignments of the matched RNA library, whose strandedness decision must be stranded, and the gene annotation.
- ``depends_on.reference``: the genome.
- ``depends_on.proteome`` (optional): a reference proteome FASTA from an ``external_data`` task (``{type: raw, tags: [proteome]}``, key ``fasta``) for the proteome similarity check of the epitopes.

pVACseq options and their inputs:

- ``pileup``: read counts of the matched RNA library at the variants (``depends_on.alignments_rna``).
- ``expression``: salmon TPMs of the matched RNA library (``depends_on.expression``).
- ``phasing``: germline variants near the somatic ones, phased with the tumor DNA (``depends_on.combined_variants``, the VEP-annotated output of ``combine_variants``, and ``depends_on.alignments_dna``, the alignments of the tumor DNA).

The matched normal library is the ``matched_normal_lib`` relationship (normal DNA of the donor); the matched RNA library is ``matched_rna_lib`` (tumor RNA of the same sample); both can be redefined in ``relationships``.
pVACfuse and pVACsplice skip tumor libraries without a matched RNA library.

``net_chop`` (with the path of a local NetChop installation) and ``netmhc_stab`` (netMHCstabpan in the container, MHC class I) post-process the epitopes of each prediction.
Without ``pvactools.container``, the task builds the container from ``docker://griffithlab/pvactools:7.1.3``.

Outputs
-------

The task publishes, below ``tasks/<task>/output/{tumor_library}/out/`` and each with a ``.md5`` file, for each enabled prediction (``pvacseq``, ``pvacfuse``, ``pvacsplice``):

- ``<prediction>/combined/{tumor_library}.Combined.filtered.tsv``, ``.all_epitopes.tsv`` and ``.all_epitopes.aggregated.tsv``: the epitopes of both MHC classes.
- with ``net_chop``: ``<prediction>/{tumor_library}.MHC_I.netchop.tsv`` and ``.MHC_II.netchop.tsv``, for the MHC classes that the HLA types and algorithms cover.
- with ``netmhc_stab``: ``<prediction>/{tumor_library}.MHC_I.netstab.tsv``.
"""

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ClassIAlgorithm, ClassIIAlgorithm
from .model import SomaticNeoepitopePrediction as SomaticNeoepitopePredictionConfigModel

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

#: The output directory and file name infix of pVACtools for each MHC class
MHC_CLASSES = {"class_i": ("MHC_Class_I", "MHC_I"), "class_ii": ("MHC_Class_II", "MHC_II")}

#: pVACtools' directory and infix for the epitopes of both classes
COMBINED = ("combined", "Combined")

#: pVACtools tables by output key; only pVACseq writes the metrics (``json``)
TABLES = {
    "all": "all_epitopes.tsv",
    "filtered": "filtered.tsv",
    "aggregated": "all_epitopes.aggregated.tsv",
    "json": "all_epitopes.aggregated.metrics.json",
}

#: The tables of both classes that the task publishes
PUBLISHED_TABLES = ("filtered", "all", "aggregated")

#: The predictions, which need a matched RNA library except pVACseq
PREDICTIONS = ("pvacseq", "pvacfuse", "pvacsplice")
NEED_RNA = ("pvacfuse", "pvacsplice")

#: The container that the task builds without ``pvactools.container``
CONTAINER_IMAGE = "docker://griffithlab/pvactools:7.1.3"
CONTAINER = "work/container/out/image.sif"

_OUT = "work/{tumor_library}/out/{tumor_library}"
_LOG = "work/{tumor_library}/log/{tumor_library}"

#: Wildcard constraints of the post-processing outputs
_PREDICTION = "{prediction,pvacseq|pvacfuse|pvacsplice}"
_MHC_CLASS = "{mhc_class,MHC_I|MHC_II}"


def _post_processed(tool: str, prediction: str = _PREDICTION, mhc_class: str = _MHC_CLASS) -> str:
    """The output of NetChop or netMHCstabpan (``tool``) for one prediction and MHC class"""
    return f"work/{{tumor_library}}/out/{prediction}/{{tumor_library}}.{mhc_class}.{tool}.tsv"


class NeoepitopeStepPart(BaseStepPart):
    """Shared helpers of the step parts"""

    entity_wildcard = "tumor_library"

    def samples(self, tumor_library: str) -> dict[str, str]:
        """The sample names in the VCF: the tumor and, if there is one, the matched normal."""
        args = {"tumor_sample": tumor_library, "tumor_library": tumor_library}
        if normal := self.parent.matched(tumor_library, "matched_normal_lib"):
            args |= {"normal_sample": normal, "normal_library": normal}
        return args

    @property
    def container(self) -> str:
        return self.config.pvactools.container or CONTAINER

    def variants(self, tumor_library: str) -> str:
        return self.parent.get_upstream_paths("variants", library_name=tumor_library).vcf


class PvacToolsStepPart(NeoepitopeStepPart):
    """Preparation shared by the predictions: container, VCF samples, HLA types"""

    name = "pvactools"
    actions = ("install", "normalize", "hla_types")

    resource_usage = {
        "install": ResourceUsage(threads=1, runtime="4h", mem="64GB"),
        "normalize": ResourceUsage(threads=1, runtime="1h", mem="4GB"),
        "hla_types": ResourceUsage(threads=1, runtime="1h", mem="4GB"),
    }
    intermediate_outputs = {"normalize": ("vcf", "vcf_tbi"), "hla_types": ("hla_types",)}

    def _get_output_files_install(self):
        return {"container": CONTAINER}

    def _get_log_file_install(self):
        return {"log": "work/container/log/image.log"}

    def _get_input_files_normalize(self, wildcards):
        return {"annotated": self.variants(wildcards.tumor_library)}

    def _get_output_files_normalize(self):
        return {"vcf": f"{_OUT}.normalized.vcf.gz", "vcf_tbi": f"{_OUT}.normalized.vcf.gz.tbi"}

    def _get_params_normalize(self, wildcards):
        return self.samples(wildcards.tumor_library)

    def _get_log_file_normalize(self):
        return {"log": f"{_LOG}.normalize.log"}

    @dictify
    def _get_input_files_hla_types(self, wildcards):
        """The HLA calls of the tumor, its matched normal and RNA library, from each HLA typing
        task that typed them."""
        tumor = wildcards.tumor_library
        libraries = [
            tumor,
            self.parent.matched(tumor, "matched_normal_lib"),
            self.parent.matched(tumor, "matched_rna_lib"),
        ]
        calls = []
        for task in self.parent.upstream_task_names("hla_types"):
            typed = set(self.parent.task_libraries(task)["library_name"])
            for library in filter(None, libraries):
                if library in typed:
                    paths = self.parent.get_upstream_paths(
                        "hla_types", upstream_task=task, library_name=library
                    )
                    calls.append(paths.calls_json)
        yield "calls", sorted(calls)

    def _get_output_files_hla_types(self):
        return {"hla_types": f"{_OUT}.hla_types.txt"}

    def _get_log_file_hla_types(self):
        return {"log": f"{_LOG}.hla_types.log"}


class PredictionStepPart(NeoepitopeStepPart):
    """A pVACtools prediction, run in the container"""

    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg = getattr(self.config.pvactools, self.name)

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=min(16, self.cfg.n_threads), runtime="24h", mem="64GB")

    published_keys = {
        "run": tuple(
            f"{key}.{COMBINED[1]}{md5}" for key in PUBLISHED_TABLES for md5 in ("", "_md5")
        )
    }

    @property
    def directory(self) -> str:
        return f"work/{{tumor_library}}/out/{self.name}"

    @property
    def intermediate_outputs(self) -> dict[str, tuple[str, ...]]:
        """The tables and sequences that only NetChop and netMHCstabpan read"""
        keys = []
        infixes = [MHC_CLASSES[name][1] for name in self.mhc_classes()]
        if self.cfg.net_chop.enabled:
            table = "all" if self.cfg.net_chop.all_epitopes else "filtered"
            keys += [f"{table}.{infix}" for infix in infixes]
            fasta = self.fasta_outputs()
            keys += [key for key in fasta if key == "fasta" or key.split(".")[1] in infixes]
        elif self.cfg.netmhc_stab.enabled and "MHC_I" in infixes:
            keys.append(f"{'all' if self.cfg.netmhc_stab.all_epitopes else 'filtered'}.MHC_I")
        return {"run": tuple(keys)}

    def fasta_outputs(self) -> dict[str, str]:
        """The epitope sequences of each MHC class, which NetChop reads"""
        return {
            f"fasta.{infix}": f"{self.directory}/{subdir}/{{tumor_library}}.fasta"
            for subdir, infix in MHC_CLASSES.values()
        }

    @dictify
    def _get_output_files_run(self):
        yield "done", f"{self.directory}/.done"
        for key, ext in TABLES.items():
            if key == "json" and self.name != "pvacseq":
                continue
            for subdir, infix in (*MHC_CLASSES.values(), COMBINED):
                path = f"{self.directory}/{subdir}/{{tumor_library}}.{infix}.{ext}"
                yield f"{key}.{infix}", path
                if infix == COMBINED[1] and key in PUBLISHED_TABLES:
                    yield f"{key}.{infix}_md5", path + ".md5"
        if self.cfg.net_chop.enabled:
            yield from self.fasta_outputs().items()

    def _get_log_file_run(self):
        return {"log": f"{_LOG}.{self.name}.log"}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "container", self.container
        yield "alleles", f"{_OUT}.hla_types.txt".format(**wildcards)
        if self.cfg.genes_of_interest:
            yield "genes", self.cfg.genes_of_interest
        if self.config.depends_on.proteome:
            yield "peptides", self.parent.get_upstream_paths("proteome").fasta

    def _get_params_run(self, wildcards):
        # The tables are written in the container's output directory, not bound one by one
        exclude = ["container", "alleles", *self._get_output_files_run()]
        exclude.remove("done")
        return {
            **self.samples(wildcards.tumor_library),
            "algorithms": self.cfg.algorithms,
            "lengths": {
                "class_i": self.cfg.class_i_epitope_length,
                "class_ii": self.cfg.class_ii_epitope_length,
            },
            "n_threads": self.cfg.n_threads,
            "exclude_bind": exclude,
            "extra_args": self.cfg.extra_args,
        }

    def mhc_classes(self) -> list[str]:
        """The MHC classes that the HLA types and the algorithms and lengths of the prediction
        cover"""
        result = []
        for name, algorithms, lengths in (
            ("class_i", ClassIAlgorithm, self.cfg.class_i_epitope_length),
            ("class_ii", ClassIIAlgorithm, self.cfg.class_ii_epitope_length),
        ):
            if name in self.parent.typed_mhc_classes and lengths:
                if any(algorithm in algorithms for algorithm in self.cfg.algorithms):
                    result.append(name)
        return result


class PvacSeqStepPart(PredictionStepPart):
    """pVACseq on somatic small variants, with RNA evidence and phased germline variants"""

    name = "pvacseq"
    actions = ("pileup", "add_expression", "phasing", "run")

    resource_usage = {
        "pileup": ResourceUsage(threads=1, runtime="4h", mem="6GB"),
        "add_expression": ResourceUsage(threads=1, runtime="4h", mem="6GB"),
        "phasing": ResourceUsage(threads=1, runtime="24h", mem="32GB"),
    }

    @property
    def intermediate_outputs(self) -> dict[str, tuple[str, ...]]:
        return {
            **super().intermediate_outputs,
            "pileup": ("vcf", "vcf_tbi"),
            "add_expression": ("vcf", "vcf_tbi"),
            "phasing": ("vcf",),
        }

    def with_expression(self, tumor_library: str) -> bool:
        """Whether pVACseq reads the VCF with RNA evidence"""
        has_rna = bool(self.parent.matched(tumor_library, "matched_rna_lib"))
        return has_rna and (self.cfg.pileup.enabled or self.cfg.expression.enabled)

    def _get_input_files_pileup(self, wildcards):
        rna = self.parent.matched(wildcards.tumor_library, "matched_rna_lib")
        alignments = self.parent.get_upstream_paths("alignments_rna", library_name=rna)
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            "loci": self.variants(wildcards.tumor_library),
            **self.parent.reference_inputs(),
        }

    def _get_output_files_pileup(self):
        return {"vcf": f"{_OUT}.pileup.vcf.gz", "vcf_tbi": f"{_OUT}.pileup.vcf.gz.tbi"}

    def _get_params_pileup(self, wildcards):
        return {"tumor_sample": wildcards.tumor_library, "extra_args": self.cfg.pileup.extra_args}

    def _get_log_file_pileup(self):
        return {"log": f"{_LOG}.pileup.log"}

    @dictify
    def _get_input_files_add_expression(self, wildcards):
        yield "annotated", f"{_OUT}.normalized.vcf.gz".format(**wildcards)
        yield "annotated_tbi", f"{_OUT}.normalized.vcf.gz.tbi".format(**wildcards)
        if self.cfg.pileup.enabled:
            yield "pileup", f"{_OUT}.pileup.vcf.gz".format(**wildcards)
            yield "pileup_tbi", f"{_OUT}.pileup.vcf.gz.tbi".format(**wildcards)
        if self.cfg.expression.enabled:
            rna = self.parent.matched(wildcards.tumor_library, "matched_rna_lib")
            expression = self.parent.get_upstream_paths("expression", library_name=rna)
            yield "gene_tpms", expression.gene_sf
            yield "transcript_tpms", expression.transcript_sf
            if table := self.cfg.expression.duplicate_transcripts_table:
                yield "duplicate_transcripts_table", table

    def _get_output_files_add_expression(self):
        return {
            "vcf": f"{_OUT}.with_expression.vcf.gz",
            "vcf_tbi": f"{_OUT}.with_expression.vcf.gz.tbi",
        }

    def _get_params_add_expression(self, wildcards):
        return {
            "format": "salmon",
            "extra_args": self.cfg.expression.extra_args,
            **self.samples(wildcards.tumor_library),
        }

    def _get_log_file_add_expression(self):
        return {"log": f"{_LOG}.add_expression.log"}

    def _get_input_files_phasing(self, wildcards):
        tumor = wildcards.tumor_library
        alignments = self.parent.get_upstream_paths("alignments_dna", library_name=tumor)
        return {
            "vcf": self.parent.get_upstream_paths("combined_variants", library_name=tumor).vcf,
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
        }

    def _get_output_files_phasing(self):
        return {"vcf": f"{_OUT}.phased.vcf.gz"}

    def _get_log_file_phasing(self):
        return {"log": f"{_LOG}.phasing.log"}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        prepared = (
            "with_expression" if self.with_expression(wildcards.tumor_library) else "normalized"
        )
        yield "vcf", f"{_OUT}.{prepared}.vcf.gz".format(**wildcards)
        yield "vcf_tbi", f"{_OUT}.{prepared}.vcf.gz.tbi".format(**wildcards)
        if self.cfg.phasing.enabled:
            yield "phased", f"{_OUT}.phased.vcf.gz".format(**wildcards)

    def _get_params_run(self, wildcards):
        params = super()._get_params_run(wildcards)
        if self.cfg.ml_predictions.enabled:
            params["extra_args"] = [
                *params["extra_args"],
                "--run-ml-predictions",
                f"--ml-threshold-accept {self.cfg.ml_predictions.accept}",
                f"--ml-threshold-reject {self.cfg.ml_predictions.reject}",
            ]
        return params


class PvacFuseStepPart(PredictionStepPart):
    """pVACfuse on the arriba fusions of the matched RNA library"""

    name = "pvacfuse"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        rna = self.parent.matched(wildcards.tumor_library, "matched_rna_lib")
        yield "fusions", self.parent.get_upstream_paths("fusions", library_name=rna).fusions

    def _get_params_run(self, wildcards):
        params = super()._get_params_run(wildcards)
        params.pop("normal_sample", None)
        params.pop("normal_library", None)
        return params


class PvacSpliceStepPart(PredictionStepPart):
    """pVACsplice on the splice junctions of the matched RNA library near somatic variants"""

    name = "pvacsplice"
    actions = ("features", "junction", "run")

    resource_usage = {
        "features": ResourceUsage(threads=1, runtime="1h", mem="4GB"),
        "junction": ResourceUsage(threads=1, runtime="4h", mem="16GB"),
    }

    @property
    def intermediate_outputs(self) -> dict[str, tuple[str, ...]]:
        return {**super().intermediate_outputs, "features": ("gtf",), "junction": ("junctions",)}

    def fasta_outputs(self):
        """pVACsplice writes the transcript sequences, which NetChop reads, once"""
        return {"fasta": f"{self.directory}/{{tumor_library}}.transcripts.fa"}

    def _get_input_files_features(self, wildcards):
        return {"features": self.parent.get_upstream_paths("features").gtf}

    def _get_output_files_features(self):
        return {"gtf": "work/features/out/features.gtf"}

    def _get_log_file_features(self):
        return {"log": "work/features/log/features.log"}

    def _get_input_files_junction(self, wildcards):
        rna = self.parent.matched(wildcards.tumor_library, "matched_rna_lib")
        alignments = self.parent.get_upstream_paths("alignments_rna", library_name=rna)
        return {
            "annotated": self.variants(wildcards.tumor_library),
            "bam": alignments.bam,
            "bai": alignments.bai,
            "strandedness": alignments.strandedness,
            "features": "work/features/out/features.gtf",
            **self.parent.reference_inputs(),
        }

    def _get_output_files_junction(self):
        return {"junctions": f"{_OUT}.junctions.tsv"}

    def _get_log_file_junction(self):
        return {"log": f"{_LOG}.junction.log"}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()
        yield "features", self.parent.get_upstream_paths("features").gtf
        yield "annotated", f"{_OUT}.normalized.vcf.gz".format(**wildcards)
        yield "annotated_tbi", f"{_OUT}.normalized.vcf.gz.tbi".format(**wildcards)
        yield "junctions", f"{_OUT}.junctions.tsv".format(**wildcards)


class PostProcessStepPart(NeoepitopeStepPart):
    """Post-processing of the epitopes of one prediction"""

    actions = ("run",)
    default_resource_usage = ResourceUsage(threads=8, runtime="6d", mem="32GB")

    def prediction(self, wildcards) -> PredictionStepPart:
        return self.parent.sub_steps[wildcards.prediction]

    def epitopes(self, wildcards, all_epitopes: bool, mhc_infix: str) -> str:
        subdir = next(d for d, infix in MHC_CLASSES.values() if infix == mhc_infix)
        table = TABLES["all" if all_epitopes else "filtered"]
        directory = self.prediction(wildcards).directory.format(**wildcards)
        return f"{directory}/{subdir}/{wildcards.tumor_library}.{mhc_infix}.{table}"


class NetChopStepPart(PostProcessStepPart):
    """Proteasomal cleavage with NetChop"""

    name = "netchop"
    published_keys = ("netchop", "netchop_md5")

    def _get_input_files_run(self, wildcards):
        prediction = self.prediction(wildcards)
        cfg = prediction.cfg.net_chop
        fasta = prediction.fasta_outputs()
        return {
            "netchop": cfg.path,
            "epitopes": self.epitopes(wildcards, cfg.all_epitopes, wildcards.mhc_class),
            "fasta": fasta.get(f"fasta.{wildcards.mhc_class}", fasta.get("fasta")).format(
                **wildcards
            ),
        }

    def _get_output_files_run(self):
        path = _post_processed("netchop")
        return {"netchop": path, "netchop_md5": path + ".md5"}

    def _get_params_run(self, wildcards):
        cfg = self.prediction(wildcards).cfg.net_chop
        return {"tool": wildcards.prediction, "method": cfg.method, "threshold": cfg.threshold}

    def _get_log_file_run(self):
        return {
            "log": "work/{tumor_library}/log/{tumor_library}.{prediction}.{mhc_class}.netchop.log"
        }


class NetStabStepPart(PostProcessStepPart):
    """Stability of the MHC class I epitopes with netMHCstabpan"""

    name = "netstab"
    published_keys = ("netstab", "netstab_md5")

    def _get_input_files_run(self, wildcards):
        cfg = self.prediction(wildcards).cfg
        if cfg.net_chop.enabled:
            path = _post_processed("netchop", wildcards.prediction, "MHC_I")
            epitopes = path.format(tumor_library=wildcards.tumor_library)
        else:
            epitopes = self.epitopes(wildcards, cfg.netmhc_stab.all_epitopes, "MHC_I")
        return {
            "container": self.container,
            "epitopes": epitopes,
            "alleles": f"{_OUT}.hla_types.txt".format(**wildcards),
        }

    def _get_output_files_run(self):
        path = _post_processed("netstab", mhc_class="MHC_I")
        return {"netstab": path, "netstab_md5": path + ".md5"}

    def _get_params_run(self, wildcards):
        cfg = self.prediction(wildcards).cfg
        return {"tool": wildcards.prediction, "lengths": cfg.class_i_epitope_length}

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.{prediction}.netstab.log"}


class SomaticNeoepitopePredictionWorkflow(BaseStep):
    """Predict neoepitopes with pVACtools"""

    name = "somatic_neoepitope_prediction"
    produces = [DataSignature(DataType.TABULAR, frozenset({"neoepitope"}))]

    config_model_class = SomaticNeoepitopePredictionConfigModel

    #: Somatic variants of tumor DNA libraries
    library_requirement = "role == 'tumor' and extraction_type == 'dna'"

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name", target="role == 'normal' and extraction_type == 'dna'"
        ),
        "matched_rna_lib": RelationshipDefinition(
            via=["donor_name", "sample_name"], target="role == 'tumor' and extraction_type == 'rna'"
        ),
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return the local paths of the combined epitopes of each enabled prediction."""
        lib = kwargs.get("library_name", "{tumor_library}")
        return {
            f"{prediction}_{key}": (
                f"work/{lib}/out/{prediction}/{COMBINED[0]}/{lib}.{COMBINED[1]}.{TABLES[key]}"
            )
            for prediction in PREDICTIONS
            if getattr(config.pvactools, prediction).enabled
            for key in PUBLISHED_TABLES
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes(
            (
                PvacToolsStepPart,
                PvacSeqStepPart,
                PvacFuseStepPart,
                PvacSpliceStepPart,
                NetChopStepPart,
                NetStabStepPart,
            )
        )
        self.check_vep_plugins()

    def check_vep_plugins(self):
        """pVACseq and pVACsplice need the VEP plugins Frameshift and Wildtype."""
        if not self.config.depends_on.variants:
            return
        vep = getattr(self.get_task_config("variants"), "vep", None)
        names = {plugin.name for plugin in getattr(vep, "plugins", None) or []}
        if missing := sorted({"Frameshift", "Wildtype"} - names):
            raise ValueError(
                f"Task {self.task_name!r}: depends_on.variants must come from VEP with the "
                f"plugins {', '.join(missing)} in vep.plugins"
            )

    @property
    def enabled_predictions(self) -> list[str]:
        return [p for p in PREDICTIONS if getattr(self.config.pvactools, p).enabled]

    def matched(self, tumor_library: str, relationship: str) -> str:
        """The library of a relationship of ``tumor_library``, or ``""``"""
        return self.library_row(tumor_library).get(relationship) or ""

    @property
    def typed_mhc_classes(self) -> set[str]:
        """The MHC classes of the HLA typing tasks; both for a task that does not declare them"""
        from snappy_pipeline.workflow_registry import WORKFLOW_REGISTRY

        result = set()
        for task in self.upstream_task_names("hla_types"):
            cls = WORKFLOW_REGISTRY[self.project.task(task).step]
            config = self.project.task_configs[task]
            classes = getattr(cls, "mhc_classes", None)
            result |= set(classes(config)) if classes else set(MHC_CLASSES)
        return result

    @listify
    def get_result_files(self):
        for prediction in self.enabled_predictions:
            part = self.sub_steps[prediction]
            for tumor in self.output_entities:
                if prediction in NEED_RNA and not self.matched(tumor, "matched_rna_lib"):
                    continue
                outputs = part._get_output_files_run()
                for key in part.published_keys["run"]:
                    yield outputs[key].replace("work/", "output/", 1).format(tumor_library=tumor)
                classes = part.mhc_classes()
                post_processed = [
                    _post_processed("netchop", prediction, MHC_CLASSES[name][1])
                    for name in classes
                    if part.cfg.net_chop.enabled
                ]
                if part.cfg.netmhc_stab.enabled and "class_i" in classes:
                    post_processed.append(_post_processed("netstab", prediction, "MHC_I"))
                for path in post_processed:
                    path = path.replace("work/", "output/", 1).format(tumor_library=tumor)
                    yield from (path, path + ".md5")
