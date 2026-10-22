"""Variant phasing with GATK or WhatsHap

Phases germline variants of pedigrees, or the combined somatic and germline variants of tumor libraries.

The ``tool`` field selects the phasing:

- ``gatk_phase_by_transmission`` runs PhaseByTransmission on the annotated variants of ``depends_on.variants``. It changes the genotypes of some variants, which is a problem when phasing *de novo* variants.
- ``gatk_read_backed_phasing`` runs ReadBackedPhasing on ``depends_on.variants`` with the alignments of the trio (``depends_on.alignments``), scattered over ``gatk_read_backed_phasing.num_jobs`` chunks of the genome. It is not fully reliable.
- ``gatk_read_backed_phasing_tumor`` runs ReadBackedPhasing on the combined somatic and germline variants of each tumor DNA library (``depends_on.combined_variants``, from ``combine_variants``) with its alignments, as pVACtools describes it for the proximal variants of neoepitope prediction.
- ``whatshap_tumor`` phases the same input with WhatsHap, a maintained read-based phasing tool, and writes the phase into the ``HP`` tag of ReadBackedPhasing (``--tag=HP``), which pVACseq reads; the genotypes stay unphased.
  WhatsHap assumes diploid data, so loss of heterozygosity, rearrangements and subclonal variants of a tumor limit it, as they limit ReadBackedPhasing.

The pedigree tools produce one result per child in a pedigree whose parents are both sequenced, the tumor tool one per tumor DNA library.

One task runs one tool.
To use both, add a ``gatk_read_backed_phasing`` task whose ``depends_on.variants`` names the ``gatk_phase_by_transmission`` task.
Only results of the GATK variant callers can be phased.
It is unclear how useful the phasing is for clinical studies.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

Here ``library_name`` is the library of the child, or the tumor library.
The output has the tags of the input plus ``phased``.
Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import os
from functools import cached_property
from typing import Any

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, select_signature

from .model import VariantPhasing as VariantPhasingConfigModel

#: Output keys and extensions of a phased VCF
_OUTPUTS = {
    "vcf": ".vcf.gz",
    "vcf_tbi": ".vcf.gz.tbi",
    "vcf_md5": ".vcf.gz.md5",
    "vcf_tbi_md5": ".vcf.gz.tbi.md5",
}

_BASE_PATH_OUT = "work/{library_name}/out/{library_name}"


class WriteTrioPedigreeStepPart(BaseStepPart):
    """Write out trio pedigree file for primary DNA sample given the index NGS library name"""

    #: Step name
    name = "write_trio_pedigree"

    #: Class available actions
    actions = ("run",)

    def _get_input_files_run(self, wildcards):
        return []

    def _get_output_files_run(self):
        return "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

    def run(self, wildcards, output):
        """Write out the pedigree information"""
        fname = self.get_output_files("run").format(**wildcards)
        os.makedirs(os.path.dirname(fname), exist_ok=True)
        libraries = self.parent.build_library_dataframe(include_background=True).set_index(
            "library_name", drop=False
        )
        child = libraries.loc[wildcards.index_ngs_library]
        family = "FAM_" + child["donor_name"]
        with open(fname, "wt") as ped_file:
            for name in (child["library_name"], child["father_name"], child["mother_name"]):
                if name == "0":
                    continue
                person = libraries.loc[name]
                sex = {"male": "1", "female": "2", "unknown": "0"}[person["sex"]]
                print(
                    "\t".join(
                        (
                            family,
                            name,
                            person["father_name"],
                            person["mother_name"],
                            sex,
                            str(person["disease_state"]),
                        )
                    ),
                    file=ped_file,
                )


class VariantPhasingBaseStep(BaseStepPart):
    """Base class of the phasing tools"""

    #: Class available actions
    actions = ("run",)

    def _yield_vcf(self, wildcards):
        """Yield the VCF to phase and its index

        Annotated variants are stored under the library of the pedigree's index, phased variants
        under the library of the child.
        """
        library_name = wildcards.library_name
        if not self.parent.variants_are_phased:
            library_name = self.parent.pedigree_index(library_name, include_background=True)
        variants = self.parent.get_upstream_paths("variants", library_name=library_name)
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}


class PhaseByTransmissionStepPart(VariantPhasingBaseStep):
    """Phasing by transmission."""

    #: Name of the step in the pipeline.
    name = "gatk_phase_by_transmission"

    published_keys = tuple(_OUTPUTS)

    @dictify
    def _get_input_files_run(self, wildcards):
        # Pedigree file required for PhaseByTransmission.
        yield (
            "ped",
            f"work/write_pedigree.{wildcards.library_name}/out/{wildcards.library_name}.ped",
        )
        yield from self._yield_vcf(wildcards)
        yield from self.parent.reference_inputs().items()

    @dictify
    def _get_output_files_run(self):
        for key, ext in _OUTPUTS.items():
            yield key, _BASE_PATH_OUT + ext

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        return {"de_novo_prior": self.config.gatk_phase_by_transmission.de_novo_prior}

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem=f"{14 * 1024}MB")


class ReadBackedPhasingStepPart(VariantPhasingBaseStep):
    """Read-backed phasing, scattered over chunks of the genome and gathered per child."""

    #: Name of the step in the pipeline.
    name = "gatk_read_backed_phasing"

    #: Class available actions
    actions = ("scatter", "run", "gather")

    extra_arguments_action = "run"

    published_keys = {"gather": tuple(_OUTPUTS)}

    temp_outputs = {"scatter": ("regions",), "run": ("vcf", "vcf_tbi")}

    resource_usage = {
        "scatter": ResourceUsage(threads=1, runtime="2m", mem="1000MB"),
        "run": ResourceUsage(threads=1, runtime="4h", mem="14336MB"),
        "gather": ResourceUsage(threads=1, runtime="4h", mem="2048MB"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.region_tpl = "work/scatter/{scatteritem}.region.bed"
        self.chunk_tpl = "work/{library_name}/chunks/{scatteritem}"

    def _scattergather(self, kind):
        return getattr(self.parent.workflow.globals.get(kind), self.name)

    def _yield_bams(self, wildcards):
        """Yield the alignments of the trio"""
        libraries = self.parent.build_library_dataframe(include_background=True)
        child = libraries[libraries["library_name"] == wildcards.library_name].iloc[0]
        if "0" in (child["father_name"], child["mother_name"]):
            return
        trio_libs = [child["library_name"], child["father_name"], child["mother_name"]]
        alignments = [
            self.parent.get_upstream_paths("alignments", library_name=lib) for lib in trio_libs
        ]
        yield "bam", [aln.bam for aln in alignments]
        yield "bai", [aln.bai for aln in alignments]

    def _get_input_files_scatter(self, wildcards):
        return {"fai": self.parent.get_upstream_paths("reference").fai}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from self._yield_bams(wildcards)
        yield from self._yield_vcf(wildcards)
        yield from self.parent.reference_inputs().items()
        yield "region", self.region_tpl.format(**wildcards)

    def _get_input_files_gather(self, wildcards):
        gather = self._scattergather("gather")
        chunk = self.chunk_tpl.replace("{scatteritem}", "{{scatteritem}}").format(**wildcards)
        return {"vcf": gather(chunk + ".vcf.gz"), "vcf_tbi": gather(chunk + ".vcf.gz.tbi")}

    def _get_params_scatter(self, wildcards):
        return {"ignore_contigs": self.config.ignore_contigs, "padding": 0}

    def _get_params_run(self, wildcards):
        return {
            "phase_quality_threshold": self.config.gatk_read_backed_phasing.phase_quality_threshold
        }

    def _get_output_files_scatter(self):
        return {"regions": self._scattergather("scatter")(self.region_tpl)}

    def _get_output_files_run(self):
        return {
            "vcf": self.chunk_tpl + ".vcf.gz",
            "vcf_tbi": self.chunk_tpl + ".vcf.gz.tbi",
        }

    @dictify
    def _get_output_files_gather(self):
        for key, ext in _OUTPUTS.items():
            yield key, _BASE_PATH_OUT + ext

    def _get_log_file_scatter(self):
        return {"log": "work/scatter/log.log"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.chunks/{scatteritem}.log"}

    def _get_log_file_gather(self):
        return super()._get_log_file_run()


class TumorPhasingStepPart(VariantPhasingBaseStep):
    """Phasing of a tumor library's combined somatic and germline variants with its alignments."""

    published_keys = tuple(_OUTPUTS)

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem="32GB")

    def _get_params_run(self, wildcards):
        return {}

    def _get_input_files_run(self, wildcards):
        library_name = wildcards.library_name
        variants = self.parent.get_upstream_paths("combined_variants", library_name=library_name)
        alignments = self.parent.get_upstream_paths("alignments", library_name=library_name)
        return {
            "vcf": variants.vcf,
            "vcf_tbi": variants.vcf_tbi,
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
        }

    @dictify
    def _get_output_files_run(self):
        for key, ext in _OUTPUTS.items():
            yield key, _BASE_PATH_OUT + ext


class TumorReadBackedPhasingStepPart(TumorPhasingStepPart):
    """GATK 3 ReadBackedPhasing of a tumor library"""

    name = "gatk_read_backed_phasing_tumor"


class TumorWhatshapStepPart(TumorPhasingStepPart):
    """WhatsHap phasing of a tumor library"""

    name = "whatshap_tumor"


_TOOL_STEP_PART = {
    "gatk_phase_by_transmission": PhaseByTransmissionStepPart,
    "gatk_read_backed_phasing": ReadBackedPhasingStepPart,
    "gatk_read_backed_phasing_tumor": TumorReadBackedPhasingStepPart,
    "whatshap_tumor": TumorWhatshapStepPart,
}

#: The tools that phase tumor libraries instead of pedigrees
_TUMOR_TOOLS = ("gatk_read_backed_phasing_tumor", "whatshap_tumor")


class VariantPhasingWorkflow(BaseStep):
    """Perform (small) variant phasing"""

    name = "variant_phasing"
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "phased"})),
        DataSignature(DataType.VARIANTS, frozenset({"combined", "phased"})),
    ]
    config_model_class = VariantPhasingConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """The variants this task reads, tagged ``phased``."""
        field = "combined_variants" if config.tool in _TUMOR_TOOLS else "variants"
        return (upstream[field].with_tags("phased"),)

    @classmethod
    def task_library_requirement(cls, config):
        """The tumor tool phases tumor DNA libraries."""
        return (
            "role == 'tumor' and extraction_type == 'dna'" if config.tool in _TUMOR_TOOLS else None
        )

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local phased-variant output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"work/{lib}/out/{lib}"
        return {"vcf": f"{prefix}.vcf.gz", "vcf_tbi": f"{prefix}.vcf.gz.tbi"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes(
            (WriteTrioPedigreeStepPart, _TOOL_STEP_PART[self.config.tool])
        )

    @cached_property
    def variants_are_phased(self) -> bool:
        """Whether ``depends_on.variants`` is a phasing task, whose results are per child."""
        dependency = self.resolve_dependency("variants")
        produced = self.project.signatures[dependency.task_name]
        return "phased" in select_signature(produced, dependency.signature).tags

    @cached_property
    def output_entities(self) -> list[str]:
        """The children whose parents are both sequenced (their DNA libraries), or the tumor
        libraries."""
        if self.config.tool in _TUMOR_TOOLS:
            return super().output_entities
        return self.children_with_both_parents()
