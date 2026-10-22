# -*- coding: utf-8 -*-
"""Somatic mutational signatures

The step tabulates the somatic variants of each tumor and deconstructs them into the known mutational signatures with the R package deconstructSigs.
The section ``deconstruct_sigs`` sets the reference signatures (``cosmic`` or ``nature2013``) and the normalization of the trinucleotide counts (``tri_counts_method``).
The step has no ``tool`` field.
It reads the bgzip-compressed and indexed VCF files of a somatic variant calling task.
The task gives one result per tumor DNA library.

Outputs
-------

The task links the result file ``{tumor_library}/out/deconstruct_sigs.{tumor_library}.tsv`` with a ``.md5`` file below ``tasks/<task>/output/`` as a symlink to ``work/``.
The table gives the fraction of variants that each signature explains.
The plot ``deconstruct_sigs.{tumor_library}.pdf`` is written to the same directory in ``work/`` and is no result file.

Downstream tasks read the key ``tsv``.
"""

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ReportOutput
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import ResourceUsage
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants

from .model import SomaticVariantSignatures as SomaticVariantSignaturesConfigModel

__author__ = "Clemens Messerschmidt"


# Default configuration variant_signatures


class SignaturesStepPart(BaseStepPart):
    """Base class for signature classes"""

    #: Class available actions
    actions = ("run",)

    entity_wildcard = "tumor_library"

    default_resource_usage = ResourceUsage(threads=2, runtime="1h", mem=f"{7 * 1024 * 2}MB")

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{self.name}.{{tumor_library}}.log"}


class TabulateVariantsStepPart(SignaturesStepPart):
    """Tabulate mutation from VCF"""

    #: Step name
    name = "tabulate_vcf"

    intermediate_outputs = {"run": ("tsv",)}

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path to input file"""
        variants: ExpectedSomaticVariants = self.parent.get_upstream_paths(
            "variants", library_name=wildcards.tumor_library
        )
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi

    def _get_output_files_run(self):
        return {"tsv": "work/{tumor_library}/out/tabulate_vcf.{tumor_library}.tsv"}

    def _get_params_run(self, wildcards):
        """The tumor library, whose sample the wrapper tabulates, with or without matched normal"""
        return {"tumor_library": wildcards.tumor_library}


class DeconstructSigsStepPart(SignaturesStepPart):
    """Use deconstructSigs R package to identify signatures from tables"""

    report_outputs = {
        "run": {"pdf": ReportOutput("deconstructSigs", caption="deconstruct_sigs.rst")}
    }

    #: Step name
    name = "deconstruct_sigs"

    result_keys = ("tsv", "tsv_md5")

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files to deconstruct signatures"""
        yield (
            "tsv",
            f"work/{wildcards.tumor_library}/out/tabulate_vcf.{wildcards.tumor_library}.tsv",
        )

    def _get_params_run(self, wildcards):
        """The options of ``whichSignatures``"""
        return self.config.deconstruct_sigs.model_dump(mode="json")

    def _get_output_files_run(self):
        base = "work/{tumor_library}/out/deconstruct_sigs.{tumor_library}"
        return {"tsv": base + ".tsv", "tsv_md5": base + ".tsv.md5", "pdf": base + ".pdf"}

    def _get_log_file(self, action):
        # RWrapper keeps the R script next to the log
        prefix = "work/{tumor_library}/log/deconstruct_sigs.{tumor_library}"
        return {"log": prefix + ".log", "script": prefix + ".R"}


class SomaticVariantSignaturesWorkflow(BaseStep):
    """Perform somatic variant signatures"""

    #: Workflow name
    name = "somatic_variant_signatures"
    produces = [DataSignature(DataType.TABULAR, frozenset({"signatures"}))]

    library_requirement = "kind == 'cancer'"

    config_model_class = SomaticVariantSignaturesConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local signature output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"tsv": f"work/{lib}/out/deconstruct_sigs.{lib}.tsv"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((TabulateVariantsStepPart, DeconstructSigsStepPart))
