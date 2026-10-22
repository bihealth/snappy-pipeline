# -*- coding: utf-8 -*-
"""Variant checking with peddy

Checks the germline variants of ``depends_on.variants`` against the pedigree.
The ``tool`` field is ``peddy``, which compares the sex, the relatedness and the heterozygosity of the samples with the pedigree.
The task produces one set of files per cohort, named after the index library.

Outputs
-------

Below ``tasks/<task>/output/``, each file with a ``.md5`` file:

- ``{library_name}/out/{library_name}.ped_check.csv``
- ``{library_name}/out/{library_name}.sex_check.csv``
- ``{library_name}/out/{library_name}.het_check.csv``
- ``{library_name}/out/{library_name}.background_pca.json``
- ``{library_name}/out/{library_name}.html``
- ``{library_name}/out/{library_name}.ped``

Downstream tasks read ``ped_check`` through ``depends_on``.
"""

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants

from .model import VariantChecking as VariantCheckingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Available tools for checking variants
VARIANT_CHECKERS = "peddy"

#: Output keys and extensions of the peddy step
_OUTPUTS = {
    "background_pca": ".background_pca.json",
    "het_check": ".het_check.csv",
    "html": ".html",
    "ped_check": ".ped_check.csv",
    "ped": ".ped",
    "sex_check": ".sex_check.csv",
}


class PeddyStepPart(BaseStepPart):
    """Compute variant statistics using peddy"""

    report_outputs = {"run": {"html": ReportOutput("peddy", caption="peddy.rst")}}

    #: Step name
    name = "peddy"

    #: Class available actions
    actions = ("run",)

    def tool_section(self):
        """peddy has no tool section: the ``extra_arguments`` are the step's own"""
        return self.config

    result_keys = tuple(f"{key}{suffix}" for key in _OUTPUTS for suffix in ("", "_md5"))
    entity_wildcard = "index_ngs_library"
    default_resource_usage = ResourceUsage(threads=1, runtime="10h", mem=f"{15 * 1024}MB")

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path to pedigree input file"""
        index = wildcards.index_ngs_library
        yield "ped", f"work/write_pedigree.{index}/out/{index}.ped"

        calling: ExpectedGermlineVariants = self.parent.get_upstream_paths(
            "variants",
            library_name=index,
        )
        yield "vcf", calling.vcf
        yield "vcf_tbi", calling.vcf_tbi

    def _get_output_files_run(self):
        prefix = "work/{index_ngs_library}/out/{index_ngs_library}"
        result = {}
        for key, ext in _OUTPUTS.items():
            result[key] = prefix + ext
            result[f"{key}_md5"] = prefix + ext + ".md5"
        return result

    def _get_params_run(self, wildcards):
        """The tool's ``extra_arguments`` only, which ``with_extra_arguments`` adds"""
        return {}

    def _get_log_file_run(self):
        return {"log": "work/{index_ngs_library}/log/snakemake.filter.log"}


class VariantCheckingWorkflow(BaseStep):
    """Perform germline variant checking"""

    #: Workflow name
    name = "variant_checking"
    produces = [DataSignature(DataType.QC, frozenset({"pedigree_check"}))]
    config_model_class = VariantCheckingConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return the ``work/`` paths of the variant-checking files for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"ped_check": f"work/{lib}/out/{lib}.ped_check.csv"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((PeddyStepPart, WritePedigreeStepPart))

        # Copy over "tools" setting from ngs_mapping/variant_calling if not set here
