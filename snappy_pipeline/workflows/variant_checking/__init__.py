# -*- coding: utf-8 -*-
"""Variant checking with peddy

The step checks the germline variants of ``depends_on.variants`` against the pedigree.
The ``tool`` field is ``peddy``, which compares the sex, the relatedness and the heterozygosity of the samples with the pedigree.
The task produces one set of files per pedigree, named after the pedigree (its ``related_group``).

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

One job per task summarizes the checks of all pedigrees in ``qc_summary/out/qc_summary.tsv``, one row per sample.
A row holds the sex check, the heterozygosity and the depth with their outlier flags, the call rate and the ancestry.
The job also renders the checks for the report below ``work/qc_summary/report/html/``.
The report shows the summary, a strip chart of the metrics and the pairs of every pedigree.
Each pair has its relatedness from the PED file and from the genotypes, shown as a table and as a scatter chart.
On the scatter chart, a swapped, unrelated or duplicated sample falls off the diagonal.
"""

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedVariants
from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart

from .model import VariantChecking as VariantCheckingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

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

    #: Step name
    name = "peddy"

    #: Class available actions
    actions = ("run",)

    wrappers = {"run": "peddy/run"}
    entity_wildcard = "related_group"
    result_keys = tuple(f"{key}{suffix}" for key in _OUTPUTS for suffix in ("", "_md5"))
    report_outputs = {"run": {"html": ReportOutput("peddy", caption="peddy.rst")}}
    default_resource_usage = ResourceUsage(threads=1, runtime="10h", mem=f"{15 * 1024}MB")

    def tool_section(self):
        """peddy has no tool section: the ``extra_arguments`` are the step's own"""
        return self.config

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path to pedigree input file"""
        index = wildcards.related_group
        yield "ped", f"work/write_pedigree.{index}/out/{index}.ped"

        calling: ExpectedVariants = self.parent.get_upstream_paths(
            "variants",
            library_name=index,
        )
        yield "vcf", calling.vcf
        yield "vcf_tbi", calling.vcf_tbi

    def _get_output_files_run(self):
        prefix = "work/{related_group}/out/{related_group}"
        result = {}
        for key, ext in _OUTPUTS.items():
            result[key] = prefix + ext
            result[f"{key}_md5"] = prefix + ext + ".md5"
        return result

    def _get_log_file_run(self):
        return {"log": "work/{related_group}/log/snakemake.filter.log"}

    def _get_params_run(self, wildcards):
        """The tool's ``extra_arguments`` only, which ``with_extra_arguments`` adds"""
        return {}


class PeddySummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the peddy checks of all pedigrees"""

    report_outputs = {
        "run": {
            "html": ReportOutput(
                "Pedigree checks", caption="peddy_summary.rst", htmlindex="index.html"
            )
        }
    }

    source = ("peddy", "run")
    source_keys = {"ped_check": "ped_check", "sex_check": "sex_check", "het_check": "het_check"}
    format = "peddy"


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
        self.register_sub_step_classes((PeddyStepPart, WritePedigreeStepPart, PeddySummaryStepPart))
