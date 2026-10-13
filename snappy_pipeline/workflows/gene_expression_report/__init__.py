# -*- coding: utf-8 -*-
"""Gene expression report for tumor RNA libraries

The step places the expression of each tumor RNA library in the context of the cohort.
It collects the expression tables of all tumor ``mRNA_seq`` libraries of the sample sheets, ranks the genes of each library, computes expression signatures, and plots the expression of selected genes against the cohort.
The task has no ``tool`` field.
It reads the expression tables from ``depends_on.expression``, a :ref:`step_gene_expression_quantification` task with ``tool: featurecounts``.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/``, for each RNA library with name ``{library_name}``:

- ``{library_name}/out/{library_name}.tsv``: the ranked expression
- ``{library_name}/out/{library_name}.pdf``: the expression signatures
- ``{library_name}/out/{library_name}.genes.pdf``: the expression distribution of genes in the cohort

Downstream tasks read the key ``tsv``.
"""

import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.gene_expression_quantification.model import ExpectedExpression

from .model import GeneExpressionReport as GeneExpressionReportConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".tsv",)


class GeneExpressionReportStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{{ngs_library}}/out/{{ngs_library}}{ext}"

    def _get_log_file(self, action):
        return "work/{ngs_library}/log/snakemake.gene_expression_quantification.log"


class GeneExpressionReportAggreateFeaturecounts(GeneExpressionReportStepPart):
    """Generate a dataframe holding each biosample and the location of the CNV
    calling results and the expression quantification files"""

    #: Step name
    name = "aggregate_counts"

    @listify
    def _get_input_files_run(self, wildcards):
        for sheet in filter(is_not_background, self.parent.sheets):
            for donor in sheet.bio_entities.values():
                for biosample in donor.bio_samples.values():
                    if biosample.extra_infos["isTumor"]:
                        for test_sample in biosample.test_samples.values():
                            for lib in test_sample.ngs_libraries.values():
                                # if there is more than one lib, cbioportal cannot use it
                                if lib.extra_infos["libraryType"] == "mRNA_seq":
                                    rna_library = lib.name
                                    expression: ExpectedExpression = self.parent.get_upstream_paths(
                                        "expression", library_name=rna_library
                                    )
                                    yield expression.tsv

    def _get_output_files_run(self):
        return {"tsv": "work/gene_exp.tsv"}


class GeneExpressionReportRankExpression(GeneExpressionReportStepPart):
    """Take normalized expression values and rank and tabulate genes by expression for MTK"""

    #: Step name
    name = "compute_ranks"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "tsv", "work/gene_exp.tsv"

    def _get_output_files_run(self):
        return {"tsv": self.base_path_out.format(ext=".tsv")}


class GeneExpressionReportComputeSignatures(GeneExpressionReportStepPart):
    """Take normalized expression values and compute set of signatures"""

    #: Step name
    name = "compute_signatures"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "tsv", self.base_path_out.format(ext=".tsv").format(**wildcards)

    def _get_output_files_run(self):
        return {"pdf": expand(self.base_path_out, ext=".pdf")}


class GeneExpressionReportPlotGeneDistribution(GeneExpressionReportStepPart):
    """Take normalized expression values and plot genes in cohort context"""

    #: Step name
    name = "plot_expression_distribution"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "tsv", self.base_path_out.format(ext=".tsv").format(**wildcards)

    def _get_output_files_run(self):
        return {"pdf": expand(self.base_path_out, ext=".genes.pdf")}


class GeneExpressionReportWorkflow(BaseStep):
    """Perform per-sample gene expression analysis"""

    #: Workflow name
    name = "gene_expression_report"
    produces = [DataSignature(DataType.TABULAR, frozenset({"expression_report"}))]

    config_model_class = GeneExpressionReportConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local gene-expression report output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"tsv": f"work/{lib}/out/{lib}.tsv"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                GeneExpressionReportAggreateFeaturecounts,
                GeneExpressionReportComputeSignatures,
                GeneExpressionReportRankExpression,
                GeneExpressionReportPlotGeneDistribution,
            )
        )
        # Inputs are resolved via get_upstream_paths() in step parts.

    @listify
    def get_result_files(self):
        name_pattern = "{library_name}"
        df = self.build_library_dataframe()
        if df.empty:
            return []
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "")
            if extraction_type.lower() == "rna":
                exts = EXT_VALUES + (".pdf", ".genes.pdf")
                yield from expand(
                    os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                    library_name=[library_name],
                    ext=exts,
                )
