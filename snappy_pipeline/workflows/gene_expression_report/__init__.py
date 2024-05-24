# -*- coding: utf-8 -*-
"""Implementation of the ``gene_expression_report`` step

"""

from collections import OrderedDict
import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

DEFAULT_CONFIG = r"""
step_config:
  gene_expression_report:
    path_gene_expression_quantification: ''  # REQUIRED
"""

#: Names of the files to create for the extension (snakemake output)
EXT_NAMES = ("tsv",)

#: Extensions of files to create as main payload
EXT_VALUES = (".tsv",)


class GeneExpressionReportStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{{tool}}.{{ngs_library}}/out/"
            "{{mapper}}.{{tool}}.{{ngs_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_log_file(self, action):
        _ = action
        return (
            "work/{{mapper}}.{tool}.{{ngs_library}}/log/"
            "snakemake.gene_expression_quantification.log"
        ).format(tool=self.__class__.name)


class GeneExpressionReportAggreateFeaturecounts(GeneExpressionReportStepPart):
    """Generate a dataframe holding each biosample and the location of the CNV
    calling results and the expression quantification files"""

    #: Step name
    name = "aggregate_counts"

    @listify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        gene_expression = self.parent.sub_workflows["gene_expression_quantification"]

        for sheet in filter(is_not_background, self.parent.sheets):
            for donor in sheet.bio_entities.values():
                for biosample in donor.bio_samples.values():
                    if biosample.extra_infos["isTumor"]:
                        for test_sample in biosample.test_samples.values():
                            for lib in test_sample.ngs_libraries.values():
                                # if there is more than one lib, cbioportal cannot use it
                                if lib.extra_infos["libraryType"] == "mRNA_seq":
                                    rna_library = lib.name
                                    exp_tpl = (
                                        "output/{mapper}.{tool}.{library_name}/out/"
                                        "{mapper}.{tool}.{library_name}.tsv"
                                    ).format(
                                        tool="featurecounts",
                                        mapper="star",
                                        library_name=rna_library,
                                    )
                                    exp_file = gene_expression(exp_tpl)
                                    yield exp_file

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "tsv", "work/gene_exp.tsv"


class GeneExpressionReportRankExpression(GeneExpressionReportStepPart):
    """Take normalized expression values and rank and tabulate genes by expression for MTK"""

    #: Step name
    name = "compute_ranks"

    @dictify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "tsv", "work/gene_exp.tsv"

    def get_output_files(self, action):
        """Return output files that sub steps must return"""
        # Validate action
        self._validate_action(action)
        return dict(zip(EXT_NAMES, expand(self.base_path_out, ext=EXT_VALUES)))


class GeneExpressionReportComputeSignatures(GeneExpressionReportStepPart):
    """Take normalized expression values and compute set of signatures"""

    #: Step name
    name = "compute_signatures"

    def get_output_files(self, action):
        """Return output files that sub steps must return"""
        # Validate action
        self._validate_action(action)
        return {"pdf": expand(self.base_path_out, ext=".pdf")}


class GeneExpressionReportPlotGeneDistribution(GeneExpressionReportStepPart):
    """Take normalized expression values and plot genes in cohort context"""

    #: Step name
    name = "plot_expression_distribution"

    def get_output_files(self, action):
        """Return output files that sub steps must return"""
        # Validate action
        self._validate_action(action)
        return {"pdf": expand(self.base_path_out, ext=".genes.pdf")}


class GeneExpressionReportWorkflow(BaseStep):
    """Perform per-sample gene expression analysis"""

    #: Workflow name
    name = "gene_expression_report"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                GeneExpressionReportAggreateFeaturecounts,
                GeneExpressionReportComputeSignatures,
                GeneExpressionReportRankExpression,
                GeneExpressionReportPlotGeneDistribution,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        if self.config["path_gene_expression_quantification"]:
            self.register_sub_workflow(
                "gene_expression_quantification", self.config["path_gene_expression_quantification"]
            )

    @listify
    def get_result_files(self):
        name_pattern = "{mapper}.{tool}.{ngs_library.name}"
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    for _test_sample in bio_sample.test_samples.values():
                        ngs_library = bio_sample.rna_ngs_library
                        if ngs_library is None:
                            break

                        exts = EXT_VALUES + (".pdf", ".genes.pdf")
                        yield from expand(
                            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                            ngs_library=ngs_library,
                            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["rna"],
                            tool="featurecounts",
                            ext=exts,
                        )
