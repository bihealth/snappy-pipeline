# -*- coding: utf-8 -*-
"""Implementation of the ``gene_expression_quantification`` step

This step takes as input alignments from ``ngs_mapping`` and performs gene
expression quantification. Further, the tool dupradar is run to help with
estimating library complexity, i.e. the number of PCR duplicates.
(This needs the bam files to be marked for duplicates, e.g. with samblaster.)

The output are tables of gene counts and gene lengths for each library for
a given annotation.

==========
Step Input
==========

The gene expression quantification step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step. It will use RNA-seq experiments only.

Additionally, salmon can be used to estimate expression directly from the FASTQs with no need
for prior mapping.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk`` and each read mapper
``mapper`` that the library has been aligned with, and the tool ``tool``, the
pipeline step will create a directory ``output/{mapper}.{tool}.{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting TSV files.

- ``{mapper}.{tool}.{lib_name}-{lib_pk}.tsv``

=====================
Default Configuration
=====================

Featurecounts needs a GTF with the gene model that will be used to count reads. By default,
strandedness is set to 0, i.e. unstranded. If you want stranded output, see the featurecounts
manual and set it appropriately (1 or 2, depending on the protocoll).

RSeQC's infer_experiment.py can be used to infer the strandedness. (You might want to run it
separately before featurecounts). This requires a 6-column bed file, which holds the positions
of transcripts and the strand they are originating from.

Salmon needs an pre-build index at the moment, the path should point to the directory, which
contains the needed files (e.g. sa.bin, rsd.bin, txpInfo.bin, etc).
Additionally, one can provide a gtf for the mapping between transcripts and genes.
"""

import os

from biomedsheets.shortcuts import GenericSampleSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    LinkOutStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

# Extensions
EXTENSIONS = {
    "featurecounts": {
        "tsv": ".tsv",
        "tsv_md5": ".tsv.md5",
        "summary": ".tsv.summary",
        "summary_md5": ".tsv.summary.md5",
    },
    "strandedness": {
        "tsv": ".tsv",
        "tsv_md5": ".tsv.md5",
        "decision": ".decision",
        "decision_md5": ".decision.md5",
    },
    "duplication": {
        "dupl_seq": ".seq.DupRate.xls",
        "dupl_seq_md5": ".seq.DupRate.xls.md5",
        "dupl_pos": ".pos.DupRate.xls",
        "dupl_pos_md5": ".pos.DupRate.xls.md5",
    },
    "dupradar": {"dupradar": ".dupradar.tsv", "dupradar_md5": ".dupradar.tsv.md5"},
    "rnaseqc": {
        "rnaseqc_metrics": ".metrics.tsv",
        "rnaseqc_metrics_md5": ".metrics.tsv.md5",
        "rnaseqc_meancov_low": ".meanCoverage_low.txt",
        "rnaseqc_meancov_low_md5": ".meanCoverage_low.txt.md5",
        "rnaseqc_meancov_medium": ".meanCoverage_medium.txt",
        "rnaseqc_meancov_medium_md5": ".meanCoverage_medium.txt.md5",
        "rnaseqc_meancov_high": ".meanCoverage_high.txt",
        "rnaseqc_meancov_high_md5": ".meanCoverage_high.txt.md5",
        "rnaseqc_meannorm_low": ".meanCoverageNorm_low.txt",
        "rnaseqc_meannorm_low_md5": ".meanCoverageNorm_low.txt.md5",
        "rnaseqc_meannorm_medium": ".meanCoverageNorm_medium.txt",
        "rnaseqc_meannorm_medium_md5": ".meanCoverageNorm_medium.txt.md5",
        "rnaseqc_meannorm_high": ".meanCoverageNorm_high.txt",
        "rnaseqc_meannorm_high_md5": ".meanCoverageNorm_high.txt.md5",
        "rnaseqc_gaplen_low": ".gapLengthHist_low.txt",
        "rnaseqc_gaplen_low_md5": ".gapLengthHist_low.txt.md5",
        "rnaseqc_gaplen_medium": ".gapLengthHist_medium.txt",
        "rnaseqc_gaplen_medium_md5": ".gapLengthHist_medium.txt.md5",
        "rnaseqc_gaplen_high": ".gapLengthHist_high.txt",
        "rnaseqc_gaplen_high_md5": ".gapLengthHist_high.txt.md5",
    },
    "stats": {"stats": ".read_alignment_report.tsv", "stats_md5": ".read_alignment_report.tsv.md5"},
    "salmon": {"transcript_sf": ".transcript.sf", "transcript_sf_md5": ".transcript.sf.md5"},
}

DEFAULT_CONFIG = r"""
step_config:
  gene_expression_quantification:
    path_link_in: ""   # OPTIONAL Override data set configuration search paths for FASTQ files
    tools:
    - strandedness
    - featurecounts
    - duplication
    - dupradar
    - rnaseqc
    - stats
    - salmon
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    strand: -1 # Use 0, 1 or 2 to force unstranded, forward or reverse strand
    featurecounts:
      path_annotation_gtf: REQUIRED  # REQUIRED
    strandedness:
      # needs column 6 with strand info, e.g. CCDS/15/GRCh37/CCDS.bed
      path_exon_bed: REQUIRED  # REQUIRED
      threshold: 0.85
    rnaseqc:
      rnaseqc_path_annotation_gtf: REQUIRED # REQUIRED
    dupradar:
      dupradar_path_annotation_gtf: REQUIRED  # REQUIRED
      num_threads: 8
    salmon:
      path_transcript_to_gene: REQUIRED  # REQUIRED
      path_index: REQUIRED # REQUIRED
      salmon_params: " --gcBias --validateMappings"
      num_threads: 16
""".lstrip()


class SalmonStepPart(BaseStepPart):
    """Gene expression quantification for raw data using salmon"""

    #: Step name
    name = "salmon"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/salmon.{{library_name}}/out/salmon.{{library_name}}{ext}"
        self.extensions = EXTENSIONS["salmon"]
        if (
            self.config["salmon"]["path_transcript_to_gene"] is not None
            and self.config["salmon"]["path_transcript_to_gene"] != ""
        ):
            self.extensions["gene_sf"] = ".gene.sf"
            self.extensions["gene_sf_md5"] = ".gene.sf.md5"
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config["path_link_in"],
        )

    @classmethod
    @dictify
    def get_input_files(cls, action):
        """Return input files"""
        assert action == "run"
        yield "done", "work/input_links/{library_name}/.done"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run"
        for k, v in self.extensions.items():
            yield k, self.base_path_out.format(ext=v)

    @dictify
    def _get_log_file(self, action):
        """Return mapping of log files."""
        assert action == "run"
        prefix = "work/salmon.{library_name}/log/salmon.{library_name}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            result = {
                "input": {
                    "reads_left": list(
                        sorted(self._collect_reads(wildcards, wildcards.library_name, ""))
                    )
                }
            }
            reads_right = list(
                sorted(self._collect_reads(wildcards, wildcards.library_name, "right-"))
            )
            if reads_right:
                result["input"]["reads_right"] = reads_right
            return result

        assert action == "run", "Unsupported actions"
        return args_function

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config["path_link_in"]:
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            time="04:00:00",  # 4 hours
            memory="2500M",
        )


class GeneExpressionQuantificationStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{tool}.{{library_name}}/out/" "{{mapper}}.{tool}.{{library_name}}{ext}"
        )

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            base_path = "output/{mapper}.{library_name}/out/" "{mapper}.{library_name}".format(
                **wildcards
            )
            return {
                "bam": ngs_mapping(base_path + ".bam"),
                "bai": ngs_mapping(base_path + ".bam.bai"),
            }

        assert action == "run", "Unsupported actions"
        return input_function

    def get_output_files(self, action):
        """Return output files that sub steps must return"""
        assert action == "run"
        return dict(
            zip(
                EXTENSIONS[self.name].keys(),
                expand(self.base_path_out, tool=[self.name], ext=EXTENSIONS[self.name].values()),
            )
        )

    @dictify
    def get_log_file(self, action):
        """Return mapping of log files."""
        assert action == "run"
        prefix = (
            "work/{{mapper}}.{tool}.{{library_name}}/log/{{mapper}}.{tool}.{{library_name}}".format(
                tool=self.__class__.name
            )
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class FeatureCountsStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification from RNA-seq using FeatureCounts"""

    #: Step name
    name = "featurecounts"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        return ResourceUsage(
            threads=2,
            time="1-00:00:00",  # 1 day
            memory="6700M",
        )


class StrandednessStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification from RNA-seq using FeatureCounts"""

    #: Step name
    name = "strandedness"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        return ResourceUsage(
            threads=2,
            time="12:00:00",  # 12 hours
            memory="6700M",
        )

    def get_strandedness_file(self, action):
        _ = action
        return expand(self.base_path_out, tool=[self.name], ext=[".decision"])


class QCStepPartDuplication(GeneExpressionQuantificationStepPart):

    #: Step name
    name = "duplication"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="3-00:00:00",  # 3 days
            memory="127G",
        )


class QCStepPartDupradar(GeneExpressionQuantificationStepPart):

    #: Step name
    name = "dupradar"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            time="4-00:00:00",  # 4 days
            memory="6700M",
        )


class QCStepPartRnaseqc(GeneExpressionQuantificationStepPart):

    #: Step name
    name = "rnaseqc"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="03:59:00",  # 3 hours and 59 minutes
            memory="16G",
        )


class QCStepPartStats(GeneExpressionQuantificationStepPart):

    #: Step name
    name = "stats"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="03:59:00",  # 3 hours and 59 minutes
            memory="4G",
        )


class GeneExpressionQuantificationWorkflow(BaseStep):
    """Perform gene expression quantification"""

    #: Workflow name
    name = "gene_expression_quantification"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

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
                StrandednessStepPart,
                FeatureCountsStepPart,
                QCStepPartDuplication,
                QCStepPartDupradar,
                QCStepPartRnaseqc,
                QCStepPartStats,
                LinkInStep,
                SalmonStepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    def get_strandedness_file(self, action):
        _ = action
        return self.sub_steps["strandedness"].get_strandedness_file("run")

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.{tool}.{ngs_library.name}"

        # Salmon special case
        salmon_name_pattern = "salmon.{ngs_library.name}"
        salmon_exts = EXTENSIONS["salmon"]
        if self.w_config["step_config"]["gene_expression_quantification"]["salmon"][
            "path_transcript_to_gene"
        ]:
            salmon_exts["gene_sf"] = ".gene.sf"
            salmon_exts["gene_sf_md5"] = ".gene.sf.md5"

        # TODO: too many ifs, use shortcut?
        # if fixed, please do the same for somatic_gene_fusion_calling
        all_fns = []
        for tool in self.config["tools"]:
            for sheet in filter(is_not_background, self.shortcut_sheets):
                for ngs_library in sheet.all_ngs_libraries:
                    extraction_type = ngs_library.test_sample.extra_infos["extractionType"]
                    if extraction_type.lower() == "rna":
                        if tool == "salmon":
                            fns = expand(
                                os.path.join(
                                    "output",
                                    salmon_name_pattern,
                                    "out",
                                    salmon_name_pattern + "{ext}",
                                ),
                                ngs_library=ngs_library,
                                ext=salmon_exts.values(),
                            )
                            all_fns.extend(fns)
                        else:
                            fns = expand(
                                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                                ngs_library=ngs_library,
                                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["rna"],
                                # tool=set(self.config['tools']),
                                tool=tool,
                                ext=EXTENSIONS[tool].values(),
                            )
                            all_fns.extend(fns)

        return all_fns

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "gene_expression_quantification", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for gene expression quantification",
        )
        self.ensure_w_config(
            (
                "step_config",
                "gene_expression_quantification",
                "featurecounts",
                "path_annotation_gtf",
            ),
            "Path to gtf file with annotations required for featurecounts",
        )
        self.ensure_w_config(
            ("step_config", "gene_expression_quantification", "strandedness", "path_exon_bed"),
            "Path to bed file with exon regions required for RSeQC",
        )
        self.ensure_w_config(
            (
                "step_config",
                "gene_expression_quantification",
                "rnaseqc",
                "rnaseqc_path_annotation_gtf",
            ),
            "Path to gtf file with annotations required for RNA-SeQC",
        )
        self.ensure_w_config(
            (
                "step_config",
                "gene_expression_quantification",
                "dupradar",
                "dupradar_path_annotation_gtf",
            ),
            "Path to gtf file with annotations required for dupradar",
        )
        self.ensure_w_config(
            ("step_config", "gene_expression_quantification", "salmon", "path_index"),
            "Path to directory containing salmon index files",
        )
