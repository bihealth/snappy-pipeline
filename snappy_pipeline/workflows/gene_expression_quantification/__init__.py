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

The tools ``featurecounts``, ``dupradar``, ``duplication``, ``rnaseqc`` and ``stats`` also need
the strandedness of each library. It comes from a separate task of this step with
``tool: strandedness``, named in ``depends_on.strandedness``:

.. code-block:: yaml

    tasks:
      - step: gene_expression_quantification
        name: strandedness
        config:
          depends_on:
            ngs_mapping: mapping
          tool: strandedness
          strandedness:
            path_exon_bed: /path/to/exons.bed
      - step: gene_expression_quantification
        name: counts
        config:
          depends_on:
            ngs_mapping: mapping
            strandedness: strandedness
          tool: featurecounts
          featurecounts:
            path_annotation_gtf: /path/to/genes.gtf

===========
Step Output
===========

For each RNA NGS library with name ``lib_name``/key ``lib_pk``, the
pipeline step will create a directory ``output/{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting TSV files.

- ``{lib_name}-{lib_pk}.tsv``

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
from typing import Any

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStepPart,
    LinkOutStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import ExpectedStrandedness
from .model import GeneExpressionQuantification as GeneExpressionQuantificationConfigModel
from .model import Salmon as SalmonConfigModel

#: Signature of the strandedness decision of a ``tool: strandedness`` task
STRANDEDNESS_SIGNATURE = DataSignature(DataType.QC, frozenset({"strandedness"}))

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
    "salmon": {
        "gene_sf": ".gene.sf",
        "gene_sf_md5": ".gene.sf.md5",
        "transcript_sf": ".transcript.sf",
        "transcript_sf_md5": ".transcript.sf.md5",
    },
}


class SalmonStepPart(BaseStepPart):
    """Gene expression quantification for raw data using salmon"""

    #: Step name
    name = "salmon"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: SalmonConfigModel = self.config.salmon
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"
        self.extensions = EXTENSIONS["salmon"]
        if self.config.salmon and self.config.salmon.path_transcript_to_gene:
            self.extensions["gene_sf"] = ".gene.sf"
            self.extensions["gene_sf_md5"] = ".gene.sf.md5"
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.parent.get_preprocessed_path(),
        )

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        if self.config.tool != self.name:
            return
        yield "done", "work/input_links/{library_name}/.done".format(**wildcards)
        yield "features", self.w_config.static_data_config.features.path
        if self.cfg and self.cfg.path_index:
            yield "indices", self.cfg.path_index

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run"
        tool = self.config.tool
        if self.name != tool:
            return {}
        for k, v in self.extensions.items():
            yield k, self.base_path_out.format(ext=v)

    @dictify
    def get_log_file(self, action):
        """Return mapping of log files."""
        assert action == "run"
        tool = self.config.tool
        if self.name != tool:
            return {}
        prefix = "work/{library_name}/log/{library_name}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper, including the input files"""
        result = {
            "input": {
                "reads_left": list(
                    sorted(self._collect_reads(wildcards, wildcards.library_name, ""))
                )
            }
        }
        reads_right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
        if reads_right:
            result["input"]["reads_right"] = reads_right
        result |= self.config.salmon.model_dump(by_alias=True)
        result["strand"] = self.config.strand
        return result

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        task_prefix = self.parent.task_path_prefix()
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.parent.get_preprocessed_path():
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            path = os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)
            path = task_prefix + path
            yield path

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            runtime="4h",  # 4 hours
            mem="32GB",
        )


class GeneExpressionQuantificationStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    #: Whether the tool reads the strandedness decision of the ``strandedness`` dependency.
    needs_strandedness = False

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"

    def _get_input_files_run(self, wildcards: Wildcards):
        """Resolve alignment inputs through the typed upstream contract broker."""
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "ngs_mapping", library_name=wildcards.library_name
        )
        input_files = {
            "bam": alignments.bam,
            "bai": alignments.bai,
        }
        if self.needs_strandedness:
            strandedness: ExpectedStrandedness = self.parent.get_upstream_paths(
                "strandedness", library_name=wildcards.library_name
            )
            input_files["decision"] = strandedness.decision
        return input_files

    def get_output_files(self, action):
        """Return output files that sub steps must return"""
        assert action == "run"
        tool = self.config.tool
        if self.name != tool:
            return {}
        return dict(
            zip(
                EXTENSIONS[self.name].keys(),
                expand(self.base_path_out, ext=EXTENSIONS[self.name].values()),
            )
        )

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"strand": self.config.strand}

    @dictify
    def get_log_file(self, action):
        """Return mapping of log files."""
        assert action == "run"
        tool = self.config.tool
        if self.name != tool:
            return {}
        prefix = "work/{library_name}/log/{library_name}"
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

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "features", self.w_config.static_data_config.features.path

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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
            runtime="1d",  # 1 day
            mem="6700MB",
        )


class StrandednessStepPart(GeneExpressionQuantificationStepPart):
    """Gene expression quantification from RNA-seq using FeatureCounts"""

    #: Step name
    name = "strandedness"

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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
            runtime="12h",  # 12 hours
            mem="6700MB",
        )

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if self.config.tool != self.name:
            return super()._get_params_run(wildcards)
        config = self.config.strandedness.model_dump(by_alias=True) | {"strand": self.config.strand}
        return {"config": config, "library_name": wildcards.library_name}


class QCStepPartDuplication(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "duplication"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            runtime="3d",  # 3 days
            mem="127GB",
        )


class QCStepPartDupradar(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "dupradar"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.tool != self.name:
            return
        yield "dupradar_path_annotation_gtf", self.config.dupradar.dupradar_path_annotation_gtf

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if self.config.tool != self.name:
            return super()._get_params_run(wildcards)
        return super()._get_params_run(wildcards) | {
            "num_threads": self.config.dupradar.num_threads,
        }

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            runtime="4d",  # 4 days
            mem="6700MB",
        )


class QCStepPartRnaseqc(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "rnaseqc"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.tool != self.name:
            return
        yield "reference", self.w_config.static_data_config.reference.path
        yield "rnaseqc_path_annotation_gtf", self.config.rnaseqc.rnaseqc_path_annotation_gtf

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            runtime="4h",  # 4 hours
            mem="16GB",
        )


class QCStepPartStats(GeneExpressionQuantificationStepPart):
    #: Step name
    name = "stats"

    needs_strandedness = True

    #: Class available actions
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            runtime="4h",  # 4 hours
            mem="4GB",
        )


class GeneExpressionQuantificationWorkflow(BaseStep):
    """Perform gene expression quantification"""

    #: Workflow name
    name = "gene_expression_quantification"

    config_model_class = GeneExpressionQuantificationConfigModel

    consumes = {DataSignature(DataType.RAW, frozenset({"rna"})): True}
    produces = [
        DataSignature(DataType.EXPRESSION, frozenset({"rna"})),
        STRANDEDNESS_SIGNATURE,
    ]

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local expression or strandedness output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        if signature is not None and signature.satisfies(STRANDEDNESS_SIGNATURE):
            return {"decision": f"output/{lib}/out/{lib}.decision"}
        if signature is not None and not signature.satisfies(DataSignature(DataType.EXPRESSION)):
            raise ValueError(
                f"GeneExpressionQuantificationWorkflow does not support signature: {signature}"
            )
        return {"tsv": f"output/{lib}/out/{lib}.tsv"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                StrandednessStepPart,
                FeatureCountsStepPart,
                QCStepPartDuplication,
                QCStepPartDupradar,
                QCStepPartRnaseqc,
                QCStepPartStats,
                LinkInStepPart,
                SalmonStepPart,
                LinkOutStepPart,
            )
        )
        # Inputs are resolved via get_upstream_paths() in step parts.

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        tool = self.config.tool
        name_pattern = "{library_name}"

        # Salmon special case
        salmon_name_pattern = "{library_name}"
        salmon_exts = EXTENSIONS["salmon"]
        if self.config.salmon and self.config.salmon.path_transcript_to_gene:
            salmon_exts["gene_sf"] = ".gene.sf"
            salmon_exts["gene_sf_md5"] = ".gene.sf.md5"

        df = self.build_library_dataframe()
        if df.empty:
            return []

        all_fns = []
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "")
            if extraction_type.lower() == "rna":
                if tool == "salmon":
                    fns = expand(
                        os.path.join(
                            "output",
                            salmon_name_pattern,
                            "out",
                            salmon_name_pattern + "{ext}",
                        ),
                        library_name=[library_name],
                        ext=salmon_exts.values(),
                    )
                    all_fns.extend(fns)
                else:
                    fns = expand(
                        os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                        library_name=[library_name],
                        ext=EXTENSIONS[tool].values(),
                    )
                    all_fns.extend(fns)

        return all_fns
