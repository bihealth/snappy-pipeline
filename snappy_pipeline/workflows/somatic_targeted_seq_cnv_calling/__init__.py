# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_target_seq_cnv_calling`` step

This step allows for the detection of CNV events for cancer samples from targeted sequenced (e.g.,
exomes or large panels).  The wrapped tools start from the aligned reads (thus off ``ngs_mapping``)
and generate CNV calls for somatic variants.

The wrapped tools implement different strategies.  Some work "reference free" and just use the
somatic BAM files for their input, some work in "matched cancer normal mode" and need the cancer
and normal BAM files, others again need both normal and cancer BAM files, and additionally a
set of non-cancer BAM files for their background.

==========
Step Input
==========

Gene somatic CNV calling for targeted sequencing starts off the aligned reads, i.e.,
``ngs_mapping``.

===========
Step Output
===========

Generally, the following links are generated to ``output/``.

.. note:: Tool-Specific Output

    As the only integrated tool is cnvkit at the moment, the output is very tailored to the result
    of this tool.  In the future, this section will contain "common" output and tool-specific
    output sub sections.

- ``{mapper}.cnvkit.{lib_name}-{lib_pk}/out/``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.bed``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.seg``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.vcf.gz``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.cnvkit.{lib_name}-{lib_pk}/report``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.diagram.pdf``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.scatter.pdf``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.heatmap.pdf``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.heatmap.chr1.pdf``
    - ...
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.scatter.chrX.pdf``
- ``{mapper}.cnvkit.{lib_name}-{lib_pk}/report``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.breaks.txt``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.genemetrics.txt``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.gender.txt``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.metrics.txt``
    - ``{mapper}.cnvkit.{lib_name}-{lib_pk}.segmetrics.txt``

For example:

::

    output/
    |-- bwa.cnvkit.P001-T1-DNA1-WES1-000007
    |   `-- out
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-000007.bed
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-000007.seg
    |       `-- bwa.cnvkit.P001-T1-DNA1-WES1-000007.vcf
    |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016
    |   `-- report
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.diagram.pdf
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.heatmap.pdf
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.scatter.pdf
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.heatmap.chr1.pdf
    |       |-- ...
    |       `-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.scatter.chrX.pdf
    |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016
    |   `-- report
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.breaks.txt
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.genemetrics.txt
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.gender.txt
    |       |-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.metrics.txt
    |       `-- bwa.cnvkit.P002-T1-DNA1-WES1-000016.segmetrics.txt
    [...]

Note that tool ``cnvetti`` doesn't follow the snappy convention above:
the tool name is followed by an underscore & the action, where the action is one of ``coverage``, ``segment`` and ``postprocess``.
For example, the output directory would contain a directory named ``bwa.cnvetti_coverage.P002-T1-DNA1-WES1-000016``.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_targeted_seq_cnv_calling.rst

=====================================
Available Somatic Targeted CNV Caller
=====================================

- ``cnvkit``

"""

from collections import OrderedDict
from itertools import chain
import os
import os.path
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the somatic_targeted_seq_cnv_calling step
DEFAULT_CONFIG = r"""
# Default configuration somatic_targeted_seq_cnv_calling
step_config:
  somatic_targeted_seq_cnv_calling:
    tools: ['cnvkit']  # REQUIRED - available: 'cnvkit', 'copywriter', 'cnvetti_on_target' and 'cnvetti_off_target'
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    cnvkit:
      path_target: REQUIRED             # Usually ../panel_of_normals/output/cnvkit.target/out/cnvkit.target.bed
      path_antitarget: REQUIRED         # Usually ../panel_of_normals/output/cnvkit.antitarget/out/cnvkit.antitarget.bed
      path_panel_of_normals: REQUIRED   # Usually ../panel_of_normals/output/{mapper}.cnvkit.create_panel/out/{mapper}.cnvkit.panel_of_normals.cnn
      plot: True                        # Generate plots (very slow)
      min_mapq: 0                       # [coverage] Mininum mapping quality score to count a read for coverage depth
      count: False                      # [coverage] Alternative couting algorithm
      gc_correction: True               # [fix] Use GC correction
      edge_correction: True             # [fix] Use edge correction
      rmask_correction: True            # [fix] Use rmask correction
      # BCBIO uses
      # seg_method: haar
      # seg_threshold: 0.0001
      # -- OR
      # seg_method: cbs
      # seg_threshold: 0.000001
      segmentation_method: cbs          # [segment] One of cbs, flasso, haar, hmm, hmm-tumor, hmm-germline, none
      segmentation_threshold: 0.000001  # [segment] Significance threshold (hmm methods: smoothing window size)
      drop_low_coverage: False          # [segment, call, genemetrics] Drop very low coverage bins
      drop_outliers: 10                 # [segment] Drop outlier bins (0 for no outlier filtering)
      smooth_cbs: True                  # [segment] Additional smoothing of CBS segmentation (WARNING- not the default value)
      center: ""                        # [call] Either one of mean, median, mode, biweight, or a constant log2 ratio value.
      filter: ampdel                    # [call] One of ampdel, cn, ci, sem (merging segments flagged with the specified filter), "" for no filtering
      calling_method: threshold         # [call] One of threshold, clonal, none
      call_thresholds: "-1.1,-0.25,0.2,0.7" # [call] Thresholds for calling integer copy number
      ploidy: 2                         # [call] Ploidy of sample cells
      purity: 0                         # [call] Estimated tumor cell fraction (0 for discarding tumor cell purity)
      gender: ""                        # [call, diagram] Specify the chromosomal sex of all given samples as male or female. Guess when missing
      male_reference: False             # [call, diagram] Create male reference
      diagram_threshold: 0.5            # [diagram] Copy number change threshold to label genes
      diagram_min_probes: 3             # [diagram] Min number of covered probes to label genes
      shift_xy: True                    # [diagram] Shift X & Y chromosomes according to sample sex
      breaks_min_probes: 1              # [breaks] Min number of covered probes for a break inside the gene
      genemetrics_min_probes: 3         # [genemetrics] Min number of covered probes to consider a gene
      genemetrics_threshold: 0.2        # [genemetrics] Min abs log2 change to consider a gene
      genemetrics_alpha: 0.05           # [genemetrics] Significance cutoff
      genemetrics_bootstrap: 100        # [genemetrics] Number of bootstraps
      segmetrics_alpha: 0.05            # [segmetrics] Significance cutoff
      segmetrics_bootstrap: 100         # [segmetrics] Number of bootstraps
      smooth_bootstrap: False           # [segmetrics] Smooth bootstrap results
    copywriter:
      path_target_regions: REQUIRED # REQUIRED
      bin_size: 20000 # TODO: make actually configurable
      plot_genes: REQUIRED  # Path to civic annotation
      genome: hg19          # Could be hg38 (consider setting prefix to 'chr' when using GRCh38.v1)
      features: EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
      prefix: ''
      nThread: 8
    cnvetti_on_target:
      path_target_regions: REQUIRED # REQUIRED
    cnvetti_off_target:
      path_target_regions: REQUIRED # REQUIRED
      window_length: 20000
"""

#: JSON key for "isCancer"
KEY_IS_CANCER = "isCancer"

#: Value for "libraryType" is whole exome sequencing
VALUE_WES = "WES"

#: Value for "libraryType" is panel sequencing
VALUE_PANEL = "Panel-seq"

#: Values for targeted sequencing
VALUES_TARGETED_SEQ = (VALUE_WES, VALUE_PANEL)

#: Standard key/extension values for BCF files
BCF_KEY_EXTS = (
    ("bcf", ".bcf"),
    ("bcf_md5", ".bcf.md5"),
    ("bcf_csi", ".bcf.csi"),
    ("bcf_csi_md5", ".bcf.csi.md5"),
)


class SomaticTargetedSeqCnvCallingStepPart(BaseStepPart):
    """Shared code for all caller classes in somatic_targeted_seq_cnv_calling"""

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.library_name]
        return pair.normal_sample.dna_ngs_library.name


class CnvettiStepPartBase(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using CNVetti; shared code.

    We group the CNVetti pipeline into three steps:

    ``coverage``
        Pompute target-wise normalized coverage of matched tumor and normal sample, the result is
        a BCF file describing the (log)fold coverage of the tumor in relation to the normal.
    ``segment``
        Perform segmentation of the (log)fold coverage.
    ``postprocess``
        Postprocessing of the segmentation, annotation with copy state and gene-wise coverage.
    """

    #: Class available actions
    actions = ("coverage", "segment", "postprocess")

    def __init__(self, parent):
        super().__init__(parent)
        # Per-action name pattern to use in paths
        self.name_pattern = "{{mapper}}.%(name)s_{action}.{{library_name}}" % {"name": self.name}

    def get_input_files(self, action):
        """Return input function for the given action.

        Actually delegates to the appropriate ``_get_input_files_{action}`` function.  The
        "coverage" action takes as input the BAI-indexed BAM files of the matched tumor/normal
        pairs, the other actions take as input the output of the previous actions.
        """
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{action}".format(action=action))()

    def _get_input_files_coverage(self):
        @dictify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shortcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            tumor_base_path = (
                "output/{mapper}.{library_name}/out/" "{mapper}.{library_name}"
            ).format(**wildcards)
            yield "normal_bam", ngs_mapping(normal_base_path + ".bam")
            yield "normal_bai", ngs_mapping(normal_base_path + ".bam.bai")
            yield "tumor_bam", ngs_mapping(tumor_base_path + ".bam")
            yield "tumor_bai", ngs_mapping(tumor_base_path + ".bam.bai")

        return input_function

    def _get_input_files_segment(self):
        @dictify
        def input_function(wildcards):
            for key, value in self._get_output_files_coverage().items():
                yield key, value.format(**wildcards)

        return input_function

    def _get_input_files_postprocess(self):
        @dictify
        def input_function(wildcards):
            for key, value in self._get_output_files_segment().items():
                yield key, value.format(**wildcards)

        return input_function

    def get_output_files(self, action):
        """Return input function for the given action.

        Actually delegates to the appropriate ``_get_input_files_{action}`` function; refer to
        documentation of the individual functions for more details.
        """
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_output_files_{action}".format(action=action))()

    @dictify
    def _get_output_files_coverage(self):
        """The "coverage" action creates a BCF file (CSI+MD5 files) with an
        entry for each target.
        """
        name_pattern = self.name_pattern.format(action="coverage")
        for key, ext in BCF_KEY_EXTS:
            yield key, os.path.join("work", name_pattern, "out", name_pattern + ext)

    @dictify
    def _get_output_files_segment(self):
        """The "segment" action creates a BCF file (CSI+MD5 files) with an entry for each target
        (infix ``.targets``) and also for each segment (infix ``.segments``).
        """
        name_pattern = self.name_pattern.format(action="segment")
        for infix in ("targets", "segments"):
            for key, ext in BCF_KEY_EXTS:
                name = "{}_{}".format(infix, key)
                yield name, os.path.join(
                    "work", name_pattern, "out", name_pattern + "." + infix + ext
                )

    @dictify
    def _get_output_files_postprocess(self):
        """The "postprocess" action creates the following text files (consumed by the export to
        cBioPortal):

        ``bins``
            The per-bin log-fold change information.
        ``segments``
            Per-segment log-fold change information.
        ``gene_call``
            Per-gene pseudo-GISTIC scores.
        ``gene_log2``
            Per-gene log-fold change information.
        """
        name_pattern = self.name_pattern.format(action="postprocess")
        for infix in ("targets", "targets_segmented", "segments", "gene_call", "gene_log2"):
            for key, ext in (("txt", ".txt"), ("md5", ".txt.md5")):
                name = "{}_{}".format(infix, key)
                yield name, os.path.join(
                    "work", name_pattern, "out", name_pattern + "_" + infix + ext
                )

    def check_config(self):
        """Check configuration"""
        if self.name not in self.config["tools"]:
            return  # skip check
        self.parent.ensure_w_config(
            ("step_config", "somatic_targeted_seq_cnv_calling", self.name, "path_target_regions"),
            "Path to target regions is missing for {}".format(self.name),
        )

    @dictify
    def _get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)
        name_pattern = self.name_pattern.format(action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, os.path.join("work", name_pattern, "log", name_pattern + ext)

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
            time="1-00:00:00",  # 1 day
            memory="7500M",
        )


class CnvettiOffTargetStepPart(CnvettiStepPartBase):
    """Perform somatic targeted CNV calling using CNVetti with off-target reads."""

    #: Step name
    name = "cnvetti_off_target"


class CnvettiOnTargetStepPart(CnvettiStepPartBase):
    """Perform somatic targeted CNV calling using CNVetti with on-target reads."""

    #: Step name
    name = "cnvetti_on_target"


def expand_id(*args):
    """Returns a dict that can be passed into expand to get identity for the given values

    ::

        >> expand_id(('foo', 'bar'))
        {'foo': ['{foo}'], 'bar': ['{bar}']}

    """
    return {key: ["{{{}}}".format(key)] for key in args}


def format_id(*args):
    """Returns a dict that can be passed into format to get identity for the given values

    ::

        >> expand_id(('foo', 'bar'))
        {'foo': '{foo}', 'bar': '{bar}'}

    """
    return {key: "{{{}}}".format(key) for key in args}


class CnvKitStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using cnvkit"""

    #: Step name
    name = "cnvkit"

    #: Class available actions
    actions = (
        "coverage",
        "fix",
        "segment",
        "call",
        "postprocess",
        "export",
        "plot",
        "report",
    )

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "plot": ResourceUsage(
            threads=1,
            time="08:00:00",  # 1 day
            memory=f"{30 * 1024}M",
        ),
        "coverage": ResourceUsage(
            threads=8,
            time="08:00:00",  # 8 hours
            memory=f"{16 * 1024}M",
        ),
        "default": ResourceUsage(
            threads=1,
            time="03:59:59",  # 1 day
            memory=f"{int(7.5 * 1024)}M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)

    def check_config(self):
        """Check configuration for cnvkit"""
        if "cnvkit" not in self.config["tools"]:
            return  # cnvkit not enabled, skip
        self.parent.ensure_w_config(
            ("step_config", "somatic_targeted_seq_cnv_calling", "cnvkit", "path_target"),
            "Path to target regions is missing for cnvkit",
        )
        self.parent.ensure_w_config(
            ("step_config", "somatic_targeted_seq_cnv_calling", "cnvkit", "path_antitarget"),
            "Path to antitarget regions is missing for cnvkit",
        )
        self.parent.ensure_w_config(
            ("step_config", "somatic_targeted_seq_cnv_calling", "cnvkit", "path_panel_of_normals"),
            "Path to panel of normals (reference) is missing for cnvkit",
        )

    def get_input_files(self, action):
        """Return input paths input function, dependent on rule"""
        # Validate action
        self._validate_action(action)
        method_mapping = {
            "coverage": self._get_input_files_coverage,
            "call": self._get_input_files_call,
            "fix": self._get_input_files_fix,
            "segment": self._get_input_files_segment,
            "postprocess": self._get_input_files_postprocess,
            "export": self._get_input_files_export,
            "plot": self._get_input_files_plot,
            "report": self._get_input_files_report,
        }
        return method_mapping[action]

    def _get_input_files_coverage(self, wildcards):
        # BAM/BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam.bai"),
        }
        return input_files

    @staticmethod
    def _get_input_files_fix(wildcards):
        tpl_base = "{mapper}.cnvkit.{library_name}"
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".{target}coverage.cnn"
        input_files = {
            "target": tpl.format(target="target", **wildcards),
            "antitarget": tpl.format(target="antitarget", **wildcards),
        }
        return input_files

    @staticmethod
    def _get_input_files_segment(wildcards):
        cnr_pattern = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cnr"
        input_files = {"cnr": cnr_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_call(wildcards):
        segment_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.segment.cns"
        )
        input_files = {"segment": segment_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_postprocess(wildcards):
        segment_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.segment.cns"
        )
        call_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.call.cns"
        )
        input_files = {
            "segment": segment_pattern.format(**wildcards),
            "call": call_pattern.format(**wildcards),
        }
        return input_files

    @staticmethod
    def _get_input_files_export(wildcards):
        cns_pattern = (
            "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.call.cns"
        )
        input_files = {"cns": cns_pattern.format(**wildcards)}
        return input_files

    @staticmethod
    def _get_input_files_plot(wildcards):
        tpl = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{ext}"
        input_files = {
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def _get_input_files_report(self, wildcards):
        tpl = "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{ext}"
        input_files = {
            "target": tpl.format(ext="targetcoverage.cnn", **wildcards),
            "antitarget": tpl.format(ext="antitargetcoverage.cnn", **wildcards),
            "cnr": tpl.format(ext="cnr", **wildcards),
            "cns": tpl.format(ext="call.cns", **wildcards),
        }
        return input_files

    def get_output_files(self, action):
        """Return output files for the given action"""
        if action == "coverage":
            return self._get_output_files_coverage()
        elif action == "fix":
            return self._get_output_files_fix()
        elif action == "segment":
            return self._get_output_files_segment()
        elif action == "call":
            return self._get_output_files_call()
        elif action == "postprocess":
            return self._get_output_files_postprocess()
        elif action == "export":
            return self._get_output_files_export()
        elif action == "plot":
            return self._get_output_files_plot()
        elif action == "report":
            return self._get_output_files_report()
        else:
            self._validate_action(action)

    @staticmethod
    def _get_output_files_coverage():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
            output_files[target + "_md5"] = output_files[target] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_fix():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cnr")
        return {"ratios": tpl, "ratios_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_segment():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".segment.cns")
        return {"segments": tpl, "segments_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_call():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".call.cns")
        return {"calls": tpl, "calls_md5": tpl + ".md5"}

    @staticmethod
    def _get_output_files_postprocess():
        name_pattern = "{mapper}.cnvkit.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cns")
        return {"final": tpl, "final_md5": tpl + ".md5"}

    @dictify
    def _get_output_files_plot(self):
        plots = (("diagram", "pdf"), ("heatmap", "pdf"), ("scatter", "png"))
        chrom_plots = (("heatmap", "pdf"), ("scatter", "png"))
        chroms = list(chain(range(1, 23), ["X", "Y"]))
        output_files = {}
        # Yield file name pairs for global plots
        tpl = (
            "work/{{mapper}}.cnvkit.{{library_name}}/report/"
            "{{mapper}}.cnvkit.{{library_name}}.{plot}.{ext}"
        )
        for plot, ext in plots:
            output_files[plot] = tpl.format(plot=plot, ext=ext)
            output_files[plot + "_md5"] = output_files[plot] + ".md5"
        # Yield file name pairs for the chromosome-wise plots
        tpl_chrom = (
            "work/{{mapper}}.cnvkit.{{library_name}}/report/"
            "{{mapper}}.cnvkit.{{library_name}}.{plot}.chr{chrom}.{ext}"
        )
        for plot, ext in chrom_plots:
            for chrom in chroms:
                key = "{plot}_chr{chrom}".format(plot=plot, chrom=chrom)
                output_files[key] = tpl_chrom.format(plot=plot, ext=ext, chrom=chrom)
                output_files[key + "_md5"] = output_files[key] + ".md5"
        return output_files

    @staticmethod
    def _get_output_files_export():
        exports = (("bed", "bed"), ("seg", "seg"), ("vcf", "vcf.gz"), ("tbi", "vcf.gz.tbi"))
        output_files = {}
        tpl = (
            "work/{{mapper}}.cnvkit.{{library_name}}/out/"
            "{{mapper}}.cnvkit.{{library_name}}.{ext}"
        )
        for export, ext in exports:
            output_files[export] = tpl.format(export=export, ext=ext)
            output_files[export + "_md5"] = output_files[export] + ".md5"
        return output_files

    @dictify
    def _get_output_files_report(self):
        reports = ("breaks", "genemetrics", "segmetrics", "sex", "metrics")
        output_files = {}
        tpl = (
            "work/{{mapper}}.cnvkit.{{library_name}}/report/"
            "{{mapper}}.cnvkit.{{library_name}}.{report}.txt"
        )
        for report in reports:
            output_files[report] = tpl.format(report=report)
            output_files[report + "_md5"] = output_files[report] + ".md5"
        return output_files

    def get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{{mapper}}.cnvkit.{{library_name}}/log/"
            "{{mapper}}.cnvkit.{action}.{{library_name}}"
        ).format(action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        log_files = {}
        for key, ext in key_ext:
            log_files[key] = prefix + ext
            log_files[key + "_md5"] = prefix + ext + ".md5"
        return log_files

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        if action == "plot":
            return self.resource_usage_dict.get("plot")
        elif action == "coverage":
            return self.resource_usage_dict.get("coverage")
        else:
            return self.resource_usage_dict.get("default")


class CopywriterStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using CopywriteR"""

    #: Step name
    name = "copywriter"

    #: Class available actions
    actions = ("prepare", "run", "call")

    #: Actions for which there are input and output methods available
    actions_w_in_out = ("run", "call")

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "prepare": ResourceUsage(
            threads=1,
            time="02:00:00",  # 2 hours
            memory="4000M",
        ),
        "run": ResourceUsage(
            threads=2,
            time="16:00:00",  # 16 hours
            memory="80000M",
        ),
        "call": ResourceUsage(
            threads=8,
            time="03:59:00",  # 3 hours and 59 minutes
            memory="8000M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.copywriter.{{library_name}}/out/"
            "{{mapper}}.copywriter.{{library_name}}{ext}"
        )

    def get_input_files(self, action):
        # Validate action
        msg = "Option available only for actions 'run' and 'call'."
        assert action in self.actions_w_in_out, msg

        def input_function_run(wildcards):
            """Helper wrapper function"""
            # Get shortcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            tumor_base_path = (
                "output/{mapper}.{library_name}/out/" "{mapper}.{library_name}"
            ).format(**wildcards)
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            }

        def input_function_call(wildcards):
            tpl = "work/{mapper}.copywriter.{library_name}/CNAprofiles/".format(**wildcards)
            exts = {
                "input": "input.Rdata",
                "segment": "segment.Rdata",
                "counts": "read_counts.txt",
                "log2": "log2_read_counts.igv",
            }
            input_files = {}
            for k, v in exts.items():
                input_files[k] = tpl + v
            return input_files

        if action == "run":
            return input_function_run
        if action == "call":
            return input_function_call

    @dictify
    def get_output_files(self, action):
        # Validate action
        msg = "Option available only for actions 'run' and 'call'."
        assert action in self.actions_w_in_out, msg

        exts = {}
        tpl = ""
        if action == "run":
            exts = {
                "input": "input.Rdata",
                "segment": "segment.Rdata",
                "counts": "read_counts.txt",
                "log2": "log2_read_counts.igv",
            }
            tpl = "work/{mapper}.copywriter.{library_name}/CNAprofiles/"
        if action == "call":
            exts = {
                "bins_txt": "bins.txt",
                "gene_call_txt": "gene_call.txt",
                "gene_log2_txt": "gene_log2.txt",
                "segments_txt": "segments.txt",
            }
            tpl = "work/{mapper}.copywriter.{library_name}/out/{mapper}.copywriter.{library_name}_"
        output_files = {}
        for k, v in exts.items():
            output_files[k] = tpl + v
            output_files[k + "_md5"] = tpl + v + ".md5"
        return output_files

    def check_config(self):
        """Check configuration"""
        if "copywriter" not in self.config["tools"]:
            return  # skip
        self.parent.ensure_w_config(
            (
                "step_config",
                "somatic_targeted_seq_cnv_calling",
                "copywriter",
                "path_target_regions",
            ),
            "Path to target regions is missing",
        )

    @dictify
    def _get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )

        if action == "prepare":
            log_file = "work/copywriter.prepare/log/snakemake.log"
            yield "log", log_file
            yield "log_md5", log_file + ".md5"
        elif action in ("call", "run"):
            tpl = (
                "work/{mapper}.copywriter.{library_name}/log/{mapper}.copywriter.{library_name}."
                + action
            )
            for key, ext in key_ext:
                yield key, tpl + ext

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return self.resource_usage_dict.get(action)


class SomaticTargetedSeqCnvCallingWorkflow(BaseStep):
    """Perform somatic targeted sequencing CNV calling"""

    #: Workflow name
    name = "somatic_targeted_seq_cnv_calling"

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
                CnvettiOffTargetStepPart,
                CnvettiOnTargetStepPart,
                CnvKitStepPart,
                CopywriterStepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the somatic targeted sequencing CNV calling step"""
        tool_actions = {
            "cnvkit": ["fix", "postprocess", "report", "export"],
            "copywriter": ("call",),
            "cnvetti_on_target": ("coverage", "segment", "postprocess"),
            "cnvetti_off_target": ("coverage", "segment", "postprocess"),
        }
        if "cnvkit" in self.config["tools"] and self.config["cnvkit"]["plot"]:
            tool_actions["cnvkit"] += ["plot"]
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                for tool in self.config["tools"]:
                    for action in tool_actions[tool]:
                        try:
                            tpls = list(self.sub_steps[tool].get_output_files(action).values())
                        except AttributeError:
                            tpls = [self.sub_steps[tool].get_output_files(action)]
                        try:
                            tpls += self.sub_steps[tool].get_log_file(action).values()
                        except AttributeError:
                            tpls += [self.sub_steps[tool].get_log_file(action)]
                        for tpl in tpls:
                            filenames = expand(
                                tpl,
                                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                                library_name=[sample_pair.tumor_sample.dna_ngs_library.name],
                            )
                            for f in filenames:
                                if ".tmp." not in f:
                                    yield f.replace("work/", "output/")

    def check_config(self):
        """Check that the necessary global configuration is present"""
        self.ensure_w_config(
            ("step_config", "somatic_targeted_seq_cnv_calling", "path_ngs_mapping"),
            (
                "Path to somatic variant calling not configured but required for "
                "targeted sequencing CNV calling"
            ),
        )
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            (
                "Path to reference FASTA file not configured but required for targeted sequencing "
                "CNV calling"
            ),
        )
