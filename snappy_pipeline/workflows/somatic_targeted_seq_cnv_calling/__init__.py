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

- ``{mapper}.cnvkit.export.{lib_name}-{lib_pk}/out/``
    - ``{mapper}.cnvkit.export.{lib_name}-{lib_pk}.bed``
    - ``{mapper}.cnvkit.export.{lib_name}-{lib_pk}.seg``
    - ``{mapper}.cnvkit.export.{lib_name}-{lib_pk}.vcf.gz``
    - ``{mapper}.cnvkit.export.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}/out``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.diagram.pdf``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.scatter.pdf``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.heatmap.pdf``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.heatmap.chr1.pdf``
    - ...
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.scatter.chrX.pdf``
- ``{mapper}.cnvkit.report.{lib_name}-{lib_pk}/out``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.breaks.txt``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.gainloss.txt``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.gender.txt``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.metrics.txt``
    - ``{mapper}.cnvkit.plot.{lib_name}-{lib_pk}.segmetrics.txt``

For example:

::

    output/
    |-- bwa.cnvkit.export.P001-T1-DNA1-WES1-000007
    |   `-- out
    |       |-- bwa.cnvkit.export.P001-T1-DNA1-WES1-000007.bed
    |       |-- bwa.cnvkit.export.P001-T1-DNA1-WES1-000007.seg
    |       `-- bwa.cnvkit.export.P001-T1-DNA1-WES1-000007.vcf
    |-- bwa.cnvkit.plot.P002-T1-DNA1-WES1-000016
    |   `-- out
    |       |-- bwa.cnvkit.plot.P002-T1-DNA1-WES1-000016.diagram.pdf
    |       |-- bwa.cnvkit.plot.P002-T1-DNA1-WES1-000016.heatmap.pdf
    |       |-- bwa.cnvkit.plot.P002-T1-DNA1-WES1-000016.scatter.pdf
    |       |-- bwa.cnvkit.plot.P002-T1-DNA1-WES1-000016.heatmap.chr1.pdf
    |       |-- ...
    |       `-- bwa.cnvkit.plot.P002-T1-DNA1-WES1-000016.scatter.chrX.pdf
    |-- bwa.cnvkit.report.P002-T1-DNA1-WES1-000016
    |   `-- out
    |       |-- bwa.cnvkit.report.P002-T1-DNA1-WES1-000016.breaks.txt
    |       |-- bwa.cnvkit.report.P002-T1-DNA1-WES1-000016.gainloss.txt
    |       |-- bwa.cnvkit.report.P002-T1-DNA1-WES1-000016.gender.txt
    |       |-- bwa.cnvkit.report.P002-T1-DNA1-WES1-000016.metrics.txt
    |       `-- bwa.cnvkit.report.P002-T1-DNA1-WES1-000016.segmetrics.txt
    [...]

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
import itertools
import os
import os.path
import sys

from biomedsheets.shortcuts import CancerCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Default configuration for the somatic_targeted_seq_cnv_calling step
DEFAULT_CONFIG = r"""
# Default configuration somatic_targeted_seq_cnv_calling
step_config:
  somatic_targeted_seq_cnv_calling:
    tools: ['cnvkit']
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    cnvkit:
      path_target_regions: REQUIRED   # REQUIRED
      seg_method: haar
      seg_threshold: 0.0001
      # BCBIO uses
      # seg_method: haar
      # seg_threshold: 0.0001
      # -- OR
      # seg_method: cbs
      # seg_threshold: 0.000001
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
    ("csi", ".bcf.csi"),
    ("csi_md5", ".bcf.csi.md5"),
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
        assert action in self.actions, "Invalid action"
        return getattr(self, "_get_input_files_{action}".format(action=action))()

    def _get_input_files_coverage(self):
        @dictify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
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
        assert action in self.actions, "Invalid action"
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
        assert action in self.actions, "Invalid action"
        name_pattern = self.name_pattern.format(action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, os.path.join("work", name_pattern, "log", name_pattern + ext)

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration with resource usage limits for
        scheduling
        """
        for action in self.actions:
            key = "somatic_targeted_seq_cnv_calling_{tool}_{action}".format(
                tool=self.name, action=action
            )
            cluster_config[key] = {"mem": 7500, "time": "24:00", "ntasks": 1}


class CnvettiOffTargetStepPart(CnvettiStepPartBase):
    """Perform somatic targeted CNV calling using CNVetti with off-target reads."""

    name = "cnvetti_off_target"


class CnvettiOnTargetStepPart(CnvettiStepPartBase):
    """Perform somatic targeted CNV calling using CNVetti with on-target reads."""

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

    name = "cnvkit"

    def __init__(self, parent):
        super().__init__(parent)

    def check_config(self):
        """Check configuration for cnvkit"""
        if "cnvkit" not in self.config["tools"]:
            return  # cnvkit not enabled, skip
        self.parent.ensure_w_config(
            ("step_config", "somatic_targeted_seq_cnv_calling", "cnvkit", "path_target_regions"),
            "Path to target regions is missing for cnvkit",
        )

    def get_input_files(self, action):
        """Return input paths input function, dependent on rule"""
        method_mapping = {
            "access": None,
            "target": self._get_input_files_target,
            "antitarget": self._get_input_files_antitarget,
            "coverage": self._get_input_files_coverage,
            "reference": self._get_input_files_reference,
            "call": self._get_input_files_call,
            "fix": self._get_input_files_fix,
            "segment": self._get_input_files_segment,
            "export": self._get_input_files_export,
            "plot": self._get_input_files_plot,
            "report": self._get_input_files_report,
        }
        assert action in method_mapping, "Unknown action"
        return method_mapping[action]

    def _get_input_files_target(self, _):
        input_files = {"access": "work/cnvkit.access/out/access.bed"}
        return input_files

    def _get_input_files_antitarget(self, _):
        input_files = {
            "access": "work/cnvkit.access/out/access.bed",
            "target": "work/cnvkit.target/out/target.bed",
        }
        return input_files

    def _get_input_files_coverage(self, wildcards):
        input_files = {
            "target": "work/cnvkit.target/out/target.bed",
            "antitarget": "work/cnvkit.antitarget/out/antitarget.bed",
        }
        # BAM/BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files["bam"] = ngs_mapping(base_path + ".bam")
        input_files["bai"] = ngs_mapping(base_path + ".bam.bai")
        return input_files

    def _get_input_files_reference(self, wildcards):
        input_files = {
            "target": "work/cnvkit.target/out/target.bed",
            "antitarget": "work/cnvkit.antitarget/out/antitarget.bed",
        }
        # BAM/BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
            normal_library=self.get_normal_lib_name(wildcards), **wildcards
        )
        input_files["bam"] = ngs_mapping(base_path + ".bam")
        input_files["bai"] = ngs_mapping(base_path + ".bam.bai")
        return input_files

    def _get_input_files_fix(self, wildcards):
        tpl_base = "{mapper}.cnvkit.{substep}.{library_name}"
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".cnn"
        input_files = {"ref": tpl.format(substep="reference", **wildcards)}
        tpl = "work/" + tpl_base + "/out/" + tpl_base + ".{target}coverage.cnn"
        for target in ("target", "antitarget"):
            input_files[target] = tpl.format(target=target, substep="coverage", **wildcards)
        return input_files

    def _get_input_files_segment(self, wildcards):
        cnr_pattern = (
            "work/{mapper}.cnvkit.fix.{library_name}/out/{mapper}.cnvkit.fix.{library_name}.cnr"
        )
        input_files = {"cnr": cnr_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_call(self, wildcards):
        segment_pattern = (
            "work/{mapper}.cnvkit.segment.{library_name}/out/"
            "{mapper}.cnvkit.segment.{library_name}.cns"
        )
        input_files = {"segment": segment_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_export(self, wildcards):
        cns_pattern = (
            "work/{mapper}.cnvkit.call.{library_name}/out/{mapper}.cnvkit.call.{library_name}.cns"
        )
        input_files = {"cns": cns_pattern.format(**wildcards)}
        return input_files

    def _get_input_files_plot(self, wildcards):
        tpl = (
            "work/{mapper}.cnvkit.{substep}.{library_name}/out/"
            "{mapper}.cnvkit.{substep}.{library_name}.{ext}"
        )
        input_files = {
            "cnr": tpl.format(substep="fix", ext="cnr", **wildcards),
            "cns": tpl.format(substep="segment", ext="cns", **wildcards),
        }
        return input_files

    def _get_input_files_report(self, wildcards):
        return self._get_input_files_plot(wildcards)

    def get_output_files(self, action):
        """Return output files for the given action"""
        method_mapping = {
            "access": self._get_output_files_access,
            "target": self._get_output_files_target,
            "antitarget": self._get_output_files_antitarget,
            "coverage": self._get_output_files_coverage,
            "reference": self._get_output_files_reference,
            "fix": self._get_output_files_fix,
            "call": self._get_output_files_call,
            "segment": self._get_output_files_segment,
            "export": self._get_output_files_export,
            "plot": self._get_output_files_plot,
            "report": self._get_output_files_report,
        }
        assert action in method_mapping, "Unknown action"
        return method_mapping[action]()

    def _get_output_files_access(self):
        return "work/cnvkit.access/out/access.bed"

    def _get_output_files_target(self):
        return "work/cnvkit.target/out/target.bed"

    def _get_output_files_antitarget(self):
        return "work/cnvkit.antitarget/out/antitarget.bed"

    def _get_output_files_coverage(self):
        name_pattern = "{mapper}.cnvkit.coverage.{library_name}"
        output_files = {}
        for target in ("target", "antitarget"):
            output_files[target] = os.path.join(
                "work", name_pattern, "out", name_pattern + ".{}coverage.cnn".format(target)
            )
        return output_files

    def _get_output_files_reference(self):
        name_pattern = "{mapper}.cnvkit.reference.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cnn")
        return tpl

    def _get_output_files_fix(self):
        name_pattern = "{mapper}.cnvkit.fix.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cnr")
        return tpl

    def _get_output_files_segment(self):
        name_pattern = "{mapper}.cnvkit.segment.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cns")
        return tpl

    def _get_output_files_call(self):
        name_pattern = "{mapper}.cnvkit.call.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".cns")
        return tpl

    @dictify
    def _get_output_files_plot(self):
        plots = ("scatter", "diagram", "heatmap")
        chroms = list(itertools.chain(range(1, 23), ["X", "Y"]))
        # Yield file name pairs for global plots
        tpl = (
            "work/{mapper}.cnvkit.plot.{library_name}/out/"
            "{mapper}.cnvkit.plot.{library_name}.{plot}.pdf"
        )
        yield from (
            (plot, tpl.format(plot=plot, **format_id("mapper", "library_name"))) for plot in plots
        )
        # Yield file name pairs for the chromosome-wise plots
        chrom_plots = ("scatter", "heatmap")
        tpl_chrom = (
            "work/{mapper}.cnvkit.plot.{library_name}/out/"
            "{mapper}.cnvkit.plot.{library_name}.{plot}.chr{chrom}.pdf"
        )
        yield from (
            (
                "{plot}_chr{chrom}".format(plot=plot, chrom=chrom),
                tpl_chrom.format(plot=plot, chrom=chrom, **format_id("mapper", "library_name")),
            )
            for plot in chrom_plots
            for chrom in chroms
        )

    def _get_output_files_export(self):
        keys = ("bed", "seg", "vcf", "tbi")
        exts = ("bed", "seg", "vcf.gz", "vcf.gz.tbi")
        name_pattern = "{mapper}.cnvkit.export.{library_name}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".{ext}")
        output_files = {}
        for key, ext in zip(keys, exts):
            output_files[key] = tpl.format(ext=ext, **format_id("mapper", "library_name"))
        return output_files

    @dictify
    def _get_output_files_report(self):
        reports = ("breaks", "gainloss", "gender", "metrics", "segmetrics")
        tpl = (
            "work/{mapper}.cnvkit.report.{library_name}/out/"
            "{mapper}.cnvkit.report.{library_name}.{report}.txt"
        )
        yield from (
            (report, tpl.format(report=report, **format_id("mapper", "library_name")))
            for report in reports
        )

    def get_log_file(self, action):
        """Return path to log file for the given action"""
        prefix = None
        if action in ("access", "target", "antitarget"):
            prefix = "work/cnvkit.{action}/log/cnvkit.{action}"
        elif action in (
            "coverage",
            "reference",
            "fix",
            "call",
            "segment",
            "export",
            "plot",
            "report",
        ):
            prefix = (
                "work/{{mapper}}.cnvkit.{action}.{{library_name}}/log/"
                "{{mapper}}.cnvkit.{action}.{{library_name}}"
            )
        else:
            raise ValueError("Unknown action {}".format(action))
        prefix = prefix.format(action=action)
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

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration with resource usage limits for
        scheduling
        """
        actions = (
            "access",
            "target",
            "antitarget",
            "coverage",
            "reference",
            "fix",
            "call",
            "segment",
            "export",
            "plot",
            "report",
        )
        for action in actions:
            key = "somatic_targeted_seq_cnv_calling_cnvkit_{}".format(action)
            if action == "plot":
                memory = 30 * 1024
            else:
                memory = int(7.5 * 1024)
            cluster_config[key] = {"mem": memory, "time": "24:00", "ntasks": 1}


class CopywriterStepPart(SomaticTargetedSeqCnvCallingStepPart):
    """Perform somatic targeted CNV calling using CopywriteR"""

    name = "copywriter"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.copywriter.{{library_name}}/out/"
            "{{mapper}}.copywriter.{{library_name}}{ext}"
        )

    def get_input_files(self, action):
        def input_function_run(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
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

        assert action in ["run", "call", "Unsupported actions"]
        if action == "run":
            return input_function_run
        if action == "call":
            return input_function_call

    @dictify
    def get_output_files(self, action):
        assert action in ["run", "call", "Unsupported actions"]
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
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )

        if action in ("prepare"):
            tpl = "work/copywriter.{action}/log/snakemake.log"
            return tpl.format(action=action)
        elif action in ("call", "run"):
            tpl = (
                "work/{mapper}.copywriter.{library_name}/log/{mapper}.copywriter.{library_name}."
                + action
            )
            for key, ext in key_ext:
                yield key, tpl + ext
        else:
            raise ValueError("Unknown action {}".format(action))

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration with resource usage limits for
        scheduling
        """
        tpl = "somatic_targeted_seq_cnv_calling_copywriter_{}"
        cluster_config[tpl.format("prepare")] = {"mem": int(4000), "time": "2:00", "ntasks": 1}
        cluster_config[tpl.format("run")] = {"mem": int(80000), "time": "16:00", "ntasks": 2}
        cluster_config[tpl.format("call")] = {"mem": int(8000), "time": "3:59:00", "ntasks": 8}


class SomaticTargetedSeqCnvCallingWorkflow(BaseStep):
    """Perform somatic targeted sequencing CNV calling"""

    name = "somatic_targeted_seq_cnv_calling"
    sheet_shortcut_class = CancerCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow,
            config,
            cluster_config,
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
            "cnvkit": ("call", "report", "export", "plot"),  # ("report", "export", "plot"),
            "copywriter": ("call",),
            "cnvetti_on_target": ("coverage", "segment", "postprocess"),
            "cnvetti_off_target": ("coverage", "segment", "postprocess"),
        }
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
        """Check that the necessary globalc onfiguration is present"""
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
