# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_cnv_checking`` step

The ``somatic_cnv_checking`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs pileups at germline locii to identify heterozygous
variants, and their B-allele frequency in the paired tumor sample. If one somatic CNV calling step
(either ``somatic_wgs_cnv_calling`` or ``somatic_targetd_seq_cnv_calling``) is present, the
log2 convarage ratio and copy number call are added to the output vcf.

==========
Step Input
==========

The somatic CNV checking step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step. Optionally, it can also use one Snakemake sub workflows for somatic CNV
calling, either``somatic_wgs_cnv_calling`` or ``somatic_targetd_seq_cnv_calling``.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk`` and each read mapper
``mapper`` that the library has been aligned with, and the CNV caller ``caller``, the
pipeline step will create a directory ``output/{mapper}.{caller}.{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.cnvkit.P001-T1-DNA1-WES1-4
    |   `-- out
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz.tbi
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz.md5
    |       `-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz.tbi.md5
    [...]

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_cnv_checking.rst

=======
Reports
=======

Currently, no reports are generated.
"""

from collections import OrderedDict
import os
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

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Default configuration for the somatic_cnv_checking schema
DEFAULT_CONFIG = r"""
# Default configuration somatic_cnv_checking
step_config:
  somatic_cnv_checking:
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    path_cnv_calling: ""              # Can use for instance ../somatic_targeted_seq_cnv_calling
    cnv_calling_tool: ""              # Can use "cnvkit"
    excluded_regions: ""              # Bed file of regions to be excluded
    max_depth: 10000                  # Max depth for pileups
    min_cov: 20                       # Minimum depth for reference and alternative alleles to consider variant
    min_baf: 0.4                      # Maximum BAF to consider variant as heterozygous (between 0 & 1/2)
"""


class SomaticCnvCheckingStepPart(BaseStepPart):
    """Base class for CNV checking sub-steps"""
        
    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Get shorcut to Snakemake sub workflow
        self.ngs_mapping = self.parent.sub_workflows["ngs_mapping"]

    def get_params(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            pair = self.tumor_ngs_library_to_sample_pair[wildcards.library]
            return {"library_name": pair.tumor_sample.dna_ngs_library.name}

        return input_function

    @dictify
    def _get_log_file(self, prefix):
        """Return dict of log files."""
        # Validate action
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class SomaticCnvCheckingPileupStepPart(SomaticCnvCheckingStepPart):
    """Perform pileups at a set of locations"""

    name = "pileup"
    actions = ("normal", "tumor")

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        def input_function_normal(wildcards):
            pair = self.tumor_ngs_library_to_sample_pair[wildcards.library]
            base_path = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                normal_library=pair.normal_sample.dna_ngs_library.name, **wildcards
            )
            return {
                "bam": self.ngs_mapping(base_path + ".bam"),
                "bai": self.ngs_mapping(base_path + ".bam.bai"),
            }

        def input_function_tumor(wildcards):
            pair = self.tumor_ngs_library_to_sample_pair[wildcards.library]
            base_path = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}".format(
                tumor_library=pair.tumor_sample.dna_ngs_library.name, **wildcards
            )
            return {
                "vcf": "work/{mapper}.{library}/out/{mapper}.{library}.normal.vcf.gz".format(
                    **wildcards
                ),
                "tbi": "work/{mapper}.{library}/out/{mapper}.{library}.normal.vcf.gz.tbi".format(
                    **wildcards
                ),
                "bam": self.ngs_mapping(base_path + ".bam"),
                "bai": self.ngs_mapping(base_path + ".bam.bai"),
            }

        if action == "normal":
            return input_function_normal
        else:
            return input_function_tumor

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        # Validate action
        self._validate_action(action)
        base_path_out = "work/{{mapper}}.{{library}}/out/{{mapper}}.{{library}}.{action}{ext}"
        return dict(zip(EXT_NAMES, expand(base_path_out, action=action, ext=EXT_VALUES)))

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self._get_log_file("work/{mapper}.{library}/log/{mapper}.{library}." + action)

    def get_resource_usage(self, action):
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="01:00:00" if action == "tumor" else "12:00:00",  # 1 hour
            memory=f"{int(3.7 * 1024 * 2)}M",
        )


class SomaticCnvCheckingCnvStepPart(SomaticCnvCheckingStepPart):
    """Merging heterozygous variants with CNV data"""

    name = "cnv"
    actions = ("run",)

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        filenames = {}
        filenames["normal"] = "work/{mapper}.{library}/out/{mapper}.{library}.normal.vcf.gz"
        filenames["normal_tbi"] = filenames["normal"] + ".tbi"
        filenames["tumor"] = "work/{mapper}.{library}/out/{mapper}.{library}.tumor.vcf.gz"
        filenames["tumor_tbi"] = filenames["tumor"] + ".tbi"
        if self.config["path_cnv_calling"]:
            tpl = "{{mapper}}.{caller}.{{library}}".format(caller=self.config["cnv_calling_tool"])
            filenames["cnv"] = os.path.join(
                self.config["path_cnv_calling"], "output", tpl, "out", tpl + ".bed.gz"
            )
            filenames["cnv_tbi"] = filenames["cnv"] + ".tbi"
        return filenames

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        if self.config["path_cnv_calling"] and self.config["cnv_calling_tool"]:
            name_pattern = "{{mapper}}.{caller}.{{library}}".format(
                caller=self.config["cnv_calling_tool"]
            )
        else:
            name_pattern = "{mapper}.{library}"
        base_path_out = "work/" + name_pattern + "/out/" + name_pattern
        return dict(zip(EXT_NAMES, [base_path_out + ext for ext in EXT_VALUES]))

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        if self.config["path_cnv_calling"] and self.config["cnv_calling_tool"]:
            name_pattern = "{{mapper}}.{caller}.{{library}}".format(
                caller=self.config["cnv_calling_tool"]
            )
        else:
            name_pattern = "{mapper}.{library}"
        return self._get_log_file(os.path.join("work", name_pattern, "log", name_pattern))


class SomaticCnvCheckingReportStepPart(SomaticCnvCheckingStepPart):
    """Plots of BAF vs CNV"""

    name = "report"
    actions = ("run",)

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{{mapper}}.{caller}.{{library}}".format(
            caller=self.config["cnv_calling_tool"]
        )
        base_path_out = "work/" + name_pattern + "/out/" + name_pattern
        return {"vcf": base_path_out + ".vcf.gz"}

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{{mapper}}.{caller}.{{library}}".format(
            caller=self.config["cnv_calling_tool"]
        )
        base_path_out = "work/" + name_pattern + "/report/" + name_pattern
        return {
            "cnv": base_path_out + ".cnv.pdf",
            "cnv_md5": base_path_out + ".cnv.pdf.md5",
            "locus": base_path_out + ".locus.pdf",
            "locus_md5": base_path_out + ".locus.pdf.md5",
        }

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{{mapper}}.{caller}.{{library}}".format(
            caller=self.config["cnv_calling_tool"]
        )
        return self._get_log_file(os.path.join("work", name_pattern, "log", name_pattern + ".report"))


class SomaticCnvCheckingWorkflow(BaseStep):
    """Perform somatic cnv checking"""

    #: Workflow name
    name = "somatic_cnv_checking"

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
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                SomaticCnvCheckingPileupStepPart,
                SomaticCnvCheckingCnvStepPart,
                SomaticCnvCheckingReportStepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        if self.config["path_cnv_calling"] and self.config["cnv_calling_tool"]:
            name_pattern = "{{mapper}}.{caller}.{{tumor_library.name}}".format(
                caller=self.config["cnv_calling_tool"]
            )
        else:
            name_pattern = "{mapper}.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            ext=[
                part + "." + ext + chksum
                for part in ("", ".normal", ".tumor")
                for ext in ("log", "conda_info.txt", "conda_list.txt")
                for chksum in ("", ".md5")
            ],
        )
        if self.config["path_cnv_calling"] and self.config["cnv_calling_tool"]:
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "report", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                ext=(".cnv.pdf", ".cnv.pdf.md5", ".locus.pdf", ".locus.pdf.md5"),
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + ".report" + "{ext}"),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                ext=[
                    "." + ext + chksum
                    for ext in ("log", "conda_info.txt", "conda_list.txt")
                    for chksum in ("", ".md5")
                ],
            )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
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
                yield from expand(
                    tpl, tumor_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "somatic_cnv_checking", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for somatic variant calling",
        )
