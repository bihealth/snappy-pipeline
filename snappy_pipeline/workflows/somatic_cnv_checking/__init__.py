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

import os
import sys
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand, Wildcards

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    SomaticTargetedSeqCnvCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_wgs_cnv_calling import SomaticWgsCnvCallingWorkflow

from .model import SomaticCnvChecking as SomaticCnvCheckingConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Default configuration for the somatic_cnv_checking schema
DEFAULT_CONFIG = SomaticCnvCheckingConfigModel.default_config_yaml_string()


class SomaticCnvCheckingStepPart(BaseStepPart):
    """Base class for CNV checking sub-steps"""

    @dictify
    def _get_log_file(self, prefix):
        """Return dict of log files."""
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

    def __init__(self, parent):
        super().__init__(parent)
        self.ngs_mapping = self.parent.sub_workflows["ngs_mapping"]

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        def input_function_normal(wildcards):
            base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(
                **wildcards
            )
            return {
                "bam": self.ngs_mapping(base_path + ".bam"),
                "bai": self.ngs_mapping(base_path + ".bam.bai"),
            }

        def input_function_tumor(wildcards):
            base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(
                **wildcards
            )
            return {
                "locii": "work/{mapper}.{normal_library}/out/{mapper}.{normal_library}.normal.vcf.gz".format(
                    normal_library=self.parent.tumor_to_normal[wildcards.library_name], **wildcards
                ),
                "locii_tbi": "work/{mapper}.{normal_library}/out/{mapper}.{normal_library}.normal.vcf.gz.tbi".format(
                    normal_library=self.parent.tumor_to_normal[wildcards.library_name], **wildcards
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
        base_path_out = (
            "work/{{mapper}}.{{library_name}}/out/{{mapper}}.{{library_name}}.{action}{ext}"
        )
        return dict(zip(EXT_NAMES, expand(base_path_out, action=action, ext=EXT_VALUES)))

    def get_args(self, **kwargs):
        def args_fn(_wildcards):
            return {
                "reference_path": self.w_config.static_data_config.reference.path,
                "min_baf": self.config.min_baf,
                "min_depth": self.config.min_depth,
                "max_depth": self.config.max_depth,
            }

        return args_fn

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self._get_log_file(
            "work/{mapper}.{library_name}/log/{mapper}.{library_name}." + action
        )

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

        def input_function(wildcards):
            normal_library = self.parent.tumor_to_normal[wildcards.library_name]
            filenames = {}
            name_pattern = "{mapper}.{normal_library}"
            tpl = os.path.join("work", name_pattern, "out", name_pattern + ".normal.vcf.gz")
            filenames["normal"] = tpl.format(normal_library=normal_library, **wildcards)
            filenames["normal_tbi"] = filenames["normal"] + ".tbi"
            name_pattern = "{mapper}.{library_name}"
            tpl = os.path.join("work", name_pattern, "out", name_pattern + ".tumor.vcf.gz")
            filenames["tumor"] = tpl.format(**wildcards)
            filenames["tumor_tbi"] = filenames["tumor"] + ".tbi"
            cnv_calling = self.parent.sub_workflows["cnv_calling"]
            base_path = "output/{mapper}.{caller}.{library_name}/out/{mapper}.{caller}.{library_name}".format(
                **wildcards
            )
            filenames["cnv"] = cnv_calling(base_path + "_dnacopy.seg")
            return filenames

        return input_function

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        name_pattern = "{mapper}.{caller}.{library_name}"
        key_ext["tsv"] = ".tsv"
        base_path_out = "work/" + name_pattern + "/out/" + name_pattern
        for key, ext in key_ext.items():
            yield (key, base_path_out + ext)
            yield (key + "_md5", base_path_out + ext + ".md5")

    def get_args(self, action: str) -> dict[str, Any]:
        # Validate action
        self._validate_action(action)
        return self.config.model_dump(by_alias=True)

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{library_name}"
        return self._get_log_file(os.path.join("work", name_pattern, "log", name_pattern))


class SomaticCnvCheckingReportStepPart(SomaticCnvCheckingStepPart):
    """Plots of BAF vs CNV"""

    name = "report"
    actions = ("run",)

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            name_pattern = "{mapper}.{caller}.{library_name}".format(**wildcards)
            base_path_out = "work/" + name_pattern + "/out/" + name_pattern
            return {"vcf": base_path_out + ".vcf.gz", "tsv": base_path_out + ".tsv"}

        return input_function

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{library_name}"
        base_path_out = "work/" + name_pattern + "/report/" + name_pattern
        return {
            "cnv": base_path_out + ".cnv.pdf",
            "cnv_md5": base_path_out + ".cnv.pdf.md5",
            "locus": base_path_out + ".locus.pdf",
            "locus_md5": base_path_out + ".locus.pdf.md5",
            "segment": base_path_out + ".segment.pdf",
            "segment_md5": base_path_out + ".segment.pdf.md5",
        }

    def get_args(self, action: str):
        # Validate action
        self._validate_action(action)

        def args_fn(wildcards: Wildcards) -> dict[str, Any]:
            return {
                "reference": self.parent.w_config.static_data_config.reference.path,
                "library_name": wildcards.library_name,
            }

        return args_fn

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{library_name}"
        return self._get_log_file(
            os.path.join("work", name_pattern, "log", name_pattern + ".report")
        )


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
            config_model_class=SomaticCnvCheckingConfigModel,
            previous_steps=(
                SomaticTargetedSeqCnvCallingWorkflow,
                SomaticWgsCnvCallingWorkflow,
                NgsMappingWorkflow,
            ),
        )
        if self.config.path_cnv_calling and self.config.cnv_assay_type:
            if self.config.cnv_assay_type == "WES":
                cnv_calling = "somatic_targeted_seq_cnv_calling"
            elif self.config.cnv_assay_type == "WES":
                cnv_calling = "somatic_wgs_cnv_calling"
            else:
                raise InvalidConfiguration(
                    "Illegal cnv_assay_type {}, must be either WES or WGS".format(
                        self.config.cnv_assay_type
                    )
                )
            self.register_sub_workflow(cnv_calling, self.config.path_cnv_calling, "cnv_calling")
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                SomaticCnvCheckingPileupStepPart,
                SomaticCnvCheckingCnvStepPart,
                SomaticCnvCheckingReportStepPart,
                LinkOutStepPart,
            )
        )
        # Assemble normal/tumor pairs
        self.tumor_to_normal = {}
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
                tumor_library = sample_pair.tumor_sample.dna_ngs_library.name
                normal_library = sample_pair.normal_sample.dna_ngs_library.name
                self.tumor_to_normal[tumor_library] = normal_library

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        # Log files from normal pileups
        name_pattern = "{mapper}.{library_name}"
        chksum = ("", ".md5")
        ext = ("log", "conda_info.txt", "conda_list.txt")
        yield from expand(
            os.path.join("output", name_pattern, "log", name_pattern + ".normal.{ext}{chksum}"),
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            library_name=set(self.tumor_to_normal.values()),
            ext=ext,
            chksum=chksum,
        )
        yield from expand(
            os.path.join("output", name_pattern, "log", name_pattern + ".tumor.{ext}{chksum}"),
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            library_name=self.tumor_to_normal.keys(),
            ext=ext,
            chksum=chksum,
        )
        # Main result: vcf & optionally segment table if CNV available
        ext = {"out": [".vcf.gz", ".vcf.gz.tbi"]}
        if self.config.path_cnv_calling:
            # CNV avaliable
            name_pattern = "{mapper}.{caller}.{library_name}"
            callers = self.w_config.step_config["somatic_targeted_seq_cnv_calling"].tools
            ext["out"] += [".tsv"]
            ext["report"] = (".cnv.pdf", ".locus.pdf", ".segment.pdf")
            ext["log"] = [
                suffix + "." + e
                for suffix in ("", ".report")
                for e in ("log", "conda_info.txt", "conda_list.txt")
            ]
        else:
            name_pattern = "{mapper}.{library_name}"
            callers = []
        for subdir, exts in ext.items():
            yield from expand(
                os.path.join("output", name_pattern, subdir, name_pattern + "{ext}{chksum}"),
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                caller=callers,
                library_name=self.tumor_to_normal.keys(),
                ext=exts,
                chksum=chksum,
            )
