# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_wgs_sv_calling`` step

The ``somatic_wgs_sv_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned NGS reads) and performs somatic SV calling on them.  The result are called SVs in VCF
format.

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk`` and each read mapper
``mapper`` that the library has been aligned with, and the variant caller ``var_caller``, the
pipeline step will create a directory ``output/{mapper}.{var_caller}.{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.manta.P001-T1-DNA1-WGS1-4
    |   `-- out
    |       |-- bwa.manta.P001-T1-DNA1-WGS1-4.vcf.gz
    |       |-- bwa.manta.P001-T1-DNA1-WGS1-4.vcf.gz.tbi
    |       |-- bwa.manta.P001-T1-DNA1-WGS1-4.vcf.gz.md5
    |       `-- bwa.manta.P001-T1-DNA1-WGS1-4.vcf.gz.tbi.md5
    [...]

Generally, these files will be unfiltered, i.e., contain low-quality variants and also variants
flagged as being non-somatic.

====================
Global Configuration
====================

- The ``static_data_config/reference/path`` has to be configured with the path to the reference
  FASTA file.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_wgs_sv_calling.rst

=============================
Available Somatic CNV Callers
=============================

The following somatic SV callers are currently available

- ``"manta"``
- ``"delly2"``

=======
Reports
=======

Currently, no reports are generated.
"""

import os
import sys
from collections import OrderedDict

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

from .model import SomaticWgsSvCalling as SomaticWgsSvCallingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available somatic variant callers
SOMATIC_VARIANT_CALLERS = ("manta", "delly2")

#: Default configuration for the somatic_wgs_sv_calling schema
DEFAULT_CONFIG = SomaticWgsSvCallingConfigModel.default_config_yaml_string()


class SomaticWgsSvCallingStepPart(BaseStepPart):
    """Base class for somatic WGS SV calling steps

    WGS SV calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    # TODO: unify with somatic (small) variant calling base class?

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{cancer_library}}/out/"
            "{{mapper}}.{var_caller}.{{cancer_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched tumor sample
        self.cancer_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.cancer_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.modules["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            cancer_base_path = (
                "output/{mapper}.{cancer_library}/out/" "{mapper}.{cancer_library}"
            ).format(**wildcards)
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(cancer_base_path + ".bam"),
                "tumor_bai": ngs_mapping(cancer_base_path + ".bam.bai"),
            }

        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.cancer_ngs_library_to_sample_pair[wildcards.cancer_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        # Validate action
        self._validate_action(action)
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, var_caller=[self.name], ext=EXT_VALUES))
        )

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return (
            "work/{{mapper}}.{var_caller}.{{cancer_library}}/log/"
            "snakemake.somatic_wgs_sv_calling.log"
        ).format(var_caller=self.__class__.name)


class MantaStepPart(SomaticWgsSvCallingStepPart):
    """Somatic WGS SV calling with Manta"""

    #: Step name
    name = "manta"

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
            threads=16,
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{int(3.75 * 1024 * 16)}M",
        )


class Delly2StepPart(BaseStepPart):
    """Somatic WGS SV identification using Delly2"""

    name = "delly2"

    #: Actions in Delly 2 workflow
    actions = ("call", "filter_normal", "genotype", "merge_genotypes", "filter_controls")

    #: Directory infixes
    dir_infixes = {
        "call": "{mapper}.delly2.call.{cancer_library}.{sv_type}",
        "filter_normal": "{mapper}.delly2.filter_normal.{cancer_library}.{sv_type}",
        "merge_calls": "{mapper}.delly2.merge_calls.{sv_type}",
        "genotype": "{mapper}.delly2.genotype.{library_name}.{sv_type}",
        "merge_genotypes": "{mapper}.delly2.merge_genotypes.{cancer_library}.{sv_type}",
        "filter_controls": "{mapper}.delly2.filter_controls.{cancer_library}.{sv_type}",
        "final_vcf": "{mapper}.delly2.{cancer_library}.{sv_type}",
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{cancer_library}}/out/"
            "{{mapper}}.{var_caller}.{{cancer_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched tumor sample
        self.cancer_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.cancer_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.cancer_ngs_library_to_sample_pair[wildcards.cancer_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        assert action in self.actions
        mapping = {
            "call": self._get_input_files_call,
            "filter_normal": self._get_input_files_filter_normal,
            "merge_calls": self._get_input_files_merge_calls,
            "genotype": self._get_input_files_genotype,
            "merge_genotypes": self._get_input_files_merge_genotypes,
            "filter_controls": self._get_input_files_filter_controls,
            "final_vcf": self._get_input_files_final_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_call(self, wildcards):
        """Return input files for "call" action: bams for matched T/N pair"""
        ngs_mapping = self.parent.modules["ngs_mapping"]
        normal_tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}{ext}"
        norm_lib = self.get_normal_lib_name(wildcards)
        for name, ext in {"normal_bam": ".bam", "normal_bai": ".bam.bai"}.items():
            yield (
                name,
                ngs_mapping(normal_tpl.format(ext=ext, normal_library=norm_lib, **wildcards)),
            )
        tumor_tpl = "output/{mapper}.{cancer_library}/out/{mapper}.{cancer_library}{ext}"
        for name, ext in {"tumor_bam": ".bam", "tumor_bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tumor_tpl.format(ext=ext, **wildcards))
        # create description of samples that went into this bcf
        infix = self.dir_infixes["call"]
        samples_file_path = os.path.join("work", infix, "out", infix + "samples.tsv").format(
            **wildcards
        )
        with open(samples_file_path, "w") as samples_file:
            samples_file.write(
                "{cancer_library}\ttumor\n{normal_library}\tcontrol\n".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )

    @listify
    def _get_input_files_filter_normal(self, wildcards):
        """Return input files for "filter_normal" action"""
        infix = self.dir_infixes["call"]
        # Samples tsv with matched T/N pair
        yield "tsv", os.path.join("work", infix, "out", infix + ".samples.tsv").format(**wildcards)
        # bcf created with call for that pair
        yield "bcf", os.path.join("work", infix, "out", infix + ".bcf").format(**wildcards)

    @listify
    def _get_input_files_merge_calls(self, wildcards):
        """Return input files for "merge_calls" action"""
        infix = self.dir_infixes["filter_normal"]
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        for pair in self._get_primary_pairs():
            yield tpl.format(library_name=pair.tumor_sample.name, **wildcards)

    @dictify
    def _get_input_files_genotype(self, wildcards):
        """Return input files for "genotype" action"""
        # Sites VCF file
        infix = self.dir_infixes["merge_calls"]
        yield "bcf", os.path.join("work", infix, "out", infix + ".bcf").format(**wildcards)
        # BAM files : we want to individually process each tumor and each normal bam
        ngs_mapping = self.parent.modules["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @listify
    def _get_input_files_merge_genotypes(self, wildcards):
        """Return input files for "merge_genotypes" action"""
        infix = self.dir_infixes["genotype"]
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        # return BCF for one tumor
        yield tpl.format(library_name=wildcards.cancer_library.name, **wildcards)

        # also create tsv with description of all the samples in this bcf
        infix = self.dir_infixes["merge_genotypes"]
        samples_file_path = os.path.join("work", infix, "out", infix + "samples.tsv")
        with open(samples_file_path, "w") as samples_file:
            # write tumor sample description
            samples_file.write("{cancer_library}\ttumor\n".format(**wildcards))

            # for all normals
            for pair in self._get_primary_pairs():
                # write sample description
                samples_file.write(
                    "{library_name}\tcontrol\n".format(
                        library_name=pair.normal_sample.dna_ngs_library.name
                    )
                )
                # return BCF
                yield tpl.format(library_name=pair.normal_sample.dna_ngs_library.name, **wildcards)

    @listify
    def _get_input_files_filter_controls(self, wildcards):
        """Return input files for "filter_controls" action"""
        infix = self.dir_infixes["merge_genotypes"]
        # Samples tsv with one tumor plus all control normals
        yield "tsv", os.path.join("work", infix, "out", infix + "samples.tsv")
        # BCF
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        yield "bcf", tpl.format(**wildcards)

    def _get_input_files_final_vcf(self, wildcards):
        """Return input files for "final_vcf" action"""
        _ = wildcards
        infix = self.dir_infixes["filter_controls"]
        yield os.path.join("work", infix, "out", infix + ".bcf")

    def _get_primary_pairs(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            yield sheet.primary_sample_pairs

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action in self.actions
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            infix = self.dir_infixes[action].replace(r",[^\.]+", "")
            if action != "final_vcf":  # generate bcf files internally
                name = name.replace("vcf", "bcf")
                ext = ext.replace("vcf.gz", "bcf")
                name = name.replace("tbi", "csi")
                ext = ext.replace(".tbi", ".csi")
            yield name, "work/" + infix + "/out/" + infix + ext
        if action in ("call", "merge_genotypes"):
            infix = self.dir_infixes[action].replace(r",[^\.]+", "")
            name = "tsv"
            ext = ".samples.tsv"
            yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action in self.actions
        infix = self.dir_infixes[action].replace(r",[^\.]+", "")
        return "work/" + infix + "/log/snakemake.log"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="6-06:00:00",  # 6 days and 6 hours
            memory=f"{int(3.75 * 1024 * 4)}M",
        )


class SomaticWgsSvCallingWorkflow(BaseStep):
    """Perform somatic variant calling"""

    #: Workflow name
    name = "somatic_wgs_sv_calling"

    #: Default biomedsheet class
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
            config_model_class=SomaticWgsSvCallingConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((Delly2StepPart, MantaStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_module("ngs_mapping", self.config.path_ngs_mapping)

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.{caller}.{cancer_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            caller=self.config.tools,
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )  # pragma: no cover
                    print(
                        msg.format(sample_pair.tumor_sample.name), file=sys.stderr
                    )  # pragma: no cover
                    continue  # pragma: no cover
                yield from expand(
                    tpl, cancer_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA file required by not available",
        )
