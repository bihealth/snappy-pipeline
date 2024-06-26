# -*- coding: utf-8 -*-
"""Implementation of the germline ``variant_phasing`` step

This step takes the result of the ``variant_annotation`` step and performs phasing of the
variants using the GATK tools.  Note that there are some issues with the GATK tools implementing
this step:

- The result of the PhaseByTransmission tool changes the genotype of some variants which is
  problematic when trying to phase *de novo* variants.
- The read backed phasing is also not 100% reliable at the moment.

Thus, the functionality of the tools is made available by this pipeline step but it is not as
fully integrated as it could because it is unclear how useful this is for clinical studies. Also,
so far only the GATK variant caller results can be phased.

Also note that this step generates one output file for each child in a pedigree where both
parents have been sequenced.

==========
Step Input
==========

The variant annotation step uses the output of the following CUBI pipeline steps:

- ``ngs_mapping``
- ``variant_annotation``

===========
Step Output
===========

For each input VCF file (i.e., for each mapper and pedigree), a directory
``output/{mapper}.{caller}.{phaser}.{index_ngs_library}/out`` will be created with the following
output files.

The ``{phaser}`` placeholder can take the values gatk_phase_by_transmission,
gatk_read_backed_phasing, and gatk_phased_both (for the latter, first phasing by transmission
and then read backed phasing is performed).

====================
Global Configuration
====================

- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_phasing.rst

=======
Reports
=======

Currently, no reports are generated.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import os
from collections import OrderedDict

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow

from .model import VariantPhasing as VariantPhasingConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Translate name in configuration to token.
CONFIG_TO_TOKEN = {
    "gatk_read_backed_phasing": "gatk_rbp",
    "gatk_phase_by_transmission": "gatk_pbt",
    "gatk_phasing_both": "gatk_pbt.gatk_rbp",
}

#: Default configuration of the wgs_sv_filtration step
DEFAULT_CONFIG = VariantPhasingConfigModel.default_config_yaml_string()


class WriteTrioPedigreeStepPart(BaseStepPart):
    """Write out trio pedigree file for primary DNA sample given the index NGS library name"""

    #: Step name
    name = "write_trio_pedigree"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to donor
        self.ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if donor.dna_ngs_library:
                        self.ngs_library_to_donor[donor.dna_ngs_library.name] = donor

    @staticmethod
    def get_output_files(action):
        assert action == "run"
        return "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

    def run(self, wildcards, output):
        """Write out the pedigree information"""
        fname = self.get_output_files("run").format(**wildcards)
        donor = self.ngs_library_to_donor[wildcards.index_ngs_library]
        family = "FAM_" + donor.name
        with open(fname, "wt") as ped_file:
            for person in (donor, donor.father, donor.mother):
                if not person:
                    continue
                name = person.dna_ngs_library.name
                father = "0"
                if person.father and person.father.dna_ngs_library:
                    father = person.father.dna_ngs_library.name
                mother = "0"
                if person.mother and person.mother.dna_ngs_library:
                    mother = person.mother.dna_ngs_library.name
                sex = {"male": "1", "female": "2", "unknown": "0"}[
                    person.extra_infos.get("sex", "unknown")
                ]
                affected = {"affected": "2", "unaffected": "1", "unknown": "0"}[
                    person.extra_infos.get("isAffected", "unknown")
                ]
                print("\t".join((family, name, father, mother, sex, affected)), file=ped_file)


class VariantPhasingBaseStep(BaseStepPart):
    """Base step for variant phasing."""

    #: The file name token.
    name_pattern = None

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        # Build mapping from ngs_library to pedigree
        self.ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if donor.dna_ngs_library:
                        self.ngs_library_to_pedigree[donor.dna_ngs_library.name] = pedigree

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out + ext

    @dictify
    def _get_log_file(self, action):
        assert action == "run"
        prefix = (
            "work/{mapper}.{caller}.jannovar_annotate_vcf.%(name)s.{index_library}/log/"
            "{mapper}.{caller}.jannovar_annotate_vcf.%(name)s.{index_library}"
        )

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, (prefix % {"name": self.name_pattern}) + ext


class PhaseByTransmissionStepPart(VariantPhasingBaseStep):
    """Phasing by transmission."""

    #: Name of the step in the pipeline.
    name = "gatk_phase_by_transmission"

    #: The file name token.
    name_pattern = "gatk_pbt"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        name_pattern = r"{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library,[^\.]+}"
        self.base_path_out = os.path.join(
            "work", name_pattern, "out", name_pattern.replace(r",[^\.]+", "")
        )

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            # Pedigree file required for PhaseByTransmission.
            yield (
                "ped",
                "work/write_pedigree.{index_library}/out/{index_library}.ped".format(**wildcards),
            )
            # Get name of real index
            real_index = self.ngs_library_to_pedigree[wildcards.index_library].index
            # Annotated variant file from variant_annotation step.
            variant_annotation = self.parent.sub_workflows["variant_annotation"]
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                input_path = (
                    "output/{mapper}.{caller}.jannovar_annotate_vcf.{real_index}/out/"
                    "{mapper}.{caller}.jannovar_annotate_vcf.{real_index}"
                ).format(real_index=real_index.dna_ngs_library.name, **wildcards)
                yield key, variant_annotation(input_path) + ext

        assert action == "run", "Unsupported actions"
        return input_function

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
            time="1-00:00:00",  # 1 day
            memory=f"{14 * 1024}M",
        )


class ReadBackedPhasingBaseStep(VariantPhasingBaseStep):
    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from library name to pedigree
        self.ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    self.ngs_library_to_donor[donor.dna_ngs_library.name] = donor

    def _yield_bams(self, wildcards):
        """Helper function used in subclass input_function"""
        donor = self.ngs_library_to_donor[wildcards.index_library]
        tpl = "output/{mapper}.{index_library}/out/{mapper}.{index_library}{ext}"
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            vals = {"mapper": wildcards.mapper, "ext": ext}
            # Note that we only perform phasing for pedigree members we have both parents, so
            # the following works.
            if (
                donor.dna_ngs_library
                and donor.father
                and donor.father.dna_ngs_library
                and donor.mother
                and donor.mother.dna_ngs_library
            ):
                files = [
                    tpl.format(index_library=donor.dna_ngs_library.name, **vals),
                    tpl.format(index_library=donor.father.dna_ngs_library.name, **vals),
                    tpl.format(index_library=donor.mother.dna_ngs_library.name, **vals),
                ]
                yield key, list(map(ngs_mapping, files))

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
        mem_mb = 8 * 1024
        return ResourceUsage(
            threads=1,
            time="1-00:00:00",  # 1 day
            memory=f"{mem_mb}M",
        )


class ReadBackedPhasingOnlyStepPart(ReadBackedPhasingBaseStep):
    """Read backed phasing (as primary and only step)"""

    #: Name of the step in the pipeline.
    name = "gatk_read_backed_phasing_only"
    #: The file name token.
    name_pattern = "gatk_rbp"

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = "{mapper}.{caller}.jannovar_annotate_vcf.gatk_rbp.{index_library}"
        self.base_path_out = os.path.join("work", name_pattern, "out", name_pattern)

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            # Get name of real index
            real_index = self.ngs_library_to_pedigree[wildcards.index_library].index
            # BAM files from ngs_mapping step.
            yield from self._yield_bams(wildcards)
            # Annotated variant file from variant_annotation step.
            variant_annotation = self.parent.sub_workflows["variant_annotation"]
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                output_path = (
                    "output/{mapper}.{caller}.jannovar_annotate_vcf.{real_index}/out/"
                    "{mapper}.{caller}.jannovar_annotate_vcf.{real_index}"
                ).format(real_index=real_index.dna_ngs_library.name, **wildcards)
                yield key, variant_annotation(output_path) + ext

        assert action == "run", "Unsupported actions"
        return input_function


class ReadBackedPhasingAlsoStepPart(ReadBackedPhasingBaseStep):
    """Read backed phasing (as step after phase by transmission)"""

    #: Name of the step in the pipeline.
    name = "gatk_read_backed_phasing_also"
    #: The file name token.
    name_pattern = "gatk_pbt.gatk_rbp"

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = "{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.gatk_rbp.{index_library}"
        self.base_path_out = os.path.join("work", name_pattern, "out", name_pattern)

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            # BAM files from ngs_mapping step.
            yield from self._yield_bams(wildcards)
            # Result of PhaseByTransmission step
            name_pattern = "{mapper}.{caller}.jannovar_annotate_vcf.gatk_pbt.{index_library}"
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                input_path = "work/" + name_pattern + "/out/" + name_pattern
                yield key, input_path.format(**wildcards) + ext

        assert action == "run", "Unsupported actions"
        return input_function


class VariantPhasingWorkflow(BaseStep):
    """Perform (small) variant phasing"""

    name = "variant_phasing"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=VariantPhasingConfigModel,
            previous_steps=(VariantAnnotationWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WriteTrioPedigreeStepPart,
                PhaseByTransmissionStepPart,
                ReadBackedPhasingOnlyStepPart,
                ReadBackedPhasingAlsoStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("variant_annotation", self.config.path_variant_annotation)
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.tools_variant_calling:
            self.config.tools_variant_calling = self.w_config.step_config["variant_calling"].tools

    @listify
    def get_result_files(self):
        """Return list of result files for the variant filtration workflow."""
        # Generate output paths without extracting individuals.
        name_pattern = "{mapper}.{caller}.jannovar_annotate_vcf.{phasing}.{index_library.name}"
        phasings = [
            token for name, token in CONFIG_TO_TOKEN.items() if name in self.config.phasings
        ]
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            phasing=phasings,
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if (
                        donor.dna_ngs_library
                        and donor.father
                        and donor.father.dna_ngs_library
                        and donor.mother
                        and donor.mother.dna_ngs_library
                    ):  # only phase if both parents present
                        yield from expand(tpl, index_library=[donor.dna_ngs_library], **kwargs)
