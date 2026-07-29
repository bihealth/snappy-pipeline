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
``output/{index_ngs_library}/out`` will be created with the following
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
from typing import Any

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
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
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

    def get_input_files(self, action: str):
        self._validate_action(action)
        return []

    @staticmethod
    def get_output_files(action):
        assert action == "run"
        return "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

    def run(self, wildcards, output):
        """Write out the pedigree information"""
        fname = self.get_output_files("run").format(**wildcards)
        os.makedirs(os.path.dirname(fname), exist_ok=True)
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
        prefix = f"work/{{index_library}}/log/{self.name_pattern}.{{index_library}}"

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext


class PhaseByTransmissionStepPart(VariantPhasingBaseStep):
    """Phasing by transmission."""

    #: Name of the step in the pipeline.
    name = "gatk_phase_by_transmission"

    #: The file name token.
    name_pattern = "gatk_pbt"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{index_library}/out/gatk_pbt.{index_library}"

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            # Pedigree file required for PhaseByTransmission.
            yield (
                "ped",
                f"work/write_pedigree.{wildcards.index_library}/out/{wildcards.index_library}.ped",
            )
            # Get name of real index
            real_index = self.ngs_library_to_pedigree[wildcards.index_library].index
            # Annotated variant file resolved via CDC broker
            upstream_vcf = self.parent.get_upstream_paths(
                self.parent.previous_step, library_name=real_index.dna_ngs_library.name
            )
            vcf = getattr(upstream_vcf, "vcf", None) or upstream_vcf["vcf"]
            vcf_tbi = getattr(upstream_vcf, "vcf_tbi", None) or upstream_vcf.get(
                "vcf_tbi", vcf + ".tbi"
            )
            yield "vcf", vcf
            yield "vcf_tbi", vcf_tbi
            yield "vcf_md5", vcf + ".md5"
            yield "vcf_tbi_md5", vcf_tbi + ".md5"
            yield "reference", self.w_config.static_data_config.reference.path

        assert action == "run", "Unsupported actions"
        return input_function

    def get_args(self, action: str) -> dict[str, Any]:
        # Validate action
        self._validate_action(action)
        return {"de_novo_prior": self.config.gatk_phase_by_transmission.de_novo_prior}

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
            runtime="1d",  # 1 day
            mem=f"{14 * 1024}MB",
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
        if not (
            donor.dna_ngs_library
            and donor.father
            and donor.father.dna_ngs_library
            and donor.mother
            and donor.mother.dna_ngs_library
        ):
            return
        trio_libs = [
            donor.dna_ngs_library.name,
            donor.father.dna_ngs_library.name,
            donor.mother.dna_ngs_library.name,
        ]
        bams = []
        bais = []
        for lib in trio_libs:
            aln: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=lib
            )
            bam = getattr(aln, "bam", None) or aln["bam"]
            bai = getattr(aln, "bai", None) or aln["bai"]
            bams.append(bam)
            bais.append(bai)
        yield "bam", bams
        yield "bai", bais

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
            runtime="1d",  # 1 day
            mem=f"{mem_mb}MB",
        )


class ReadBackedPhasingOnlyStepPart(ReadBackedPhasingBaseStep):
    """Read backed phasing (as primary and only step)"""

    #: Name of the step in the pipeline.
    name = "gatk_read_backed_phasing_only"
    #: The file name token.
    name_pattern = "gatk_rbp"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{index_library}/out/gatk_rbp.{index_library}"

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            real_index = self.ngs_library_to_pedigree[wildcards.index_library].index
            # BAM files from ngs_mapping step.
            yield from self._yield_bams(wildcards)
            # Annotated variant file resolved via CDC broker
            upstream_vcf = self.parent.get_upstream_paths(
                self.parent.previous_step, library_name=real_index.dna_ngs_library.name
            )
            vcf = getattr(upstream_vcf, "vcf", None) or upstream_vcf["vcf"]
            vcf_tbi = getattr(upstream_vcf, "vcf_tbi", None) or upstream_vcf.get(
                "vcf_tbi", vcf + ".tbi"
            )
            yield "vcf", vcf
            yield "vcf_tbi", vcf_tbi
            yield "vcf_md5", vcf + ".md5"
            yield "vcf_tbi_md5", vcf_tbi + ".md5"

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
        self.base_path_out = "work/{index_library}/out/gatk_pbt.gatk_rbp.{index_library}"

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            # BAM files from ngs_mapping step.
            yield from self._yield_bams(wildcards)
            # Result of PhaseByTransmission step
            infix = "work/{index_library}/out/gatk_pbt.{index_library}"
            base_in = infix.format(**wildcards)
            yield "vcf", base_in + ".vcf.gz"
            yield "vcf_tbi", base_in + ".vcf.gz.tbi"
            yield "vcf_md5", base_in + ".vcf.gz.md5"
            yield "vcf_tbi_md5", base_in + ".vcf.gz.tbi.md5"

        assert action == "run", "Unsupported actions"
        return input_function


class VariantPhasingWorkflow(BaseStep):
    """Perform (small) variant phasing"""

    name = "variant_phasing"
    consumes = {DataSignature(DataType.VARIANTS, frozenset({"germline"})): True}
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "phased"}))]
    config_model_class = VariantPhasingConfigModel
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local phased-variant output paths for downstream consumers."""
        cls.require_signature(signature)
        if "phasing" not in kwargs:
            raise ValueError(
                "Parameter 'phasing' is required when requesting output paths from "
                "'variant_phasing' (e.g., phasing='gatk_pbt.gatk_rbp', "
                "phasing='gatk_pbt', or phasing='gatk_rbp')."
            )
        lib = kwargs.get("library_name", "{library_name}")
        phasing = kwargs["phasing"]
        prefix = f"output/{lib}/out/{phasing}.{lib}"
        return {"vcf": f"{prefix}.vcf.gz", "vcf_tbi": f"{prefix}.vcf.gz.tbi"}

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str | None = None,
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=(VariantAnnotationWorkflow, NgsMappingWorkflow),
            task_name=task_name,
            **kwargs,
        )

        for prev in ("variant_annotation", "variant_calling"):
            if getattr(self.config.depends_on, prev, None):
                self.previous_step = prev
                break
        else:
            self.previous_step = "variant_annotation"

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

    @listify
    def get_result_files(self):
        """Return list of result files for the variant phasing workflow."""
        phasings = [
            token for name, token in CONFIG_TO_TOKEN.items() if name in self.config.phasings
        ]
        for phasing_token in phasings:
            yield from self._yield_result_files(
                f"output/{{index_library.name}}/out/{phasing_token}.{{index_library.name}}{{ext}}",
                ext=EXT_VALUES,
            )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
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
                        yield from expand(
                            tpl,
                            index_library=[donor.dna_ngs_library],
                            **kwargs,
                        )
