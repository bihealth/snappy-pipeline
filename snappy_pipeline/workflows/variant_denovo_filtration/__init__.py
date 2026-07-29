# -*- coding: utf-8 -*-
"""Implementation of the ``variant_denovo_filtration`` step.

This step implements filtration of variants to *de novo* variants.  This step was introduced for
the "Ionizing Radiation" study in ca. 2016 and the aim here is to get a set of high-confidence
*de novo* sequence variants (both SNVs and indels, although the latter turned out to be less
reliable).  Further, if the variants are phased, assigning to paternal or maternal allele can
be attempted.  This allows to study paternal age effects.

Note that in contrast to ``variant_calling`` and ``variant_annotation`` but in consistency with
``variant_phasing``, the central individual here are children and not the index of pedigrees.

==========
Step Input
==========

The step reads in the variant call files from one of the following steps:

- ``variant_calling``
- ``variant_annotation``
- ``variant_phasing``

Of course, assignment to parental allele can only be performed on phased variants.  Further, only
filtering annotated variants is really useful as one wants to excludes variants in problematic
genomic regions.

===========
Step Output
===========

For all children with both parents present, variant *de novo* annotation will be attempted on
the primary DNA NGS library of that child.  The name of this library will be used as the
identification token in the output file and file name.
For each pedigree, the following files will be generated:

- ``de_novos.{lib_name}.vcf.gz.tbi``
- ``de_novos.{lib_name}.vcf.gz``
- ``de_novos.{lib_name}.vcf.gz.md5``
- ``de_novos.{lib_name}.vcf.gz.tbi.md5``
- ``de_novos_hard.{lib_name}.vcf.gz``
- ``de_novos_hard.{lib_name}.vcf.gz.tbi``
- ``de_novos_hard.{lib_name}.vcf.gz.md5``
- ``de_novos_hard.{lib_name}.vcf.gz.tbi.md5``
- ``de_novos_hard.{lib_name}.summary.txt``
- ``de_novos_hard.{lib_name}.summary.txt.md5``

The the ``annotation`` and ``phasing`` will only be persent when the input is read from the
``variant_annotation`` or ``variant_phasing`` steps, respectively.

====================
Global Configuration
====================

No global configuration is in use.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_denovo_filtration.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import itertools
import os
from collections import OrderedDict
from typing import Any

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow

from .model import VariantDenovoFiltration as VariantDenovoFiltrationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Default configuration for the variant_denovo_filtration step
DEFAULT_CONFIG = VariantDenovoFiltrationConfigModel.default_config_yaml_string()


class FilterDeNovosBaseStepPart(BaseStepPart):
    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        #: Name of the previous step and token
        self.previous_step = self.parent.previous_step
        self.prev_token = self.parent.prev_token
        #: Mapping from ngs_library to pedigree, only used when previous input is not
        #: variant_phasing.
        self.ngs_library_to_pedigree = OrderedDict()
        self.ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if donor.dna_ngs_library:
                        self.ngs_library_to_pedigree[donor.dna_ngs_library.name] = pedigree
                        self.ngs_library_to_donor[donor.dna_ngs_library.name] = donor

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
            runtime="2h",  # 2 hours
            mem=f"{2 * 1024}MB",
        )


class FilterDeNovosStepPart(FilterDeNovosBaseStepPart):
    """Step for soft-filtering variants (annotations and adding soft-filters)."""

    #: Step name
    name = "filter_denovo"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{index_library}/out/{index_library}.soft"
        self.path_log = "work/{index_library}/log/filter_denovo.{index_library}.log"

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            ped_entry = self.ngs_library_to_pedigree.get(wildcards.index_library)
            if not ped_entry:
                return {}
            real_index = ped_entry.index
            real_path = f"work/write_pedigree.{real_index.dna_ngs_library.name}/out/{real_index.dna_ngs_library.name}.ped"
            yield "ped", real_path

            _aln: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=wildcards.index_library
            )
            yield "bam", getattr(_aln, "bam", None)
            yield "bai", getattr(_aln, "bai", None)

            extra_kwargs = {}
            if self.previous_step == "variant_phasing":
                phasing_cfg = self.parent.get_task_config("variant_phasing")
                phasings = getattr(phasing_cfg, "phasings", ["gatk_phasing_both"])
                token_map = {
                    "gatk_read_backed_phasing": "gatk_rbp",
                    "gatk_phase_by_transmission": "gatk_pbt",
                    "gatk_phasing_both": "gatk_pbt.gatk_rbp",
                }
                extra_kwargs["phasing"] = token_map.get(phasings[0], "gatk_pbt.gatk_rbp")

            upstream_vcf = self.parent.get_upstream_paths(
                self.previous_step,
                library_name=real_index.dna_ngs_library.name,
                **extra_kwargs,
            )
            yield "vcf", getattr(upstream_vcf, "vcf", None)
            yield "vcf_tbi", getattr(upstream_vcf, "vcf_tbi", None)

        return input_function

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out + ext

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.path_log

    def get_args(self, action: str):
        self._validate_action(action)

        def args_fn(wildcards: Wildcards) -> dict[str, Any]:
            return {
                "besenbacher": self.config.params_besenbacher.model_dump(by_alias=True),
                "index_library": wildcards.index_library,
            }

        return args_fn

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


class FilterDeNovosHardStepPart(FilterDeNovosBaseStepPart):
    """Step for hard-filtering variants."""

    #: Step name
    name = "filter_denovo_hard"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/{index_library}/out/{index_library}.soft"
        self.base_path_out = "work/{index_library}/out/{index_library}"
        self.path_log = "work/{index_library}/log/filter_denovo_hard.{index_library}.log"

    @dictify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "vcf", self.base_path_in + ".vcf.gz"
        yield "vcf_tbi", self.base_path_in + ".vcf.gz.tbi"

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out + ext
        yield "summary", self.base_path_out + ".summary.txt"
        yield "summary_md5", self.base_path_out + ".summary.txt.md5"

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.path_log

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        def args_function(wildcards: Wildcards) -> dict[str, Any]:
            donor = self.ngs_library_to_donor[wildcards.index_library]
            return {
                "index_library": wildcards.index_library,
                "father": donor.father.dna_ngs_library.name,
                "mother": donor.mother.dna_ngs_library.name,
                "bad_regions_expression": self.config.bad_regions_expression,
            }

        return args_function


class SummarizeCountsStepPart(FilterDeNovosBaseStepPart):
    """Summarizing counts."""

    #: Step name
    name = "summarize_counts"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        self.name_pattern = "summarize_counts"
        self.base_path_out = os.path.join("work", self.name_pattern, "out", self.name_pattern)
        self.path_log = os.path.join("work", self.name_pattern, "log", self.name_pattern + ".log")

    @listify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        name_pattern = "%sde_novos_hard.{index_library.name}" % (self.prev_token,)
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if not donor.dna_ngs_library:
                        continue
                    elif not donor.father or not donor.father.dna_ngs_library:
                        continue
                    elif not donor.mother or not donor.mother.dna_ngs_library:
                        continue
                    else:
                        yield from expand(
                            os.path.join("work", name_pattern, "out", name_pattern + "{ext}"),
                            index_library=[donor.dna_ngs_library],
                            ext=(".summary.txt",),
                        )

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "txt", self.base_path_out + ".txt"
        yield "txt_md5", self.base_path_out + ".txt.md5"

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.path_log


class CollectMsdnStepPart(FilterDeNovosBaseStepPart):
    """Step part for collecting the MSDN."""

    #: Step name
    name = "collect_msdn"

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        result = {"gatk3_hc": [], "gatk_ug": []}
        name_pattern = "%sde_novos_hard.{index_library}" % (self.prev_token,)
        tpl = "work/" + name_pattern + "/out/" + name_pattern + ".summary.txt"
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if not donor.dna_ngs_library:
                        continue
                    elif not donor.father or not donor.father.dna_ngs_library:
                        continue
                    elif not donor.mother or not donor.mother.dna_ngs_library:
                        continue
                    else:
                        for caller in result.keys():
                            result[caller].append(
                                tpl.format(
                                    mapper="",
                                    caller=caller,
                                    index_library=donor.dna_ngs_library.name,
                                )
                            )
        return result

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "txt", "work/multisite_de_novo/out/multisite_de_novo.txt"
        yield "txt_md5", "work/multisite_de_novo/out/multisite_de_novo.txt.md5"

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return "work/multisite_de_novo/log/multisite_de_novo.log"


class SummarizeDeNovoCountsStepPart(FilterDeNovosBaseStepPart):
    """Step part for creating summary counts."""

    #: Step name
    name = "summarize_counts"

    @listify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        name_pattern = "%sde_novos_hard.{index_library}" % (self.prev_token,)
        tpl = "work/" + name_pattern + "/out/" + name_pattern + ".summary.txt"
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if not donor.dna_ngs_library:
                        continue
                    elif not donor.father or not donor.father.dna_ngs_library:
                        continue
                    elif not donor.mother or not donor.mother.dna_ngs_library:
                        continue
                    else:
                        yield tpl.format(
                            index_library=donor.dna_ngs_library.name,
                        )

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "txt", "work/denovo_count_summary/out/denovo_count_summary.txt"
        yield (
            "txt_md5",
            ("work/denovo_count_summary/out/denovo_count_summary.txt.md5"),
        )

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return "work/denovo_count_summary/log/denovo_count_summary.log"


class VariantDeNovoFiltrationWorkflow(BaseStep):
    """Perform (small) variant de novo filtration"""

    #: Workflow name
    name = "variant_denovo_filtration"
    consumes = {DataSignature(DataType.VARIANTS, frozenset({"germline"})): True}
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "denovo"}))]
    config_model_class = VariantDenovoFiltrationConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local de-novo filtration output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"output/{lib}/out/{lib}"
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
            previous_steps=(VariantPhasingWorkflow, VariantAnnotationWorkflow, NgsMappingWorkflow),
            task_name=task_name,
            **kwargs,
        )
        for prev in ("variant_phasing", "variant_annotation", "variant_calling"):
            if getattr(self.config.depends_on, prev, None):
                self.previous_step = prev
                break
        else:
            raise Exception("No previous step given!")  # pragma: no cover
        #: Name token for input
        self.prev_token = {
            "variant_phasing": "jannovar_annotate_vcf.gatk_pbt.gatk_rbp.",
            "variant_annotation": "jannovar_annotate_vcf.",
            "variant_calling": "",
        }[self.previous_step]
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                (WritePedigreeStepPart, (True, True)),
                FilterDeNovosStepPart,
                FilterDeNovosHardStepPart,
                CollectMsdnStepPart,
                SummarizeDeNovoCountsStepPart,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the variant de novo filtration workflow."""
        # Hard-filtered results
        ext_values = list(itertools.chain(EXT_VALUES, (".summary.txt", ".summary.txt.md5")))
        yield from self._yield_result_files(
            "output/{index_library.name}/out/{index_library.name}{ext}",
            ext=ext_values,
        )
        # Summarise counts
        yield from expand(
            "output/denovo_count_summary/out/denovo_count_summary{ext}",
            ext=(".txt", ".txt.md5"),
        )
        # Collect MSDN statistics
        if self.get_task_config("variant_denovo_filtration").collect_msdn:
            yield from expand(
                "output/multisite_de_novo/out/multisite_de_novo{ext}",
                ext=(".txt", ".txt.md5"),
            )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                any_trio = False
                for donor in pedigree.donors:  # skip if no trio in pedigree
                    if (
                        donor.dna_ngs_library
                        and donor.father
                        and donor.father.dna_ngs_library
                        and donor.mother
                        and donor.mother.dna_ngs_library
                    ):
                        any_trio = True
                        break
                if any_trio:
                    for donor in pedigree.donors:
                        if not donor.dna_ngs_library:
                            continue
                        elif not donor.father or not donor.father.dna_ngs_library:
                            continue
                        elif not donor.mother or not donor.mother.dna_ngs_library:
                            continue
                        else:
                            yield from expand(
                                tpl,
                                index_library=[donor.dna_ngs_library],
                                **kwargs,
                            )

    def check_config(self):
        pass
