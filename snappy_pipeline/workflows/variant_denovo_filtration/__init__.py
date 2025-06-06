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
identification token in the output file and file name.  For each read mapper, variant caller,
and pedigree, the following files will be generated:

- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos.{lib_name}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos.{lib_name}.vcf.gz``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos.{lib_name}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos.{lib_name}.vcf.gz.tbi.md5``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos_hard.{lib_name}.vcf.gz``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos_hard.{lib_name}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos_hard.{lib_name}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos_hard.{lib_name}.vcf.gz.tbi.md5``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos_hard.{lib_name}.summary.txt``
- ``{mapper}.{var_caller}.{annotation}.{phasing}.de_novos_hard.{lib_name}.summary.txt.md5``

The the ``annotation`` and ``phasing`` will only be persent when the input is read from the
``variant_annotation`` or ``variant_phasing`` steps, respectively.

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.gatk3_hc.de_novos.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.gatk3_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk3_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz.md5
    |       |-- bwa.gatk3_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.gatk3_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz.tbi.md5
    |       |-- bwa.gatk3_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk3_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz.md5
    |       |-- bwa.gatk3_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.gatk3_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz.tbi.md5
    |       |-- bwa.gatk3_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk3_hc.de_novos_hard.P001-N1-DNA1-WES1.summary.txt
    |       `-- bwa.gatk3_hc.de_novos_hard.P001-N1-DNA1-WES1.summary.txt.md5
    [...]

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

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
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
            time="02:00:00",  # 2 hours
            memory=f"{2 * 1024}M",
        )


class FilterDeNovosStepPart(FilterDeNovosBaseStepPart):
    """Step for soft-filtering variants (annotations and adding soft-filters)."""

    #: Step name
    name = "filter_denovo"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        self.name_pattern = r"{mapper}.{caller}.%sde_novos.{index_library,[^\.]+}" % (
            self.prev_token,
        )
        self.base_path_out = os.path.join(
            "work", self.name_pattern, "out", self.name_pattern.replace(r",[^\.]+", "")
        )
        self.path_log = os.path.join(
            "work", self.name_pattern, "log", self.name_pattern.replace(r",[^\.]+", "") + ".log"
        )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            # Get name of real index, used when input is not variant_phasing
            real_index = self.ngs_library_to_pedigree[wildcards.index_library].index
            # Pedigree file required for PhaseByTransmission.
            real_path = "work/write_pedigree.{real_index}/out/{real_index}.ped".format(
                real_index=real_index.dna_ngs_library.name, **wildcards
            )
            yield "ped", real_path
            # BAM and BAI file of the offspring
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            path_bam = ("output/{mapper}.{index_library}/out/{mapper}.{index_library}.bam").format(
                **wildcards
            )
            yield "bam", ngs_mapping(path_bam)
            yield "bai", ngs_mapping(path_bam + ".bai")
            # Input file comes from previous step.
            prev_step = self.parent.sub_workflows[self.previous_step]
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                name_pattern = self.name_pattern.replace(r",[^\.]+", "").replace("de_novos.", "")
                if self.previous_step != "variant_phasing":
                    name_pattern = name_pattern.replace("{index_library}", "{real_index}")
                input_path = ("output/" + name_pattern + "/out/" + name_pattern).format(
                    real_index=real_index.dna_ngs_library.name, **wildcards
                )
                yield key, prev_step(input_path) + ext

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


class FilterDeNovosHardStepPart(FilterDeNovosBaseStepPart):
    """Step for hard-filtering variants."""

    #: Step name
    name = "filter_denovo_hard"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        self.name_pattern = r"{mapper}.{caller}.%sde_novos_hard.{index_library,[^\.]+}" % (
            self.prev_token,
        )
        self.base_path_out = os.path.join(
            "work", self.name_pattern, "out", self.name_pattern.replace(r",[^\.]+", "")
        )
        self.base_path_in = self.base_path_out.replace("de_novos_hard", "de_novos")
        self.path_log = os.path.join(
            "work", self.name_pattern, "log", self.name_pattern.replace(r",[^\.]+", "") + ".log"
        )

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

        def args_function(wildcards):
            donor = self.ngs_library_to_donor[wildcards.index_library]
            return {
                "father": donor.father.dna_ngs_library.name,
                "mother": donor.mother.dna_ngs_library.name,
            }

        return args_function


class SummarizeCountsStepPart(FilterDeNovosBaseStepPart):
    """Summarizing counts."""

    #: Step name
    name = "summarize_counts"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        self.name_pattern = "{mapper}.{caller}.summarize_counts"
        self.base_path_out = os.path.join("work", self.name_pattern, "out", self.name_pattern)
        self.path_log = os.path.join("work", self.name_pattern, "log", self.name_pattern + ".log")

    @listify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        name_pattern = "{{mapper}}.{{caller}}.%sde_novos_hard.{index_library.name}" % (
            self.prev_token,
        )
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
        name_pattern = "{mapper}.{caller}.%sde_novos_hard.{index_library}" % (self.prev_token,)
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
                                    mapper="{mapper}",
                                    caller=caller,
                                    index_library=donor.dna_ngs_library.name,
                                )
                            )
        return result

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "txt", "work/{mapper}.multisite_de_novo/out/{mapper}.multisite_de_novo.txt"
        yield "txt_md5", "work/{mapper}.multisite_de_novo/out/{mapper}.multisite_de_novo.txt.md5"

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return "work/{mapper}.multisite_de_novo/log/{mapper}.multisite_de_novo.log"


class SummarizeDeNovoCountsStepPart(FilterDeNovosBaseStepPart):
    """Step part for creating summary counts."""

    #: Step name
    name = "summarize_counts"

    @listify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        name_pattern = "{mapper}.{caller}.%sde_novos_hard.{index_library}" % (self.prev_token,)
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
                        for caller in self.config.tools_variant_calling:
                            yield tpl.format(
                                mapper="{mapper}",
                                caller=caller,
                                index_library=donor.dna_ngs_library.name,
                            )

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "txt", "work/{mapper}.denovo_count_summary/out/{mapper}.denovo_count_summary.txt"
        yield (
            "txt_md5",
            ("work/{mapper}.denovo_count_summary/out/{mapper}.denovo_count_summary.txt.md5"),
        )

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return "work/{mapper}.denovo_count_summary/log/{mapper}.denovo_count_summary.log"


class VariantDeNovoFiltrationWorkflow(BaseStep):
    """Perform (small) variant de novo filtration"""

    #: Workflow name
    name = "variant_denovo_filtration"

    #: Default biomed sheet class
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
            config_model_class=VariantDenovoFiltrationConfigModel,
            previous_steps=(VariantPhasingWorkflow, VariantAnnotationWorkflow, NgsMappingWorkflow),
        )
        # Register sub workflows
        for prev in ("variant_phasing", "variant_annotation", "variant_calling"):
            if cfg := self.config.get(f"path_{prev}"):
                self.previous_step = prev
                self.register_sub_workflow(prev, cfg)
                break
        else:
            raise Exception("No path to previous step given!")  # pragma: no cover
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
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
        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.tools_variant_calling:
            self.config.tools_variant_calling = self.w_config.step_config["variant_calling"].tools

    @listify
    def get_result_files(self):
        """Return list of result files for the variant de novo filtration workflow."""
        # Hard-filtered results
        name_pattern = "{mapper}.{caller}.%sde_novos_hard.{index_library.name}" % (self.prev_token,)
        ext_values = list(itertools.chain(EXT_VALUES, (".summary.txt", ".summary.txt.md5")))
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            ext=ext_values,
        )
        # Summarise counts
        yield from expand(
            "output/{mapper}.denovo_count_summary/out/{mapper}.denovo_count_summary{ext}",
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            ext=(".txt", ".txt.md5"),
        )
        # Collect MSDN statistics
        if self.w_config.step_config["variant_denovo_filtration"].collect_msdn:
            yield from expand(
                "output/{mapper}.multisite_de_novo/out/{mapper}.multisite_de_novo{ext}",
                mapper=self.config.tools_ngs_mapping,
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
                            yield from expand(tpl, index_library=[donor.dna_ngs_library], **kwargs)

    def check_config(self):
        if not self.config.tools_ngs_mapping:
            self.ensure_w_config(
                ("step_config", "ngs_mapping", "tools", "dna"),
                "Either define tools_ngs_mapping or provide a configuration for ngs_mapping",
            )
        if not self.config.tools_variant_calling:
            self.ensure_w_config(
                ("step_config", "variant_calling", "tools"),
                "Either define tools_variant_calling or provide a configuration for variant_calling",
            )
