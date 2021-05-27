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
    +-- bwa.gatk_hc.de_novos.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.gatk_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz.md5
    |       |-- bwa.gatk_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.gatk_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz.tbi.md5
    |       |-- bwa.gatk_hc.de_novos.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz.md5
    |       |-- bwa.gatk_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.gatk_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz.tbi.md5
    |       |-- bwa.gatk_hc.de_novos_hard.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.gatk_hc.de_novos_hard.P001-N1-DNA1-WES1.summary.txt
    |       `-- bwa.gatk_hc.de_novos_hard.P001-N1-DNA1-WES1.summary.txt.md5
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

from collections import OrderedDict
import itertools
import os

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Default configuration for the variant_denovo_filtration step
DEFAULT_CONFIG = r"""
step_config:
  variant_denovo_filtration:
    # One of the following must be given!
    path_variant_phasing: ''
    path_variant_annotation: ''
    path_variant_calling: ''
    path_ngs_mapping: ../ngs_mapping
    tools_ngs_mapping: null          # defaults to ngs_mapping tool
    tools_variant_calling: null      # defaults to variant_annotation tool
    info_key_reliable_regions: []    # optional INFO keys with reliable regions
    info_key_unreliable_regions: []  # optional INFO keys with unreliable regions
    params_besenbacher:              # parameters for Besenbacher quality filter
      min_gq: 50
      min_dp: 10
      max_dp: 120
      min_ab: 0.20
      max_ab: 0.9
      max_ad2: 1
    bad_region_expressions: []
    # e.g.,
    # - 'UCSC_CRG_MAPABILITY36 == 1'
    # - 'UCSC_SIMPLE_REPEAT == 1'
    collect_msdn: True               # whether or not to collect MSDN (requires GATK HC+UG)
"""


class FilterDeNovosBaseStepPart(BaseStepPart):
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

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_denovo_filtration_{}_run".format(self.name)] = {
            "mem": 2 * 1024,
            "time": "02:00",
            "ntasks": 1,
        }


class FilterDeNovosStepPart(FilterDeNovosBaseStepPart):
    """Step for soft-filtering variants (annotations and adding soft-filters)."""

    name = "filter_denovo"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        self.name_pattern = "{mapper}.{caller}.%sde_novos.{index_library,[^\.]+}" % (
            self.prev_token,
        )
        self.base_path_out = os.path.join(
            "work", self.name_pattern, "out", self.name_pattern.replace(r",[^\.]+", "")
        )
        self.path_log = os.path.join(
            "work", self.name_pattern, "log", self.name_pattern.replace(r",[^\.]+", "") + ".log"
        )

    def get_input_files(self, action):
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
            path_bam = (
                "output/{mapper}.{index_library}/out/{mapper}." "{index_library}.bam"
            ).format(**wildcards)
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

        assert action == "run", "Unsupported action"
        return input_function

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out + ext

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_denovo_filtration_{}_run".format(self.name)] = {
            "mem": 14 * 1024,
            "time": "24:00",
            "ntasks": 1,
        }


class FilterDeNovosHardStepPart(FilterDeNovosBaseStepPart):
    """Step for hard-filtering variants."""

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
        yield "vcf", self.base_path_in + ".vcf.gz"
        yield "tbi", self.base_path_in + ".vcf.gz.tbi"

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out + ext
        yield "summary", self.base_path_out + ".summary.txt"
        yield "summary_md5", self.base_path_out + ".summary.txt.md5"

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log

    def get_args(self, action):
        def args_function(wildcards):
            donor = self.ngs_library_to_donor[wildcards.index_library]
            return {
                "father": donor.father.dna_ngs_library.name,
                "mother": donor.mother.dna_ngs_library.name,
            }

        assert action == "run"
        return args_function


class SummarizeCountsStepPart(FilterDeNovosBaseStepPart):
    """Summarizing counts."""

    name = "summarize_counts"

    def __init__(self, parent):
        super().__init__(parent)
        # Output and log paths
        self.name_pattern = "{mapper}.{caller}.summarize_counts"
        self.base_path_out = os.path.join("work", self.name_pattern, "out", self.name_pattern)
        self.path_log = os.path.join("work", self.name_pattern, "log", self.name_pattern + ".log")

    @listify
    def get_input_files(self, action):
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
        assert action == "run"
        yield "txt", self.base_path_out + ".txt"
        yield "txt_md5", self.base_path_out + ".txt.md5"

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log


class CollectMsdnStepPart(FilterDeNovosBaseStepPart):
    """Step part for collecting the MSDN."""

    name = "collect_msdn"

    def get_input_files(self, action):
        result = {"gatk_hc": [], "gatk_ug": []}
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
        assert action == "run"
        yield "txt", "work/{mapper}.multisite_de_novo/out/{mapper}.multisite_de_novo.txt"
        yield "txt_md5", "work/{mapper}.multisite_de_novo/out/{mapper}.multisite_de_novo.txt.md5"

    def get_log_file(self, action):
        assert action == "run"
        return "work/{mapper}.multisite_de_novo/log/{mapper}.multisite_de_novo.log"


class SummarizeDeNovoCountsStepPart(FilterDeNovosBaseStepPart):
    """Step part for creating summary counts."""

    name = "summarize_counts"

    @listify
    def get_input_files(self, action):
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
                        for caller in self.config["tools_variant_calling"]:
                            yield tpl.format(
                                mapper="{mapper}",
                                caller=caller,
                                index_library=donor.dna_ngs_library.name,
                            )

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        yield "txt", "work/{mapper}.denovo_count_summary/out/{mapper}.denovo_count_summary.txt"
        yield "txt_md5", (
            "work/{mapper}.denovo_count_summary/out/{mapper}.denovo_count_summary.txt.md5"
        )

    def get_log_file(self, action):
        assert action == "run"
        return "work/{mapper}.denovo_count_summary/log/{mapper}.denovo_count_summary.log"


class VariantDeNovoFiltrationWorkflow(BaseStep):
    """Perform (small) variant de novo filtration"""

    name = "variant_denovo_filtration"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
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
            (VariantPhasingWorkflow, VariantAnnotationWorkflow, NgsMappingWorkflow),
        )
        # Register sub workflows
        for prev in ("variant_phasing", "variant_annotation", "variant_calling"):
            if self.config["path_%s" % prev]:
                self.previous_step = prev
                self.register_sub_workflow(prev, self.config["path_%s" % prev])
                break
        else:
            raise Exception("No path to previous step given!")  # pragma: no cover
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
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
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_variant_calling"]:
            self.config["tools_variant_calling"] = self.w_config["step_config"]["variant_calling"][
                "tools"
            ]

    @listify
    def get_result_files(self):
        """Return list of result files for the variant de novo filtration workflow."""
        # Hard-filtered results
        name_pattern = "{mapper}.{caller}.%sde_novos_hard.{index_library.name}" % (self.prev_token,)
        ext_values = list(itertools.chain(EXT_VALUES, (".summary.txt", ".summary.txt.md5")))
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
            ext=ext_values,
        )
        # Summarise counts
        yield from expand(
            "output/{mapper}.denovo_count_summary/out/{mapper}.denovo_count_summary{ext}",
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
            ext=(".txt", ".txt.md5"),
        )
        # Collect MSDN statistics
        if self.w_config["step_config"]["variant_denovo_filtration"]["collect_msdn"]:
            yield from expand(
                "output/{mapper}.multisite_de_novo/out/{mapper}.multisite_de_novo{ext}",
                mapper=self.config["tools_ngs_mapping"],
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
        """Check that the path to the variant annotation step is present."""
        self.ensure_w_config(
            ("step_config", "variant_denovo_filtration", "path_ngs_mapping"),
            ("Path to ngs_mapping not configured but required for variant_denovo_filtration"),
        )
