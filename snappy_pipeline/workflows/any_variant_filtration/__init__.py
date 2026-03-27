# -*- coding: utf-8 -*-
"""Generic variant filtration step

The assumption is that variant filtration is similar for all kind of variants,
germline or somatic.

This implementation attempts to filter a set of variants stored in vcf format,
using a generic sample sheet (_i.e._ without restriction on tumor/normal pairs
or the availability of pedigree).

Unfortunately, I haven't found a way to design a completely generic implementation.
This is because during the step initialisation, the pipeline validates the configuration
of previous steps, and there is no common configuration of the germline & somatic variant
calling steps. So there is no configration possible of a generic variant calling step.
Therefore, which variant calling step has been run must be known at the time of the
initialisation of the variant filtration step. So a workaround has been implemented,
and the user must set a ``variant_origin`` parameter in order to define which type of
variant calling step has been run. It is a poor solution, because it prevents the
filtration of both somatic & germline variants in the same pipelie run, until the
configuration of steps is replaced by a list. However, it is possible to annotate a
combined vcf with both germline & somatic variants, which is needed by the neo-epitope
step, for example. But of course, the file naming is another problem, which is dealt
with by the ``combined_phasing`` step for example.

==========
Step Input
==========

For each library, a vcf following the naming pattern ``{mapper}.{caller}.{library}.vcf.gz``
in the folder ``{path_variant}/output/{mapper}.{caller}.{library}/out``, where
``path_variant`` is defined by the user in the configuration.

Typically, ``{caller}`` is ``mutect2`` for somatic variants, ``gatk_hc`` for germline
ones, and ``combined`` after the ``combined_phasing`` step.

If the vcf has been annotated before filtration, the naming pattern becomes
``{mapper}.{caller}.{annotator}.{library}`` for the folder and the file name.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name`` and each read mapper
``mapper`` that the library has been aligned with, and the variant caller ``caller``, the
pipeline step will create a directory ``output/{mapper}.{caller}.filtered.{lib_name}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

Two ``vcf`` files are produced:

- ``{mapper}.{caller}.filtered.{lib_name}.vcf.gz`` which contains only the variants that have passed all filters, or that were protected, and
- ``{mapper}.{caller}.filtered.{lib_name}.full.vcf.gz`` which contains all variants, with the reason for rejection in the ``FILTER`` column.

Of course, when annotation has been done prior to filtration, an ``.{annotator}``
will be inserted between the caller and ``.filtered``.

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.mutect2.vep.filtered.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.vcf.gz
    |       |-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.vcf.gz.md5
    |       `-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.vcf.gz.tbi.md5
    |       |-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.full.vcf.gz
    |       |-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.full.vcf.gz.tbi
    |       |-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.full.vcf.gz.md5
    |       `-- bwa.mutect2.vep.filtered.P001-T1-DNA1-WES1.full.vcf.gz.tbi.md5
    [...]

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_filtration.rst

=======
Filters
=======

The following filters are implemented:

- ``dkfz``: uses orientiation biases to remove sequencing & PCR artifacts.
  The current implementation doesn't allow any parametrisation of this filter.
  This filter will add ``bSeq`` or ``bPcr`` to the FILTER column of rejected variants.
- ``ebfilter``: Bayesian statistical model to score variants.
  Variants with a score lower than ``ebfilter_threshold`` are rejected.
  The scoring algorithm can be parameterised from the coniguration.
  This filter will add ``ebfilter_<n>`` to the FILTER column of rejected variants.
  The implementation of this filter is old, and was tuned against a specific dataset.
  In that sense, its use is discouraged (it should be replaced by ``varlociraptor``)
- ``bcftools``: flexible filter based on `bcftools expressions <https://samtools.github.io/bcftools/bcftools.html#expressions>`_.
  The expression can be designed to ``include`` or ``exclude`` variants.
  This filter will add ``bcftools_<n>`` to the FILTER column of rejected variants.
- ``regions``: filter to exclude variants outside of user's defined regions.
  Typically used to reject variants outside of coding regions.
  This filter will add ``regions_<n>`` to the FILTER column of rejected variants.
- ``protected``: anti-filter to avoid variants in protected regions to be otherwise filtered out.
  This filter "whitelists" variants in specific regions. This is valuable to protect
  known drivers against being filtered out, even if there is little experimental support for them.
  This filter will add ``PROTECTED`` to the FILTER column of rejected variants.

In the above description, ``<n>`` is here the sequence number of the filter in the filter list.

The filters can be used or not, and can be used multiple times. For example, it is possible to
use the ``bcftools`` filter to reject differentially potential FFPE artifacts. The filter list would then be:

.. code-block:: yaml

  filter_list:
  - dkfz: {}
  - bcftools:
    exclude: "AD[1:0]+AD[1:1]<50 | AD[1:1]<5 | AD[1:1]/(AD[1:0]+AD[1:1])<0.05"
  - bcftools:
    exclude: "((REF='C' & ALT='T') | (REF='G' & ALT='A')) & AD[1:1]/(AD[1:0]+AD[1:1])<0.10"
  - protected:
    path_bed: hotspots_locii.bed

This list of filters would apply the DKFZBiasFilter, reject all variants with depth lower than 50, less than 5 reads supporting the alternative allele, or with a variant allele fraction below 5%.
It would also reject all C-to-T and G-to-A variants with a VAF lower than 10%, because they might be FFPE artifacts.
All variants overlapping with hotspots locii would be protected against filtration.
"""

import os
from typing import Any

from snakemake.io import expand, Wildcards
from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.samplesheet import filter_table_by_modality, sample_sheets
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from .model import AnyVariantFiltration as AnyVariantFiltrationConfigModel


__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Extensions of log files (TODO: should go in a generic location)
LOG_EXT_VALUES = (
    ".log",
    ".log.md5",
    ".conda_info.txt",
    ".conda_info.txt.md5",
    ".conda_list.txt",
    ".conda_list.txt.md5",
    ".merged.tar.gz",
    ".merged.tar.gz.md5",
)

#: Default configuration for the any_variant_calling step
DEFAULT_CONFIG = AnyVariantFiltrationConfigModel.default_config_yaml_string()


class AnyVariantFiltrationStepPart(BaseStepPart):
    """Shared code for all tools in somatic_variant_filtration"""

    def __init__(self, parent):
        super().__init__(parent)
        self.config = parent.config
        self.name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            self.name_pattern += ".{annotator}"
        self.name_pattern += ".{tumor_library}"


class OneFilterStepPart(AnyVariantFiltrationStepPart):
    """Performs one filtration step using checkpoints rather than rules"""

    #: Step name
    name = "one_filter"

    #: Class available actions
    actions = ("run",)

    #: Default filtration resource usage (should be light)
    resource_usage = {"run": ResourceUsage(threads=1, time="02:00:00", memory=f"{8 * 1024}M")}

    def get_input_files(self, action):
        """Return path to input or previous filter vcf file"""
        # Validate action
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards):
        filter_nb = int(wildcards["filter_nb"])
        name_pattern = self.name_pattern.format(**wildcards)
        if filter_nb > 1:
            prev = list(self.config.filter_list[filter_nb - 2].keys())[0]
            n = filter_nb - 1
            yield (
                "vcf",
                os.path.join("work", name_pattern, "out", name_pattern + f".{prev}_{n}.vcf.gz"),
            )
        else:
            variant = self.parent.sub_workflows["variant"]
            base_path = os.path.join("output", name_pattern, "out", name_pattern)
            yield "vcf", variant(base_path.format(**wildcards) + ".vcf.gz")

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        prefix = os.path.join(
            "work",
            self.name_pattern,
            "out",
            self.name_pattern + "." + self.filter_name + "_{filter_nb}",
        )
        key_ext = {
            "vcf": ".vcf.gz",
            "vcf_tbi": ".vcf.gz.tbi",
            "vcf_md5": ".vcf.gz.md5",
            "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in key_ext.items():
            yield key, prefix + ext

    @dictify
    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield (
                key,
                os.path.join(
                    "work",
                    self.name_pattern,
                    "log",
                    self.name_pattern + "." + self.filter_name + "_{filter_nb}" + ext,
                ),
            )
            yield (
                key + "_md5",
                os.path.join(
                    "work",
                    self.name_pattern,
                    "log",
                    self.name_pattern + "." + self.filter_name + "_{filter_nb}" + ext + ".md5",
                ),
            )

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        return self._get_args

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
        filter_nb = int(wildcards["filter_nb"])
        params = self.config.filter_list[filter_nb - 1][self.filter_name].model_dump(by_alias=True)
        params["filter_name"] = "{}_{}".format(self.filter_name, wildcards["filter_nb"])
        return params


class OneFilterWithBamStepPart(OneFilterStepPart):
    @dictify
    def _get_input_files_run(self, wildcards):
        parent = super(OneFilterWithBamStepPart, self)._get_input_files_run
        yield from parent(wildcards).items()

        yield "reference", self.w_config.static_data_config.reference.path

        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        name_pattern = "{mapper}.{tumor_library}".format(**wildcards)
        base_path = os.path.join("output", name_pattern, "out", name_pattern)
        yield "bam", ngs_mapping(base_path + ".bam")


class OneFilterDkfzStepPart(OneFilterWithBamStepPart):
    name = "one_dkfz"
    filter_name = "dkfz"
    resource_usage = {"run": ResourceUsage(threads=1, time="12:00:00", memory=f"{3 * 1024}M")}


class OneFilterBcftoolsStepPart(OneFilterStepPart):
    name = "one_bcftools"
    filter_name = "bcftools"


class OneFilterRegionsStepPart(OneFilterStepPart):
    name = "one_regions"
    filter_name = "regions"


class OneFilterProtectedStepPart(OneFilterStepPart):
    name = "one_protected"
    filter_name = "protected"


class LastFilterStepPart(AnyVariantFiltrationStepPart):
    """Mark last filter as final output"""

    #: Step name
    name = "last_filter"

    #: Class available actions
    actions = ("run",)

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        filter_names = [list(filter_name.keys())[0] for filter_name in self.config.filter_list]
        filter_nb = len(self.config.filter_list)
        filter_name = filter_names[filter_nb - 1]
        vcf = os.path.join(
            "work",
            self.name_pattern,
            "out",
            self.name_pattern + f".{filter_name}_{filter_nb}.vcf.gz",
        )
        prefix = os.path.join("work", self.name_pattern, "log", self.name_pattern)
        logs = [
            prefix + "." + filter_name + "_" + str(filter_nb + 1) + "." + e + m
            for filter_nb, filter_name in enumerate(filter_names)
            for e in ("log", "conda_list.txt", "conda_info.txt")
            for m in ("", ".md5")
        ]
        return {"vcf": vcf, "logs": logs}

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        name_pattern += ".filtered.{tumor_library}"
        vcf = os.path.join("work", name_pattern, "out", name_pattern)
        merged_log = os.path.join("work", name_pattern, "log", name_pattern + ".merged.tar.gz")
        return {
            "vcf": vcf + ".vcf.gz",
            "vcf_tbi": vcf + ".vcf.gz.tbi",
            "vcf_md5": vcf + ".vcf.gz.md5",
            "vcf_tbi_md5": vcf + ".vcf.gz.tbi.md5",
            "full": vcf + ".full.vcf.gz",
            "full_tbi": vcf + ".full.vcf.gz.tbi",
            "full_md5": vcf + ".full.vcf.gz.md5",
            "full_tbi_md5": vcf + ".full.vcf.gz.tbi.md5",
            "log": merged_log,
            "log_md5": merged_log + ".md5",
        }

    @dictify
    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        name_pattern += ".filtered.{tumor_library}"
        tpl = os.path.join("work", name_pattern, "log", name_pattern)
        return {
            "log": tpl + ".log",
            "log_md5": tpl + ".log.md5",
            "conda_list": tpl + ".conda_list.txt",
            "conda_list_md5": tpl + ".conda_list.txt.md5",
            "conda_info": tpl + ".conda_info.txt",
            "conda_info_md5": tpl + ".conda_info.txt.md5",
        }


class AnyVariantFiltrationWorkflow(BaseStep):
    """Perform somatic variant filtration"""

    #: Workflow name
    name = "any_variant_filtration"
    variant_origin = VariantOrigin.ANY
    model_class = AnyVariantFiltrationConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        # Ugly hack to allow exchanging the order of somatic_variant_annotation &
        # somatic_variant_filtration steps.
        # The import of the other workflow must be dependent on the config:
        # if in the somatic_variant_filtration config, has_annotation is True,
        # then the filtration step will include the annotation workflow as a
        # previous step.
        # THIS IMPLIES THAT NO FILTRATION OCCURED BEFORE ANNOTATION
        # This protects against circular import of workflows.
        #
        # This must be done before initialisation of the workflow.
        if config["step_config"][self.name]["has_annotation"]:
            from snappy_pipeline.workflows.any_variant_annotation import (
                AnyVariantAnnotationWorkflow,
            )

            previous_step = AnyVariantAnnotationWorkflow
        else:
            if self.variant_origin == VariantOrigin.SOMATIC:
                from snappy_pipeline.workflows.somatic_variant_calling import (
                    SomaticVariantCallingWorkflow,
                )

                previous_step = SomaticVariantCallingWorkflow
            elif self.variant_origin == VariantOrigin.GERMLINE:
                from snappy_pipeline.workflows.germline_variant_calling import (
                    GermlineVariantCallingWorkflow,
                )

                previous_step = GermlineVariantCallingWorkflow
            else:
                from snappy_pipeline.workflows.any_variant_calling import AnyVariantCallingWorkflow

                previous_step = AnyVariantCallingWorkflow
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=self.model_class,
            previous_steps=(previous_step, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                OneFilterDkfzStepPart,
                OneFilterBcftoolsStepPart,
                OneFilterRegionsStepPart,
                OneFilterProtectedStepPart,
                LastFilterStepPart,
                LinkOutStepPart,
            )
        )

        # Might be used by some filter...
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

        # Register sub-workflows and copy tools by variant origin
        if self.config.variant_origin == VariantOrigin.SOMATIC:
            if self.config.has_annotation:
                self.register_sub_workflow(
                    "somatic_variant_annotation", self.config.path_variant, "variant"
                )
            else:
                self.register_sub_workflow(
                    "somatic_variant_calling", self.config.path_variant, "variant"
                )
            if not self.config.tools_variant_calling:
                self.config.tools_variant_calling = self.w_config.step_config[
                    "somatic_variant_calling"
                ].tools
        elif self.config.variant_origin == VariantOrigin.GERMLINE:
            if self.config.has_annotation:
                self.register_sub_workflow(
                    "germline_variant_annotation", self.config.path_variant, "variant"
                )
            else:
                self.register_sub_workflow(
                    "germline_variant_calling", self.config.path_variant, "variant"
                )
            if not self.config.tools_variant_calling:
                self.config.tools_variant_calling = self.w_config.step_config[
                    "germline_variant_calling"
                ].tools
        else:
            if self.config.has_annotation:
                self.register_sub_workflow(
                    "any_variant_annotation", self.config.path_variant, "variant"
                )
            else:
                self.register_sub_workflow(
                    "any_variant_calling", self.config.path_variant, "variant"
                )
        assert self.config.tools_variant_calling, "No tool configured for variant calling"

        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"]["tools"]["dna"]

        if not self.config.has_annotation:
            self.config.tools_variant_annotation = []

        self.table = filter_table_by_modality(sample_sheets(self.sheets), modality="dna")
        assert self.table.shape[1] > 0, "No valid samples"

    @listify
    def _get_result_only_files(self):
        """Returns the vcf files (vcf, vcf.tbi & checksums)"""
        if self.config.has_annotation:
            name_pattern = "{mapper}.{caller}.{annotator}.filtered.{library}"
        else:
            name_pattern = "{mapper}.{caller}.filtered.{library}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            annotator=self.config.tools_variant_annotation,
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + ".full{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            annotator=self.config.tools_variant_annotation,
            ext=EXT_VALUES,
        )

    @listify
    def _get_log_only_files(self):
        """Returns the log files"""
        if self.config.has_annotation:
            name_pattern = "{mapper}.{caller}.{annotator}.filtered.{library}"
        else:
            name_pattern = "{mapper}.{caller}.filtered.{library}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            annotator=self.config.tools_variant_annotation,
            ext=LOG_EXT_VALUES,
        )

    def get_result_files(self):
        return self._get_result_only_files() + self._get_log_only_files()

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        match self.config.variant_origin:
            case VariantOrigin.SOMATIC:
                libraries = self.table[self.table["isTumor"]]
            case VariantOrigin.GERMLINE:
                libraries = self.table[~self.table["isTumor"]]
            case _:
                libraries = self.table
        for library in libraries["ngs_library"]:
            yield from expand(tpl, library=[library], **kwargs)
