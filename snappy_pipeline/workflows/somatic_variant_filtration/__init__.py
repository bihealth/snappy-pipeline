# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_filtration`` step

The current implementation supports two filtration schema:

- the *legacy* schema, now deprecated, always runs the `DKFZBiasFilter <https://github.com/DKFZ-ODCF/DKFZBiasFilter>`_ &
  `EBFilter <https://doi.org/10.1093/nar/gkt126>`_, and produces files for all combinations of available filters.
- the *new* schema focuses on flexibility, allows any combination of filters, and returns a single fitlered file for each sample.

The *new* schema is used when the configuration option ``filter_list`` is not empty.
The following document describes only this *new* schema.

==========
Step Input
==========

The step requires ``vcf`` files from either the ``somatic_variant_calling`` or ``somatic_variant_annotation`` steps.
In the former case, the configuration option ``has_annotation`` must be set to ``False``.

In both cases, it will use the regular output ``vcf`` file, not ``*.full.vcf.gz``.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name`` and each read mapper
``mapper`` that the library has been aligned with, and the variant caller ``var_caller``, the
pipeline step will create a directory ``output/{mapper}.{var_caller}.{annotator}.{lib_name}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

Two ``vcf`` files are produced:

- ``{mapper}.{var_caller}.{annotator}.{lib_name}.vcf.gz`` which contains only the variants that have passed all filters, or that were protected, and
- ``{mapper}.{var_caller}.{annotator}.{lib_name}.full.vcf.gz`` which contains all variants, with the reason for rejection in the ``FILTER`` column.

When the ``somatic_variant_annotation`` step has been omitted, and the filtration is done directly from the output of the ``somatic_variant_calling`` step,
then the output files are stored in the ``output/{mapper}.{var_caller}.{lib_name}/out`` directory, under the names ``{mapper}.{var_caller}.{lib_name}.vcf.gz`` &
``{mapper}.{var_caller}.{lib_name}.full.vcf.gz``

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
  - ebfilter:
    ebfilter_threshold: 2.4
  - bcftools:
    exclude: "AD[1:0]+AD[1:1]<50 | AD[1:1]<5 | AD[1:1]/(AD[1:0]+AD[1:1])<0.05"
  - bcftools:
    exclude: "((REF='C' & ALT='T') | (REF='G' & ALT='A')) & AD[1:1]/(AD[1:0]+AD[1:1])<0.10"
  - protected:
    path_bed: hotspots_locii.bed

This list of filters would apply the DKFZBiasFilter, the EBFilter, reject all variants with depth lower than 50, less than 5 reads supporting the alternative allele, or with a variant allele fraction below 5%.
It would also reject all C-to-T and G-to-A variants with a VAF lower than 10%, because they might be FFPE artifacts.
All variants overlapping with hotspots locii would be protected against filtration.

Note that the parallelisation of ``ebfilter`` has been removed, even though this operation can be slow when there are many variants (from WGS data for example).
"""

from collections import OrderedDict
import os
import random
import sys

from snakemake.io import expand

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.somatic_variant_annotation import (
    ANNOTATION_TOOLS,
    SomaticVariantAnnotationWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)

from .model import SomaticVariantFiltration as SomaticVariantFiltrationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = SomaticVariantFiltrationConfigModel.default_config_yaml_string()


class SomaticVariantFiltrationStepPart(BaseStepPart):
    """Shared code for all tools in somatic_variant_filtration"""

    def __init__(self, parent):
        super().__init__(parent)
        self.config = parent.config
        self.name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            self.name_pattern += ".{annotator}"
        self.name_pattern += ".{tumor_library}"
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor
        # Build mapping from tumor library to normal library
        self.tumor_to_normal_library = OrderedDict()
        for tumor_library, normal_sample in self.tumor_ngs_library_to_sample_pair.items():
            for test_sample in normal_sample.normal_sample.bio_sample.test_samples.values():
                for ngs_library in test_sample.ngs_libraries.values():
                    if tumor_library not in self.tumor_to_normal_library:
                        self.tumor_to_normal_library[tumor_library] = (
                            test_sample.name + "-" + ngs_library.secondary_id
                        )

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair.get(wildcards.tumor_library, None)
        return pair.normal_sample.dna_ngs_library.name if pair else None

    def get_params(self, action):
        """Return arguments to pass down."""
        _ = action

        def params_function(wildcards):
            if wildcards.tumor_library not in self.donors:
                return {
                    "tumor_library": wildcards.tumor_library,
                    "normal_library": self.get_normal_lib_name(wildcards),
                }
            else:
                return {}

        return params_function


class OneFilterStepPart(SomaticVariantFiltrationStepPart):
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

        @dictify
        def input_function(wildcards):
            filter_nb = int(wildcards["filter_nb"])
            if filter_nb > 1:
                prev = list(self.config.filter_list[filter_nb - 2].keys())[0]
                n = filter_nb - 1
                yield (
                    "vcf",
                    os.path.join(
                        "work", self.name_pattern, "out", self.name_pattern + f".{prev}_{n}.vcf.gz"
                    ),
                )
            else:
                yield (
                    "vcf",
                    os.path.join(
                        self.config.path_somatic_variant,
                        "output",
                        self.name_pattern,
                        "out",
                        self.name_pattern + ".vcf.gz",
                    ),
                )

        return input_function

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

    def get_params(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            return {"filter_name": "{}_{}".format(self.filter_name, wildcards["filter_nb"])}

        return input_function


class OneFilterWithBamStepPart(OneFilterStepPart):
    def get_input_files(self, action):
        """Return path to input or previous filter vcf file & normal/tumor bams"""
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            parent = super(OneFilterWithBamStepPart, self).get_input_files(action)
            yield from parent(wildcards).items()

            yield (
                "bam",
                os.path.join(
                    self.config.path_ngs_mapping,
                    "output",
                    "{mapper}.{tumor_library}",
                    "out",
                    "{mapper}.{tumor_library}.bam",
                ),
            )
            if normal_library := self.tumor_to_normal_library.get(wildcards["tumor_library"], None):
                yield (
                    "normal",
                    os.path.join(
                        self.config.path_ngs_mapping,
                        "output",
                        f"{{mapper}}.{normal_library}",
                        "out",
                        f"{{mapper}}.{normal_library}.bam",
                    ),
                )

        return input_function


class OneFilterDkfzStepPart(OneFilterWithBamStepPart):
    name = "one_dkfz"
    filter_name = "dkfz"
    resource_usage = {"run": ResourceUsage(threads=1, time="12:00:00", memory=f"{3 * 1024}M")}


class OneFilterEbfilterStepPart(OneFilterWithBamStepPart):
    name = "one_ebfilter"
    filter_name = "ebfilter"
    resource_usage = {"run": ResourceUsage(threads=1, time="24:00:00", memory=f"{2 * 1024}M")}

    def get_input_files(self, action):
        """Return path to input or previous filter vcf file & normal/tumor bams"""
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            parent = super(OneFilterEbfilterStepPart, self).get_input_files(action)
            yield from parent(wildcards).items()

            yield "txt", self._get_output_files_write_panel()["txt"].format(**wildcards)

        return input_function

    @dictify
    def _get_output_files_write_panel(self):
        yield (
            "txt",
            (
                "work/{mapper}.eb_filter.panel_of_normals/out/{mapper}.eb_filter."
                "panel_of_normals.txt"
            ),
        )

    def get_params(self, action):
        """Return add EBFilter parameters to parameters"""
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            parent = super(OneFilterEbfilterStepPart, self).get_params(action)
            parameters = parent(wildcards)
            filter_nb = int(wildcards["filter_nb"])
            ebfilter_config = self.config.filter_list[filter_nb - 1][self.filter_name]
            parameters.update(ebfilter_config)
            parameters["has_annotation"] = self.config.has_annotation
            return parameters

        return input_function


class OneFilterBcftoolsStepPart(OneFilterStepPart):
    name = "one_bcftools"
    filter_name = "bcftools"

    def get_params(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            parent = super(OneFilterBcftoolsStepPart, self).get_params(action)
            parameters = parent(wildcards)
            filter_nb = int(wildcards["filter_nb"])
            filter = self.config.filter_list[filter_nb - 1][self.filter_name]
            keywords = filter.keywords()
            parameters.update(keywords)
            return parameters

        return input_function


class OneFilterRegionsStepPart(OneFilterStepPart):
    name = "one_regions"
    filter_name = "regions"

    def get_params(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            parent = super(OneFilterRegionsStepPart, self).get_params(action)
            parameters = parent(wildcards)
            filter_nb = int(wildcards["filter_nb"])
            filter = self.config.filter_list[filter_nb - 1][self.filter_name]
            keywords = filter.keywords()
            parameters.update(keywords)
            return parameters

        return input_function


class OneFilterProtectedStepPart(OneFilterStepPart):
    name = "one_protected"
    filter_name = "protected"

    def get_params(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            parent = super(OneFilterProtectedStepPart, self).get_params(action)
            parameters = parent(wildcards)
            filter_nb = int(wildcards["filter_nb"])
            filter = self.config.filter_list[filter_nb - 1][self.filter_name]
            keywords = filter.keywords()
            parameters.update(keywords)
            return parameters

        return input_function


class LastFilterStepPart(SomaticVariantFiltrationStepPart):
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


class DkfzBiasFilterStepPart(SomaticVariantFiltrationStepPart):
    """Flag variants with the DKFZ bias filter"""

    #: Step name
    name = "dkfz_bias_filter"

    #: Class available actions
    actions = ("run",)

    @dictify
    def get_input_files(self, action):
        """Return path to jannovar-annotated vcf input file"""
        # Validate action
        self._validate_action(action)

        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"

        # VCF file and index
        tpl = f"output/{name_pattern}." "{tumor_library}/" f"out/{name_pattern}." "{tumor_library}"
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        somatic_variant = self.parent.sub_workflows["somatic_variant"]
        for key, ext in key_ext.items():
            yield key, somatic_variant(tpl + ext)
        # BAM file and index
        tpl = "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}"
        key_ext = {"bam": ".bam", "bai": ".bam.bai"}
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for key, ext in key_ext.items():
            yield key, ngs_mapping(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)

        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"

        prefix = (
            rf"work/{name_pattern}.{self.name}."
            r"{tumor_library,[^\.]+}/"
            rf"out/{name_pattern}.{self.name}."
            r"{tumor_library}"
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
    def _get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)

        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        name_pattern += f".{self.name}"
        prefix = os.path.join(
            "work",
            name_pattern + r".{tumor_library,[^\.]+}",
            "log",
            name_pattern + ".{tumor_library}",
        )

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

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
            time="3-00:00:00",  # 3 days
            memory=f"{3 * 1024}M",
        )


class EbFilterStepPart(SomaticVariantFiltrationStepPart):
    """Flag variants with EBFilter"""

    #: Step name
    name = "eb_filter"

    #: Class available actions
    actions = ("run", "write_panel")

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_run(self, wildcards):
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        # VCF file and index
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        tpl = (
            f"work/{name_pattern}."
            "dkfz_bias_filter.{tumor_library}/"
            f"out/{name_pattern}."
            "dkfz_bias_filter."
            "{tumor_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, tpl.format(**wildcards) + ext
        # BAM file and index
        tpl = r"output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}"
        key_ext = {"bam": ".bam", "bai": ".bam.bai"}
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for key, ext in key_ext.items():
            yield key, ngs_mapping(tpl.format(**wildcards) + ext)
        # Panel of normals TXT file
        yield "txt", self._get_output_files_write_panel()["txt"].format(**wildcards)

    def _get_input_files_write_panel(self, wildcards):
        bam_paths = self._get_panel_of_normal_bams(wildcards)
        return {"bam": bam_paths, "bai": [p + ".bai" for p in bam_paths]}

    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    def get_params(self, action):
        """Return EBFilter parameters from the config"""
        # Validate action
        self._validate_action(action)
        parameters = self.config.eb_filter
        parameters.update(self.config.filter_sets.dkfz_and_ebfilter)
        parameters["has_annotation"] = self.config.has_annotation
        return parameters

    @dictify
    def _get_output_files_run(self):
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        # VCF file and index
        prefix = (
            rf"work/{name_pattern}."
            r"dkfz_bias_filter.eb_filter.{tumor_library,[^\.]+}/"
            rf"out/{name_pattern}."
            r"dkfz_bias_filter.eb_filter.{tumor_library}"
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
    def _get_output_files_write_panel(self):
        yield (
            "txt",
            (
                "work/{mapper}.eb_filter.panel_of_normals/out/{mapper}.eb_filter."
                "panel_of_normals.txt"
            ),
        )

    @dictify
    def _get_log_file(self, action):
        """Return path to log file for the given action"""
        # Validate action
        self._validate_action(action)

        if action == "write_panel":
            # No log returned by EbFilterStepPart.write_panel_of_normals_file
            return {}
        else:
            name_pattern = "{mapper}.{var_caller}"
            if self.config.has_annotation:
                name_pattern += ".{annotator}"
            name_pattern += ".dkfz_bias_filter.eb_filter"
            prefix = os.path.join(
                "work",
                name_pattern + r".{tumor_library,[^\.]+}",
                "log",
                name_pattern + ".{tumor_library}",
            )
            key_ext = (
                ("log", ".log"),
                ("conda_info", ".conda_info.txt"),
                ("conda_list", ".conda_list.txt"),
            )
            for key, ext in key_ext:
                yield key, prefix + ext

    def write_panel_of_normals_file(self, wildcards):
        """Write out file with paths to panels-of-normal"""
        output_path = self.get_output_files("write_panel")["txt"].format(**wildcards)
        with open(output_path, "wt") as outf:
            for bam_path in self._get_panel_of_normal_bams(wildcards):
                print(bam_path, file=outf)

    @listify
    def _get_panel_of_normal_bams(self, wildcards):
        """Return list of "panel of normal" BAM files."""
        libraries = []
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    if not bio_sample.extra_infos["isTumor"]:
                        libraries.append(bio_sample.dna_ngs_library.name)
        libraries.sort()
        random.seed(self.config.eb_filter.shuffle_seed)
        lib_count = self.config.eb_filter.panel_of_normals_size
        random.shuffle(libraries)
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}"
        for library in libraries[:lib_count]:
            yield ngs_mapping(tpl.format(normal_library=library, **wildcards) + ".bam")

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
            time="04:00:00",  # 4 hours
            memory=f"{8 * 1024}M",
        )


class ApplyFiltersStepPart(SomaticVariantFiltrationStepPart):
    """Apply the configured filters."""

    #: Step name
    name = "apply_filters"

    #: Class available actions
    actions = ("run",)

    #: Default filtration resource usage (should be light)
    resource_usage = {"run": ResourceUsage(threads=1, time="02:00:00", memory=f"{8 * 1024}M")}

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        name_pattern += ".dkfz_bias_filter.eb_filter"
        self.base_path_out = os.path.join(
            "work",
            name_pattern + r".{tumor_library,[^\.]+}.{filter_set,[^\.]+}",
            "out",
            name_pattern + ".{tumor_library}.{filter_set}" + "{ext}",
        )
        self.path_log = os.path.join(
            "work",
            name_pattern + r".{tumor_library,[^\.]+}.{filter_set,[^\.]+}",
            "log",
            name_pattern + ".{tumor_library}.{filter_set}" + ".log",
        )

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        def args_function(wildcards):
            result = {
                "normal_sample": self.get_normal_lib_name(wildcards),
                "tumor_sample": wildcards.tumor_library,
            }
            return result

        return args_function

    @dictify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        name_pattern += ".dkfz_bias_filter.eb_filter.{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, tpl + ext

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{ext}", ext)

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.path_log


class FilterToExonsStepPart(SomaticVariantFiltrationStepPart):
    """Apply the configured filters."""

    #: Step name
    name = "filter_to_exons"

    #: Class available actions
    actions = ("run",)

    #: Default filtration resource usage (should be light)
    resource_usage = {"run": ResourceUsage(threads=1, time="02:00:00", memory=f"{8 * 1024}M")}

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = "{mapper}.{var_caller}"
        if self.config.has_annotation:
            name_pattern += ".{annotator}"
        name_pattern += ".dkfz_bias_filter.eb_filter"
        self.base_path_out = os.path.join(
            "work",
            name_pattern + r".{tumor_library,[^\.]+}.{filter_set,[^\.]+}.{exon_list}",
            "out",
            name_pattern + ".{tumor_library}.{filter_set}.{exon_list}" + "{ext}",
        )
        self.path_log = os.path.join(
            "work",
            name_pattern + r".{tumor_library,[^\.]+}.{filter_set,[^\.]+}.{exon_list}",
            "log",
            name_pattern + ".{tumor_library}.{filter_set}.{exon_list}" + ".log",
        )
        self.base_path_in_ = os.path.join(
            "work",
            name_pattern + ".{tumor_library}.{filter_set}",
            "out",
            name_pattern + ".{tumor_library}.{filter_set}" + "{ext}",
        )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                yield (
                    key,
                    self.base_path_in_.format(
                        tumor_library=wildcards.tumor_library,
                        mapper=wildcards.mapper,
                        var_caller=wildcards.var_caller,
                        annotator=wildcards.get("annotator", ""),
                        filter_set=wildcards.filter_set,
                        exon_list=wildcards.exon_list,
                        ext=ext,
                    ),
                )

        return input_function

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{ext}", ext)

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.path_log


class SomaticVariantFiltrationWorkflow(BaseStep):
    """Perform somatic variant filtration"""

    #: Workflow name
    name = "somatic_variant_filtration"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=False)
    }

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
            config_model_class=SomaticVariantFiltrationConfigModel,
            previous_steps=(
                SomaticVariantAnnotationWorkflow,
                SomaticVariantCallingWorkflow,
                NgsMappingWorkflow,
            ),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                DkfzBiasFilterStepPart,
                EbFilterStepPart,
                OneFilterDkfzStepPart,
                OneFilterEbfilterStepPart,
                OneFilterBcftoolsStepPart,
                OneFilterRegionsStepPart,
                OneFilterProtectedStepPart,
                LastFilterStepPart,
                ApplyFiltersStepPart,
                FilterToExonsStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow(
            (
                "somatic_variant_annotation"
                if self.config.has_annotation
                else "somatic_variant_calling"
            ),
            self.config.path_somatic_variant,
            "somatic_variant",
        )
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.tools_somatic_variant_calling:
            self.config.tools_somatic_variant_calling = self.w_config.step_config[
                "somatic_variant_calling"
            ].tools
        if not self.config.tools_somatic_variant_annotation:
            self.config.tools_somatic_variant_annotation = self.w_config.step_config[
                "somatic_variant_annotation"
            ].tools

    @listify
    def get_result_files(self):
        """Return list of result files
        Process all primary DNA libraries and perform pairwise calling for tumor/normal pairs
        """
        mappers = set(self.config.tools_ngs_mapping) & set(
            self.w_config.step_config["ngs_mapping"].tools.dna
        )
        callers = set(self.config.tools_somatic_variant_calling) & set(
            SOMATIC_VARIANT_CALLERS_MATCHED
        )
        if self.config.has_annotation:
            annotators = set(self.config.tools_somatic_variant_annotation) & set(ANNOTATION_TOOLS)
        else:
            annotators = []

        log_ext = [e + m for e in ("log", "conda_list.txt", "conda_info.txt") for m in ("", ".md5")]

        if self.config.filtration_schema == "list":
            name_pattern = "{mapper}.{caller}"
            if self.config.has_annotation:
                name_pattern += ".{annotator}"
            name_pattern += ".filtered.{tumor_library}"

            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                ext=[f + e for f in ("", ".full") for e in EXT_VALUES],
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + ".{ext}"),
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                ext=log_ext,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + ".merged.tar.gz{ext}"),
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                ext=("", ".md5"),
            )
        else:
            filter_sets = {"no_filter"}
            if self.config.filter_sets:
                filter_sets |= self.config.filter_sets.keys()
            exon_lists = {"genome_wide"}
            if self.config.exon_lists:
                exon_lists |= self.config.exon_lists.keys()

            name_pattern = "{mapper}.{caller}"
            if self.config.has_annotation:
                name_pattern += ".{annotator}"
            name_pattern += ".dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}.{exon_list}"

            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                filter_set=filter_sets,
                exon_list=exon_lists,
                ext=EXT_VALUES,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + ".{ext}"),
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                filter_set=filter_sets,
                exon_list=exon_lists,
                ext=("log",),
            )

        # TODO: filtration for joint calling not implemented yet

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for bio_entity in sheet.sheet.bio_entities.values():
                for bio_sample in bio_entity.bio_samples.values():
                    if not bio_sample.extra_infos.get("isTumor", False):
                        continue
                    for test_sample in bio_sample.test_samples.values():
                        extraction_type = test_sample.extra_infos.get("extractionType", "unknown")
                        if extraction_type.lower() != "dna":
                            if extraction_type == "unknown":
                                msg = "INFO: sample {} has missing extraction type, ignored"
                                print(msg.format(test_sample.name), file=sys.stderr)
                            continue
                        for ngs_library in test_sample.ngs_libraries.values():
                            yield from expand(tpl, tumor_library=[ngs_library.name], **kwargs)
