# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_annotation`` step

The ``somatic_variant_annotation`` step takes as the input the results of the
``somatic_variant_calling`` step (bgzip-ed and indexed VCF files) and performs annotation of the
somatic variants.  The result are annotated versions of the somatic variant VCF files (again
bgzip-ed and indexed VCF files).

==========
Step Input
==========

The somatic variant annotation step uses Snakemake sub workflows for using the result of the
``somatic_variant_calling`` step. It can now also use the output from the ``somatic_variant_filtration``
step.

The main assumption is that each VCF file contains the two matched normal and tumor samples.

===========
Step Output
===========

Users can annotate all genes & transcripts overlapping with the variant locus, or
they can select one representative gene and transcript for annotation.
In the latter case, the output vcf file will only contain one annotation per variant, while
in the former case, there might be over 100 annotations for each variant.

The ordering of features driving the representative annotation choice is under user control.
The default order is:

1. ``biotype``: protein coding genes come first, it is unclear what is the order for other types of genes
2. ``mane``: the `MANE transcript <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ is selected before other transcripts
3. ``appris``: the `APPRIS principal isoform <https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii89/6701991>`_ is selected before alternates
4. ``tsl``: `Transcript Support Level <http://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html>`_ values in increasing order
5. ``ccds``: Transcripts with `CCDS <https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi>`_ ids are selected before those without
6. ``canonical``: ENSEMBL canonical transcripts are selected before the others
7. ``rank``: VEP internal ranking is used
8. ``length``: longer transcripts are preferred to shorter ones

This order is (hopefully) suitable for cBioPortal export, as well defined transcripts from protein-coding genes are selected when possible.
However, it is recommended to check the full annotation for variants in or nearby disease-relevant genes.

All annotators generate a vcf with one annotation per transcript, and some annotators
(only ENSEMBL's Variant Effect Predictor in the current implementation) can also produce another
output containing all annotations.
The single annotation vcf is named ``<mapper>.<caller>.<annotator>.vcf.gz`` and
the full annotation output is named ``<mapper>.<caller>.<annotator>.full.vcf.gz``

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_annotation.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import os
import sys
from collections import OrderedDict

from snakemake.io import expand

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS,
    SomaticVariantCallingWorkflow,
)

from .model import SomaticVariantAnnotation as SomaticVariantAnnotationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Names of the annotator tools
ANNOTATION_TOOLS = ("vep",)

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = SomaticVariantAnnotationConfigModel.default_config_yaml_string()


class AnnotateSomaticVcfStepPart(BaseStepPart):
    """Annotate VCF file from somatic calls

    .. note:

        The ``tumor_library`` wildcard can actually be the name of a donor!
    """

    #: Only creates vcf with one annotation per variant
    has_full = False

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancre sample
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

    @staticmethod
    def _name_template(
        config: SomaticVariantAnnotationConfigModel, annotator: str | None = None
    ) -> str:
        tpl = "{mapper}.{var_caller}"
        if annotator:
            tpl += f".{annotator}"
        if config.is_filtered:
            tpl += ".filtered.{tumor_library}"
        else:
            tpl += ".{tumor_library}"
        return tpl

    @dictify
    def get_input_files(self, action):
        """Return path to somatic vcf input file"""
        # Validate action
        self._validate_action(action)
        tpl = self._name_template(self.config)
        tpl = os.path.join("output", tpl, "out", tpl)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["somatic_variant"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        tpl = self._name_template(self.config, annotator=self.annotator)
        prefix = os.path.join("work", tpl, "out", tpl)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        if self.has_full:
            key_ext["full"] = ".full.vcf.gz"
            key_ext["full_tbi"] = ".full.vcf.gz.tbi"
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        """Return mapping of log files."""
        # Validate action
        self._validate_action(action)
        tpl = self._name_template(self.config, annotator=self.annotator)
        prefix = os.path.join("work", tpl, "log", tpl)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext


class VepAnnotateSomaticVcfStepPart(AnnotateSomaticVcfStepPart):
    """Annotate VCF file from somatic calling using ENSEMBL's VEP"""

    #: Step name
    name = "vep"

    #: Annotator name to construct output paths
    annotator = "vep"

    #: Class available actions
    actions = ("run",)

    #: Also creates vcf with all annotations
    has_full = True

    #: Allowed keywords for pick order
    PICK_ORDER = (
        "biotype",
        "mane_select",
        "mane_plus_clinical",
        "appris",
        "tsl",
        "ccds",
        "canonical",
        "rank",
        "length",
    )

    @dictify
    def get_input_files(self, action: str):
        input_files = super().get_input_files(action)
        for k, v in input_files.items():
            yield k, v
        yield "reference", self.w_config.static_data_config.reference.path

    def get_args(self, action):
        """Return arguments to pass down."""
        self._validate_action(action)
        return {"config": self.config.get(self.name).model_dump(by_alias=True)}

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.vep.num_threads,
            time="24:00:00",  # 24 hours
            memory=f"{16 * 1024 * 1}M",
        )


class SomaticVariantAnnotationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    name = "somatic_variant_annotation"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        # Ugly hack to allow exchanging the order of somatic_variant_annotation &
        # somatic_variant_filtration steps.
        # The import of the other workflow must be dependent on the config:
        # if in the somatic_variant_annotation config, the filtration_schema is not unfiltered,
        # then the annotation step will include the filtration workflow as a
        # previous step.
        # THIS IMPLIES THAT NO ANNOTATION OCCURED BEFORE FILTRATION
        # This protects against circular import of workflows.
        #
        # This must be done before initialisation of the workflow.
        previous_steps = [SomaticVariantCallingWorkflow, NgsMappingWorkflow]
        if config["step_config"]["somatic_variant_annotation"].get("is_filtered", False):
            from snappy_pipeline.workflows.somatic_variant_filtration import (
                SomaticVariantFiltrationWorkflow,
            )

            previous_steps.insert(0, SomaticVariantFiltrationWorkflow)
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticVariantAnnotationConfigModel,
            previous_steps=previous_steps,
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((VepAnnotateSomaticVcfStepPart, LinkOutStepPart))
        # Register sub workflows
        if self.config.is_filtered:
            self.register_sub_workflow(
                "somatic_variant_filtration", self.config.path_somatic_variant, "somatic_variant"
            )
        else:
            self.register_sub_workflow(
                "somatic_variant_calling", self.config.path_somatic_variant, "somatic_variant"
            )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.tools_somatic_variant_calling:
            self.config.tools_somatic_variant_calling = self.w_config.step_config[
                "somatic_variant_calling"
            ].tools

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        annotators = set(self.config.tools) & set(ANNOTATION_TOOLS)
        callers = set(self.config.tools_somatic_variant_calling)
        name_pattern = AnnotateSomaticVcfStepPart._name_template(
            self.config, annotator="{annotator}"
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            var_caller=callers & set(SOMATIC_VARIANT_CALLERS),
            annotator=annotators,
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            var_caller=callers & set(SOMATIC_VARIANT_CALLERS),
            annotator=annotators,
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )
        # Annotators with full output
        full = list(
            filter(
                lambda x: self.sub_steps[x].has_full,
                set(self.config.tools) & set(ANNOTATION_TOOLS),
            ),
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + ".full{ext}"),
            mapper=self.config.tools_ngs_mapping,
            var_caller=callers & set(SOMATIC_VARIANT_CALLERS),
            annotator=full,
            ext=EXT_VALUES,
        )
        # joint calling
        name_pattern = AnnotateSomaticVcfStepPart._name_template(
            self.config, annotator="{annotator}"
        )

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

    def _yield_result_files_joint(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the joint somatic variant callers such as
        "Bcftools joint".
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for donor in sheet.donors:
                yield from expand(tpl, tumor_library=[donor], **kwargs)
