# -*- coding: utf-8 -*-
"""Generic variant annotation step

The assumption is that variant annotation is similar for all kind of variants,
germline or somatic.

This implementation attepmts to provide annotation for a set of variants stored in vcf
format, using a generic sample sheet (_i.e._ without restriction on tumor/normal pairs
or the availability of pedigree).

Unfortunately, I haven't found a way to design a completely generic implementation.
This is because during the step initialisation, the pipeline validates the configuration
of previous steps, and there is no common configuration of the germline & somatic variant
calling steps. So there is no configration possible of a generic variant calling step.
Therefore, which variant calling step has been run must be known at the time of the
initialisation of the variant annotation step. So a workaround has been implemented,
and the user must set a ``variant_origin`` parameter in order to define which type of
variant calling step has been run. It is a poor solution, because it prevents the
annotation of both somatic & germline variants in the same pipelie run, until the
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

If the vcf has been filtered before annotation, the naming pattern becomes
``{mapper}.{caller}.filtered.{library}`` for the folder and the file name.

===========
Step Output
===========

Annotations depend on the choice of transcript overlapping with the variant.
Annotators can either select a single transcript for their annotation, or annotate
all transcripts overlapping with the variant.
When a tool can do both, two output vcf files are produced, the single annotation
named ``<mapper>.<caller>.<annotator>.vcf.gz`` and the full annotation named
``<mapper>.<caller>.<annotator>.full.vcf.gz``.
Some tools put the single annotation feature selection under user control.

---
VEP
---

The default order for single annotation feature selection is:

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

VEP plugins can be used, the user can either specify the path to the plugin (when already stored),
or a URL to download it. In the latter case, the downloaded plugin will be in
``work/vep_plugins/out{plugin.name}.pm``, in the former case, a symlink will be created in that
directory (vep requires that all plugins are stored in the same directory).
The legacy argument ``plugins_dir`` is ignored.

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

from snakemake.io import expand

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.samplesheet import filter_table_by_modality, sample_sheets
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from .model import AnyVariantAnnotation as AnyVariantAnnotationConfigModel
from snappy_pipeline.models.annotation import Vep as VepConfigModel

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
)

#: Names of the annotator tools
ANNOTATION_TOOLS = ("vep",)

#: Default configuration for the variant_calling step
DEFAULT_CONFIG = AnyVariantAnnotationConfigModel.default_config_yaml_string()


class AnnotateVcfStepPart(BaseStepPart):
    """Annotate VCF file from somatic or germline or any calls"""

    #: Only creates vcf with one annotation per variant
    has_full = False

    @staticmethod
    def _name_template(
        config: AnyVariantAnnotationConfigModel, annotator: str | None = None
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
        """Return path to vcf input file"""
        # Validate action
        self._validate_action(action)
        tpl = self._name_template(self.config)
        tpl = os.path.join("output", tpl, "out", tpl)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["variant"]
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
    def get_log_file(self, action):
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
            yield key + "_md5", prefix + ext + ".md5"


class VepAnnotateVcfStepPart(AnnotateVcfStepPart):
    """Annotate VCF file from variant calling using ENSEMBL's VEP"""

    #: Step name
    name = "vep"

    #: Annotator name to construct output paths
    annotator = "vep"

    #: Class available actions
    actions = ("plugins", "run")

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
        cfg: VepConfigModel = self.config.get("vep")
        yield "plugins", [f"work/vep_plugins/out/{plugin.name}.pm" for plugin in cfg.plugins]

    def get_output_files(self, action: str):
        self._validate_action(action)
        if action == "plugins":
            cfg: VepConfigModel = self.config.get("vep")
            return [f"work/vep_plugins/out/{plugin.name}.pm" for plugin in cfg.plugins]
        return super().get_output_files(action)

    def get_log_file(self, action: str):
        self._validate_action(action)
        return super().get_log_file(action)

    def get_args(self, action):
        """Return arguments to pass down."""
        self._validate_action(action)
        if action == "plugins":
            return self.config.get(self.name).get("plugins", [])
        vep_config = dict(self.config.get(self.name).model_dump(by_alias=True))
        vep_config["plugins"] = [plugin["name"] for plugin in vep_config["plugins"]]
        vep_config["plugins_dir"] = "work/vep_plugins/out"
        return {"config": vep_config}

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


class AnyVariantAnnotationWorkflow(BaseStep):
    """Perform variant annotation"""

    name = "any_variant_annotation"
    variant_origin = VariantOrigin.ANY
    model_class = AnyVariantAnnotationConfigModel

    sheet_shortcut_class = GenericSampleSheet

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
        if config["step_config"][self.name]["is_filtered"]:
            from snappy_pipeline.workflows.any_variant_filtration import (
                AnyVariantFiltrationWorkflow,
            )

            previous_step = AnyVariantFiltrationWorkflow
        else:
            if self.variant_origin == VariantOrigin.SOMATIC:
                from snappy_pipeline.workflows.somatic_variant_calling import (
                    SomaticVariantCallingWorkflow,
                )

                previous_step = SomaticVariantCallingWorkflow
            elif config["step_config"][self.name].get("variant_origin", VariantOrigin.GERMLINE):
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
            previous_steps=(previous_step,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((VepAnnotateVcfStepPart, LinkOutStepPart))

        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"]["tools"]["dna"]

        # Register sub-workflows and copy tools by variant origin
        if self.config.variant_origin == VariantOrigin.SOMATIC:
            if self.config.is_filtered:
                self.register_sub_workflow(
                    "somatic_variant_filtration", self.config.path_variant, "variant"
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
            if self.config.is_filtered:
                self.register_sub_workflow(
                    "germline_variant_filtration", self.config.path_variant, "variant"
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
            if self.config.is_filtered:
                self.register_sub_workflow(
                    "any_variant_filtration", self.config.path_variant, "variant"
                )
            else:
                self.register_sub_workflow(
                    "any_variant_calling", self.config.path_variant, "variant"
                )
        assert self.config.tools_variant_calling, "No tool configured for variant calling"

        self.table = filter_table_by_modality(sample_sheets(self.sheets), modality="dna")
        assert self.table.shape[1] > 0, "No valid samples"

    @listify
    def _get_result_only_files(self):
        """Returns the vcf files (vcf, vcf.tbi & checksums)"""
        if self.config.is_filtered:
            name_pattern = "{mapper}.{caller}.{annotator}.filtered.{library}"
        else:
            name_pattern = "{mapper}.{caller}.{annotator}.{library}"
        for tool in self.config.tools:
            yield from self._yield_result_files(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=self.config.tools_ngs_mapping,
                caller=self.config.tools_variant_calling,
                annotator=[tool],
                ext=EXT_VALUES,
            )
            if self.sub_steps[tool].has_full:
                yield from self._yield_result_files(
                    os.path.join("output", name_pattern, "out", name_pattern + ".full{ext}"),
                    mapper=self.config.tools_ngs_mapping,
                    caller=self.config.tools_variant_calling,
                    annotator=[tool],
                    ext=EXT_VALUES,
                )

    @listify
    def _get_log_only_files(self):
        """Returns the log files"""
        if self.config.is_filtered:
            name_pattern = "{mapper}.{caller}.{annotator}.filtered.{library}"
        else:
            name_pattern = "{mapper}.{caller}.{annotator}.{library}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            annotator=self.config.tools,
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
