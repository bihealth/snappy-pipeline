# -*- coding: utf-8 -*-
"""Implementation of the ``igv_session_generation`` step

This step takes as the input the output of the following steps and generates an IGV session XML
file that displays the results as genome tracks:

- ``ngs_mapping``
- ``variant_annotation`` or ``variant_calling``

==========
Step Input
==========

The IGV session generation step takes as the input of the following CUBI pipeline steps:

- ``ngs_mapping``
- ``variant_annotation`` or ``variant_calling``

===========
Step Output
===========

.. note: TODO

====================
Global Configuration
====================

.. note: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_igv_session_generation.rst

=======
Reports
=======

Currently, no reports are generated.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import os

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from jinja2 import Environment, FileSystemLoader
from snakemake import shell
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_phasing import VariantPhasingWorkflow
from .model import IgvSessionGeneration as IgvSessionGenerationConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".igv_session.xml", ".igv_session.xml.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("xml", "xml_md5")

#: Default configuration of the wgs_sv_filtration step
DEFAULT_CONFIG = IgvSessionGenerationConfigModel.default_config_yaml_string()


class WriteIgvSessionFileStepPart(BaseStepPart):
    """Write out the IGV session XML file."""

    #: Step name
    name = "write_igv_session_file"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.previous_step = self.parent.previous_step
        self.prev_token = self.parent.prev_token
        self.ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    self.ngs_library_to_pedigree[donor.dna_ngs_library.name] = pedigree

    def _get_path_bam(self, wildcards, donor):
        # TODO: This cannot be correct. For each set of donors it will return the same index bam.
        # TODO: For instance, given pedigree (P001, P002, P003) it will return three time the
        # TODO: same value: 'NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam'
        _ = donor
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        return ngs_mapping(
            "output/{mapper}.{index_library}/out/{mapper}.{index_library}.bam".format(**wildcards)
        )

    def _get_path_vcf(self, wildcards, real_index):
        prev_step = self.parent.sub_workflows[self.previous_step]
        name_pattern = "{mapper}.{caller}{prev_token}.{real_index_library}"
        input_path = ("output/" + name_pattern + "/out/" + name_pattern).format(
            prev_token=self.prev_token,
            real_index_library=real_index.dna_ngs_library.name,
            **wildcards
        )
        return prev_step(input_path + ".vcf.gz")

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            # Get name of real index, used when input is not variant_phasing
            pedigree = self.ngs_library_to_pedigree[wildcards.index_library]
            # BAM and BAI file of the offspring
            yield "bam", [self._get_path_bam(wildcards, donor) for donor in pedigree.donors]
            # Input file comes from previous step.
            yield "vcf", self._get_path_vcf(wildcards, pedigree.index)

        return input_function

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{index_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + "%s")
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, tpl % ext

    @staticmethod
    def render_from_template(directory, template_name, **kwargs):
        loader = FileSystemLoader(directory)
        env = Environment(loader=loader)
        template = env.get_template(template_name)
        return template.render(**kwargs)

    def run(self, wildcards, output):
        """Write out the IGV session file information"""
        # Get information for writing out the template.
        pedigree = self.ngs_library_to_pedigree[wildcards.index_library]
        template_vars = {
            "bams": [
                {
                    "label": donor.name,
                    "path": os.path.relpath(self._get_path_bam(wildcards, donor), output.xml),
                }
                for donor in pedigree.donors
            ],
            "vcfs": [
                {
                    "label": "Variant Calls",
                    "path": os.path.relpath(
                        self._get_path_vcf(wildcards, pedigree.index), output.xml
                    ),
                }
            ],
        }
        # Write output file.
        with open(output.xml, "wt") as session_file:
            xml = self.__class__.render_from_template(
                os.path.dirname(__file__), "igv_session.xml", **template_vars
            )
            print(xml, file=session_file)
        # Write MD5 sum file.
        shell(
            r"""
            pushd $(dirname {output.xml})
            md5sum $(basename {output.xml}) >$(basename {output.xml}).md5
            """.format(
                output=output
            )
        )


class IgvSessionGenerationWorkflow(BaseStep):
    """Perform IGV session generation"""

    #: Workflow name
    name = "igv_session_generation"

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
            config_model_class=IgvSessionGenerationConfigModel,
            previous_steps=(VariantPhasingWorkflow, VariantAnnotationWorkflow, NgsMappingWorkflow),
        )
        # Register sub workflows
        for prev in ("variant_phasing", "variant_annotation", "variant_calling"):
            if self.config["path_%s" % prev]:
                self.previous_step = prev
                self.register_sub_workflow(prev, self.config["path_%s" % prev])
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
        self.register_sub_step_classes((WriteIgvSessionFileStepPart, LinkOutStepPart))
        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.tools_variant_calling:
            self.config.tools_variant_calling = self.w_config.step_config["variant_calling"].tools

    @listify
    def get_result_files(self):
        """Return list of result files for the workflow."""
        # Hard-filtered results
        name_pattern = "{mapper}.{caller}%s.{index_library.name}" % (self.prev_token,)
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        # TODO: currently only families with children are supported
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if not donor.dna_ngs_library:
                        continue
                    elif not donor.father or not donor.father.dna_ngs_library:
                        continue
                    elif not donor.mother or not donor.mother.dna_ngs_library:
                        continue
                    else:
                        yield from expand(tpl, index_library=[donor.dna_ngs_library], **kwargs)
