# -*- coding: utf-8 -*-
"""Implementation of the ``roh_calling`` step

The roh_calling step allows for calling runs of homozygosity from germline VCF files.

Currently, it only employs ``bcftools roh`` with simple parameters but in the future more
involved approaches could be added.

==========
Stability
==========

The step itself is quite new.  Except for default parameter changes in the ROH calling, no large
changes are expected.

==========
Step Input
==========

The ``roh_calling`` step uses Snakemake sub workflows for using the result of the
``variant_calling`` step.

===========
Step Output
===========

Run-of-homozygosity will be performed on the input VCF files, on all samples at the same time.
The only tool for this currently implemented is ``bcftools roh``.  The result of this is one
gzip-compressed text files containing the assignment at the variant loci and the resulting
segmentation.

- ``{mapper}.{var_caller}.{roh_caller}.{ngs_library_name}.txt.gz``
- ``{mapper}.{var_caller}.{roh_caller}.{ngs_library_name}.txt.gz.tbi``

This file is then further processed to generate one (tabixed) BED file for each individual in the
pedigree listing the ROH variants.

- ``{mapper}.{var_caller}.{roh_caller}.{ngs_library_name}.vcf.gz``
- ``{mapper}.{var_caller}.{roh_caller}.{ngs_library_name}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{roh_caller}.{ngs_library_name}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{roh_caller}.{ngs_library_name}.vcf.gz.tbi.md5``

For example, it could look like this:

::

    output/
    +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1
    |   `-- out
    │       +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz
    │       +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.md5
    │       +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.tbi
    │       +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1.bed.gz.tbi.md5
    │       +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz
    │       +-- bwa.gatk_ug.bcftools_roh.P001-N1-DNA1-WGS1.txt.gz.md5
    │       +-- bwa.gatk_ug.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz
    │       +-- bwa.gatk_ug.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.md5
    │       +-- bwa.gatk_ug.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi
    │       +-- bwa.gatk_ug.bcftools_roh.P002-N1-DNA1-WGS1.bed.gz.tbi.md5
    │       +-- bwa.gatk_ug.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz
    │       +-- bwa.gatk_ug.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.md5
    │       +-- bwa.gatk_ug.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.tbi
    │       `-- bwa.gatk_ug.bcftools_roh.P003-N1-DNA1-WGS1.bed.gz.tbi.md5
    [...]

====================
Global Configuration
====================

No global configuration is used.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_roh_calling.rst

=====================
Available ROH Callers
=====================

At the moment, only ``bcftools_roh`` is available as a caller.

=======
Reports
=======

No reports are generated.
"""

from collections import OrderedDict
import os  # noqa: F401
import sys

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
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the roh_calling step
DEFAULT_CONFIG = r"""
# Default configuration roh_calling
step_config:
  roh_calling:
    path_variant_calling: ../variant_calling
    tools_ngs_mapping:
      - bwa
    tools_variant_calling:
      - gatk_hc
    tools: [bcftools_roh] # REQUIRED, available: 'bcftools_roh'.
    bcftools_roh:
      path_targets: null
      path_af_file: null
      ignore_homref: false
      skip_indels: false
      rec_rate: 1e-8
"""


#: Key => extension mapping for BED file.
BED_EXTENSIONS = {
    "bed": ".bed.gz",
    "vcf_tbi": ".bed.gz.tbi",
    "bed_md5": ".bed.gz.md5",
    "tbi_vcf_md5": ".bed.gz.tbi.md5",
}


class BcftoolsRohStepPart(BaseStepPart):
    """Perform ROH calling using ``bcftools roh``."""

    name = "bcftools_roh"

    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{name}.{{library_name}}/out/{{mapper}}.{name}.{{library_name}}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return input function for bcftools_roh rule"""
        assert action in self.actions
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path to input VCF file"""
        tpl = (
            "output/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["variant_calling"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl.format(**wildcards) + ext)

    def get_output_files(self, action):
        """Return output function for bcftools_roh rule"""
        assert action in self.actions
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_run(self):
        """Return output files for ROH calling"""
        path = (
            "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}"
        )
        yield "txt", path + ".regions.txt.gz"
        yield "txt_md5", path + ".regions.txt.gz.md5"

    def get_log_file(self, action):
        assert action in self.actions
        return (
            "work/{mapper}.{var_caller}.bcftools_roh.{index_ngs_library}/"
            "log/snakemake.bcftools_roh.log"
        )

    def get_resource_usage(self, action):
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
        return ResourceUsage(
            threads=1,
            time="00:10:00",
            memory="4G",
        )


class RohCallingWorkflow(BaseStep):
    """Perform run-of-homozygosity calling workflow"""

    name = "roh_calling"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            (VariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((BcftoolsRohStepPart, LinkOutStepPart))
        # Register sub workflows
        self.register_sub_workflow("variant_calling", self.config["path_variant_calling"])
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
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        yield from self._yield_result_files(
            (
                "output/{mapper}.{caller}.bcftools_roh.{index_library}/out/"
                "{mapper}.{caller}.bcftools_roh.{donor_library}{ext}"
            ),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:  # pragma: no cover
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                elif not pedigree.index.dna_ngs_library:  # pragma: no cover
                    msg = "INFO: pedigree index without DNA NGS library (names: {})"
                    print(
                        msg.format(  # pragma: no cover
                            list(sorted(d.name for d in pedigree.donors))
                        ),
                        file=sys.stderr,
                    )
                    continue  # pragma: no cover
                # Text file only for the index
                yield from expand(
                    tpl,
                    ext=[".regions.txt.gz", ".regions.txt.gz.md5"],
                    donor_library=[pedigree.index.dna_ngs_library.name],
                    index_library=[pedigree.index.dna_ngs_library.name],
                    **kwargs,
                )
