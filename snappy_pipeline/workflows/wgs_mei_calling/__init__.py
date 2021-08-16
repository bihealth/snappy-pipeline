# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_mei_calling`` step

The wgs_mei_calling step allows for the calling of mobile element insertions from WGS data.

==========
Step Input
==========

The MEI calling step uses Snakemake sub workflows for using the result of the ``ngs_mapping``
step.

===========
Step Output
===========

For each read mapper, MEI calling tool, and pedigree a directory
``output/{mapper}.{caller}.{index_dna_library}/out`` will be created.  This directory will contain
the following links to files with the MEI calling results in VCF format and supporting files.

- ``{mapper}.{caller}.{index_dna_library}.vcf.gz``
- ``{mapper}.{caller}.{index_dna_library}.vcf.gz.md5``
- ``{mapper}.{caller}.{index_dna_library}.vcf.gz.tbi``
- ``{mapper}.{caller}.{index_dna_library}.vcf.gz.tbi.md5``

At the moment, the files contain the calls of the whole study, with the samples for each pedigree
as the first entries.

.. note:

    In the future, the ME insertion events not observed in the pedigree will probably be removed
    and an additional "whole patient cohort" VCF file will be provided.

For example, the output might look as follows

::

    output/
    |-- bwa.melt.151388-DNA1-DNA1-WGS1-000004
    |   `-- out
    |       |-- bwa.melt.151388-DNA1-DNA1-WGS1-000004.vcf.gz
    |       |-- bwa.melt.151388-DNA1-DNA1-WGS1-000004.vcf.gz.md5
    |       |-- bwa.melt.151388-DNA1-DNA1-WGS1-000004.vcf.gz.tbi
    |       `-- bwa.melt.151388-DNA1-DNA1-WGS1-000004.vcf.gz.tbi.md5
    [...]

====================
Global Configuration
====================

- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_wgs_mei_calling.rst

===========
Limitations
===========

Currently only works with ``hs37d5`` reference
"""

from collections import OrderedDict
import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand, touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Available somatic variant callers
MEI_CALLERS = ("melt",)

#: Actions used for Melt
MELT_ACTIONS = (
    "preprocess",
    "indiv_analysis",
    "group_analysis",
    "genotype",
    "make_vcf",
    "merge_vcf",
    "reorder_vcf",
)

#: Default configuration for the wgs_mei_calling step
DEFAULT_CONFIG = r"""
step_config:
  wgs_mei_calling:
    tools:
    - melt
    path_ngs_mapping: ../ngs_mapping
    melt:
      me_refs_infix: 1KGP_Hg19
      me_types:
      - ALU
      - LINE1
      - SVA
      genes_file: add_bed_files/1KGP_Hg19/hg19.genes.bed  # adjust, e.g., Hg38/Hg38.genes.bed
"""


class MeltStepPart(BaseStepPart):
    """MEI calling using Melt

    We perform cohort-wide genotyping of the mobile elements.  This requires a somewhat complicated
    workflow.  The methods input/output/etc. delegate to appropriate privated functions, thus
    this class is somewhat lengthy.  In the, however, it's not that complex.
    """

    name = "melt"

    def __init__(self, parent):
        super().__init__(parent)
        #: All individual's primary NGS libraries
        self.all_dna_ngs_libraries = []
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    self.all_dna_ngs_libraries.append(donor.dna_ngs_library.name)
        #: Linking NGS libraries to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        assert action in MELT_ACTIONS, "Invalid action: {}".format(action)
        mapping = {
            "preprocess": self._get_input_files_preprocess,
            "indiv_analysis": self._get_input_files_indiv_analysis,
            "group_analysis": self._get_input_files_group_analysis,
            "genotype": self._get_input_files_genotype,
            "make_vcf": self._get_input_files_make_vcf,
            "merge_vcf": self._get_input_files_merge_vcf,
            "reorder_vcf": self._get_input_files_reorder_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_preprocess(self, wildcards):
        # Get shorcut to NGS mapping sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_input_files_indiv_analysis(self, wildcards):
        tpl = "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}{ext}"
        yield "orig_bam", tpl.format(ext=".bam", **wildcards)
        yield "orig_bai", tpl.format(ext=".bam.bai", **wildcards)
        yield "disc_bam", tpl.format(ext=".bam.disc", **wildcards)
        yield "disc_bai", tpl.format(ext=".bam.disc.bai", **wildcards)

    @listify
    def _get_input_files_group_analysis(self, wildcards):
        for library_name in self.all_dna_ngs_libraries:
            yield "work/{mapper}.melt.indiv_analysis.{me_type}/out/.done.{library_name}".format(
                library_name=library_name, **wildcards
            )

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield "done", "work/{mapper}.melt.group_analysis.{me_type}/out/.done".format(**wildcards)
        yield "bam", "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam".format(
            **wildcards
        )

    @listify
    def _get_input_files_make_vcf(self, wildcards):
        yield "work/{mapper}.melt.group_analysis.{me_type}/out/.done".format(**wildcards)
        # Probably, best create a create-list intermediate target
        for library_name in self.all_dna_ngs_libraries:
            yield "work/{mapper}.melt.genotype.{me_type}/out/.done.{library_name}".format(
                library_name=library_name, **wildcards
            )

    @listify
    def _get_input_files_merge_vcf(self, wildcards):
        for me_type in self.config["melt"]["me_types"]:
            yield "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz".format(
                me_type=me_type, **wildcards
            )

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        yield "vcf", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz".format(
            **wildcards
        )
        yield "tbi", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi".format(
            **wildcards
        )

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    def get_output_files(self, action):
        assert action in MELT_ACTIONS
        mapping = {
            "preprocess": self._get_output_files_preprocess,
            "indiv_analysis": self._get_output_files_indiv_analysis,
            "group_analysis": self._get_output_files_group_analysis,
            "genotype": self._get_output_files_genotype,
            "make_vcf": self._get_output_files_make_vcf,
            "merge_vcf": self._get_output_files_merge_vcf,
            "reorder_vcf": self._get_output_files_reorder_vcf,
        }
        return mapping[action]()

    @dictify
    def _get_output_files_preprocess(self):
        # Note that mapper is not part of the output BAM file as MELT infers sample file from BAM
        # file name instead of using sample name from BAM header.
        tpl = "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}%s"
        yield "orig_bam", tpl % ".bam"
        yield "orig_bai", tpl % ".bam.bai"
        yield "disc_bam", tpl % ".bam.disc"
        yield "disc_bai", tpl % ".bam.disc.bai"
        yield "disc_fq", tpl % ".bam.fq"

    @dictify
    def _get_output_files_indiv_analysis(self):
        yield "done", touch("work/{mapper}.melt.indiv_analysis.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_output_files_group_analysis(self):
        yield "done", touch("work/{mapper}.melt.group_analysis.{me_type}/out/.done")

    @dictify
    def _get_output_files_genotype(self):
        yield "done", touch("work/{mapper}.melt.genotype.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_output_files_make_vcf(self):
        yield "list_txt", "work/{mapper}.melt.genotype.{me_type}/out/list.txt"
        yield "done", touch("work/{mapper}.melt.make_vcf.{me_type}/out/.done")
        yield "vcf", "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz"
        yield "tbi", "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz.tbi"

    @dictify
    def _get_output_files_merge_vcf(self):
        yield "vcf", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz"
        yield "tbi", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi"
        yield "vcf_md5", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.md5"
        yield "tbi_md5", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi.md5"

    @dictify
    def _get_output_files_reorder_vcf(self):
        tpl = "work/{mapper}.melt.{index_library_name}/out/{mapper}.melt.{index_library_name}.%s%s"
        key_ext = {"vcf": "vcf.gz", "tbi": "vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, tpl % (ext, "")
            yield key + "_md5", tpl % (ext, ".md5")

    def get_log_file(self, action):
        assert action in MELT_ACTIONS
        if action == "preprocess":
            return "work/{mapper}.melt.preprocess.{library_name}/log/snakemake.wgs_mei_calling.log"
        elif action in ("indiv_analysis", "genotype"):
            return (
                "work/{{mapper}}.melt.{action}.{{me_type}}/log/"
                "snakemake.wgs_mei_calling.{{library_name}}.log"
            ).format(action=action)
        elif action in ("indiv_analysis", "group_analysis", "genotype", "make_vcf"):
            return (
                "work/{{mapper}}.melt.{action}.{{me_type}}/log/" "snakemake.wgs_mei_calling.log"
            ).format(action=action)
        elif action == "merge_vcf":
            return "work/{mapper}.melt.merge_vcf/log/snakemake.wgs_mei_calling.log"
        elif action == "reorder_vcf":
            return (
                "work/{mapper}.melt.reorder_vcf.{index_library_name}/log/"
                "snakemake.wgs_mei_calling.log"
            )
        else:
            assert False, "Invalid action"

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration for running Melt"""
        for action in MELT_ACTIONS:
            cluster_config["wgs_mei_calling_melt_{}".format(action)] = {
                "mem": int(3.75 * 1024 * 6),
                "time": "128:00",
                "ntasks": 6,
            }


class WgsMeiCallingWorkflow(BaseStep):
    """Perform germline WGS MEI calling"""

    name = "wgs_mei_calling"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((MeltStepPart, LinkOutStepPart))
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_pattern = "{mapper}.{caller}.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=self.config["tools"],
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "wgs_mei_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for MEI calling",
        )
