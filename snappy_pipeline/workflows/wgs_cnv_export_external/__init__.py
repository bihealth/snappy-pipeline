# -*- coding: utf-8 -*-
"""VarFish export of externally generated WGS CNV calls

The step annotates externally generated CNV VCF files with ``varfish-annotator-cli annotate-sv`` and creates the files that VarFish Server imports.
The VCFs come from an ``external_data`` task.
The step has no ``tool`` field.

The task gives one result per pedigree, named after the primary DNA library of its index.
With ``merge_vcf_flag``, VCFs of a pedigree that are still separate are merged first.
``varfish_server_compatibility`` builds output for VarFish Server v1.2 and early versions of v2.

This step is experimental.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``varfish_annotated.{library_name}/out/varfish_annotated.{library_name}.gts.tsv.gz``
- ``varfish_annotated.{library_name}/out/varfish_annotated.{library_name}.feature-effects.tsv.gz``
- ``varfish_annotated.{library_name}/out/varfish_annotated.{library_name}.db-infos.tsv.gz``

Here ``library_name`` is the index library of the pedigree.
Downstream tasks read the keys ``gts`` and ``db_infos``.
"""

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeSampleNameStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import WgsCnvExportExternal as WgsCnvExportExternalConfigModel

#: Infixes to use for file name generation
INFIXES = ("gts", "feature-effects", "db-infos")

#: VCF key to file extensions
KEY_EXT = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}

#: Default configuration for the wgs_cnv_export_external step


class VarfishAnnotatorExternalStepPart(BaseStepPart):
    """Annotate externally generated VCF file using "varfish-annotator annotate-sv"."""

    #: Step name
    name = "varfish_annotator_external"

    #: Class available actions
    actions = ("annotate", "merge_vcf")

    intermediate_outputs = {"merge_vcf": ("vcf", "vcf_tbi")}

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Define mapper+caller tag
        self.mapper_caller_tag = self._get_mapper_caller_tag()

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        if self.config.merge_vcf_flag:
            yield "vcf", sorted(self._collect_vcfs(wildcards))
        else:
            yield "vcf", [self._vcf(wildcards.index_ngs_library)]

    @dictify
    def _get_input_files_annotate(self, wildcards):
        """"""
        # Pedigree
        tpl = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        yield "ped", tpl.format(**wildcards)
        # Reference
        yield "reference", self.parent.get_upstream_paths("reference").fasta
        yield "db", self.config.path_db
        yield "refseq_ser", self.config.path_refseq_ser
        yield "ensembl_ser", self.config.path_ensembl_ser
        # VCF
        tpl = (
            f"work/{self.mapper_caller_tag}{{index_ngs_library}}/out/"
            f"{self.mapper_caller_tag}{{index_ngs_library}}"
        )
        for key, ext in KEY_EXT.items():
            yield key, tpl.format(**wildcards) + ext

    @dictify
    def _get_output_files_merge_vcf(self):
        tpl = (
            f"work/{self.mapper_caller_tag}{{index_ngs_library}}/out/"
            f"{self.mapper_caller_tag}{{index_ngs_library}}"
        )
        for key, ext in KEY_EXT.items():
            yield key, tpl + ext

    @dictify
    def _get_output_files_annotate(self):
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
        )
        for infix in INFIXES:
            key = infix.replace("-", "_")
            yield key, prefix + f".{infix}.tsv.gz"
            yield key + "_md5", prefix + f".{infix}.tsv.gz.md5"

    def get_result_files(self) -> list[str]:
        """The files of the export, for the index of each pedigree"""
        indexes = [
            pedigree.index.dna_ngs_library.name
            for sheet in filter(is_not_background, self.parent.shortcut_sheets)
            for pedigree in sheet.cohort.pedigrees
            if pedigree.index and pedigree.index.dna_ngs_library
        ]
        return [
            path
            for tpl in self.get_output_files("annotate").values()
            for path in expand(tpl.replace("work/", "output/", 1), index_ngs_library=indexes)
        ]

    def _get_log_file_merge_vcf(self):
        return {
            "log": (
                f"work/{self.mapper_caller_tag}{{index_ngs_library}}/log/"
                f"{self.mapper_caller_tag}{{index_ngs_library}}.merge_vcf.log"
            )
        }

    def _get_log_file_annotate(self):
        return {
            "log": "work/varfish_annotated.{index_ngs_library}/log/"
            "varfish_annotated.{index_ngs_library}.log"
        }

    resource_usage = {
        "annotate": ResourceUsage(threads=2, runtime="4d", mem=f"{7 * 1024 * 2}MB"),
    }
    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem=f"{7 * 1024 * 2}MB")

    def _get_params_merge_vcf(self, wildcards):
        result = {
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
            "merge_option": self.config.merge_option,
            "gvcf_option": False,
        }
        return result

    def _get_params_annotate(self, wildcards):
        varfish_server_compatibility_flag = self.config.varfish_server_compatibility
        return {
            "step_name": "wgs_cnv_export_external",
            "varfish_server_compatibility": varfish_server_compatibility_flag,
            "config": self.config.model_dump(
                by_alias=True, exclude={"path_db", "path_refseq_ser", "path_ensembl_ser"}
            ),
        }

    def _vcf(self, library_name):
        return self.parent.get_upstream_paths("variants", library_name=library_name).vcf

    def _collect_vcfs(self, wildcards):
        """Yield the paths of the VCFs of the pedigree's DNA libraries"""
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
            yield self._vcf(donor.dna_ngs_library.name)

    def _collect_sample_ids(self, wildcards):
        """Yield sample ids in pedigree"""
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
            yield donor.name

    def _get_mapper_caller_tag(self):
        """Get mapper and caller tag

        :return: Returns tag to be used to name intermediate and final files. Tag based on
        information provided in configuration. Output examples: 'bwa.delly2.', 'dragen.', or ''.
        """
        mapper = self.config.tool_ngs_mapping
        caller = self.config.tool_wgs_cnv_calling
        if mapper and caller:
            return ""
        elif mapper or caller:
            mapper = mapper or ""
            caller = caller or ""
            return ""
        else:
            return ""


class WgsCnvExportExternalWorkflow(BaseStep):
    """Perform germline WGS CNV export external"""

    #: Workflow name
    name = "wgs_cnv_export_external"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"external"}))]
    config_model_class = WgsCnvExportExternalConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local external WGS CNV export output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"work/varfish_annotated.{lib}/out/varfish_annotated.{lib}"
        return {
            "gts": f"{prefix}.gts.tsv.gz",
            "db_infos": f"{prefix}.db-infos.tsv.gz",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Load external data search information
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeSampleNameStepPart,
                VarfishAnnotatorExternalStepPart,
            )
        )
