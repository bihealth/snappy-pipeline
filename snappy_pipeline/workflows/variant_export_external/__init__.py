# -*- coding: utf-8 -*-
"""Implementation of the ``variant_export_external`` step.

The ``variant_export_external`` takes externally generated VCF files and uses
``varfish-annotator-cli annotate`` command to create files fit for import into VarFish
Server. If externally generated BAM files are also available, it can optionally generate QC reports.

==========
Stability
==========

This step is considered experimental, use it at your own discretion.

==========
Step Input
==========

The variant export external step uses VCFs externally generated as input.

===========
Step Output
===========

For all samples, the workflow will execute ``varfish-annotator-cli``, pedigrees with multiple VCFs
will be merged based on `id`. The name of the index's primary DNA NGS library will be used as an
identification token in the output file. The following is an example of the files that
will be generated if BAM files are also available:


::
    output/
    +-- varfish_annotated.P001-N1-DNA1-WGS1
    |   `-- out
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.bam-qc.tsv.gz
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.bam-qc.tsv.gz.md5
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz.md5
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz.md5
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz
    |   |   +-- varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz.md5
    |
    [...]

====================
Global Configuration
====================

Not applicable.

=====================
Default Configuration
=====================

.. include:: DEFAULT_CONFIG_variant_export_external.rst

==================
Parallel Execution
==================

Parallel execution is not performed currently.
"""

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeSampleNameStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import TargetCovReportStepPart

from .model import VariantExportExternal as VariantExportExternalConfigModel

#: Default configuration for the somatic_variant_calling step


class BamReportsExternalStepPart(TargetCovReportStepPart):
    """Build target coverage report and QC report for external BAM files"""

    #: Step name
    name = "bam_reports"

    #: Class available actions
    actions = ("bam_qc", "collect", "run")

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def _alignments(self, library_name):
        return self.parent.get_upstream_paths("alignments", library_name=library_name)

    @dictify
    def _get_input_files_bam_qc(self, wildcards):
        alignments = self._alignments(wildcards.library_name)
        yield "bam", alignments.bam
        yield "bai", alignments.bai

    @dictify
    def _get_input_files_run(self, wildcards):
        alignments = self._alignments(wildcards.library_name)
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        if self.config.target_coverage_report.path_targets_bed:
            yield "targets_bed", self.config.target_coverage_report.path_targets_bed

    @listify
    def _get_input_files_collect(self, wildcards):
        _ = wildcards
        mapper_lib = "{library_name}"
        yield f"work/{mapper_lib}/report/alfred_qc/{mapper_lib}.alfred.json.gz"

    @dictify
    def _get_output_files_collect(self):
        yield "report", "work/target_cov_report/out/target_coverage_summary.txt"

    @dictify
    def _get_output_files_bam_qc(self):
        for report in ("bamstats", "flagstats", "idxstats"):
            report_path = f"work/{{library_name}}/report/bam_qc/{{library_name}}.bam.{report}.txt"
            yield report, report_path

    def _get_log_file_run(self):
        return "work/{library_name}/log/snakemake.target_coverage.log"

    def _get_log_file_collect(self):
        return "work/target_cov_report/log/snakemake.target_coverage.log"

    def _get_log_file_bam_qc(self):
        return {"log": "work/{library_name}/log/{library_name}.bam_qc.log"}

    def get_result_files(self):
        """The reports are not published."""
        return []

    def _collect_bam_files(self, wildcards):
        """Yield the BAM file of the library."""
        yield self._alignments(wildcards.library_name).bam


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    #: Step name
    name = "varfish_annotator_external"

    #: Class available actions
    actions = ("annotate", "bam_qc", "gvcf_to_vcf", "merge_vcf")

    intermediate_outputs = {"gvcf_to_vcf": ("vcf", "vcf_tbi"), "merge_vcf": ("vcf", "vcf_tbi")}

    def get_result_files(self) -> list[str]:
        """The files of the export, for the index of each pedigree"""
        indexes = [
            pedigree.index.dna_ngs_library.name
            for sheet in filter(is_not_background, self.parent.shortcut_sheets)
            for pedigree in sheet.cohort.pedigrees
            if pedigree.index and pedigree.index.dna_ngs_library
        ]
        keys = ["gts", "gts_md5", "db_infos", "db_infos_md5"]
        outputs = self.get_output_files("annotate")
        if self.config.bam_available_flag:
            outputs |= self.get_output_files("bam_qc")
            keys += ["bam_qc", "bam_qc_md5"]
        return [
            path
            for key in keys
            for path in expand(
                outputs[key].replace("work/", "output/", 1), index_ngs_library=indexes
            )
        ]

    #: VCF key to file extensions
    vcf_key_ext_dict = {
        "vcf": ".vcf.gz",
        "vcf_tbi": ".vcf.gz.tbi",
    }

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # External tool prefix - included in merged VCF name
        self.external_tool_prefix = self._get_external_source_prefix()

    def _get_external_source_prefix(self):
        """Get external source prefix

        :return: Returns external tool prefix if any provided, example: "dragen.". Otherwise,
        returns empty string.
        """
        if self.config.external_tool:
            return self.config.external_tool.lower() + "."
        return ""

    def _vcf(self, library_name):
        return self.parent.get_upstream_paths("variants", library_name=library_name).vcf

    @dictify
    def _get_input_files_gvcf_to_vcf(self, wildcards):
        yield "vcf", sorted(self._collect_gvcf(wildcards))
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        if self.config.merge_vcf_flag:
            yield "vcf", sorted(self._collect_vcfs(wildcards))
        else:
            yield "vcf", [self._vcf(wildcards.index_ngs_library)]
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    @dictify
    def _get_input_files_annotate(self, wildcards):
        # Pedigree
        yield (
            "ped",
            "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped".format(
                **wildcards
            ),
        )
        # Reference
        yield "reference", self.parent.get_upstream_paths("reference").fasta
        yield "db", self.config.path_db
        yield "refseq_ser", self.config.path_refseq_ser
        yield "ensembl_ser", self.config.path_ensembl_ser
        if self.config.path_exon_bed:
            yield "exon_bed", self.config.path_exon_bed
        # VCF
        tpl = (
            f"work/{self.external_tool_prefix}{{index_ngs_library}}/out/"
            f"{self.external_tool_prefix}{{index_ngs_library}}"
        )
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, tpl.format(**wildcards) + ext

    def _get_input_files_bam_qc(self, wildcards):
        # Get names of primary libraries of the selected pedigree.  The pedigree is selected
        # by the primary DNA NGS library of the index.
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "cov_qc": []}
        for donor in pedigree.donors:
            _mapper = self.external_tool_prefix[:-1]  # strip trailing dot
            library_name = donor.dna_ngs_library.name
            if not donor.dna_ngs_library:
                continue
            tpl = f"work/{library_name}/report/bam_qc/{library_name}.bam.%s.txt"
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(tpl % key)
            if donor.dna_ngs_library.name not in self.parent.ngs_library_list:
                continue
            path = f"work/{library_name}/report/alfred_qc/{library_name}.alfred.json.gz"
            result["cov_qc"].append(path)

        return result

    @dictify
    def _get_output_files_gvcf_to_vcf(self):
        tpl = (
            f"work/{self.external_tool_prefix}{{index_ngs_library}}/out/"
            f"{self.external_tool_prefix}{{index_ngs_library}}"
        )
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, tpl + ext

    _get_output_files_merge_vcf = _get_output_files_gvcf_to_vcf

    @dictify
    def _get_output_files_annotate(self):
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
        )
        paths_work = {}
        for infix in ("gts", "db-infos"):
            key = infix.replace("-", "_")
            paths_work[key] = prefix + f".{infix}.tsv.gz"
            paths_work[key + "_md5"] = prefix + f".{infix}.tsv.gz.md5"
        paths_work["ped"] = f"{prefix}.ped"
        yield from paths_work.items()

    @dictify
    def _get_output_files_bam_qc(self):
        key = "bam_qc"
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
        )
        yield key, prefix + ".bam-qc.tsv.gz"
        yield key + "_md5", prefix + ".bam-qc.tsv.gz.md5"

    def _get_log_file(self, action):
        if action in ("gvcf_to_vcf", "merge_vcf"):
            prefix = (
                f"work/{self.external_tool_prefix}{{index_ngs_library}}/"
                f"log/{self.external_tool_prefix}{{index_ngs_library}}.{action}"
            )
        else:
            prefix = (
                f"work/varfish_annotated.{{index_ngs_library}}/log/"
                f"varfish_annotated.{action}.{{index_ngs_library}}"
            )
        return {"log": prefix + ".log"}

    resource_usage = {
        "annotate": ResourceUsage(threads=2, runtime="4d", mem=f"{7 * 1024 * 2}MB"),
    }
    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem=f"{7 * 1024 * 2}MB")

    def _get_params_gvcf_to_vcf(self, wildcards):
        result = {
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
        }
        return result

    def _get_params_merge_vcf(self, wildcards):
        result = {
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
            "merge_option": self.config.merge_option,
        }
        return result

    def _get_params_annotate(self, wildcards):
        return {
            "config": self.config.model_dump(
                by_alias=True,
                exclude={"path_db", "path_refseq_ser", "path_ensembl_ser", "path_exon_bed"},
            )
        }

    def _get_params_bam_qc(self, wildcards):
        """Get parameters for wrapper ``variant_annotator/bam_qc``

        Creates dictionary that links library name to identifier that should be used in output file.
        The wrapper will derive the library name from the input file name, for analysis using
        externally generated data, the values will be the sample name as provided by the external
        source (sample name). For snappy-based analysis, both keys and values will be the
        library name.

        Dictionary expected structure:
        {
            "P001-N1-DNA1-WGS1": "P001",
            "P002-N1-DNA1-WGS1": "P002",
            "P003-N1-DNA1-WGS1": "P003",
        }

        :return: Dictionary linking library name to identifier that should be used in output file.
        Key: library name; Value: identifier to be used in file.
        """
        library_name_to_file_identifier = {}
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                library_name_to_file_identifier[donor.dna_ngs_library.name] = donor.secondary_id
        return library_name_to_file_identifier

    def _collect_gvcf(self, wildcards):
        """Yield the path of the index's gVCF"""
        yield self._vcf(wildcards.index_ngs_library)

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


class VariantExportExternalWorkflow(BaseStep):
    """Perform germline variant export for externally generated VCFs"""

    #: Workflow name
    name = "variant_export_external"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"external"}))]
    config_model_class = VariantExportExternalConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local external-export output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"work/varfish_annotated.{lib}/out/varfish_annotated.{lib}"
        return {
            "gts": f"{prefix}.gts.tsv.gz",
            "db_infos": f"{prefix}.db-infos.tsv.gz",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeSampleNameStepPart,
                VarfishAnnotatorAnnotateStepPart,
                BamReportsExternalStepPart,
            )
        )
        self.ngs_library_list = self._build_ngs_library_list()

    @listify
    def _build_ngs_library_list(self):
        for sheet in self.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    for test_sample in bio_sample.test_samples.values():
                        for library in test_sample.ngs_libraries.values():
                            yield library.name
