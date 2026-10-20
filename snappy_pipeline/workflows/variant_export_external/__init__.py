# -*- coding: utf-8 -*-
"""VarFish export of externally generated germline variants

The step annotates externally generated (g)VCF files with ``varfish-annotator-cli annotate`` and creates the files that VarFish Server imports.
The VCFs and, optionally, the BAM files come from ``external_data`` tasks.
The step has no ``tool`` field.

The task gives one result per pedigree, named after the primary DNA library of its index.
With ``merge_vcf_flag``, VCFs of a pedigree that are still separate are merged first.
With ``gvcf_option``, the inputs are converted from genomic VCFs.
With ``bam_available_flag``, the step also writes a BAM QC table from the BAM files in ``depends_on.alignments``.
The target regions of each library's kit come from ``depends_on.kits``, so every library must be a targeted library (``WES``, ``Panel_seq``) whose kit that task has.

This step is experimental.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``varfish_annotated.{library_name}/out/varfish_annotated.{library_name}.gts.tsv.gz``
- ``varfish_annotated.{library_name}/out/varfish_annotated.{library_name}.db-infos.tsv.gz``
- ``varfish_annotated.{library_name}/out/varfish_annotated.{library_name}.bam-qc.tsv.gz``, only with ``bam_available_flag``.

Here ``library_name`` is the index library of the pedigree.
Downstream tasks read the keys ``gts`` and ``db_infos``.
"""

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
        yield "targets_bed", self.parent.library_kit_paths("kits", wildcards.library_name).targets

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
        indexes = list(self.parent.pedigrees())
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
        yield from self.parent.reference_inputs().items()

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        if self.config.merge_vcf_flag:
            yield "vcf", sorted(self._collect_vcfs(wildcards))
        else:
            yield "vcf", [self._vcf(wildcards.index_ngs_library)]
        yield from self.parent.reference_inputs().items()

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
        yield from self.parent.reference_inputs().items()
        yield "db", self.config.path_db
        yield "refseq_ser", self.config.path_refseq_ser
        yield "ensembl_ser", self.config.path_ensembl_ser
        if self.config.depends_on.exons:
            yield "exon_bed", self.parent.get_upstream_paths("exons").bed
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
        members = self.parent.pedigree_members(wildcards.index_ngs_library, include_background=True)
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "cov_qc": []}
        for library_name in members["library_name"]:
            tpl = f"work/{library_name}/report/bam_qc/{library_name}.bam.%s.txt"
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(tpl % key)
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
                exclude={"path_db", "path_refseq_ser", "path_ensembl_ser"},
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
        members = self.parent.pedigree_members(wildcards.index_ngs_library, include_background=True)
        return {
            library_name: secondary_id
            for library_name, secondary_id in zip(members["library_name"], members["secondary_id"])
        }

    def _collect_gvcf(self, wildcards):
        """Yield the path of the index's gVCF"""
        yield self._vcf(wildcards.index_ngs_library)

    def _collect_vcfs(self, wildcards):
        """Yield the paths of the VCFs of the pedigree's DNA libraries"""
        members = self.parent.pedigree_members(wildcards.index_ngs_library, include_background=True)
        for library_name in members["library_name"]:
            yield self._vcf(library_name)

    def _collect_sample_ids(self, wildcards):
        """Yield sample ids in pedigree"""
        members = self.parent.pedigree_members(wildcards.index_ngs_library, include_background=True)
        yield from members["donor_name"]


class VariantExportExternalWorkflow(BaseStep):
    """Perform germline variant export for externally generated VCFs"""

    #: Workflow name
    name = "variant_export_external"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"external"}))]
    config_model_class = VariantExportExternalConfigModel

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
        if self.config.bam_available_flag:
            self.check_library_kits(
                "kits", self.primary_dna_libraries(include_background=True)["library_name"]
            )
