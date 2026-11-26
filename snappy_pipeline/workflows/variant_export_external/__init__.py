# -*- coding: utf-8 -*-
"""VarFish export of externally generated germline variants

The step annotates externally generated (g)VCF files with ``varfish-annotator-cli annotate`` and creates the files that VarFish Server imports.
The VCFs and the optional BAM files come from ``external_data`` tasks.
The step has no ``tool`` field.

The task gives one result per pedigree, named after the pedigree (its ``family_id``).
With ``merge_vcf_flag``, the step first merges the VCFs of a pedigree that are still separate.
With ``gvcf_option`` (the default), the inputs are genomic VCFs, which the task converts to VCFs before the annotation; without it, the task merges the VCFs.
With ``bam_available_flag``, the step also writes a BAM QC table from the BAM files in ``depends_on.alignments``.
The target regions of each library's kit come from ``depends_on.kits``.
Every library must therefore be a targeted library (``WES``, ``Panel_seq``) whose kit is in that task.

This step is experimental.

With ``depends_on.exons``, a BED file of exons, the export keeps only the small variants near the exons.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``varfish_annotated.{related_group}/out/varfish_annotated.{related_group}.gts.tsv.gz``
- ``varfish_annotated.{related_group}/out/varfish_annotated.{related_group}.db-infos.tsv.gz``
- ``varfish_annotated.{related_group}/out/varfish_annotated.{related_group}.bam-qc.tsv.gz``, only with ``bam_available_flag``.

Here the group is the name of the pedigree, whose index library holds the external files.
Downstream tasks read the keys ``gts`` and ``db_infos``.
"""

from snakemake.io import expand

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeSampleNameStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import VariantExportExternal as VariantExportExternalConfigModel


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    #: Step name
    name = "varfish_annotator_external"

    #: Class available actions
    actions = ("annotate", "bam_qc", "gvcf_to_vcf", "merge_vcf")

    @property
    def wrappers(self) -> dict[str, str]:
        """The VCF comes from ``gvcf_to_vcf`` with ``gvcf_option``, else from ``merge_vcf``"""
        vcf = (
            {"gvcf_to_vcf": "bcftools/gvcf_to_vcf"}
            if self.config.gvcf_option
            else {"merge_vcf": "bcftools/merge_snv_vcf"}
        )
        return {
            "bam_qc": "varfish_annotator/bam_qc",
            "annotate": "varfish_annotator/annotate",
            **vcf,
        }

    intermediate_outputs = {"gvcf_to_vcf": ("vcf", "vcf_tbi"), "merge_vcf": ("vcf", "vcf_tbi")}
    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem=f"{7 * 1024 * 2}MB")

    resource_usage = {
        "annotate": ResourceUsage(threads=2, runtime="4d", mem=f"{7 * 1024 * 2}MB"),
    }

    #: VCF key to file extensions
    vcf_key_ext_dict = {
        "vcf": ".vcf.gz",
        "vcf_tbi": ".vcf.gz.tbi",
    }

    #: The (merged) VCF of a pedigree, without the extension
    vcf_tpl = "work/{related_group}/out/{related_group}"

    def _get_log_file(self, action):
        if action in ("gvcf_to_vcf", "merge_vcf"):
            prefix = f"work/{{related_group}}/log/{{related_group}}.{action}"
        else:
            prefix = (
                f"work/varfish_annotated.{{related_group}}/log/"
                f"varfish_annotated.{action}.{{related_group}}"
            )
        return {"log": prefix + ".log"}

    def get_result_files(self) -> list[str]:
        """The files of the export, for each pedigree whose index has a DNA library"""
        indexes = [p for p in self.parent.pedigrees() if self.parent.pedigree_index(p)]
        keys = ["gts", "gts_md5", "db_infos", "db_infos_md5"]
        outputs = self.get_output_files("annotate")
        if self.config.bam_available_flag:
            outputs |= self.get_output_files("bam_qc")
            keys += ["bam_qc", "bam_qc_md5"]
        return [path for key in keys for path in expand(outputs[key], related_group=indexes)]

    @dictify
    def _get_input_files_annotate(self, wildcards):
        # Pedigree
        yield (
            "ped",
            "work/write_pedigree.{related_group}/out/{related_group}.ped".format(**wildcards),
        )
        # Reference
        yield from self.parent.reference_inputs().items()
        yield "db", self.config.path_db
        yield "refseq_ser", self.config.path_refseq_ser
        yield "ensembl_ser", self.config.path_ensembl_ser
        if self.config.depends_on.exons:
            yield "exon_bed", self.parent.get_upstream_paths("exons").bed
        # VCF
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, self.vcf_tpl.format(**wildcards) + ext

    @dictify
    def _get_output_files_annotate(self):
        prefix = "work/varfish_annotated.{related_group}/out/varfish_annotated.{related_group}"
        paths_work = {}
        for infix in ("gts", "db-infos"):
            key = infix.replace("-", "_")
            paths_work[key] = prefix + f".{infix}.tsv.gz"
            paths_work[key + "_md5"] = prefix + f".{infix}.tsv.gz.md5"
        paths_work["ped"] = f"{prefix}.ped"
        yield from paths_work.items()

    def _get_params_annotate(self, wildcards):
        return {
            "config": self.config.model_dump(
                by_alias=True,
                exclude={"path_db", "path_refseq_ser", "path_ensembl_ser"},
            )
        }

    def _get_input_files_bam_qc(self, wildcards):
        # Get names of primary libraries of the selected pedigree.  The pedigree is selected
        # by the primary DNA NGS library of the index.
        members = self.parent.pedigree_members(wildcards.related_group, include_background=True)
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "alfred_qc": []}
        for library_name in members["library_name"]:
            tpl = f"work/{library_name}/report/bam_qc/{library_name}.bam.%s.txt"
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(tpl % key)
            path = f"work/{library_name}/out/{library_name}.target_coverage.json.gz"
            result["alfred_qc"].append(path)

        return result

    @dictify
    def _get_output_files_bam_qc(self):
        key = "bam_qc"
        prefix = "work/varfish_annotated.{related_group}/out/varfish_annotated.{related_group}"
        yield key, prefix + ".bam-qc.tsv.gz"
        yield key + "_md5", prefix + ".bam-qc.tsv.gz.md5"

    def _get_params_bam_qc(self, wildcards):
        """The names of the pedigree members in the table, the sample names of the external
        source, in the order of the ``bamstats``, ``flagstats`` and ``idxstats`` inputs"""
        members = self.parent.pedigree_members(wildcards.related_group, include_background=True)
        return {"sample_names": members["secondary_id"].tolist()}

    @dictify
    def _get_input_files_gvcf_to_vcf(self, wildcards):
        yield "vcf", sorted(self._collect_gvcf(wildcards))
        yield from self.parent.reference_inputs().items()

    @dictify
    def _get_output_files_gvcf_to_vcf(self):
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, self.vcf_tpl + ext

    def _get_params_gvcf_to_vcf(self, wildcards):
        result = {
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
        }
        return result

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        if self.config.merge_vcf_flag:
            yield "vcf", sorted(self._collect_vcfs(wildcards))
        else:
            yield "vcf", [self._vcf(self.parent.pedigree_index(wildcards.related_group))]
        yield from self.parent.reference_inputs().items()

    _get_output_files_merge_vcf = _get_output_files_gvcf_to_vcf

    def _get_params_merge_vcf(self, wildcards):
        result = {
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
            "merge_option": self.config.merge_option,
        }
        return result

    def _vcf(self, library_name):
        return self.parent.get_upstream_paths("variants", library_name=library_name).vcf

    def _collect_gvcf(self, wildcards):
        """Yield the path of the index's gVCF"""
        yield self._vcf(self.parent.pedigree_index(wildcards.related_group))

    def _collect_vcfs(self, wildcards):
        """Yield the paths of the VCFs of the pedigree's DNA libraries"""
        members = self.parent.pedigree_members(wildcards.related_group, include_background=True)
        for library_name in members["library_name"]:
            yield self._vcf(library_name)

    def _collect_sample_ids(self, wildcards):
        """Yield sample ids in pedigree"""
        members = self.parent.pedigree_members(wildcards.related_group, include_background=True)
        yield from members["donor_name"]


class BamReportsExternalStepPart(BaseStepPart):
    """Build target coverage report and QC report for external BAM files"""

    #: Step name
    name = "bam_reports"

    #: Class available actions
    actions = ("bam_qc", "run")

    wrappers = {"run": "alfred/qc", "bam_qc": "samtools/bam_qc_report"}
    default_resource_usage = ResourceUsage(threads=2, runtime="4h", mem="20GB")

    def get_result_files(self):
        """The reports are no result files."""
        return []

    @dictify
    def _get_input_files_bam_qc(self, wildcards):
        alignments = self._alignments(wildcards.library_name)
        yield "bam", alignments.bam
        yield "bai", alignments.bai

    @dictify
    def _get_output_files_bam_qc(self):
        for report in ("bamstats", "flagstats", "idxstats"):
            report_path = f"work/{{library_name}}/report/bam_qc/{{library_name}}.bam.{report}.txt"
            yield report, report_path

    def _get_log_file_bam_qc(self):
        return {"log": "work/{library_name}/log/{library_name}.bam_qc.log"}

    @dictify
    def _get_input_files_run(self, wildcards):
        alignments = self._alignments(wildcards.library_name)
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield from self.parent.reference_inputs().items()
        yield "targets_bed", self.parent.library_kit_paths("kits", wildcards.library_name).targets

    def _get_output_files_run(self):
        """The alfred report, without md5: it is no result file."""
        return {"json": "work/{library_name}/out/{library_name}.target_coverage.json.gz"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.target_cov_report.log"}

    def _get_params_run(self, wildcards):
        return {"ignore_read_groups": self.config.ignore_read_groups}

    def _alignments(self, library_name):
        return self.parent.get_upstream_paths("alignments", library_name=library_name)


class VariantExportExternalWorkflow(BaseStep):
    """Perform germline variant export for externally generated VCFs"""

    #: Workflow name
    name = "variant_export_external"

    produces = [DataSignature(DataType.EXPORTS, frozenset({"external"}))]
    config_model_class = VariantExportExternalConfigModel

    #: The pedigrees include the members of background libraries
    uses_background = True

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
