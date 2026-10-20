# -*- coding: utf-8 -*-
"""Somatic CNV checking with B-allele frequencies

The step finds heterozygous germline variants in the matched normal and measures their B-allele frequency in the tumor.
It then merges these frequencies with the segments of a somatic CNV caller.
This shows whether the copy number calls agree with the allelic imbalance in the tumor.
The variants in the regions of ``depends_on.excluded_regions`` are left out, if it is set.

The step has no ``tool`` field.
It gives one result per tumor DNA library, paired with its matched normal.
It publishes files only when ``depends_on.copy_number`` names a CNV calling task.
Without it, the pileups stay in ``work/`` as intermediate files.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``{tumor_library}/out/{tumor_library}.vcf.gz`` with ``.tbi``: the heterozygous variants with their B-allele frequency, and the log2 ratio and copy number call of the segment.
- ``{tumor_library}/out/{tumor_library}.tsv``
- ``{tumor_library}/report/{tumor_library}.cnv.pdf``, ``{tumor_library}.locus.pdf`` and ``{tumor_library}.segment.pdf``: plots of the B-allele frequency against the copy number.

Downstream tasks read the key ``vcf``.
"""

import os
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import SomaticCnvChecking as SomaticCnvCheckingConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


class SomaticCnvCheckingPileupStepPart(BaseStepPart):
    """Perform pileups at a set of locations"""

    name = "pileup"
    actions = ("normal", "tumor")

    intermediate_outputs = {"normal": ("vcf", "vcf_tbi"), "tumor": ("vcf", "vcf_tbi")}

    resource_usage = {
        "normal": ResourceUsage(threads=2, runtime="12h", mem=f"{int(3.7 * 1024 * 2)}MB"),
        "tumor": ResourceUsage(threads=2, runtime="1h", mem=f"{int(3.7 * 1024 * 2)}MB"),
    }

    def _get_input_files_normal(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.tumor_library
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
        }

    def _get_input_files_tumor(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.tumor_library
        )
        normal_lib = self.parent._get_normal_lib(wildcards.tumor_library)
        return {
            "locii": "work/{normal_library}/out/{normal_library}.normal.vcf.gz".format(
                normal_library=normal_lib, **wildcards
            ),
            "locii_tbi": "work/{normal_library}/out/{normal_library}.normal.vcf.gz.tbi".format(
                normal_library=normal_lib, **wildcards
            ),
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
        }

    @staticmethod
    def _pileup_outputs(action):
        """Return the pileup VCF + TBI file (intermediate, so without md5 files)"""
        base = f"work/{{tumor_library}}/out/{{tumor_library}}.{action}.vcf.gz"
        return {"vcf": base, "vcf_tbi": base + ".tbi"}

    def _get_output_files_normal(self):
        return self._pileup_outputs("normal")

    def _get_output_files_tumor(self):
        return self._pileup_outputs("tumor")

    def _get_params_normal(self, wildcards):
        return {
            "min_baf": self.config.min_baf,
            "min_depth": self.config.min_depth,
            "max_depth": self.config.max_depth,
        }

    _get_params_tumor = _get_params_normal

    def _get_log_file(self, action):
        return {"log": f"work/{{tumor_library}}/log/{{tumor_library}}.{action}.log"}


class SomaticCnvCheckingCnvStepPart(BaseStepPart):
    """Merging heterozygous variants with CNV data"""

    name = "cnv"
    actions = ("run",)

    entity_wildcard = "tumor_library"

    published_keys = ("vcf", "vcf_md5", "tbi", "tbi_md5", "tsv", "tsv_md5")

    def _get_input_files_run(self, wildcards):
        normal_library = self.parent._get_normal_lib(wildcards.tumor_library)
        filenames = {}
        name_pattern = "{normal_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".normal.vcf.gz")
        filenames["normal"] = tpl.format(normal_library=normal_library, **wildcards)
        filenames["normal_tbi"] = filenames["normal"] + ".tbi"
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "out", name_pattern + ".tumor.vcf.gz")
        filenames["tumor"] = tpl.format(**wildcards)
        filenames["tumor_tbi"] = filenames["tumor"] + ".tbi"
        filenames["cnv"] = self.parent.get_upstream_paths(
            "copy_number", library_name=wildcards.tumor_library
        ).dnacopy_seg
        if self.config.depends_on.excluded_regions:
            filenames["excluded_regions"] = self.parent.get_upstream_paths("excluded_regions").bed
        return filenames

    @dictify
    def _get_output_files_run(self):
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi", "tsv": ".tsv"}
        for key, ext in key_ext.items():
            yield key, "work/{tumor_library}/out/{tumor_library}" + ext
            yield key + "_md5", "work/{tumor_library}/out/{tumor_library}" + ext + ".md5"

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.log"}


class SomaticCnvCheckingReportStepPart(BaseStepPart):
    """Plots of BAF vs CNV"""

    name = "report"
    actions = ("run",)

    entity_wildcard = "tumor_library"

    published_keys = (
        "cnv",
        "cnv_md5",
        "locus",
        "locus_md5",
        "segment",
        "segment_md5",
    )

    def _get_input_files_run(self, wildcards):
        name_pattern = "{tumor_library}".format(**wildcards)
        base_path_out = "work/" + name_pattern + "/out/" + name_pattern
        return {
            "vcf": base_path_out + ".vcf.gz",
            "tsv": base_path_out + ".tsv",
            **self.parent.reference_inputs(),
        }

    def _get_output_files_run(self):
        base_path_out = "work/{tumor_library}/report/{tumor_library}"
        return {
            "cnv": base_path_out + ".cnv.pdf",
            "cnv_md5": base_path_out + ".cnv.pdf.md5",
            "locus": base_path_out + ".locus.pdf",
            "locus_md5": base_path_out + ".locus.pdf.md5",
            "segment": base_path_out + ".segment.pdf",
            "segment_md5": base_path_out + ".segment.pdf.md5",
        }

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "tumor_library": wildcards.tumor_library,
        }

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.report.log"}


class SomaticCnvCheckingWorkflow(BaseStep):
    """Perform somatic cnv checking"""

    #: Workflow name
    name = "somatic_cnv_checking"
    produces = [DataSignature(DataType.QC, frozenset({"cnv_check"}))]

    config_model_class = SomaticCnvCheckingConfigModel

    library_requirement = "kind == 'cancer'"

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local CNV checking output paths for downstream consumers."""
        lib = kwargs.get("tumor_library", "{tumor_library}")
        return {"vcf": f"work/{lib}/out/{lib}.vcf.gz"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.has_cnv_calling = bool(self.config.depends_on.copy_number)
        # Register sub step classes so the sub steps are available
        sub_steps = [SomaticCnvCheckingPileupStepPart]
        if self.has_cnv_calling:
            sub_steps += [SomaticCnvCheckingCnvStepPart, SomaticCnvCheckingReportStepPart]
        self.register_sub_step_classes(tuple(sub_steps))

    def _get_normal_lib(self, tumor_library: str) -> str | None:
        """Resolve matched normal library for a tumor library via the DataFrame."""
        df = self.build_library_dataframe()
        tumor_df = df[df["library_name"] == tumor_library]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None
