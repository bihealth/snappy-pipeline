# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_cnv_checking`` step

The ``somatic_cnv_checking`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs pileups at germline locii to identify heterozygous
variants, and their B-allele frequency in the paired tumor sample. If one somatic CNV calling step
(either ``somatic_wgs_cnv_calling`` or ``somatic_targetd_seq_cnv_calling``) is present, the
log2 convarage ratio and copy number call are added to the output vcf.

==========
Step Input
==========

The somatic CNV checking step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step. Optionally, it can also use one Snakemake sub workflows for somatic CNV
calling, either``somatic_wgs_cnv_calling`` or ``somatic_targetd_seq_cnv_calling``.

===========
Step Output
===========

For each tumor DNA NGS library with name ``lib_name``/key ``lib_pk``
the pipeline step will create a directory ``output/{lib_name}-{lib_pk}/out``
with symlinks of the following names to the resulting VCF, TBI, and MD5 files.

- ``{lib_name}-{lib_pk}.vcf.gz``
- ``{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.cnvkit.P001-T1-DNA1-WES1-4
    |   `-- out
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz.tbi
    |       |-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz.md5
    |       `-- bwa.cnvkit.P001-T1-DNA1-WES1-4.vcf.gz.tbi.md5
    [...]

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_cnv_checking.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import os
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
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
            "reference": self.parent.get_upstream_paths("reference").fasta,
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
            "reference": self.parent.get_upstream_paths("reference").fasta,
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
        if self.config.excluded_regions:
            filenames["excluded_regions"] = self.config.excluded_regions
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
            "reference": self.parent.get_upstream_paths("reference").fasta,
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

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

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
