# -*- coding: utf-8 -*-
"""Purity and ploidy estimation for somatic samples

The step estimates the tumor cell content (purity) and the ploidy of tumor samples from the B-allele frequencies and the coverage of the tumor and its matched normal.
The ``tool`` field has the value ``ascat``, which runs ASCAT.
The B-allele frequencies are taken at the loci of the task in ``depends_on.snp_loci``.
The task gives one result per tumor DNA library with a matched normal.
It leaves out the other tumor libraries, with a warning.

Whole exome data is not supported yet.
The coverage of the exome tumor and normal would have to come from a CNV caller that writes bins, and no such caller exists.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``{tumor_library}/out/{tumor_library}.ploidy.txt``: the estimated ploidy.
- ``{tumor_library}/out/{tumor_library}.aberrantcellfraction.txt``: the estimated purity.
- ``{tumor_library}/out/{tumor_library}.goodness_of_fit.txt``
- ``{tumor_library}/out/{tumor_library}.segments.txt`` and ``{tumor_library}/out/{tumor_library}.segments_raw.txt``: the fitted and the raw copy number segments.
- ``{tumor_library}/out/{tumor_library}.tumour.png``, ``.germline.png`` and ``.ASPCF.png``: the raw tumor and germline data, and the segmentation.

Downstream tasks read the keys ``ploidy``, ``purity``, ``goodness_of_fit``, ``segments`` and ``segments_raw``.
The B-allele frequency and coverage files per library are intermediate files.
ASCAT also writes the profile plots (``.ASCATprofile.png``, ``.rawprofile.png`` and ``.sunrise.png``) when it finds a solution.
They are not kept.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    MATCHED_NORMAL_REQUIREMENT,
    ReportOutput,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import ResourceUsage

from .model import SomaticPurityPloidyEstimate as SomaticPurityPloidyEstimateConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Tools for estimating purity and ploidy.
PURITY_PLOIDY_TOOLS = "ascat"


class AscatStepPart(BaseStepPart):
    """Estimation of purity and ploidy using ASCAT.

    Notes
    -----

    - Although we name the virtual probes "SNP${num}", they are not guaranteed
      to
    """

    report_outputs = {
        "run_ascat": {
            key: ReportOutput("ASCAT", caption="ascat_plots.rst")
            for key in ("plot_tumour", "plot_germline", "plot_aspcf")
        }
    }

    # Name of the step.
    name = "ascat"

    # The actions for generating BAF/CNV files for the tumor and/nor normal
    # sample, and finally to run the ASCAT pipeline.
    actions = (
        "baf_tumor",
        "baf_normal",
        "cnv_tumor",
        "cnv_normal",
        "run_ascat",
    )

    intermediate_outputs = {
        action: ("txt",) for action in ("baf_tumor", "baf_normal", "cnv_tumor", "cnv_normal")
    }

    entity_wildcard = "tumor_library"

    #: Files of ``run_ascat``: output key and suffix
    result_suffixes = {
        "goodness_of_fit": "goodness_of_fit.txt",
        "ploidy": "ploidy.txt",
        "purity": "aberrantcellfraction.txt",
        "segments": "segments.txt",
        "segments_raw": "segments_raw.txt",
        "plot_tumour": "tumour.png",
        "plot_germline": "germline.png",
        "plot_aspcf": "ASPCF.png",
    }

    result_keys = {"run_ascat": tuple(k + md5 for k in result_suffixes for md5 in ("", "_md5"))}

    default_resource_usage = ResourceUsage(threads=8, runtime="2d", mem=f"{10 * 1024 * 8}MB")

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.tumor_library]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def _get_input_files_baf_tumor(self, wildcards):
        """Return input files for generating BAF file for the tumor."""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.tumor_library
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
            "b_af_loci": self.parent.get_upstream_paths("snp_loci").bed,
        }

    def _get_input_files_baf_normal(self, wildcards):
        """Return input files for generating BAF file for the normal."""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.normal_library
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
            "b_af_loci": self.parent.get_upstream_paths("snp_loci").bed,
        }

    def _get_input_files_cnv_tumor(self, wildcards):
        """Return input files for generating BAF file for the tumor."""
        return self._get_input_files_baf_tumor(wildcards)

    def _get_input_files_cnv_normal(self, wildcards):
        """Return input files for generating CNV file for the normal."""
        return self._get_input_files_baf_normal(wildcards)

    @dictify
    def _get_input_files_run_ascat(self, wildcards):
        """Return input files for actually running ASCAT."""
        result = {
            "baf_tumor": "work/{tumor_library}/out/{tumor_library}.baf_tumor.txt",
            "baf_normal": "work/{normal_library}/out/{normal_library}.baf_normal.txt",
            "cnv_tumor": "work/{tumor_library}/out/{tumor_library}.cnv_tumor.txt",
            "cnv_normal": "work/{normal_library}/out/{normal_library}.cnv_normal.txt",
        }
        normal_library = self.get_normal_lib_name(wildcards)
        for key, value in result.items():
            yield key, value.format(normal_library=normal_library, **wildcards)

    @staticmethod
    def _get_output_files_baf_tumor():
        """Return output files for generating BAF file for the tumor."""
        return {"txt": "work/{tumor_library}/out/{tumor_library}.baf_tumor.txt"}

    @staticmethod
    def _get_output_files_baf_normal():
        """Return output files for generating BAF file for the normal."""
        return {"txt": "work/{normal_library}/out/{normal_library}.baf_normal.txt"}

    @staticmethod
    def _get_output_files_cnv_tumor():
        """Return output files for generating BAF file for the tumor."""
        return {"txt": "work/{tumor_library}/out/{tumor_library}.cnv_tumor.txt"}

    @staticmethod
    def _get_output_files_cnv_normal():
        """Return output files for generating CNV file for the normal."""
        return {"txt": "work/{normal_library}/out/{normal_library}.cnv_normal.txt"}

    @dictify
    def _get_output_files_run_ascat(self):
        """Return output files for actually running ASCAT."""
        for key, suffix in self.result_suffixes.items():
            path = f"work/{{tumor_library}}/out/{{tumor_library}}.{suffix}"
            yield key, path
            yield key + "_md5", path + ".md5"

    def _get_params_run_ascat(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"tumor_library": wildcards.tumor_library}

    def _get_log_file(self, action):
        """Return path to log file"""
        log_dict = {
            "baf_tumor": "work/{tumor_library}/log/{tumor_library}.baf_tumor.log",
            "baf_normal": "work/{normal_library}/log/{normal_library}.baf_normal.log",
            "cnv_tumor": "work/{tumor_library}/log/{tumor_library}.cnv_tumor.log",
            "cnv_normal": "work/{normal_library}/log/{normal_library}.cnv_normal.log",
            "run_ascat": "work/{tumor_library}/log/{tumor_library}.run.log",
        }
        return {"log": log_dict[action]}


class SomaticPurityPloidyEstimateWorkflow(BaseStep):
    """Perform purity and ploidy estimation"""

    #: Workflow name
    name = "somatic_purity_ploidy_estimate"
    produces = [DataSignature(DataType.TABULAR, frozenset({"purity_ploidy"}))]

    config_model_class = SomaticPurityPloidyEstimateConfigModel

    library_requirement = MATCHED_NORMAL_REQUIREMENT

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local purity/ploidy output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {
            key: f"work/{lib}/out/{lib}.{suffix}"
            for key, suffix in AscatStepPart.result_suffixes.items()
            if not key.startswith("plot_")
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((AscatStepPart,))
