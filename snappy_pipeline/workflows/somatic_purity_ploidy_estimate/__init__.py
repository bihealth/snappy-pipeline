# -*- coding: utf-8 -*-
"""Purity and ploidy estimation for somatic samples

The step estimates the tumor cell content (purity) and the ploidy of tumor samples from the B-allele frequencies and the coverage of the tumor and its matched normal.
The ``tool`` field has the value ``ascat``, which runs ASCAT.
The B-allele frequencies are taken at the loci of ``ascat.b_af_loci``.
The task gives one result per tumor DNA library.

Whole exome data is not supported yet.
The coverage of the exome tumor and normal would have to come from a CNV caller that writes bins, and no such caller exists.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``{tumor_library}/out/{tumor_library}.ploidy.txt``: the estimated ploidy.
- ``{tumor_library}/out/{tumor_library}.goodness_of_fit.txt``
- ``{tumor_library}/out/{tumor_library}.segments.txt`` and ``{tumor_library}/out/{tumor_library}.segments_raw.txt``: the fitted and the raw copy number segments.

Downstream tasks read the key ``done``, the marker file ``{tumor_library}/out/.done`` in ``work/``.
The B-allele frequency and coverage files per library are intermediate files.
The purity (``{tumor_library}.aberrantcellfraction.txt``), the ASCAT plots and other ASCAT results are written to the same directory in ``work/``.
They are not published.
"""

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart
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

    # Name of the step.
    name = "ascat"

    # The actions for generating BAF/CNV files for the tumor and/nor normal
    # sample, and finally to run the ASCAT pipeline.
    actions = (
        "baf_tumor",
        "baf_normal",
        "cnv_tumor",
        "cnv_normal",
        "cnv_tumor_wes",
        "cnv_normal_wes",
        "run_ascat",
    )

    intermediate_outputs = {
        action: ("txt",) for action in ("baf_tumor", "baf_normal", "cnv_tumor", "cnv_normal")
    }

    entity_wildcard = "tumor_library"

    published_keys = {
        "run_ascat": (
            "goodness_of_fit",
            "goodness_of_fit_md5",
            "ploidy",
            "ploidy_md5",
            "segments",
            "segments_md5",
            "segments_raw",
            "segments_raw_md5",
        )
    }

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
            "reference": self.parent.get_upstream_paths("reference").fasta,
            "b_af_loci": self.config.ascat.b_af_loci,
        }

    def _get_input_files_baf_normal(self, wildcards):
        """Return input files for generating BAF file for the normal."""
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.normal_library
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            "reference": self.parent.get_upstream_paths("reference").fasta,
            "b_af_loci": self.config.ascat.b_af_loci,
        }

    def _get_input_files_cnv_tumor(self, wildcards):
        """Return input files for generating BAF file for the tumor."""
        return self._get_input_files_baf_tumor(wildcards)

    def _get_input_files_cnv_normal(self, wildcards):
        """Return input files for generating CNV file for the normal."""
        return self._get_input_files_baf_normal(wildcards)

    def _get_input_files_cnv_tumor_wes(self, wildcards):
        """Return input files for generating CNV file from copywriter for tumor."""
        raise NotImplementedError(
            "Copywriter-based WES CNV input is not available: no somatic_targeted_seq_cnv_calling "
            "tool writes copywriter bins."
        )

    def _get_input_files_cnv_normal_wes(self, wildcards):
        """Return input files for generating CNV file from copywriter for normal."""
        raise NotImplementedError(
            "Copywriter-based WES CNV input is not available: no somatic_targeted_seq_cnv_calling "
            "tool writes copywriter bins."
        )

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
        yield "done", touch("work/{tumor_library}/out/.done")
        infixes = ("goodness_of_fit", "ploidy", "segments", "segments_raw")
        for infix in infixes:
            path = ("work/{tumor_library}/out/{tumor_library}.%s.txt") % infix
            yield infix, path
            yield infix + "_md5", path + ".md5"

    def _get_params_cnv_tumor(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"tumor_library": wildcards.tumor_library}

    def _get_params_run_ascat(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"tumor_library": wildcards.tumor_library}

    def _get_log_file(self, action):
        """Return path to log file"""
        # TODO: implement log option for actions `cnv_tumor_wes` and `cnv_normal_wes`.
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
        """Return local purity/ploidy output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"work/{lib}/out/.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((AscatStepPart,))
