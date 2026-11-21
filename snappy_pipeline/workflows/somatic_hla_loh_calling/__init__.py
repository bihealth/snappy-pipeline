# -*- coding: utf-8 -*-
"""Somatic HLA loss of heterozygosity calling

The step detects loss of heterozygosity at the HLA locus in cancer samples (whole genomes, exomes or large panels) with LOHHLA.
LOHHLA compares the aligned reads of the tumor and its matched normal, using the HLA types of the normal from the ``hla_typing`` task in ``depends_on.hla_types``.

The step has no ``tool`` field.
It gives one result per tumor DNA library, paired with its matched normal.
It leaves out tumor libraries without a matched normal, with a warning.
The tumor purity and ploidy estimates of LOHHLA are not used.

Outputs
-------

The task links the result directory ``{tumor_library}/out/`` below ``tasks/<task>/output/`` as a symlink to ``work/``.
It holds the files that LOHHLA writes, and it is the key ``out`` for downstream tasks.
"""

from snakemake.io import directory

from snappy_pipeline.models import ToolModel
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ReportOutput,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.hla_typing.model import ExpectedHlaTyping

from .model import SomaticHlaLohCalling as SomaticHlaLohCallingConfigModel

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"


class LohhlaStepPart(BaseStepPart):
    """Perform LOHHLA analysis of tumor/normal WES/WGS pairs."""

    #: Step name
    name = "lohhla"

    #: Class available actions
    actions = ("run",)
    wrappers = {"run": "lohhla"}

    entity_wildcard = "tumor_library"

    result_keys = ("out",)

    def tool_section(self) -> ToolModel:
        """The options are on the step model"""
        return self.config

    def _get_input_files_run(self, wildcards):
        tumor = wildcards.tumor_library
        hla_typing: ExpectedHlaTyping = self.parent.get_upstream_paths(
            "hla_types", library_name=self.parent.matched_normal(tumor)
        )
        return {
            **self.parent.tumor_normal_alignments(tumor),
            "hla": hla_typing.txt,
            "hla_dat": self.config.path_hla_dat,
            "hla_fasta": self.config.path_hla_fasta,
            "picard_dir": self.config.path_picard_dir,
        }

    def _get_params_run(self, wildcards):
        return {
            "tumor_library": wildcards.tumor_library,
            "min_coverage_filter": self.config.min_coverage_filter,
        }

    #: The plots of each HLA gene in LOHHLA's ``Figures`` directory; it writes none when all
    #: HLA alleles are homozygous
    report_outputs = {
        "run": {
            "out": ReportOutput("LOHHLA", caption="lohhla.rst", patterns=("Figures/{plot}.pdf",))
        }
    }

    def _get_output_files_run(self):
        # LOHHLA writes many files that the step does not name
        return {"out": directory("work/{tumor_library}/out")}

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.log"}


class SomaticHlaLohCallingWorkflow(BaseStep):
    """Perform somatic hla loh calling"""

    #: Workflow name
    name = "somatic_hla_loh_calling"
    produces = [DataSignature(DataType.TABULAR, frozenset({"hla_loh"}))]

    config_model_class = SomaticHlaLohCallingConfigModel

    library_requirement = MATCHED_NORMAL_REQUIREMENT

    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HLA LOH output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"out": f"work/{lib}/out"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((LohhlaStepPart,))
