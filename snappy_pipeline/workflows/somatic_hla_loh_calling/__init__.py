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

The task publishes ``{tumor_library}/out/{tumor_library}.done`` below ``tasks/<task>/output/`` as a symlink to ``work/``.
This marker file is the key ``done`` for downstream tasks.
LOHHLA writes its result files into the same directory in ``work/``.
They are not published.
"""

from snappy_pipeline.models import RelationshipDefinition, ToolModel
from snappy_pipeline.workflows.abstract import MATCHED_NORMAL_REQUIREMENT, BaseStep, BaseStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.hla_typing.model import ExpectedHlaTyping
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import SomaticHlaLohCalling as SomaticHlaLohCallingConfigModel

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"


class LohhlaStepPart(BaseStepPart):
    """Perform LOHHLA analysis of tumor/normal WES/WGS pairs."""

    #: Step name
    name = "lohhla"

    #: Class available actions
    actions = ("run",)

    entity_wildcard = "tumor_library"

    published_keys = ("done",)

    def tool_section(self) -> ToolModel:
        """The options are on the step model"""
        return self.config

    def _get_input_files_run(self, wildcards):
        normal_lib = self.get_normal_lib_name(wildcards)
        tumor_lib = wildcards.tumor_library
        normal: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=normal_lib
        )
        tumor: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=tumor_lib
        )
        hla_typing: ExpectedHlaTyping = self.parent.get_upstream_paths(
            "hla_types", library_name=normal_lib
        )
        return {
            "normal_bam": normal.bam,
            "normal_bai": normal.bai,
            "tumor_bam": tumor.bam,
            "tumor_bai": tumor.bai,
            "hla": hla_typing.txt,
            "hla_dat": self.config.path_hla_dat,
            "hla_fasta": self.config.path_hla_fasta,
            "picard_dir": self.config.path_picard_dir,
        }

    def _get_params_run(self, wildcards):
        return {"tumor_library": wildcards.tumor_library}

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.tumor_library]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def _get_output_files_run(self):
        return {"done": "work/{tumor_library}/out/{tumor_library}.done"}

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.log"}


class SomaticHlaLohCallingWorkflow(BaseStep):
    """Perform somatic hla loh calling"""

    #: Workflow name
    name = "somatic_hla_loh_calling"
    produces = [DataSignature(DataType.TABULAR, frozenset({"hla_loh"}))]

    config_model_class = SomaticHlaLohCallingConfigModel

    library_requirement = MATCHED_NORMAL_REQUIREMENT

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HLA LOH output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"work/{lib}/out/{lib}.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((LohhlaStepPart,))
