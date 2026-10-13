# -*- coding: utf-8 -*-
"""Somatic microsatellite instability calling

The step detects microsatellite instability in cancer samples (whole genomes, exomes or large panels).
The ``tool`` field has the value ``mantis_msi2``, which runs `MANTIS2 <https://github.com/nh13/MANTIS2>`_ on the aligned reads of each tumor and its matched normal.
MANTIS is not maintained anymore, and MANTIS2 is the only supported version.
The task gives one result per tumor DNA library.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``{tumor_library}/out/{tumor_library}.results.txt``
- ``{tumor_library}/out/{tumor_library}.results.txt.status``

Downstream tasks read the key ``results``.
"""

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import SomaticMsiCalling as SomaticMsiCallingConfigModel

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

#: Output keys and extensions of the MANTIS2 results
_OUTPUTS = {
    "result": ".results.txt",
    "status": ".results.txt.status",
    "result_md5": ".results.txt.md5",
    "status_md5": ".results.txt.status.md5",
}


class Mantis2StepPart(BaseStepPart):
    """Perform somatic microsatellite instability with MANTIS_msi2"""

    name = "mantis_msi2"

    actions = ("run",)

    entity_wildcard = "tumor_library"

    published_keys = tuple(_OUTPUTS)

    # Single thread because of https://github.com/OSU-SRLab/MANTIS/issues/57
    default_resource_usage = ResourceUsage(threads=1, runtime="24h", mem=f"{30 * 1024 * 3}MB")

    def _get_input_files_run(self, wildcards):
        # Get names of primary libraries of the selected cancer bio sample and the
        # corresponding primary normal sample
        normal_lib = self.get_normal_lib_name(wildcards)
        tumor_lib = wildcards.tumor_library
        normal: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=normal_lib
        )
        tumor: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=tumor_lib
        )
        return {
            "normal_bam": normal.bam,
            "normal_bai": normal.bai,
            "tumor_bam": tumor.bam,
            "tumor_bai": tumor.bai,
            "reference": self.parent.get_upstream_paths("reference").fasta,
            "loci_bed": self.config.loci_bed,
        }

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.tumor_library]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def _get_output_files_run(self):
        return {
            key: f"work/{{tumor_library}}/out/{{tumor_library}}{ext}"
            for key, ext in _OUTPUTS.items()
        }

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.log"}


class SomaticMsiCallingWorkflow(BaseStep):
    """Perform somatic microsatellite instability analysis"""

    #: Step name
    name = "somatic_msi_calling"
    produces = [DataSignature(DataType.TABULAR, frozenset({"msi"}))]

    config_model_class = SomaticMsiCallingConfigModel

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
        """Return local MSI calling output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"results": f"work/{lib}/out/{lib}.results.txt"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((Mantis2StepPart,))
