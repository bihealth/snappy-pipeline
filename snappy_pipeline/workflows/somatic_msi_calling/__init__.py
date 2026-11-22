# -*- coding: utf-8 -*-
"""Somatic microsatellite instability calling

The step detects microsatellite instability in cancer samples (whole genomes, exomes or large panels).
The ``tool`` field has the value ``mantis_msi2``, which runs `MANTIS2 <https://github.com/nh13/MANTIS2>`_ on the aligned reads of each tumor and its matched normal.
MANTIS is not maintained anymore, and MANTIS2 is the only supported version.
The fields ``min_read_length``, ``min_read_quality``, ``min_locus_quality``, ``min_locus_coverage`` and ``min_repeat_reads`` set the MANTIS thresholds.
Their defaults are the values that the MANTIS authors recommend for exomes.
The microsatellite loci come from the task in ``depends_on.microsatellites``.
The task gives one result per tumor DNA library with a matched normal.
It leaves out the other tumor libraries, with a warning.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.

- ``{tumor_library}/out/{tumor_library}.results.txt``
- ``{tumor_library}/out/{tumor_library}.results.txt.status``

Downstream tasks read the key ``results``.
"""

from snappy_pipeline.models import ToolModel
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import SomaticMsiCalling as SomaticMsiCallingConfigModel

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

#: The path of the MANTIS2 results without suffix
_PREFIX = "work/{library}/out/{library}"

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
    wrappers = {"run": "mantis/mantis_msi2/run"}
    entity_wildcard = "tumor_library"
    result_keys = tuple(_OUTPUTS)

    # Single thread because of https://github.com/OSU-SRLab/MANTIS/issues/57
    default_resource_usage = ResourceUsage(threads=1, runtime="24h", mem=f"{30 * 1024 * 3}MB")

    def tool_section(self) -> ToolModel:
        """The options are on the step model"""
        return self.config

    def _get_input_files_run(self, wildcards):
        return {
            **self.parent.tumor_normal_alignments(wildcards.tumor_library),
            **self.parent.reference_inputs(),
            "loci_bed": self.parent.get_upstream_paths("microsatellites").bed,
        }

    def _get_output_files_run(self):
        prefix = _PREFIX.format(library="{tumor_library}")
        return {key: prefix + ext for key, ext in _OUTPUTS.items()}

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/{tumor_library}.log"}

    def _get_params_run(self, wildcards):
        return {
            name: getattr(self.config, name)
            for name in (
                "min_read_length",
                "min_read_quality",
                "min_locus_quality",
                "min_locus_coverage",
                "min_repeat_reads",
            )
        }


class SomaticMsiCallingWorkflow(BaseStep):
    """Perform somatic microsatellite instability analysis"""

    #: Step name
    name = "somatic_msi_calling"

    produces = [DataSignature(DataType.TABULAR, frozenset({"msi"}))]
    config_model_class = SomaticMsiCallingConfigModel
    library_requirement = MATCHED_NORMAL_REQUIREMENT
    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local MSI calling output paths for downstream consumers."""
        prefix = _PREFIX.format(library=kwargs.get("library_name", "{library_name}"))
        return {"results": prefix + _OUTPUTS["result"]}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((Mantis2StepPart,))
