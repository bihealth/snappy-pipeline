# -*- coding: utf-8 -*-
"""Homologous recombination deficiency score

The step computes the scarHRD score of tumor samples.
The score combines loss of heterozygosity, large-scale transitions and telomeric imbalances.
It is a proxy for the homologous recombination deficiency scores of SNP arrays, as described by `Sztupinszki et al. <https://doi.org/10.1038/s41523-018-0066-6>`_.
The ``tool`` field has the value ``scarHRD``.

The task gives one result per tumor DNA library with a matched normal, as Sequenza does.
It leaves out the other tumor libraries, with a warning.
It reads the ``seqz`` file of a ``somatic_targeted_seq_cnv_calling`` task with ``tool: sequenza`` (``depends_on.copy_number``).
Other CNV callers are rejected.
scarHRD repeats some Sequenza steps with its own arguments.
This is why the Sequenza output is required.
Decoupling scarHRD from Sequenza would need work, because scarHRD is poorly documented.

scarHRD is not part of CRAN, Bioconductor or bioconda.
It also needs versions of `sequenza <https://github.com/oicr-gsi/sequenza/tree/master>`_ and `copynumber <https://github.com/aroneklund/copynumber>`_ that are no longer on CRAN or Bioconductor.
The task installs these R packages in a job of its own.

Outputs
-------

The task publishes ``{library_name}/out/{library_name}.hrd.json`` with a ``.md5`` file below ``tasks/<task>/output/`` as a symlink to ``work/``.
Downstream tasks read the key ``json``.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import HomologousRecombinationDeficiency as HomologousRecombinationDeficiencyConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


class ScarHRDStepPart(BaseStepPart):
    """Computes homologous recombination deficiency score with scarHRD"""

    #: Step name
    name = "scarHRD"

    #: Class available actions
    actions = (
        "install",
        "run",
    )

    published_keys = {"run": ("hrd", "hrd_md5")}

    resource_usage = {"run": ResourceUsage(threads=1, mem="32GB", runtime="24h")}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "done", "work/R_packages/out/scarHRD.done"
        copy_number = self.parent.get_upstream_paths(
            "copy_number", library_name=wildcards.library_name
        )
        yield "seqz", copy_number.seqz
        yield from self.parent.reference_inputs().items()

    def _get_output_files_install(self):
        return {"done": "work/R_packages/out/scarHRD.done"}

    def _get_output_files_run(self):
        return {
            "hrd": "work/{library_name}/out/{library_name}.hrd.json",
            "hrd_md5": "work/{library_name}/out/{library_name}.hrd.json.md5",
        }

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return self.config.scarHRD.model_dump(by_alias=True)

    def _get_log_file_install(self):
        return {"log": "work/R_packages/log/scarHRD.log"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}


class HomologousRecombinationDeficiencyWorkflow(BaseStep):
    """Compute Homologous Recombination Deficiency score"""

    #: Step name
    name = "homologous_recombination_deficiency"
    produces = [DataSignature(DataType.TABULAR, frozenset({"hrd"}))]

    config_model_class = HomologousRecombinationDeficiencyConfigModel

    library_requirement = MATCHED_NORMAL_REQUIREMENT

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HRD output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"json": f"work/{lib}/out/{lib}.hrd.json"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((ScarHRDStepPart,))

    def check_config(self):
        """Check that the upstream CNV calling tool is supported.

        ``scarHRD`` can currently only consume ``sequenza`` copy number output.
        """
        tool = self.get_task_config("copy_number").tool
        if tool != "sequenza":
            raise InvalidConfiguration(
                f"Tool '{tool}' of the depends_on.copy_number task "
                f"'{self.project.dependencies[self.task_name]['copy_number']}' is not supported by "
                "'homologous_recombination_deficiency'; expected 'sequenza'"
            )
