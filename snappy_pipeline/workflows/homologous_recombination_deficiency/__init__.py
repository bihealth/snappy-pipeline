# -*- coding: utf-8 -*-
"""Homologous recombination deficiency score

The step computes the scarHRD score of tumor samples.
The score combines loss of heterozygosity, large-scale transitions and telomeric imbalances.
It is a proxy for the homologous recombination deficiency scores of SNP arrays, as described by `Sztupinszki et al. <https://doi.org/10.1038/s41523-018-0066-6>`_.
The only value of the ``tool`` field is ``scarHRD``.

The task writes one result per tumor DNA library with a matched normal library, as Sequenza does.
It leaves out the other tumor libraries and warns about them.
It reads the ``seqz`` file of the task in ``depends_on.copy_number``, a ``somatic_targeted_seq_cnv_calling`` task with ``tool: sequenza``.
Other CNV callers are rejected.
scarHRD repeats some Sequenza steps with its own arguments, so it needs the Sequenza output.
Decoupling scarHRD from Sequenza would need work, because scarHRD is poorly documented.

scarHRD is not part of CRAN, Bioconductor or bioconda.
It also needs versions of `sequenza <https://bitbucket.org/sequenzatools/sequenza>`_ and `copynumber <https://github.com/aroneklund/copynumber>`_ that are no longer on CRAN or Bioconductor.
When Snakemake creates the conda environment of the wrapper, a post-deploy script installs these R packages at fixed commits.

Outputs
-------

Below ``tasks/<task>/output/``, the task links the result file ``{library_name}/out/{library_name}.hrd.json`` and its ``.md5`` file as symlinks to ``work/``.
Downstream tasks read the key ``json``.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
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
    actions = ("run",)

    wrappers = {"run": "scarHRD/run"}
    result_keys = {"run": ("hrd", "hrd_md5")}
    resource_usage = {"run": ResourceUsage(threads=1, mem="32GB", runtime="24h")}

    @dictify
    def _get_input_files_run(self, wildcards):
        copy_number = self.parent.get_upstream_paths(
            "copy_number", library_name=wildcards.library_name
        )
        yield "seqz", copy_number.seqz

    def _get_output_files_run(self):
        return {
            "hrd": "work/{library_name}/out/{library_name}.hrd.json",
            "hrd_md5": "work/{library_name}/out/{library_name}.hrd.json.md5",
        }

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return self.config.scarHRD.model_dump(by_alias=True)


class HomologousRecombinationDeficiencyWorkflow(BaseStep):
    """Compute Homologous Recombination Deficiency score"""

    #: Step name
    name = "homologous_recombination_deficiency"

    produces = [DataSignature(DataType.TABULAR, frozenset({"hrd"}))]
    config_model_class = HomologousRecombinationDeficiencyConfigModel
    library_requirement = MATCHED_NORMAL_REQUIREMENT
    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HRD output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"json": f"work/{lib}/out/{lib}.hrd.json"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((ScarHRDStepPart,))
