"""Germline structural variant calling for targeted sequencing data

The step calls structural variants in exomes and panels from the aligned reads, using depth of coverage, read pair and split read signals.
The ``tool`` field selects the caller:

- ``delly2`` calls the variants per library, merges the sites of a pedigree and genotypes them in every member.
  The regions in ``depends_on.excluded_regions`` are optional and excluded.
- ``manta`` calls the variants of all members of a pedigree together.
- ``gcnv`` runs GATK gCNV in case mode, with the model of the task in ``depends_on.gcnv_model``.
  That task is a :ref:`step_gcnv_model` task with ``intervals: targeted``.
  Its model per library kit and its intervals are the ones of the calls.
- ``melt`` calls mobile element insertions, with the genes of the task in ``depends_on.gene_regions``.

The task gives one result per pedigree, named after its index library.
``skip_libraries`` of a tool skips the whole pedigree of a library.
The step gives no results if the mapper of the alignments is not a DNA mapper.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
All tools publish the same files.

- ``{library_name}/out/{library_name}.vcf.gz``
- ``{library_name}/out/{library_name}.vcf.gz.tbi``

Here ``library_name`` is the index library of the pedigree.
Downstream tasks read the key ``vcf``.
"""

from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart

from .model import SvCallingTargeted as SvCallingTargetedConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class SvCallingTargetedWorkflow(BaseStep):
    """Perform germline targeted sequencing CNV calling"""

    #: Workflow name
    name = "sv_calling_targeted"
    #: The pedigrees of background data sets are called, too
    uses_background = True

    config_model_class = SvCallingTargetedConfigModel
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"}))]

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        selected_tool = self.config.tool
        match selected_tool:
            case Tool.gcnv:
                selected_sub_step = RunGcnvStepPart
            case Tool.delly2:
                selected_sub_step = Delly2StepPart
            case Tool.manta:
                selected_sub_step = MantaStepPart
            case Tool.melt:
                selected_sub_step = MeltStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {selected_tool}")
        # Register only the selected tool's step part.
        self.register_sub_step_classes((WritePedigreeStepPart, selected_sub_step))

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local targeted SV output paths for downstream consumers."""
        library_name = kwargs.get("library_name", "{library_name}")
        return {"vcf": f"work/{library_name}/out/{library_name}.vcf.gz"}
