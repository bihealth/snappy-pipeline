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

The task gives one result per pedigree, named after the pedigree (its ``related_group``).
``skip_libraries`` of a tool skips the whole pedigree of a library.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``, each with a ``.md5`` file.
All tools have the same result files.

- ``{library_name}/out/{library_name}.vcf.gz``
- ``{library_name}/out/{library_name}.vcf.gz.tbi``

Here ``library_name`` is the name of the pedigree.
Downstream tasks read the key ``vcf``.
"""

from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart
from snappy_pipeline.workflows.common.sv_calling import SvCallingWorkflow

from .model import SvCallingTargeted as SvCallingTargetedConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class SvCallingTargetedWorkflow(SvCallingWorkflow):
    """Perform germline targeted sequencing CNV calling"""

    #: Workflow name
    name = "sv_calling_targeted"

    config_model_class = SvCallingTargetedConfigModel

    tool_parts = {
        Tool.gcnv: RunGcnvStepPart,
        Tool.delly2: Delly2StepPart,
        Tool.manta: MantaStepPart,
        Tool.melt: MeltStepPart,
    }

    #: Exomes and panels: the callers take the target regions into account
    targeted = True
