"""Germline structural variant calling for targeted sequencing data

The step calls structural variants in exomes and panels from the aligned reads, using depth of coverage, read pair and split read signals.
The ``tool`` field selects the caller:

- ``delly2`` calls the variants per library, merges the sites of a pedigree and genotypes them in every member.
- ``manta`` calls the variants of all members of a pedigree together.
- ``gcnv`` runs GATK gCNV in case mode, with the models in ``gcnv.precomputed_model_paths``.
  The models come from the step :ref:`step_helper_gcnv_model_targeted`.
  ``gcnv.path_target_interval_list_mapping`` maps the library kit to its target intervals.
- ``melt`` calls mobile element insertions.

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

from collections import Counter

from biomedsheets.shortcuts import GermlineCaseSheet

from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_common import build_ngs_library_to_kit
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart

from .model import SvCallingTargeted as SvCallingTargetedConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Minimum number of samples per kit to apply gCNV calling criteria to be analyzed
GCNV_MIN_KIT_SAMPLES = 10

#: Default configuration for the sv_calling_targeted step


class GcnvTargetedStepPart(RunGcnvStepPart):
    """Targeted seq. CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self.parent.ngs_library_to_kit


class SvCallingTargetedWorkflow(BaseStep):
    """Perform germline targeted sequencing CNV calling"""

    #: Workflow name
    name = "sv_calling_targeted"
    config_model_class = SvCallingTargetedConfigModel
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"}))]

    sheet_shortcut_class = GermlineCaseSheet

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        selected_tool = self.config.tool
        # gCNV-specific shortcuts must be initialized BEFORE registering sub-step classes
        # so that GcnvTargetedStepPart.__init__ can access them.
        if selected_tool == Tool.gcnv:
            self.ngs_library_to_kit = build_ngs_library_to_kit(
                self.primary_dna_libraries(include_background=True),
                self.get_task_config("sv_calling_targeted").gcnv.path_target_interval_list_mapping,
            )
            self.library_kit_counts_dict = dict(Counter(self.ngs_library_to_kit.values()))
        else:
            self.ngs_library_to_kit = {}
            self.library_kit_counts_dict = {}

        match selected_tool:
            case Tool.gcnv:
                selected_sub_step = GcnvTargetedStepPart
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

    def get_library_count(self, library_kit):
        """Get library count.

        :param library_kit: Library kit name.
        :type library_kit: str

        :return: Returns number of samples with inputted library kit. If library name not defined,
        it returns zero.
        """
        return self.library_kit_counts_dict.get(library_kit, 0)
