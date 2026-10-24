"""Processing of the alignments of a mapping task

Processes the BAM files of the task named in ``depends_on.alignments``, an ``ngs_mapping`` task or another ``alignment_processing`` task.
The ``tool`` field selects the processing, and one task applies one.
Chain several tasks for a sequence, for example duplicate marking before BQSR.
Downstream tasks read the BAM files of the last task of the chain, through the keys ``bam`` and ``bai`` as from a mapping task.

``gatk4_bqsr``
    Base quality score recalibration (BQSR) with GATK ``BaseRecalibrator`` and ``ApplyBQSR``.
    ``BaseRecalibrator`` models the errors of the base qualities per read group, by quality score, cycle and sequence context, and ignores the known variant sites of ``depends_on.common_variants`` (for example ``small_exac_common_3.vcf`` of the GATK resource bundle) and reads marked as duplicates.
    ``ApplyBQSR`` writes the BAM file with the recalibrated base qualities.

The task processes the DNA libraries of its upstream task, those of background data sets included, unless ``library_selection`` narrows them.
Its output has the input's tags plus ``recalibrated`` (``gatk4_bqsr``).
With ``delete_intermediates: [bam]``, the BAM file of the task is temporary, for a later task that processes it further (see :ref:`usage`).

Outputs
-------

The task links the following result files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.
Each file has an ``.md5`` file.

- ``{library_name}/out/{library_name}.bam`` and ``{library_name}/out/{library_name}.bam.bai``
- ``gatk4_bqsr``: ``{library_name}/report/{library_name}.recal_table.txt``, the recalibration table

Downstream tasks read the keys ``bam`` and ``bai``.
"""

from snappy_pipeline.models import Intermediate
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ReportOutput, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import AlignmentProcessing as AlignmentProcessingConfigModel
from .model import Tool

#: The tag that each tool adds to the tags of its input
TOOL_TAGS = {Tool.gatk4_bqsr: "recalibrated"}

#: The processed BAM file of a library
BAM = "work/{library_name}/out/{library_name}.bam"


class AlignmentProcessingStepPart(BaseStepPart):
    """The BAM file of each library, processed by one tool"""

    actions = ("run",)

    #: The BAM, which a later task may process further
    bam_outputs = {"run": ("bam", "bam_bai")}

    def _get_input_files_run(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        return {"bam": alignments.bam, "bai": alignments.bai}

    def _get_output_files_run(self):
        outputs = {"bam": BAM, "bam_bai": BAM + ".bai"}
        if Intermediate.bam not in self.config.delete_intermediates:
            outputs |= {"bam_md5": BAM + ".md5", "bam_bai_md5": BAM + ".bai.md5"}
        return outputs


class Gatk4BqsrStepPart(AlignmentProcessingStepPart):
    """Base quality score recalibration with GATK"""

    name = "gatk4_bqsr"

    result_keys = ("bam", "bam_bai", "bam_md5", "bam_bai_md5", "recal_table", "recal_table_md5")

    report_outputs = {"run": {"recal_table": ReportOutput("BQSR", caption="recal_table.rst")}}

    default_resource_usage = ResourceUsage(threads=1, runtime="24h", mem="16GB")

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()
        yield from self.parent.vcf_inputs("common_variants").items()

    def _get_output_files_run(self):
        table = "work/{library_name}/report/{library_name}.recal_table.txt"
        return super()._get_output_files_run() | {
            "recal_table": table,
            "recal_table_md5": table + ".md5",
        }

    def _get_params_run(self, wildcards):
        return {}  # get_params adds the extra_arguments

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.recalibration.log"}


class AlignmentProcessingWorkflow(BaseStep):
    """Process the alignments of a mapping task"""

    name = "alignment_processing"

    #: The DNA libraries, background included
    library_requirement = "extraction_type == 'dna'"
    library_default = ""
    background_results = True

    produces = [DataSignature(DataType.ALIGNMENTS, frozenset(TOOL_TAGS.values()))]

    config_model_class = AlignmentProcessingConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """The alignments this task reads, with the tag of its tool."""
        return (upstream["alignments"].with_tags(TOOL_TAGS[config.tool]),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        lib = kwargs.get("library_name", "{library_name}")
        bam = BAM.format(library_name=lib)
        return {"bam": bam, "bai": bam + ".bai"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        tool_to_class = {Tool.gatk4_bqsr: Gatk4BqsrStepPart}
        self.register_sub_step_classes((tool_to_class[self.config.tool],))
