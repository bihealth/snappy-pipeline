"""Mobile element insertion calling for targeted sequencing

The step finds mobile element insertions (MEI) in the alignments from ``depends_on.alignments`` (:ref:`step_ngs_mapping`).
The only value of ``tool`` is ``scramble``, which runs Scramble.
The step writes one VCF per library with a DNA library.
It is experimental.

Outputs
-------

The task publishes ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi`` and ``.md5`` files below ``tasks/<task>/output/``.
The paths are the same as below ``work/``.
The key ``vcf`` is the one downstream tasks read through ``depends_on``.
"""

import os

from snakemake.io import expand

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import TargetedSeqMeiCalling as TargetedSeqMeiCallingConfigModel

#: Extensions of files to create as main payload.
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")


#: Default configuration for the targeted_seq_mei_calling step.


class ScrambleStepPart(BaseStepPart):
    """Mobile element insertion detection using GeneDx::scramble"""

    #: Step name
    name = "scramble"

    #: Class available actions
    actions = ("cluster", "analysis")

    intermediate_outputs = {"cluster": ("txt",)}

    default_resource_usage = ResourceUsage(threads=1, runtime="6h", mem=f"{8 * 1024}MB")

    def __init__(self, parent):
        super().__init__(parent)

    def _get_log_file(self, action):
        return f"work/{{library_name}}/log/{{library_name}}.{action}.log"

    @listify
    def _get_input_files_cluster(self, wildcards):
        """Yield BAM files based on subworkflow `ngs_mapping` results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield alignments.bam

    @dictify
    def _get_input_files_analysis(self, wildcards):
        """Yield input files' pattern for rule `annotate` - based on scramble call results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        name_pattern = "{library_name}"
        base_name_out = "work/{name_pattern}/out/{name_pattern}.cluster.{ext}".format(
            name_pattern=name_pattern, ext="txt"
        )
        yield "cluster", base_name_out.format(**wildcards)
        yield "reference_genome", self.config.scramble.blast_ref
        if self.config.scramble.mei_refs:
            yield "mei_refs", self.config.scramble.mei_refs

    def _get_output_files_cluster(self):
        """Return output files' patterns for scramble cluster call."""
        return {"txt": "work/{library_name}/out/{library_name}.cluster.txt"}

    @dictify
    def _get_output_files_analysis(self):
        """Yield output files' patterns for scramble call."""
        ext_dict = {
            "txt": ".MEIs.txt",
            "vcf": ".vcf",
            "vcf_gz": ".vcf.gz",
            "vcf_gz_md5": ".vcf.gz.md5",
            "vcf_tbi": ".vcf.gz.tbi",
            "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in ext_dict.items():
            yield key, f"work/{{library_name}}/out/{{library_name}}{ext}"

    def _get_params_analysis(self, wildcards):
        """Get parameters.

        :param wildcards: Snakemake rule wildcards (unused).
        :type wildcards: snakemake.io.Wildcards

        :return: Returns parameters required to run analysis part of scramble.

        :raises InvalidConfiguration: if information provided in configuration isn't enough to run
        the analysis.
        """
        blast_ref_path = self.config.scramble.blast_ref
        try:
            if not os.path.isfile(blast_ref_path):
                raise InvalidConfiguration(
                    f"Provided path to reference genome ('blast_ref') "
                    f"is not a real file: {blast_ref_path}"
                )
        except TypeError as e:
            raise TypeError("Path to reference genome ('blast_ref') cannot be empty.") from e
        params = {
            "n_cluster": self.config.scramble.n_cluster,
            "mei_score": self.config.scramble.mei_score,
            "indel_score": self.config.scramble.indel_score,
            "mei_polya_frac": self.config.scramble.mei_polya_frac,
        }
        return params


class MeiWorkflow(BaseStep):
    """Perform germline mobile element insertion detection."""

    #: Workflow name
    name = "targeted_seq_mei_calling"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "mei"}))]
    config_model_class = TargetedSeqMeiCallingConfigModel

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((ScrambleStepPart,))

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local MEI calling output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"vcf": f"work/{lib}/out/{lib}.vcf.gz"}

    @listify
    def get_result_files(self):
        """Get list of result files.

        :return: Return list of result files for the germline mobile element insertion
        detection workflow.
        """
        # Initialise variable
        name_pattern = "{library_name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for library_name in self.primary_dna_libraries()["library_name"]:
            yield from expand(tpl, library_name=[library_name], **kwargs)
