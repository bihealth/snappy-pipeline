# -*- coding: utf-8 -*-
"""Variant annotation with VEP or Mehari

Annotates the VCF files of the task named in ``depends_on.variants``.
The input can be somatic or germline.
The ``tool`` field selects the annotator: ``vep`` runs Ensembl VEP, ``mehari`` runs Mehari.
One task runs one tool, so a second annotation needs a second task.

The tool's database comes from ``depends_on.database``, an ``external_data`` task of type ``database`` tagged with the tool.
Its ``files`` give the keys below.

- ``vep``: ``cache`` (the VEP cache directory).
- ``mehari``: ``transcript_db`` (one path or a list), and optionally ``frequency_db`` and ``clinvar_db``.

The task produces one annotated VCF per input library.
Its output has the input's tags plus ``annotated``.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ExpectedVariantVcf
from .model import VariantAnnotation as VariantAnnotationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Output keys and extensions of the annotated VCF
_OUTPUTS = {
    "vcf": ".vcf.gz",
    "vcf_tbi": ".vcf.gz.tbi",
    "vcf_md5": ".vcf.gz.md5",
    "vcf_tbi_md5": ".vcf.gz.tbi.md5",
}


class VariantAnnotationStepPart(BaseStepPart):
    """Base class for concrete variant annotation tool step parts."""

    actions = ("run",)

    published_keys = tuple(_OUTPUTS)

    @dictify
    def _get_input_files_run(self, wildcards):
        lib = wildcards.library_name
        variant = ExpectedVariantVcf.model_validate(
            self.parent.get_upstream_paths("variants", library_name=lib)
        )
        yield "vcf", variant.vcf
        yield "vcf_tbi", variant.vcf_tbi
        yield from self.parent.reference_inputs().items()

    def _get_output_files_run(self):
        return {
            key: f"work/{{library_name}}/out/{{library_name}}{ext}" for key, ext in _OUTPUTS.items()
        }

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards):
        return {"config": getattr(self.config, self.name).model_dump(by_alias=True)}


class VepStepPart(VariantAnnotationStepPart):
    name = "vep"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "cache", self.parent.get_upstream_paths("database").cache

    def _get_params_run(self, wildcards):
        return {"config": self.config.vep.model_dump(by_alias=True)}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        num_threads = self.config.vep.num_threads
        return ResourceUsage(threads=num_threads, runtime="1d", mem=f"{2 * num_threads}GB")


class MehariStepPart(VariantAnnotationStepPart):
    name = "mehari"

    def _get_params_run(self, wildcards):
        """Pass the options; the files (reference, databases) are inputs."""
        return {"config": self.config.mehari.model_dump(by_alias=True)}

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        database = self.parent.get_upstream_paths("database")
        yield "transcripts", database.transcript_db
        if database.frequency_db:
            yield "frequencies", database.frequency_db
        if database.clinvar_db:
            yield "clinvar", database.clinvar_db

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.mehari.threads, runtime="20m", mem="8GB")


_TOOL_STEP_PART = {
    "vep": VepStepPart,
    "mehari": MehariStepPart,
}


class VariantAnnotationWorkflow(BaseStep):
    """Annotate germline or somatic variant calls with a single selected tool."""

    name = "variant_annotation"
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel", "annotated"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "annotated"})),
    ]

    @classmethod
    def task_produces(cls, config, upstream):
        """The variants this task reads, tagged ``annotated``."""
        return (upstream["variants"].with_tags("annotated"),)

    config_model_class = VariantAnnotationConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local annotated VCF output paths for somatic or germline variants."""
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "vcf": f"work/{lib}/out/{lib}.vcf.gz",
            "vcf_tbi": f"work/{lib}/out/{lib}.vcf.gz.tbi",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        tool_cls = _TOOL_STEP_PART[self.config.tool]
        self.register_sub_step_classes((tool_cls,))
