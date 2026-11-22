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
  Mehari 0.51 reads only transcript databases built by mehari 0.47 or newer (mehari-data-tx 0.14 or newer); the VarFish export wrappers run the same mehari version, so one database task serves both.

VEP annotates one consequence per variant (``--pick`` with ``vep.pick_order``); with ``vep.pick: false`` it annotates all transcripts.
A task that needs all transcripts, such as neoepitope prediction, reads such an annotation task.

VEP runs the plugins in ``vep.plugins``, each a name and the path of its ``.pm`` file; neoepitope prediction needs ``Frameshift`` and ``Wildtype`` from pVACtools.

The task produces one annotated VCF per input library.
Its output has the input's tags plus ``annotated`` and the tool (``vep`` or ``mehari``), so a consumer can require an annotator.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

VEP also writes its summary, ``{library_name}/report/{library_name}.summary.html``, which ``snappy report`` shows.
It is not a result file, and ``vep.output_options: [no_stats]`` switches it off.
VEP's warnings go to the log.

- with ``calls_report``: ``{library_name}/report/calls/{library_name}.calls.tsv``, the calls of the VCF as a table (one row per call), and the same table as a datavzrd site below ``report/calls/html/``, which ``snappy report`` shows.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.common.calls_report import CallsReportStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import Tool
from .model import VariantAnnotation as VariantAnnotationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: The annotated VCF of an entity, without ``.vcf.gz``
_OUT_PREFIX = "work/{library_name}/out/{library_name}"


class VariantAnnotationStepPart(BaseStepPart):
    """Base class for concrete variant annotation tool step parts."""

    actions = ("run",)
    result_keys = tuple(vcf_outputs(_OUT_PREFIX))

    @dictify
    def _get_input_files_run(self, wildcards):
        variant = self.parent.get_upstream_paths("variants", library_name=wildcards.library_name)
        yield "vcf", variant.vcf
        yield "vcf_tbi", variant.vcf_tbi
        yield from self.parent.reference_inputs().items()

    def _get_output_files_run(self):
        return vcf_outputs(_OUT_PREFIX)

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}


class VepStepPart(VariantAnnotationStepPart):
    name = "vep"
    wrappers = {"run": "vep/run"}
    report_outputs = {"run": {"summary": ReportOutput("VEP", caption="vep_summary.rst")}}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        num_threads = self.config.vep.num_threads
        return ResourceUsage(threads=num_threads, runtime="1d", mem=f"{2 * num_threads}GB")

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "cache", self.parent.get_upstream_paths("database").cache
        if self.config.vep.plugins:
            yield "plugins", [plugin.path for plugin in self.config.vep.plugins]

    def _get_output_files_run(self):
        """The annotated VCF and, unless ``output_options`` has ``no_stats``, VEP's summary"""
        outputs = super()._get_output_files_run()
        if "no_stats" not in self.config.vep.output_options:
            outputs["summary"] = "work/{library_name}/report/{library_name}.summary.html"
        return outputs

    def _get_params_run(self, wildcards):
        config = self.config.vep.model_dump(by_alias=True, exclude={"extra_arguments"})
        config["plugins"] = [plugin.name for plugin in self.config.vep.plugins]
        return {"config": config}


class MehariStepPart(VariantAnnotationStepPart):
    name = "mehari"
    wrappers = {"run": "mehari/run"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.mehari.threads, runtime="20m", mem="8GB")

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        database = self.parent.get_upstream_paths("database")
        yield "transcripts", database.transcript_db
        if database.frequency_db:
            yield "frequencies", database.frequency_db
        if database.clinvar_db:
            yield "clinvar", database.clinvar_db

    def _get_params_run(self, wildcards):
        """Pass the options; the files (reference, databases) are inputs."""
        return {"config": self.config.mehari.model_dump(by_alias=True, exclude={"extra_arguments"})}


class VariantAnnotationWorkflow(BaseStep):
    """Annotate germline or somatic variant calls with a single selected tool."""

    name = "variant_annotation"

    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel", "annotated"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "annotated"})),
    ]

    config_model_class = VariantAnnotationConfigModel
    tool_parts = {Tool.vep: VepStepPart, Tool.mehari: MehariStepPart}

    @classmethod
    def task_produces(cls, config, upstream):
        """The variants this task reads, tagged ``annotated``, with the annotator,
        so that a consumer can require one (neoepitope prediction needs VEP),
        and with the genome build of the annotation (``GRCh37``, ``GRCh38``),
        which the calls report reads for its build-dependent linkouts."""
        assembly = getattr(config, str(config.tool)).assembly
        return (upstream["variants"].with_tags("annotated", str(config.tool), assembly),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local annotated VCF output paths for somatic or germline variants."""
        lib = kwargs.get("library_name", "{library_name}")
        paths = vcf_outputs(_OUT_PREFIX.format(library_name=lib))
        return {key: paths[key] for key in ("vcf", "vcf_tbi")}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        if self.config.calls_report:
            self.register_sub_step_classes((CallsReportStepPart,))
