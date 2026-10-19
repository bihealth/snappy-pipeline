"""Repeat expansion analysis with ExpansionHunter

The step runs the tool selected by ``tool`` on the alignments from ``depends_on.alignments`` (:ref:`step_ngs_mapping`).
The only value is ``expansionhunter``, which calls repeat expansions against the key ``catalog`` of ``depends_on.database`` and annotates them with its key ``annotation``.
The database is an ``external_data`` task with ``{type: database, tags: [expansionhunter]}``.
The step writes one result set per library with a DNA library.
It is experimental.

Outputs
-------

The task publishes these files below ``tasks/<task>/output/``, each with a ``.md5`` file.
The paths are the same as below ``work/``.

- ``{library_name}/out/{library_name}.vcf``, the repeat expansion calls
- ``annotated.{library_name}/out/annotated.{library_name}.json``, the calls with the annotation

The key ``vcf`` is the one downstream tasks read through ``depends_on``.
"""

import os
from functools import cached_property

from biomedsheets.shortcuts import GermlineCaseSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.repeat_expansion.annotate_expansionhunter import (
    AnnotateExpansionHunter,
)

from .model import RepeatExpansion as RepeatExpansionConfigModel

#: Extensions of files to create as main payload - JSON.
EXT_JSON = (".json", ".json.md5")
#: Extensions of files to create as main payload - VCF.
EXT_VCF = (".vcf", ".vcf.md5")
#: Default configuration for the repeat_expansion step.


class ExpansionHunterStepPart(BaseStepPart):
    """Repeat expansion analysis with Illumina::ExpansionHunter"""

    #: Step name
    name = "expansionhunter"

    #: Valid actions
    actions = ("run", "annotate")

    intermediate_outputs = {"run": ("json",)}

    @cached_property
    def library_name_to_sex(self) -> dict[str, str]:
        """Return the sex of each library's donor, including background data sets."""
        df = self.parent.build_library_dataframe(include_background=True)
        return dict(zip(df["library_name"], df["sex"]))

    def _get_input_files_run(self, wildcards):
        """Yield BAM files based on subworkflow `ngs_mapping` results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            "reference": self.parent.get_upstream_paths("reference").fasta,
            "repeat_catalog": self.parent.get_upstream_paths("database").catalog,
        }

    @dictify
    def _get_input_files_annotate(self, wildcards):
        """Yield input files' pattern for rule `annotate` - based on ExpansionHunter call results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        name_pattern = wildcards.library_name
        yield (
            "json",
            "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                name_pattern=name_pattern, ext="json"
            ),
        )
        yield "annotation", self.parent.get_upstream_paths("database").annotation

    @dictify
    def _get_output_files_run(self):
        """Yield output files' patterns for rule `run` - ExpansionHunter call."""
        for key, ext in {"json": "json", "vcf": "vcf", "vcf_md5": "vcf.md5"}.items():
            yield key, f"work/{{library_name}}/out/{{library_name}}.{ext}"

    @dictify
    def _get_output_files_annotate(self):
        """Yield output files' patterns for rule `annotate`."""
        for key, ext in {"json": "json", "json_md5": "json.md5"}.items():
            yield key, f"work/annotated.{{library_name}}/out/annotated.{{library_name}}.{ext}"

    def _get_log_file_run(self):
        """Return the log file pattern for rule `run` - ExpansionHunter call."""
        return "work/{library_name}/log/{library_name}.log"

    def _get_params_run(self, wildcards):
        """Get donor's sex.

        :param wildcards: Snakemake wildcards associated with rule (unused).
        :type wildcards: snakemake.io.Wildcards

        :return: Returns donor's sex as found in sample sheet: 'female' or 'male', otherwise None.
        """
        sex = self.library_name_to_sex[wildcards.library_name]
        return {"sex": sex if sex in ("female", "male") else None}

    def annotate_results(self, _wildcards, sm_input, sm_output):
        """Annotate/Explain ExpansionHunter results.

        :param _wildcards: Snakemake wildcards associated with rule (unused).
        :type _wildcards: snakemake.io.Wildcards

        :param sm_input: Snakemake input associated with rule.
        :type sm_input: snakemake.io.Namedlist

        :param sm_output: Snakemake output associated with rule.
        :type sm_output: snakemake.io.Namedlist
        """
        # Absolute path from input and output
        input_path = os.path.join(os.getcwd(), sm_input.json)
        output_path = os.path.join(os.getcwd(), sm_output.json)
        # Annotate
        AnnotateExpansionHunter(
            eh_json=input_path,
            annotation_json=sm_input.annotation,
            output_path=output_path,
        ).run()


class RepeatExpansionWorkflow(BaseStep):
    """Perform germline repeat expansion analysis."""

    #: Workflow name
    name = "repeat_expansion"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "repeats"}))]

    config_model_class = RepeatExpansionConfigModel

    #: Sample sheet shortcut class
    sheet_shortcut_class = GermlineCaseSheet

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((ExpansionHunterStepPart,))

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local repeat expansion output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"vcf": f"work/{lib}/out/{lib}.vcf"}

    @listify
    def get_result_files(self):
        """Return list of result files for the germline repeat expansion analysis workflow."""
        # Yield the JSON annotated results files
        name_pattern = "annotated.{library_name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_JSON,
        )
        # Yield the VCF results files
        name_pattern = "{library_name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_VCF,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for library_name in self.primary_dna_libraries()["library_name"]:
            yield from expand(tpl, library_name=[library_name], **kwargs)
