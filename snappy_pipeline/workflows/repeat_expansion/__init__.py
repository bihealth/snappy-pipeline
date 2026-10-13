"""Repeat expansion analysis with ExpansionHunter

The step runs the tool selected by ``tool`` on the alignments from ``depends_on.alignments`` (:ref:`step_ngs_mapping`).
The only value is ``expansionhunter``, which calls repeat expansions against ``repeat_catalog`` and annotates them with ``repeat_annotation``.
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
from collections import OrderedDict

from biomedsheets.shortcuts import KEY_SEX, GermlineCaseSheet, is_not_background
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

    def __init__(self, *args, **kwargs):
        """Constructor."""
        super().__init__(*args, **kwargs)
        #: Build shortcut from library name to sex
        self.library_name_to_sex = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_sex.update(self._library_name_to_sex(sheet))

    @staticmethod
    def _library_name_to_sex(sheet):
        """Library name to sex.

        :param sheet: Sample sheet.
        :type sheet: biomedsheets.shortcuts.GermlineCaseSheet

        :return: Yields (library name, sex).
        """
        for donor in sheet.donors:
            sex = donor.extra_infos.get(KEY_SEX)
            for bio_sample in donor.bio_samples.values():
                for test_sample in bio_sample.test_samples.values():
                    for ngs_library in test_sample.ngs_libraries.values():
                        yield ngs_library.name, sex

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
            "repeat_catalog": self.config.repeat_catalog,
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
        yield "annotation", self.config.repeat_annotation

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

        :return: Returns donor's sex as found in sample sheet: 'female', 'male' or 'unknown'.
        """
        return {"sex": self.library_name_to_sex[wildcards.library_name]}

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
    def _all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the germline repeat expansion analysis workflow."""
        # Yield the JSON annotated results files
        name_pattern = "annotated.{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_JSON,
        )
        # Yield the VCF results files
        name_pattern = "{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_VCF,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for donor in self._all_donors(include_background=False):
            if donor.dna_ngs_library:  # ignores samples without DNA library
                yield from expand(tpl, donor=[donor], **kwargs)
