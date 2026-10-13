"""IGV session generation

The step writes an IGV session XML file per pedigree.
The session shows the alignments of all pedigree members and the variant calls of the index as genome tracks.
The task has no ``tool`` field.

``depends_on.alignments`` names the ``ngs_mapping`` task.
The variants come from exactly one of ``depends_on.variants`` (:ref:`step_variant_calling`), ``depends_on.annotated_variants`` (:ref:`step_variant_annotation`) or ``depends_on.phased_variants`` (:ref:`step_variant_phasing`).
The session refers to the published files in ``output/`` of these tasks by relative path.
Only pedigrees whose index has both parents with a DNA library get a session.

Outputs
-------

The task publishes ``{library_name}/out/{library_name}.igv_session.xml`` and its ``.md5`` file below ``tasks/<task>/output/``.
The path is the same as below ``work/``, and ``{library_name}`` is the DNA library of the index.
The key ``xml`` is the one downstream tasks read through ``depends_on``.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import os
import re

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from jinja2 import Environment, FileSystemLoader
from snakemake import shell
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import IgvSessionGeneration as IgvSessionGenerationConfigModel

#: Output keys and extensions of the session file
_OUTPUTS = {"xml": ".igv_session.xml", "xml_md5": ".igv_session.xml.md5"}

#: Default configuration of the wgs_sv_filtration step


class WriteIgvSessionFileStepPart(BaseStepPart):
    """Write out the IGV session XML file."""

    #: Step name
    name = "write_igv_session_file"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.previous_step = self.parent.previous_step
        self.prev_token = self.parent.prev_token
        self.ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    self.ngs_library_to_pedigree[donor.dna_ngs_library.name] = pedigree

    def _get_path_bam(self, wildcards, donor):
        # TODO: This cannot be correct. For each set of donors it will return the same index bam.
        # TODO: For instance, given pedigree (P001, P002, P003) it will return three runtime the
        # TODO: same value: '../ngs_mapping/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam'
        _ = donor
        return self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.index_library
        ).bam

    def _get_path_vcf(self, wildcards, real_index):
        kwargs = {"library_name": real_index.dna_ngs_library.name}
        if self.previous_step == "phased_variants":
            kwargs["phasing"] = "gatk_pbt.gatk_rbp"
        return self.parent.get_upstream_paths(self.previous_step, **kwargs).vcf

    @staticmethod
    def _published(path):
        """Return the published ``output/`` path of an upstream ``work/`` path.

        The session file references the published files, not the ``work/`` files other tasks read.
        """
        return re.sub(r"(^|/)work/", r"\1output/", path, count=1)

    @dictify
    def _get_input_files_run(self, wildcards):
        # Get name of real index, used when input is not variant_phasing
        pedigree = self.ngs_library_to_pedigree[wildcards.index_library]
        # BAM and BAI file of the offspring
        yield "bam", [self._get_path_bam(wildcards, donor) for donor in pedigree.donors]
        # Input file comes from previous step.
        yield "vcf", self._get_path_vcf(wildcards, pedigree.index)

    def _get_output_files_run(self):
        return {
            key: f"work/{{index_library}}/out/{{index_library}}{ext}"
            for key, ext in _OUTPUTS.items()
        }

    @staticmethod
    def render_from_template(directory, template_name, **kwargs):
        loader = FileSystemLoader(directory)
        env = Environment(loader=loader)
        template = env.get_template(template_name)
        return template.render(**kwargs)

    def run(self, wildcards, output):
        """Write out the IGV session file information"""
        # Get information for writing out the template.
        pedigree = self.ngs_library_to_pedigree[wildcards.index_library]
        template_vars = {
            "bams": [
                {
                    "label": donor.name,
                    "path": os.path.relpath(
                        self._published(self._get_path_bam(wildcards, donor)), output.xml
                    ),
                }
                for donor in pedigree.donors
            ],
            "vcfs": [
                {
                    "label": "Variant Calls",
                    "path": os.path.relpath(
                        self._published(self._get_path_vcf(wildcards, pedigree.index)),
                        output.xml,
                    ),
                }
            ],
        }
        # Write output file.
        with open(output.xml, "wt") as session_file:
            xml = self.__class__.render_from_template(
                os.path.dirname(__file__), "igv_session.xml", **template_vars
            )
            print(xml, file=session_file)
        # Write MD5 sum file.
        shell(
            r"""
            pushd $(dirname {output.xml})
            md5sum $(basename {output.xml}) >$(basename {output.xml}).md5
            """.format(output=output)
        )


class IgvSessionGenerationWorkflow(BaseStep):
    """Perform IGV session generation"""

    #: Workflow name
    name = "igv_session_generation"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"igv"}))]

    config_model_class = IgvSessionGenerationConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local IGV session output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        token = kwargs.get("token", "")
        if token:
            token = token + "."
        prefix = f"work/{token}{lib}/out/{token}{lib}"
        return {"xml": f"{prefix}.igv_session.xml"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        for prev in ("phased_variants", "annotated_variants", "variants"):
            if getattr(self.config.depends_on, prev, None):
                self.previous_step = prev
                break
        else:
            raise Exception("No previous step given!")  # pragma: no cover
        #: Name token for input
        self.prev_token = {
            "phased_variants": "jannovar_annotate_vcf.gatk_pbt.gatk_rbp.",
            "annotated_variants": "jannovar_annotate_vcf.",
            "variants": "",
        }[self.previous_step]
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((WriteIgvSessionFileStepPart,))
        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here

    @listify
    def get_result_files(self):
        """Return list of result files for the workflow."""
        # Hard-filtered results
        name_pattern = "%s.{index_library.name}" % (self.prev_token,)
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=_OUTPUTS.values(),
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        # TODO: currently only families with children are supported
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                for donor in pedigree.donors:
                    if not donor.dna_ngs_library:
                        continue
                    elif not donor.father or not donor.father.dna_ngs_library:
                        continue
                    elif not donor.mother or not donor.mother.dna_ngs_library:
                        continue
                    else:
                        yield from expand(
                            tpl,
                            index_library=[donor.dna_ngs_library],
                            **kwargs,
                        )
