# -*- coding: utf-8 -*-
"""
Creates a set of variants (germline AND somatic) suitable to infer purity & somatic CNVs.

Different CNV callers require different annotations for the variants, in order to classify them
into somatic or germline.

``cnvkit`` `needs <https://cnvkit.readthedocs.io/en/stable/fileformats.html#vcf>`_ ``INFO/SOMATIC``
for somatic variants.
It also suggests adding a ``PEDIGREE`` header line, but it doesn't seem to follow vcf format,
and therefore it is not recommended.
The ``somatic_cnv_calling`` implementation of ``cnvkit`` doesn't need this header to work.

``PureCN`` doesn't give a formal requirement of the ``vcf``, but the following fields are listed in its
`documentation <https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/PureCN.pdf>`_:
``FORMAT/AD`` (support for the reference & alternate alleles in the normal & tumor samples),
``INFO/DB`` & ``INFO/POP_AF`` to flag a variant as germline, with its observed population frequency, and
``INFO/SOMATIC`` to flag the variant as somatic.
The documentation also mentions ``BQ`` flag (for base quality), but ``bcftools mpileup`` doesn't emit
this flag (anymore?).

``THetA2`` also requires a file with SNP loci. It is unclear whether it needs only germline, or
germline & somatic locations.
In any case, the format is NOT ``vcf``, and requires a separate re-formatting step
(included in the ``cnvkit`` implementation)

In summary, to create a ``vcf`` that can be used for ``cnvkit``, ``PureCN`` & ``THetA2``,
the following configuration is suggested::

    somatic_variants_for_cnv:
      tools: [bcftools]
      bcftools:
        annotate:
        - path_annotation: /path/to/GATK_af_only_gnomAD.vcf.gz
          columns:
          - INFO/POP_AF:=INFO/AF
          mark_sites: +INFO/DB

The ``INFO/SOMATIC`` flag is automatically inserted by the tool.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variants_for_cnv.rst

"""

import os

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet
from biomedsheets.io_tsv import EXTRACTION_TYPE_DNA
from snakemake.io import Wildcards, InputFiles, pipe

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.samplesheet import sample_sheets, tumor_to_normal_mapping
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.germline_snvs import GermlineSnvsWorkflow
from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow

from .model import SomaticVariantsForCnv as SomaticVariantsForCnvModel
from .model import BcfTools as BcfToolsModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticVariantsForCnvModel.default_config_yaml_string()


class BcftoolsStepPart(BaseStepPart):
    """Simple-minded germline variants pileups (meant for somatic CNV BAF estimation)"""

    # Name of the step.
    name = "bcftools"

    actions = ("merge", "pileup", "annotate", "flag", "remove_unseen")

    resource_usage = {
        "merge": ResourceUsage(threads=1, time="2:00:00", memory=f"{4 * 1024 * 1}M"),
        "pileup": ResourceUsage(threads=2, time="24:00:00", memory=f"{4 * 1024 * 1}M"),
        "annotate": ResourceUsage(threads=1, time="2:00:00", memory=f"{4 * 1024 * 1}M"),
        "flag": ResourceUsage(threads=1, time="2:00:00", memory=f"{4 * 1024 * 1}M"),
        "remove_unseen": ResourceUsage(threads=1, time="2:00:00", memory=f"{4 * 1024 * 1}M"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: BcfToolsModel = getattr(self.config, self.name, {})

    def get_input_files(self, action: str):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def get_output_files(self, action: str):
        base_out = os.path.join("work", "{mapper}.bcftools.{library_name}", "out")

        if action == "flag":
            tpl = os.path.join(base_out, "{mapper}.flag.{library_name}.vcf")
            yield "vcf", pipe(tpl)
        else:
            if action == "merge":
                tpl = os.path.join(base_out, "{mapper}.merge.{library_name}.bed.gz")
                output_files = {"bed": tpl, "bed_tbi": tpl + ".tbi"}
            elif action == "pileup":
                tpl = os.path.join(base_out, "{mapper}.pileup.{library_name}.vcf.gz")
                output_files = {"vcf": tpl, "vcf_tbi": tpl + ".tbi"}
            elif action == "annotate":
                tpl = os.path.join(base_out, "{mapper}.annotate_{n}.{library_name}.vcf.gz")
                output_files = {"vcf": tpl, "vcf_tbi": tpl + ".tbi"}
            elif action == "remove_unseen":
                tpl = os.path.join(base_out, "{mapper}.bcftools.{library_name}.vcf.gz")
                output_files = {"vcf": tpl, "vcf_tbi": tpl + ".tbi"}
            else:
                actions_str = ", ".join(self.actions)
                error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
                raise UnsupportedActionException(error_message)

            for k, v in output_files.items():
                yield k, v
                yield k + "_md5", v + ".md5"

    def get_result_files(self) -> list[str]:
        tpl = "work/{mapper}.bcftools.{library_name}/out/{mapper}.bcftools.{library_name}."
        results = [tpl + "vcf.gz", tpl + "vcf.gz.tbi"]
        return results

    @dictify
    def get_log_file(self, action: str):
        if action in ("merge", "pileup", "flag"):
            tpl = action
        elif action == "annotate":
            tpl = "annotate_{n}"
        elif action == "remove_unseen":
            tpl = "bcftools"
        else:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)

        base_log = "work/{mapper}.bcftools.{library_name}/log/{mapper}." + tpl + ".{library_name}."
        for k, ext in (
            ("conda_info", "conda_info.txt"),
            ("conda_list", "conda_list.txt"),
            ("log", "log"),
            ("script", "sh"),
        ):
            yield k, base_log + ext
            yield k + "_md5", base_log + ext + ".md5"

    @staticmethod
    def _collapsed_arg_value(k: str, v: Any):
        if v is None:
            return None
        k = k.replace("_", "-")
        if isinstance(v, bool):
            if v:
                return f"--{k}"
            else:
                return None
        if isinstance(v, str):
            return f"--{k} '{v}'"
        if isinstance(v, list):
            return "--{k} '{l}'".format(k=k, l=",".join(map(str, v)))
        return f"--{k} {v}"

    @staticmethod
    def collapse_args(args: dict[str, Any]) -> list[str]:
        """
        Creates list of arguments for the ``bcftools`` wrapper.

        The wrapper expects a list of options as ``--option value``
        or simply ``--option`` for flags.

        Many options contain dash signs (``-``), which are incompatible with pydantic models.
        The conversion from underscore to dash characters is done in here.
        This is because, unfortunately, the validation & serialization facilities offered by pydantic
        cannot be used to normalize the option names.
        (There seems to be a conflict between alias generators in the step part models, and
        the complete configuration validation in the abstract initialisation method of step parts)
        """
        collapsed_args = []
        extra_args: dict[str, Any] = args.pop("extra_args")
        for k, v in args.items():
            if k not in extra_args:
                x = BcftoolsStepPart._collapsed_arg_value(k, v)
                if x is not None:
                    collapsed_args.append(x)
        for k, v in extra_args.items():
            x = BcftoolsStepPart._collapsed_arg_value(k, v)
            if x is not None:
                collapsed_args.append(x)
        return sorted(collapsed_args)

    def get_args(self, action: str):
        # Validate action
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_input_files_merge(self, wildcards: Wildcards) -> dict[str, str]:
        germline_snvs = self.parent.sub_workflows["germline_snvs"]
        tpl = "{mapper}.{tool}.{library_name}".format(
            mapper=wildcards.mapper,
            tool=self.config.tool_germline_snvs,
            library_name=self.parent.normals[wildcards.library_name],
        )
        germline = germline_snvs(os.path.join("output", tpl, "out", tpl + ".vcf.gz"))

        somatic_variant_calling = self.parent.sub_workflows["somatic_variant_calling"]
        tpl = "{mapper}.{tool}.{library_name}".format(
            mapper=wildcards.mapper,
            tool=self.config.tool_somatic_variant_calling,
            library_name=wildcards.library_name,
        )
        somatic = somatic_variant_calling(os.path.join("output", tpl, "out", tpl + ".vcf.gz"))

        input_files = {"germline": germline, "somatic": somatic}

        if len(self.config.ignore_chroms) + len(self.cfg.ignore_chroms) > 0:
            input_files["regions"] = "work/bcftools/out/regions.bed.gz".format(**wildcards)

        return input_files

    def _get_args_merge(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {}

    def _get_input_files_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {
            "fasta_ref": self.w_config.static_data_config.reference.path,
            "regions_file": "work/{mapper}.bcftools.{library_name}/out/{mapper}.merge.{library_name}.bed.gz".format(
                **wildcards
            ),
        }

        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{{library_name}}/out/{mapper}.{{library_name}}.bam".format(
            **wildcards
        )
        input_files["bams"] = [
            ngs_mapping(base_path.format(library_name=self.parent.normals[wildcards.library_name])),
            ngs_mapping(base_path.format(library_name=wildcards.library_name)),
        ]
        input_files["bais"] = [bam + ".bai" for bam in input_files["bams"]]

        return input_files

    def _get_args_pileup(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        params = self.cfg.pileup.model_dump(by_alias=True)
        for k in ("skip_all_set", "skip_any_set", "skip_all_unset", "skip_any_unset"):
            if len(params[k]) == 0:
                del params[k]
        return {"index": True, "extra_args": BcftoolsStepPart.collapse_args(params)}

    def _get_input_files_annotate(self, wildcards: Wildcards) -> dict[str, str]:
        index = int(wildcards.n)
        if index == 0:
            label = "pileup"
        else:
            label = "annotate_{}".format(index - 1)
        return {
            "vcf": "work/{mapper}.bcftools.{library_name}/out/{mapper}.{label}.{library_name}.vcf.gz".format(
                label=label, **wildcards
            ),
            "annotations": self.cfg.annotate[int(wildcards.n)].path_annotation,
        }

    def _get_args_annotate(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        index = int(wildcards.n)
        params = self.cfg.annotate[index].model_dump(by_alias=True)
        del params["path_annotation"]
        if params["header_line"] is None:
            del params["header_line"]
        return {"extra_args": BcftoolsStepPart.collapse_args(params), "index": True}

    def _get_input_files_flag(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "work/{mapper}.bcftools.{library_name}/out/{mapper}.{substep}.{library_name}.vcf.gz"
        last = len(self.cfg.annotate)
        if last == 0:
            vcf = tpl.format(substep="pileup", **wildcards)
        else:
            vcf = tpl.format(substep=f"annotate_{last - 1}", **wildcards)
        somatic_variant_calling = self.parent.sub_workflows["somatic_variant_calling"]
        base_out = "output/{mapper}.{caller}.{library_name}/out/{mapper}.{caller}.{library_name}.vcf.gz".format(
            mapper=wildcards.mapper,
            caller=self.config.tool_somatic_variant_calling,
            library_name=wildcards.library_name,
        )
        return {"vcf": vcf, "annotations": somatic_variant_calling(base_out)}

    def _get_args_flag(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {"extra_args": ["--mark-sites '+SOMATIC'", "--regions-overlap 2"]}

    def _get_input_files_remove_unseen(self, wildcards: Wildcards) -> dict[str, str]:
        return {
            "vcf": "work/{mapper}.bcftools.{library_name}/out/{mapper}.flag.{library_name}.vcf".format(
                **wildcards
            )
        }

    def _get_args_remove_unseen(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {"extra_args": ["--trim-alt-alleles", "--trim-unseen-allele"], "index": True}


class SomaticVariantsForCnvWorkflow(BaseStep):
    """Extract germline variants for BAF"""

    #: Workflow name
    name = "somatic_variants_for_cnv"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticVariantsForCnvModel,
            previous_steps=(
                NgsMappingWorkflow,
                GermlineSnvsWorkflow,
                SomaticVariantCallingWorkflow,
            ),
        )
        self.table = sample_sheets(self.sheets)
        self.table = self.table[self.table["extractionType"] == EXTRACTION_TYPE_DNA]
        self.register_sub_step_classes((BcftoolsStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        self.register_sub_workflow("germline_snvs", self.config.path_germline_snvs)
        self.register_sub_workflow(
            "somatic_variant_calling", self.config.path_somatic_variant_calling
        )

        self.normals = tumor_to_normal_mapping(self.table)

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        if self.table.shape[0] == 0:
            self.logger.warning("No sample to process")
        for sub_step in self.sub_steps.values():
            tool = sub_step.name
            tool_config = getattr(self.config, tool, None)
            if tool_config is None:
                continue
            for ngs_library in self.table[self.table["isTumor"]].index.to_list():
                files_in_work = sub_step.get_result_files()
                for file_in_work in files_in_work:
                    file_in_output = file_in_work.format(
                        mapper=self.config.tool_ngs_mapping, library_name=ngs_library
                    ).replace("work/", "output/", 1)
                    yield file_in_output
                    yield file_in_output + ".md5"
