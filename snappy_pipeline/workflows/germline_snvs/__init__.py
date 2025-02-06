# -*- coding: utf-8 -*-
"""
Computes strongly supported germline variants for CNV estimation from B-allele fractions (BAF).
This step produces the input for the ``somatic_variants_for_cnv`` step, which outputs a vcf file
containing the depth, vaf & annotations at germline & somatic variants positions.
This can then be used by CNV callers such as ``PureCN``, ``cnvkit`` & others.

The current implementation is very simple, based on bcftools only.
The algorithm is as follows. It can be appiled to any DNA sample, but it is meant only for normal samples:

1. Identify variants using bcftools mpileup, on the whole genome, or on user-supplied regions
2. Call the genotype for all these variants
3. Filter variants for strong support
   (by default, depth>=50, suppport for alt between 0.82 & 1.22 times the support for ref)
4. Optionally annotate variants with gnomAD or dbSNP or others.

When the aim is to compute BAF for CNV, annotations are not required at this point.

.. note::

    Multiple annotations can be added. Each of them will trigger an independent run of the wrapper.
    The output & log files will be numbered, and the last job will take the final annotation output
    and just remove unseen alleles.

----------------------
Implementation details
----------------------

Except for the annotations, the pileup, call & filter rules are linked through named pipes.
(This is because ``bcftools annotate`` doesn't accept input from ``stdin``)
It has an impact on the resources, as they should all share the same time resources.
The number of threads is also affected, the total requirement must count all sub-steps.

.. note::

    When there are no annotations, the whole procedure could be connected through named pipes.
    It is not the case at the moment, but it is on the future developments roadmap.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_germline_snvs.rst

"""

import os

from typing import Any

from biomedsheets.shortcuts import GenericSampleSheet
from biomedsheets.io_tsv import EXTRACTION_TYPE_DNA
from snakemake.io import Wildcards, InputFiles, pipe

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.samplesheet import sample_sheets
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage

from .model import GermlineSnvs as GermlineSnvsModel
from .model import BcfTools as BcfToolsModel
from .model import BcftoolsCaller

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = GermlineSnvsModel.default_config_yaml_string()


class BcftoolsStepPart(BaseStepPart):
    """Simple-minded variants discovery (meant for somatic CNV BAF estimation)"""

    # Name of the step.
    name = "bcftools"

    actions = ("pileup", "call", "filter", "annotate", "remove_unseen")

    resource_usage = {
        "pileup": ResourceUsage(threads=2, time="24:00:00", memory=f"{4 * 1024 * 1}M"),
        "call": ResourceUsage(threads=2, time="24:00:00", memory=f"{4 * 1024 * 1}M"),
        "filter": ResourceUsage(threads=1, time="24:00:00", memory=f"{4 * 1024 * 1}M"),
        "annotate": ResourceUsage(threads=1, time="2:00:00", memory=f"{4 * 1024 * 1}M"),
        "last": ResourceUsage(threads=1, time="2:00:00", memory=f"{4 * 1024 * 1}M"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: BcfToolsModel = getattr(self.config, self.name, {})

    def get_input_files(self, action: str):
        """
        Return input files functions

        The ``pileup`` step takes the sample bam file (& reference) as input,
        the ``call`` & ``filter`` steps have a vcf named pipe as sole input
        (produced by the former sub-step).
        Unfortunately, ``bcftools annotate`` cannot use piped input, so it must be
        taken out of the generic input.
        The last step is to remove unseen alleles from the final result.
        """
        # Validate action
        self._validate_action(action)

        if action == "pileup":
            return self._get_input_files_pileup
        if action == "annotate":
            return self._get_input_files_annotate
        if action == "remove_unseen":
            return self._get_input_files_remove_unseen

        i = self.actions.index(action)

        def input_function(wildcards: Wildcards):
            previous = BcftoolsStepPart.actions[i - 1]
            return {
                "vcf": "work/{mapper}.bcftools.{library_name}/out/{mapper}.{action}.{library_name}.vcf".format(
                    mapper=wildcards.mapper,
                    library_name=wildcards.library_name,
                    action=previous,
                )
            }

        return input_function

    @dictify
    def get_output_files(self, action: str):
        # Validate action
        self._validate_action(action)
        # Pipes output are uncompressed vcfs (except for last sub-step)
        if action in ("pileup", "call"):
            yield (
                "vcf",
                pipe(
                    f"work/{{mapper}}.bcftools.{{library_name}}/out/{{mapper}}.{action}.{{library_name}}.vcf"
                ),
            )
        else:
            match action:
                case "filter":
                    tpl = "{mapper}.filter.{library_name}."
                case "annotate":
                    tpl = "{mapper}.annotate_{n}.{library_name}."
                case "remove_unseen":
                    tpl = "{mapper}.bcftools.{library_name}."
            tpl = os.path.join("work", "{mapper}.bcftools.{library_name}", "out", tpl)
            output_files = {"vcf": tpl + "vcf.gz", "vcf_tbi": tpl + "vcf.gz.tbi"}
            for k, v in output_files.items():
                yield k, v
                yield k + "_md5", v + ".md5"

    def get_result_files(self) -> list[str]:
        tpl = "work/{mapper}.bcftools.{library_name}/out/{mapper}.bcftools.{library_name}."
        results = [tpl + "vcf.gz", tpl + "vcf.gz.tbi"]
        return results

    @dictify
    def get_log_file(self, action: str):
        # Validate action
        self._validate_action(action)
        if action == "annotate":
            action = "annotate_{n}"
        if action == "remove_unseen":
            action = "bcftools"
        tpl = f"work/{{mapper}}.bcftools.{{library_name}}/log/{{mapper}}.{action}.{{library_name}}."
        for ext in ("conda_info.txt", "conda_list.txt", "log"):
            k = ext.replace(".txt", "")
            v = tpl + ext
            yield k, v
            yield k + "_md5", v + ".md5"
        k = "script"
        v = tpl + "sh"
        yield k, v
        yield k + "_md5", v + ".md5"

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

    def _get_input_files_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {
            "bams": [ngs_mapping(base_path + ".bam")],
            "bais": [ngs_mapping(base_path + ".bam.bai")],
            "fasta_ref": self.w_config.static_data_config.reference.path,
        }
        if len(self.config.ignore_chroms) + len(self.cfg.ignore_chroms) > 0:
            input_files["regions_file"] = "work/bcftools/out/regions.bed.gz".format(**wildcards)
        return input_files

    def _get_args_pileup(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        params = self.cfg.pileup.model_dump(by_alias=True)
        for k in ("skip_all_set", "skip_any_set", "skip_all_unset", "skip_any_unset"):
            if len(params[k]) == 0:
                del params[k]
        if self.cfg.skip_indels:
            params["skip_indels"] = True
        return {"extra_args": BcftoolsStepPart.collapse_args(params)}

    def _get_args_call(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        params = self.cfg.call.model_dump(by_alias=True)
        if params["caller"] == BcftoolsCaller.CONSENSUS:
            params["consensus_caller"] = True
        else:
            params["multiallelic_caller"] = True
        del params["caller"]
        if self.cfg.skip_indels:
            params["skip_variants"] = "indels"
        return {"extra_args": BcftoolsStepPart.collapse_args(params)}

    def _get_args_filter(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        params = self.cfg.filter.model_dump(by_alias=True)

        if "include" not in params["extra_args"] and "exclude" not in params["extra_args"]:
            # Default filtration:
            # Depth > min_depth & rmin*#ref < #alt < rmax*#ref & heterozygous & genotype_quality > min_qual
            depth = f"FORMAT/AD[0:0] + FORMAT/AD[0:1] >= {self.cfg.filter.min_depth}"
            r_min = self.cfg.filter.min_baf / (1 - self.cfg.filter.min_baf)
            r_max = 1 / r_min
            baf = f"({r_min} * FORMAT/AD[0:0] <= FORMAT/AD[0:1]) & (FORMAT/AD[0:1] <= {r_max} * FORMAT/AD[0:0])"
            genotype = '(FORMAT/GT[0] == "0/1") | (FORMAT/GT[0] == "0|1")'
            qual = f"QUAL >= {self.cfg.filter.min_qual}"
            include = f"({depth}) & ({baf}) & ({genotype}) & ({qual})"
            if self.cfg.skip_indels:
                include += ' & (TYPE = "snp")'
            params["include"] = include

        for k in ("min_depth", "min_baf", "min_qual"):
            del params[k]

        return {"extra_args": BcftoolsStepPart.collapse_args(params), "index": True}

    def _get_input_files_annotate(self, wildcards: Wildcards) -> dict[str, str]:
        index = int(wildcards.n)
        if index == 0:
            label = "filter"
        else:
            label = "annotate_{}".format(index - 1)
        return {
            "vcf": "work/{mapper}.bcftools.{library_name}/out/{mapper}.{label}.{library_name}.vcf.gz".format(
                label=label, **wildcards
            ),
            "annotations": self.cfg.annotate[index].path_annotation,
        }

    def _get_args_annotate(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        index = int(wildcards.n)
        params = self.cfg.annotate[index].model_dump(by_alias=True)
        del params["path_annotation"]
        return {"extra_args": BcftoolsStepPart.collapse_args(params), "index": True}

    def _get_input_files_remove_unseen(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "work/{mapper}.bcftools.{library_name}/out/{mapper}.{substep}.{library_name}.vcf.gz"
        last = len(self.cfg.annotate)
        if last == 0:
            vcf = tpl.format(substep="filter", **wildcards)
        else:
            vcf = tpl.format(substep=f"annotate_{last - 1}", **wildcards)
        return {"vcf": vcf}

    def _get_args_remove_unseen(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {"index": True, "extra_args": ["--trim-alt-alleles", "--trim-unseen-allele"]}


class GermlineSnvsWorkflow(BaseStep):
    """Extract germline variants for BAF"""

    #: Workflow name
    name = "germline_snvs"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

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
            config_model_class=GermlineSnvsModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        self.table = sample_sheets(self.sheets)
        self.table = self.table[self.table["extractionType"] == EXTRACTION_TYPE_DNA]
        self.register_sub_step_classes((BcftoolsStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        if self.table.shape[0] == 0:
            self.logger.warning("No sample to process")
        else:
            for sub_step in self.sub_steps.values():
                tool = sub_step.name
                tool_config = getattr(self.config, tool, None)
                if tool_config is None:
                    continue
                for ngs_library in self.table.index.to_list():
                    files_in_work = sub_step.get_result_files()
                    for file_in_work in files_in_work:
                        file_in_output = file_in_work.format(
                            mapper=self.config.tool_ngs_mapping, library_name=ngs_library
                        ).replace("work/", "output/", 1)
                        yield file_in_output
                        yield file_in_output + ".md5"
