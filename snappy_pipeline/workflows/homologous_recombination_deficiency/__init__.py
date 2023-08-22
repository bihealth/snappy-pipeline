# -*- coding: utf-8 -*-
"""Implementation of the ``homologous_recombination_deficiency`` step

This step allows for the computation of the scarHRD score, which is a composite of
loss-of-heterozygocity, large-scale transitions and telomeric imbalances.
The score is a proxy for holomogous recombination deficiencies scores produced from SNP arrays.
The software is described in `Sztupinszki et al.<https://doi.org/10.1038/s41523-018-0066-6>`,
but it is not part of CRAN, Bioconductor or bioconda (currently).
The implementation also relies on versions of `sequenza<https://github.com/oicr-gsi/sequenza/tree/master>` &
`copynumber<https://github.com/aroneklund/copynumber>` which are not part of CRAN or Bioconductor anymore.
The most recent version of sequanza is downloaded from anaconda (for the R scripts),
and from bioconda for the python utilities. Note that anaconda lists the r-sequenza package in the
bioconda directory, while this package is not found when searching `bioconda<https://bioconda.github.io/index.html>`.
copynumber is obtained from a fork of the official Bioconductor deprecated package.

==========
Step Input
==========

``homologous_recombination_deficiency`` starts off the aligned reads, i.e. ``ngs_mapping``.
Both the normal & tumor samples are required to generate the score.

===========
Step Output
===========

Generally, the following links are generated to ``output/``.

.. note:: Tool-Specific Output

    As the only integrated tool is scarHRD at the moment, the output is very tailored to the result
    of this tool.  In the future, this section might contain "common" output and tool-specific
    output sub sections.

- ``{mapper}.scarHRD.{lib_name}-{lib_pk}/out/``
    - ``{mapper}.scarHRD.{lib_name}-{lib_pk}.seqz.gz``
    - ``{mapper}.scarHRD.{lib_name}-{lib_pk}.json``

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_homologous_recombination_deficiency.rst

=====================================
Available HRD tools
=====================================

- ``scarHRD``

"""

from collections import OrderedDict
import sys

from biomedsheets.shortcuts import CancerCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the homologous recombination deficiency step
DEFAULT_CONFIG = r"""
# Default configuration homologous_recombination_deficiency
step_config:
  homologous_recombination_deficiency:
    tools: ['scarHRD']  # REQUIRED - available: 'mantis'
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    scarHRD:
      genome_name: "grch37"  # Must be either "grch37", "grch38" or "mouse"
      chr_prefix: False
      length: 50             # Wiggle track for GC reference file
"""


class ScarHRDStepPart(BaseStepPart):
    """Computes homologous recombination deficiency score with scarHRD"""

    #: Step name
    name = "scarHRD"

    #: Class available actions
    actions = (
        "install",
        "gcreference",
        "run",
    )

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.scarHRD.{{library_name}}/out/{{mapper}}.scarHRD.{{library_name}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.library_name]
        return pair.normal_sample.dna_ngs_library.name

    def get_input_files(self, action):
        def input_function_run(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            tumor_base_path = (
                "output/{mapper}.{library_name}/out/" "{mapper}.{library_name}"
            ).format(**wildcards)
            return {
                "lib_path": "work/R_packages/out/.done",
                "gc": "work/static_data/out/{genome_name}_{length}.wig.gz".format(
                    genome_name=self.config["scarHRD"]["genome_name"],
                    length=self.config["scarHRD"]["length"],
                ),
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            }

        if action == "install":
            return None
        elif action == "gcreference":
            return None
        elif action == "run":
            return input_function_run
        else:
            raise UnsupportedActionException(
                "Action '{action}' is not supported. Valid options: {valid}".format(
                    action=action, valid=", ".join(self.actions)
                )
            )

    def get_output_files(self, action):
        if action == "install":
            return {"lib_path": "work/R_packages/out/.done"}
        elif action == "gcreference":
            return {
                "gc": "work/static_data/out/{genome_name}_{length}.wig.gz".format(
                    genome_name=self.config["scarHRD"]["genome_name"],
                    length=self.config["scarHRD"]["length"],
                )
            }
        elif action == "run":
            return {
                "sequenza": "work/{mapper}.scarHRD.{library_name}/out/{mapper}.scarHRD.{library_name}.seqz.gz",
                "scarHRD": "work/{mapper}.scarHRD.{library_name}/out/{mapper}.scarHRD.{library_name}.json",
            }
        else:
            raise UnsupportedActionException(
                "Action '{action}' is not supported. Valid options: {valid}".format(
                    action=action, valid=", ".join(self.actions)
                )
            )

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        if action == "install":
            prefix = "work/R_packages/log/R_packages"
        elif action == "gcreference":
            prefix = "work/static_data/log/{genome_name}_{length}".format(
                genome_name=self.config["scarHRD"]["genome_name"],
                length=self.config["scarHRD"]["length"],
            )
        elif action == "run":
            prefix = "work/{mapper}.scarHRD.{library_name}/log/{mapper}.scarHRD.{library_name}"
        else:
            raise UnsupportedActionException(
                "Action '{action}' is not supported. Valid options: {valid}".format(
                    action=action, valid=", ".join(self.actions)
                )
            )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        if action == "install" or action == "gcreference":
            return ResourceUsage(
                threads=1,
                time="02:00:00",  # 2 hours
                memory="4096M",
                partition="short",
            )
        elif action == "run":
            return ResourceUsage(
                threads=2,
                time="48:00:00",  # 2 hours
                memory="32G",
            )
        else:
            raise UnsupportedActionException(
                "Action '{action}' is not supported. Valid options: {valid}".format(
                    action=action, valid=", ".join(self.actions)
                )
            )


class HomologousRecombinationDeficiencyWorkflow(BaseStep):
    """Compute Homologous Recombination Deficiency score"""

    #: Step name
    name = "homologous_recombination_deficiency"

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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((ScarHRDStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the somatic targeted sequencing CNV calling step"""
        tool_actions = {"scarHRD": ("run",)}
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                for tool in self.config["tools"]:
                    for action in tool_actions[tool]:
                        try:
                            tpls = self.sub_steps[tool].get_output_files(action).values()
                        except AttributeError:
                            tpls = self.sub_steps[tool].get_output_files(action)
                        for tpl in tpls:
                            filenames = expand(
                                tpl,
                                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                                library_name=[sample_pair.tumor_sample.dna_ngs_library.name],
                            )
                            for f in filenames:
                                if ".tmp." not in f:
                                    yield f.replace("work/", "output/")

    def check_config(self):
        """Check that the necessary globalc onfiguration is present"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA file not configured but required",
        )
