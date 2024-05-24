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

TODO: The implementation needs to uncouple ``scarHRD`` from ``sequenza``.
It should be possible as ``scarHRD`` could use the output from other CNV callers.
However, it is poorly documented, and will require significant work.

==========
Step Input
==========

``homologous_recombination_deficiency`` starts off the ``sequenza`` copy number output.
Some steps performed by ``sequenza`` are repeated by ``scanHRD``, as ``scanHRD`` uses a specific
set of arguments. This is the reason why ``sequenza`` is required by the step.

===========
Step Output
===========

Generally, the following links are generated to ``output/``.

.. note:: Tool-Specific Output

    As the only integrated tool is scarHRD at the moment, the output is very tailored to the result
    of this tool.  In the future, this section might contain "common" output and tool-specific
    output sub sections.

- ``{mapper}.{caller}.scarHRD.{lib_name}-{lib_pk}/out/``
    - ``{mapper}.{caller}.scarHRD.{lib_name}-{lib_pk}.seqz.gz``
    - ``{mapper}.{caller}.scarHRD.{lib_name}-{lib_pk}.json``

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
from snappy_pipeline.workflows.somatic_targeted_seq_cnv_calling import (
    SomaticTargetedSeqCnvCallingWorkflow,
)
from .model import HomologousRecombinationDeficiency as HomologousRecombinationDeficiencyConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the homologous recombination deficiency step
DEFAULT_CONFIG = HomologousRecombinationDeficiencyConfigModel.default_config_yaml_string()


class ScarHRDStepPart(BaseStepPart):
    """Computes homologous recombination deficiency score with scarHRD"""

    #: Step name
    name = "scarHRD"

    #: Class available actions
    actions = (
        "install",
        "run",
    )

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        self._validate_action(action)

        return self._get_input_files_run

    @dictify
    def _get_input_files_run(self, wildcards):
        self.cnv_calling = self.parent.sub_workflows["cnv_calling"]
        base_name = f"{wildcards.mapper}.{wildcards.caller}.{wildcards.library_name}"
        yield "done", "work/R_packages/out/scarHRD.done"
        yield "seqz", self.cnv_calling(f"output/{base_name}/out/{base_name}.seqz.gz")

    def get_output_files(self, action):
        if action == "install":
            return {"done": "work/R_packages/out/scarHRD.done"}
        elif action == "run":
            return {
                "scarHRD": "work/{mapper}.{caller}.scarHRD.{library_name}/out/{mapper}.{caller}.scarHRD.{library_name}.json",
                "scarHRD_md5": "work/{mapper}.{caller}.scarHRD.{library_name}/out/{mapper}.{caller}.scarHRD.{library_name}.json.md5",
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
            prefix = "work/R_packages/log/scarHRD"
        elif action == "run":
            prefix = "work/{mapper}.{caller}.scarHRD.{library_name}/log/{mapper}.{caller}.scarHRD.{library_name}"
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
        self._validate_action(action)
        if action == "run":
            return ResourceUsage(
                threads=1,
                memory="32G",
                time="24:00:00",
            )
        else:
            return super().get_resource_usage(action)


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
            config_model_class=HomologousRecombinationDeficiencyConfigModel,
            previous_steps=(SomaticTargetedSeqCnvCallingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((ScarHRDStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow(
            "somatic_targeted_seq_cnv_calling", self.config["path_cnv_calling"], "cnv_calling"
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the homologous recombination deficiency step"""
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
                        tpls = list(tpls)
                        tpls += list(self.sub_steps[tool].get_log_file(action).values())
                        for tpl in tpls:
                            filenames = expand(
                                tpl,
                                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                                caller=["sequenza"],
                                library_name=[sample_pair.tumor_sample.dna_ngs_library.name],
                            )
                            for f in filenames:
                                if ".tmp." not in f and not f.endswith(".done"):
                                    yield f.replace("work/", "output/")

    def check_config(self):
        """Check that the necessary globalc onfiguration is present"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA file not configured but required",
        )
        assert (
            "sequenza" in self.w_config["step_config"]["somatic_targeted_seq_cnv_calling"]["tools"]
        )
