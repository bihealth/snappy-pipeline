# -*- coding: utf-8 -*-
"""Implementation of the ``reference`` step

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_adapter_trimming.rst

"""

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStepPart, BaseStep, LinkInStep
from .model import ReferenceModel as ReferenceConfigModel

#: Default configuration for the reference
DEFAULT_CONFIG = ReferenceConfigModel.default_config_yaml_string()


class ReferenceStepPart(BaseStepPart):
    """Reference retrieval common features"""

    #: Step name
    name = ""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{source}.{{library_name}}"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        # Validate action
        self._validate_action(action)
        return (("out_done", self.base_path_out.format(source=self.name) + "/out/.done"),)

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)
        _ = action
        prefix = "work/{source}/log/{source}.{{reference_name}}".format(source=self.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        yield (
            "done",
            "work/{source}.{{reference_name}}/log/.done".format(source=self.name),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            return {}

        # Validate action
        self._validate_action(action)
        return args_function


class ReferenceWorkflow(BaseStep):
    """Automatically retrieve reference data"""

    #: Step name
    name = "reference"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, config_model_class=ReferenceConfigModel)
        self.register_sub_step_classes((LinkInStep,))

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    @listify
    def get_result_files(self):
        """Return list of result files for the reference workflow"""
        tpls = ("output/{source}/{reference_name}/out/.done",)
        for name, reference in self.config["references"]:
            for tpl in tpls:
                yield tpl.format(source=reference.source, reference_name=name)
