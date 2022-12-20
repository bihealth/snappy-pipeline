# -*- coding: utf-8 -*-
"""Implementation of the gCNV CASE mode run methods.
"""

import glob
import json
import os
import re

from snakemake.io import touch

from snappy_pipeline.base import InvalidConfiguration, UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.gcnv.gcnv_common import GcnvStepPart


class RunGcnvStepPart(GcnvStepPart):
    """Class with methods to run GATK4 gCNV"""

    #: Class available actions
    actions = (
        "preprocess_intervals",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "merge_cohort_vcfs",
        "extract_ped",
    )

    def __init__(self, parent):
        super().__init__(parent)
        # Validate configuration, precomputed models must be present
        self.validate_request()

    def validate_request(self):
        """Validate request.

        Checks if the request can be performed given the provided information and parameters:
        Given a path to a model provided in the config (``precomputed_model_paths``), it will
        check that the path were provided for each available library kit and that the paths
        contain the necessary files for analysis.

        :raises InvalidConfiguration: if information provided in configuration isn't enough to run
        the analysis.
        """
        if "gcnv" not in self.config["tools"]:
            return

        # Get precomputed models from configurations
        path_to_models = self.config["gcnv"]["precomputed_model_paths"]

        # No model provided
        if not path_to_models:
            msg_tpl = "Precomputed model paths must be configured (key: 'precomputed_model_paths')."
            raise InvalidConfiguration(msg_tpl)
        else:
            # Validate configuration - check if only expected keys are present
            self.validate_precomputed_model_paths_config(config=path_to_models)

            # Check model directories content
            for model in path_to_models:
                # Validate ploidy-model
                ploidy_path = model.get("contig_ploidy")
                if not self.validate_ploidy_model_directory(path=ploidy_path):
                    msg_tpl = (
                        "Provided path either not a directory or "
                        "does not contain all required contig-ploidy model files: {0}"
                    )
                    raise InvalidConfiguration(msg_tpl.format(str(model)))
                # Find all model directories
                model_pattern = model.get("model_pattern")
                path_list = self.get_model_dir_list(pattern=model_pattern)
                if len(path_list) == 0:
                    msg_tpl = "Could not find any directory path with the provided pattern: '{0}' "
                    raise InvalidConfiguration(msg_tpl.format(model_pattern))
                # Validate directory
                for path in path_list:
                    if not self.validate_call_model_directory(path=path):
                        msg_tpl = (
                            "Provided path either not a directory or "
                            "does not contain all required call model files: {0}"
                        )
                        raise InvalidConfiguration(msg_tpl.format(str(model)))

    def validate_precomputed_model_paths_config(self, config):
        """Validate precomputed model config.

        Evaluates if provided configuration has the following format:

        precomputed_model_paths:
          - library: "Agilent SureSelect Human All Exon V6"
            contig_ploidy": /path/to/ploidy-model
            model_pattern: /path/to/model_*

        :param config: List of precomputed model configuration dictionary.
        :type config: list

        :raises InvalidConfiguration: if configuration not as expected for
        ``precomputed_model_paths`` list.
        """
        # Initialise variables
        expected_keys = ("library", "model_pattern", "contig_ploidy")
        expected_format = (
            '{\n    "library": "Agilent SureSelect Human All Exon V6"\n'
            '    "contig_ploidy": /path/to/ploidy-model\n'
            '    "model_pattern": "/path/to/model_*"\n}'
        )
        # Test
        for model in config:
            # Test keys
            n_keys_pass = len(model) == 3
            keys_pass = all(key in expected_keys for key in model)
            # Test values
            values_pass = all(isinstance(value, str) for value in model.values())
            # Validate
            if not (n_keys_pass and keys_pass and values_pass):
                msg_tpl = (
                    "Provided configuration not as expected...\n"
                    "Expected:\n{e_}\nObserved:\n{o_}\n"
                )
                pretty_model = self._pretty_print_config(config=model)
                raise InvalidConfiguration(msg_tpl.format(e_=expected_format, o_=pretty_model))

    @staticmethod
    def _pretty_print_config(config):
        """Pretty format configuration.

        :param config: Configuration dictionary to be formatted.
        :type config: OrderedDict

        :return: Configuration as a nicely formatted string.
        """
        return str(json.dumps(config, sort_keys=False, indent=4))

    @staticmethod
    def validate_ploidy_model_directory(path):
        """Validate gCNV ploidy-model directory.

        :param path: Path to gCNV ploidy-model directory.
        :return: Returns True if path is a directory and contains all required files; otherwise,
        False.
        """
        # Model required files
        model_files = [
            "contig_ploidy_prior.tsv",
            "gcnvkernel_version.json",
            "interval_list.tsv",
            "mu_mean_bias_j_lowerbound__.tsv",
            "mu_psi_j_log__.tsv",
            "ploidy_config.json",
            "std_mean_bias_j_lowerbound__.tsv",
            "std_psi_j_log__.tsv",
        ]
        # Check if path is a directory
        if not os.path.isdir(path):
            return False
        # Check if directory contains required model files
        if not all(os.path.isfile(os.path.join(path, m_file)) for m_file in model_files):
            return False
        return True

    @staticmethod
    def validate_call_model_directory(path):
        """Validate gCNV call-model directory.

        :param path: Path to gCNV call-model directory. Files created using COHORT mode.
        :return: Returns True if path is a directory and contains all required files; otherwise,
        False.
        """
        # Model required files
        model_files = [
            "calling_config.json",
            "gcnvkernel_version.json",
            "log_q_tau_tk.tsv",
            "mu_ard_u_log__.tsv",
            "mu_psi_t_log__.tsv",
            "std_ard_u_log__.tsv",
            "std_psi_t_log__.tsv",
            "denoising_config.json",
            "interval_list.tsv",
            "mu_W_tu.tsv",
            "mu_log_mean_bias_t.tsv",
            "std_W_tu.tsv",
            "std_log_mean_bias_t.tsv",
        ]
        # Check if path is a directory
        if not os.path.isdir(path):
            return False
        # Check if directory contains required model files
        if not all(os.path.isfile(os.path.join(path, m_file)) for m_file in model_files):
            return False
        return True

    @staticmethod
    def get_model_dir_list(pattern):
        """Get model directories.

        :param pattern: Pattern of model directory paths. Expects models to be based on scattered
        step: ``gatk IntervalListTools``.
        :type pattern: str

        :return: Returns list with all directories that match the inputted pattern.
        """
        return [path_ for path_ in glob.glob(pattern) if os.path.isdir(path_)]

    def get_params(self, action):
        """
        :param action: Action (i.e., step) in the workflow. Currently available for:
        'ploidy-model', 'model', and 'postgermline_models'.
        :type action: str

        :return: Returns input function for gCNV rule based on inputted action.

        :raises UnsupportedActionException: if invalid action.
        """
        # Actions with parameters
        valid_actions = ("model", "ploidy_model", "postgermline_models")
        # Validate inputted action
        if action not in valid_actions:
            error_message = "Action '{action}' is not supported. Valid options: {options}.".format(
                action=action, options=", ".join(valid_actions)
            )
            raise UnsupportedActionException(error_message)

        # Return requested function
        return getattr(self, "_get_params_{}".format(action))

    def _get_params_ploidy_model(self, wildcards):
        """Get ploidy-model parameters.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'library_kit'
        (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns ploidy-model parameters dictionary if analysis type is 'case_mode';
        otherwise, returns empty dictionary. Step: Calling autosomal and allosomal contig ploidy
        with `DetermineGermlineContigPloidy`.
        """
        path = "__no_ploidy_model_for_library_in_config__"
        for model in self.config["gcnv"]["precomputed_model_paths"]:
            # Adjust library kit name from config to wildcard
            library_to_wildcard = model.get("library").strip().replace(" ", "_")
            if library_to_wildcard == wildcards.library_kit:
                path = model.get("contig_ploidy")
        return {"model": path}

    def _get_params_model(self, wildcards):
        """Get model parameters.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'library_kit'
        (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns model parameters dictionary if analysis type is 'case_mode';
        otherwise, returns empty dictionary. Step: Calling copy number variants with
        `GermlineCNVCaller`.
        """
        path = "__no_model_for_library_in_config__"
        for model in self.config["gcnv"]["precomputed_model_paths"]:
            # Adjust library kit name from config to wildcard
            library_to_wildcard = model.get("library").strip().replace(" ", "_")
            if library_to_wildcard == wildcards.library_kit:
                pattern = model.get("model_pattern")
                model_dir_dict = self._get_model_dir_to_dict(pattern)
                path = model_dir_dict.get(wildcards.shard)
        return {"model": path}

    def _get_params_postgermline_models(self, wildcards):
        """Get post germline model parameters.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'library_name'
        (e.g., 'P001-N1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards

        :return: Returns model parameters dictionary if analysis type is 'case_mode';
        otherwise, returns empty dictionary. Step: consolidating the scattered
        `GermlineCNVCaller` results, performs segmentation and calls copy number states with
        `PostprocessGermlineCNVCalls `.
        """
        paths = ["__no_model_available_for_library__"]
        for model in self.config["gcnv"]["precomputed_model_paths"]:
            # Adjust library kit name from config to wildcard
            library_to_wildcard = model.get("library").strip().replace(" ", "_")
            # Get library kit associated with library name
            library_kit = self.ngs_library_to_kit[wildcards.library_name]
            if library_to_wildcard == library_kit:
                pattern = model.get("model_pattern")
                model_dir_dict = self._get_model_dir_to_dict(pattern)
                paths = list(model_dir_dict.values())
        return {"model": paths}

    def _get_model_dir_to_dict(self, pattern):
        """Get model directories dictionary.

        :param pattern: Pattern of model directory paths. Expects models to be based on scattered
        step: ``gatk IntervalListTools``.
        :type pattern: str

        :return: Returns dictionary with model directories' information.
        Key: shard (e.g., '01_of_42'); Value: directory path (e.g., '/path/to/model_01_of_42').
        """
        # Initialise variables
        out_dict = {}
        default_key = "001"
        re_pattern = pattern.replace("*", "(.*)")
        # Get all model directories
        path_list = self.get_model_dir_list(pattern)
        # Populate dictionary
        # Assumption: if no group, a single path was provided instead of a pattern.
        for path in path_list:
            try:
                key = re.search(re_pattern, path).group(1)
                out_dict[key] = path
            except IndexError:
                out_dict[default_key] = path
        return out_dict

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards):
        """Yield input files for ``contig_ploidy`` rule  in CASE MODE using precomputed model.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                tsvs.append(
                    "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                        name_pattern=name_pattern, ext=ext
                    )
                )
        yield ext, tsvs

    @dictify
    def _get_input_files_call_cnvs(self, wildcards):
        """Yield input files for ``call_cnvs`` in CASE mode.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        # Initialise variables
        tsv_ext = "tsv"
        tsv_path_pattern = "{mapper}.gcnv_coverage.{library_name}"
        ploidy_ext = "ploidy"
        ploidy_path_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}"

        # Yield coverage tsv files for all library associated with kit
        coverage_files = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                path_pattern = tsv_path_pattern.format(mapper=wildcards.mapper, library_name=lib)
                coverage_files.append(
                    "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                        name_pattern=path_pattern, ext=tsv_ext
                    )
                )
        yield tsv_ext, coverage_files

        # Yield ploidy files
        path_pattern = ploidy_path_pattern.format(**wildcards)
        yield ploidy_ext, "work/{name_pattern}/out/{name_pattern}/.done".format(
            name_pattern=path_pattern
        )

    @dictify
    def _get_input_files_post_germline_calls(self, wildcards):
        """Yield input files for ``post_germline_calls`` in CASE mode.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_name' (e.g., 'P001-N1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards
        """
        # Get library kit associated with library
        library_kit = self.ngs_library_to_kit.get(wildcards.library_name)

        # Get shards - based on scattered step
        model_dir_dict = None
        for model in self.config["gcnv"]["precomputed_model_paths"]:
            # Adjust library name to wildcard
            library_to_wildcard = model.get("library").strip().replace(" ", "_")
            if library_to_wildcard == library_kit:
                pattern = model.get("model_pattern")
                model_dir_dict = self._get_model_dir_to_dict(pattern)
                break

        # Validate model
        if not model_dir_dict:
            msg_error = f"Configuration doesn't contain model for library '{library_kit}'"
            raise InvalidConfiguration(msg_error)

        # Yield cnv calls output
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield "calls", [
            "work/{name_pattern}.{shard}/out/{name_pattern}.{shard}/.done".format(
                name_pattern=name_pattern, shard=shard
            )
            for shard in model_dir_dict
        ]

        # Yield contig-ploidy output
        ext = "ploidy"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield ext, "work/{name_pattern}/out/{name_pattern}/.done".format(name_pattern=name_pattern)

    @listify
    def _get_input_files_merge_cohort_vcfs(self, wildcards):
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = "{mapper}.gcnv_post_germline_calls.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                yield "work/{name_pattern}/out/{name_pattern}.vcf.gz".format(
                    name_pattern=name_pattern
                )

    @dictify
    def _get_input_files_extract_ped(self, wildcards):
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        name_pattern = "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        for key, ext in (("vcf", ".vcf.gz"), ("tbi", ".vcf.gz.tbi")):
            yield key, "work/{name_pattern}/out/{name_pattern}{ext}".format(
                name_pattern=name_pattern, ext=ext
            )

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    @staticmethod
    @dictify
    def _get_output_files_contig_ploidy():
        """Yield dictionary with output files for ``contig_ploidy`` rule in CASE MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}"
        yield ext, touch(
            "work/{name_pattern}/out/{name_pattern}/.{ext}".format(
                name_pattern=name_pattern, ext=ext
            )
        )

    @staticmethod
    @dictify
    def _get_output_files_call_cnvs():
        """Yield dictionary with output files for ``call_cnvs`` rule in CASE MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}"
        yield ext, touch(
            "work/{name_pattern}/out/{name_pattern}/.{ext}".format(
                name_pattern=name_pattern, ext=ext
            )
        )

    @staticmethod
    @dictify
    def _get_output_files_post_germline_calls():
        name_pattern = "{mapper}.gcnv_post_germline_calls.{library_name}"
        pairs = {"ratio_tsv": ".ratio.tsv", "itv_vcf": ".interval.vcf.gz", "seg_vcf": ".vcf.gz"}
        for key, ext in pairs.items():
            yield key, touch(
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                )
            )

    @staticmethod
    @dictify
    def _get_output_files_merge_cohort_vcfs():
        name_pattern = "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}"
        pairs = {
            "vcf": ".vcf.gz",
            "vcf_md5": ".vcf.gz.md5",
            "tbi": ".vcf.gz.tbi",
            "tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in pairs.items():
            yield key, "work/{name_pattern}/out/{name_pattern}{ext}".format(
                name_pattern=name_pattern, ext=ext
            )

    @staticmethod
    @dictify
    def _get_output_files_extract_ped():
        name_pattern = "{mapper}.gcnv.{library_name}"
        kvs = (
            ("vcf", ".vcf.gz"),
            ("vcf_md5", ".vcf.gz.md5"),
            ("tbi", ".vcf.gz.tbi"),
            ("tbi_md5", ".vcf.gz.tbi.md5"),
        )
        for key, suffix in kvs:
            yield key, "work/{name_pattern}/out/{name_pattern}{suffix}".format(
                name_pattern=name_pattern, suffix=suffix
            )