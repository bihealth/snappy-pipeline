# -*- coding: utf-8 -*-
"""Implementation of the gCNV CASE mode run methods.
"""

import glob
from itertools import chain
import json
import os
import re
import warnings

from snakemake.io import Wildcards, expand, touch

from snappy_pipeline.base import InvalidConfiguration, UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.gcnv.gcnv_common import (
    GcnvCommonStepPart,
    InconsistentLibraryKitsWarning,
)

#: Predefined file name keys and extensions for the overall results files.
RESULT_EXTENSIONS = {
    "vcf": ".vcf.gz",
    "vcf_md5": ".vcf.gz.md5",
    "vcf_tbi": ".vcf.gz.tbi",
    "vcf_tbi_md5": ".vcf.gz.tbi.md5",
}

#: Predefined file name keys and extensions for log files.
LOG_EXTENSIONS = {
    "conda_info": ".conda_info.txt",
    "conda_info_md5": ".conda_info.txt.md5",
    "conda_list": ".conda_list.txt",
    "conda_list_md5": ".conda_list.txt.md5",
    "wrapper": ".wrapper.py",
    "wrapper_md5": ".wrapper.py.md5",
    "env_yaml": ".environment.yaml",
    "env_yaml_md5": ".environment.yaml.md5",
    "log": ".log",
    "log_md5": ".log.md5",
}


def get_model_dir_list(pattern):
    """Get model directories.

    :param pattern: Pattern of model directory paths. Expects models to be based on scattered
    step: ``gatk IntervalListTools``.
    :type pattern: str

    :return: Returns list with all directories that match the inputted pattern.
    """
    return [path_ for path_ in glob.glob(pattern) if os.path.isdir(path_)]


def get_model_dir_to_dict(pattern):
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
    path_list = get_model_dir_list(pattern)
    # Populate dictionary
    # Assumption: if no group, a single path was provided instead of a pattern.
    for path in path_list:
        try:
            key = re.search(re_pattern, path).group(1)
            out_dict[key] = path
        except IndexError:
            out_dict[default_key] = path
    return out_dict


class ValidationMixin:
    """Mixin that provides validation methods for input model"""

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
                path_list = get_model_dir_list(pattern=model_pattern)
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
                pretty_model = self._pretty_print_config(config=model)
                msg = (
                    "Provided configuration not as expected...\n"
                    f"\nn_keys_pass={n_keys_pass}, keys_pass={keys_pass}, values_pass={values_pass}\n"
                    f"Expected:\n{expected_format}\nObserved:\n{pretty_model}\n"
                )
                raise InvalidConfiguration(msg)

    def _pretty_print_config(self, config):
        """Pretty format configuration.

        :param config: Configuration dictionary to be formatted.
        :type config: OrderedDict

        :return: Configuration as a nicely formatted string.
        """
        return str(json.dumps(config, sort_keys=False, indent=4))

    def validate_ploidy_model_directory(self, path):
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

    def validate_call_model_directory(self, path):
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


class ContigPloidyMixin:
    """Methods for ``contig_ploidy``."""

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards: Wildcards):
        """Yield input files for ``contig_ploidy`` rule  in CASE MODE using precomputed model.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        """
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = f"{wildcards.mapper}.gcnv_coverage.{lib}"
                tsvs.append(f"work/{name_pattern}/out/{name_pattern}.{ext}")
        yield ext, tsvs

    @dictify
    def _get_output_files_contig_ploidy(self):
        """Yield dictionary with output files for ``contig_ploidy`` rule in CASE MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}"
        yield ext, touch(f"work/{name_pattern}/out/{name_pattern}/.{ext}")

    def _get_params_contig_ploidy(self, wildcards: Wildcards):
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


class CallCnvsMixin:
    """Methods for ``call_cnvs``."""

    @dictify
    def _get_input_files_call_cnvs(self, wildcards: Wildcards):
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
                name_pattern = tsv_path_pattern.format(mapper=wildcards.mapper, library_name=lib)
                coverage_files.append(f"work/{name_pattern}/out/{name_pattern}.{tsv_ext}")
        yield tsv_ext, coverage_files

        # Yield ploidy files
        name_pattern = ploidy_path_pattern.format(**wildcards)
        yield ploidy_ext, f"work/{name_pattern}/out/{name_pattern}/.done"

    @dictify
    def _get_output_files_call_cnvs(self):
        """Yield dictionary with output files for ``call_cnvs`` rule in CASE MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}"
        yield ext, touch(f"work/{name_pattern}/out/{name_pattern}/.{ext}")

    def _get_params_call_cnvs(self, wildcards):
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
                model_dir_dict = get_model_dir_to_dict(pattern)
                path = model_dir_dict.get(wildcards.shard)
        return {"model": path}


class PostGermlineCallsMixin:
    """Methods for ``post_germline_calls``."""

    @dictify
    def _get_output_files_post_germline_calls(self):
        name_pattern = "{mapper}.gcnv_post_germline_calls.{library_name}"
        extensions = {
            "ratio_tsv": ".ratio.tsv",
            "itv_vcf": ".interval.vcf.gz",
            "seg_vcf": ".vcf.gz",
        }
        for key, ext in extensions.items():
            yield key, touch(f"work/{name_pattern}/out/{name_pattern}{ext}")

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
                model_dir_dict = get_model_dir_to_dict(pattern)
                break

        # Validate model
        if not model_dir_dict:
            msg_error = f"Configuration doesn't contain model for library '{library_kit}'"
            raise InvalidConfiguration(msg_error)

        # Yield cnv calls output
        name_pattern = f"{wildcards.mapper}.gcnv_call_cnvs.{library_kit}"
        yield "calls", [
            f"work/{name_pattern}.{shard}/out/{name_pattern}.{shard}/.done"
            for shard in model_dir_dict
        ]

        # Yield contig-ploidy output
        ext = "ploidy"
        name_pattern = f"{wildcards.mapper}.gcnv_contig_ploidy.{library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}/.done"

    def _get_params_post_germline_calls(self, wildcards):
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
                model_dir_dict = get_model_dir_to_dict(pattern)
                paths = list(model_dir_dict.values())
        return {"model": paths}


class JointGermlineCnvSegmentationMixin:
    """Methods for ``post_germline_calls``."""

    @dictify
    def _get_output_files_joint_germline_cnv_segmentation(self):
        name_pattern = "{mapper}.gcnv.{library_name}"
        for key, suffix in RESULT_EXTENSIONS.items():
            yield key, f"work/{name_pattern}/out/{name_pattern}{suffix}"

    @dictify
    def _get_log_file_joint_germline_cnv_segmentation(self):
        """Return log file **pattern** for the step ``joint_germline_cnv_segmentation``."""
        name_pattern = "{mapper}.gcnv.{library_name}"
        for key, ext in LOG_EXTENSIONS.items():
            yield key, f"work/{name_pattern}/log/{name_pattern}.joint_germline_segmentation{ext}"

    @dictify
    def _get_input_files_joint_germline_cnv_segmentation(self, wildcards):
        # Yield list of paths to input VCF files
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        ped_ngs_library_names = [
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        ]
        vcfs = []
        for library_name in sorted(ped_ngs_library_names):
            name_pattern = f"{wildcards.mapper}.gcnv_post_germline_calls.{library_name}"
            vcfs.append(f"work/{name_pattern}/out/{name_pattern}.vcf.gz")
        yield "vcf", vcfs
        # Yield path to interval list file
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        name_pattern = f"gcnv_preprocess_intervals.{library_kit}"
        yield "interval_list", f"work/{name_pattern}/out/{name_pattern}.interval_list"
        # Yield path to pedigree file
        name_pattern = f"write_pedigree.{wildcards.library_name}"
        yield "ped", f"work/{name_pattern}/out/{wildcards.library_name}.ped"


class RunGcnvStepPart(
    ValidationMixin,
    ContigPloidyMixin,
    CallCnvsMixin,
    PostGermlineCallsMixin,
    JointGermlineCnvSegmentationMixin,
    GcnvCommonStepPart,
):
    """Class with methods to run GATK4 gCNV calling in CASE mode"""

    #: Class available actions
    actions = (
        "preprocess_intervals",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "joint_germline_cnv_segmentation",
    )

    def __init__(self, parent):
        super().__init__(parent)
        # Validate configuration, precomputed models must be present
        self.validate_request()

    def get_params(self, action: str):
        """
        :param action: Action (i.e., step) in the workflow. Currently available for:
        'ploidy-model', 'model', and 'post_germline_calls'.

        :return: Returns input function for gCNV rule based on inputted action.

        :raises UnsupportedActionException: if invalid action.
        """
        # Actions with parameters
        valid_actions = ("call_cnvs", "contig_ploidy", "post_germline_calls")
        # Validate inputted action
        if action not in valid_actions:
            error_message = f"Action '{action}' is not supported. Valid options: {valid_actions}."
            raise UnsupportedActionException(error_message)

        # Return requested function
        return getattr(self, f"_get_params_{action}")

    @listify
    def get_result_files(self):
        """Return list of **concrete** paths to result files for the given configuration and sample sheets.

        The function will skip pedigrees where samples have inconsistent library kits and print a warning.
        """

        def path_work_to_output(work_path):
            """Helper to convert a work to an output path"""
            return re.sub(r"^work/", "output/", work_path)

        # Get list with all result path template strings (output and log files).  This is done using the
        # functions generating the output patterns of the joint germline CNV segmentation step.
        result_path_tpls = list(
            chain(
                self._get_output_files_joint_germline_cnv_segmentation().values(),  # output files
                self._get_log_file_joint_germline_cnv_segmentation().values(),  # log files
            )
        )

        # Iterate over all pedigrees.  Check library kit consistency for each pedigree.  In case of inconsistent
        # library kits within one pedigree, raise a warning and skip this pedigree.
        for index_library_name, pedigree in self.index_ngs_library_to_pedigree.items():
            # Obtain all library names
            library_names = [
                donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
            ]
            library_kits = [self.ngs_library_to_kit[library_name] for library_name in library_names]
            if len(set(library_kits)) != 1:
                names_kits = list(zip(library_names, library_kits))
                msg = (
                    "Found inconsistent library kits (more than one kit!) for pedigree with index "
                    f"{index_library_name}.  The library name/kit pairs are {names_kits}.  This pedigree "
                    "will be SKIPPED in CNV call generation."
                )
                warnings.warn(InconsistentLibraryKitsWarning(msg))
                continue
            else:
                # The library kits are consistent in the pedigree.  Yield all concrete output paths, replacing
                # prefix "work/" by "output/".
                for path_tpl in result_path_tpls:
                    yield from map(
                        path_work_to_output,
                        expand(
                            path_tpl,
                            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                            library_name=[index_library_name],
                        ),
                    )
