# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_gene_fusion_calling`` step

The somatic_gene_fusion calling step allows for the detection of gene fusions from RNA-seq data
in cancer.  The wrapped tools start at the raw RNA-seq reads and generate filtered lists of
predicted gene fusions.

==========
Step Input
==========

Gene fusion calling starts at the raw RNA-seq reads.  Thus, the input is very similar to one of
:ref:`ngs_mapping step <step_ngs_mapping>`.

See :ref:`ngs_mapping_step_input` for more information.

===========
Step Output
===========

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_gene_fusion_calling.rst

=============================
Available Gene Fusion Callers
=============================

- ``fusioncatcher``

"""

import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import touch

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    LinkOutStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: HLA typing tools
GENE_FUSION_CALLERS = ("fusioncatcher", "jaffa", "pizzly", "hera", "star_fusion")

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  somatic_gene_fusion_calling:
    path_link_in: ""  # OPTIONAL Override data set configuration search paths for FASTQ files
    tools: ['fusioncatcher', 'jaffa', 'arriba']
    fusioncatcher:
      data_dir: REQUIRED   # REQUIRED
      configuration: null  # optional
      num_threads: 16
    pizzly:
      kallisto_index: REQUIRED    # REQUIRED
      transcripts_fasta: REQUIRED # REQUIRED
      annotations_gtf: REQUIRED       # REQUIRED
      kmer_size: 31
    hera:
      path_index: REQUIRED   # REQUIRED
      path_genome: REQUIRED  # REQUIRED
    star_fusion:
      path_ctat_resource_lib: REQUIRED
    defuse:
      path_dataset_directory: REQUIRED
    arriba:
      path_index: REQUIRED       # REQUIRED  STAR path index (preferably 2.7.10 or later)
      features: REQUIRED         # REQUIRED  Gene features (for ex. ENCODE or ENSEMBL) in gtf format
      blacklist: ""              # optional (provided in the arriba distribution, see /fast/work/groups/cubi/projects/biotools/static_data/app_support/arriba/v2.3.0)
      known_fusions: ""          # optional
      tags: ""                   # optional (can be set to the same path as known_fusions)
      structural_variants: ""    # optional
      protein_domains: ""        # optional
      num_threads: 8
      trim_adapters: false
      num_threads_trimming: 2
      star_parameters:
      - " --outFilterMultimapNmax 50"
      - " --peOverlapNbasesMin 10"
      - " --alignSplicedMateMapLminOverLmate 0.5"
      - " --alignSJstitchMismatchNmax 5 -1 5 5"
      - " --chimSegmentMin 10"
      - " --chimOutType WithinBAM HardClip"
      - " --chimJunctionOverhangMin 10"
      - " --chimScoreDropMax 30"
      - " --chimScoreJunctionNonGTAG 0"
      - " --chimScoreSeparation 1"
      - " --chimSegmentReadGapMax 3"
      - " --chimMultimapNmax 50"
""".lstrip()


class SomaticGeneFusionCallingStepPart(BaseStepPart):
    """Base class for somatic gene fusion calling"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{name}.{{library_name}}/out/.done".format(name=self.name)
        # Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config["path_link_in"],
        )

    @dictify
    def get_input_files(self, action):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        yield "done", os.path.join(self.base_path_in, ".done")

    @dictify
    def get_output_files(self, action):
        """Return output files that all read mapping sub steps must return (BAM + BAI file)"""
        # Validate action
        self._validate_action(action)
        yield "done", touch(self.base_path_out)

    def get_log_file(self, action):
        """Return path to log file"""
        # Validate action
        self._validate_action(action)
        return "work/{name}.{{library_name}}/log/snakemake.gene_fusion_calling.log".format(
            name=self.name
        )

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        _ = library_name
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config["path_link_in"]:
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)


class FusioncatcherStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Fusioncatcher"""

    #: Step name
    name = "fusioncatcher"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            tumor = flatten(
                zip(
                    sorted(self._collect_reads(wildcards, wildcards.library_name, "")),
                    sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")),
                )
            )
            tumor = list(map(os.path.abspath, tumor))
            return {"normal": [], "tumor": tumor}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{7500 * 4}M",
        )


class JaffaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using JAFFA"""

    #: Step name
    name = "jaffa"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{40 * 1024 * 4}M",
        )


class PizzlyStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Kallisto+Pizzly"""

    #: Step name
    name = "pizzly"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{20 * 1024 * 4}M",
        )


class StarFusionStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using STAR-Fusion"""

    #: Step name
    name = "star_fusion"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{30 * 1024 * 4}M",
        )


class DefuseStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Defuse"""

    #: Step name
    name = "defuse"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            time="5-00:00:00",  # 5 days
            memory=f"{10 * 1024 * 8}M",
        )


class HeraStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Hera"""

    #: Step name
    name = "hera"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            time="5-00:00:00",  # 5 days
            memory=f"{20 * 1024 * 8}M",
        )


class ArribaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using arriba"""

    #: Step name
    name = "arriba"

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute STAR & arriba are present in the
        configuration. It further checks that the provided index has all the expected file
        extensions. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if self.__class__.name not in self.config["tools"]:
            return  # arriba not run, don't check configuration  # pragma: no cover

        # Check required configuration settings present
        self.parent.ensure_w_config(
            config_keys=("step_config", "somatic_gene_fusion_calling", "arriba", "path_index"),
            msg="Path to STAR indices is required",
        )
        self.parent.ensure_w_config(
            config_keys=("step_config", "somatic_gene_fusion_calling", "arriba", "features"),
            msg="Path to genomic features gtf file is required",
        )

        # Check that the path to the STAR index is valid.
        for fn in ("Genome", "SA", "SAindex"):
            expected_path = self.config["arriba"]["path_index"] + "/" + fn
            if not os.path.exists(expected_path):  # pragma: no cover
                tpl = "Expected STAR indices input path {expected_path} does not exist!".format(
                    expected_path=expected_path
                )
                raise InvalidConfiguration(tpl)

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        base_path_out = "work/{name}.{{library_name}}/out/{name}.{{library_name}}.{ext}"
        key_ext = (
            ("fusions", "fusions.tsv"),
            ("discarded", "discarded_fusions.tsv.gz"),
        )
        for key, ext in key_ext:
            yield key, base_path_out.format(name=self.name, ext=ext)
            yield key + "_md5", base_path_out.format(name=self.name, ext=ext) + ".md5"
        yield "done", "work/arriba.{library_name}/out/.done"

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        _ = action
        prefix = "work/{name}.{{library_name}}/log/{name}.{{library_name}}".format(name=self.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action):
        """Get resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["arriba"]["num_threads"], time="24:00:00", memory=f"{64 * 1024}M"
        )  # 1 day


class SomaticGeneFusionCallingWorkflow(BaseStep):
    """Perform somatic gene fusion calling"""

    #: Workflow name
    name = "somatic_gene_fusion_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(workflow, config, config_lookup_paths, config_paths, workdir)
        self.register_sub_step_classes(
            (
                FusioncatcherStepPart,
                JaffaStepPart,
                PizzlyStepPart,
                HeraStepPart,
                StarFusionStepPart,
                DefuseStepPart,
                ArribaStepPart,
                LinkInStep,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        # Convert sheet parsing into method
        library_names_list = self._get_all_rna_ngs_libraries()
        # Get results
        name_pattern = "{fusion_caller}.{ngs_library.name}"
        for fusion_caller in self.config["tools"]:
            for ngs_library in library_names_list:
                # Constant to all callers
                name_pattern_value = name_pattern.format(
                    fusion_caller=fusion_caller, ngs_library=ngs_library
                )
                yield os.path.join("output", name_pattern_value, "out", ".done")
                # Caller specific stuff...
                if fusion_caller == "arriba":
                    yield from self._yield_arriba_files(ngs_library)
                else:
                    yield os.path.join(
                        "output", name_pattern_value, "log", "snakemake.gene_fusion_calling.log"
                    )

    def _yield_arribe_files(self, ngs_library):
        tpl = "output/arriba.{library_name}/out/arriba.{library_name}.{ext}"
        for ext in ("fusions.tsv", "discarded_fusions.tsv.gz"):
            yield tpl.format(library_name=ngs_library.name, ext=ext)
            yield tpl.format(library_name=ngs_library.name, ext=ext + ".md5")
        tpl = "output/arriba.{library_name}/log/arriba.{library_name}.{ext}"
        for ext in ("log", "conda_list", "conda_info"):
            yield tpl.format(library_name=ngs_library.name, ext=ext)
            yield tpl.format(library_name=ngs_library.name, ext=ext + ".md5")

    def check_config(self):
        """Check that the required configurations are present."""
        # TODO: implement check for REQUIRED configurations.
