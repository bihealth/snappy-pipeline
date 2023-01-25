"""Implementation of the ``ngs_mapping`` step

The ngs_mapping step allows for the alignment of NGS data with standard read mappers, such as BWA
for DNA data and STAR for RNA data.  Also, it provides functionality to compute post-alignment
statistics, such as the coverage of target (e.g., exome or panel) regions.

There is a distinction made between "normal" DNA reads (short reads from Illumina) and "long"
DNA reads, such as PacBio/Oxford Nanopore.  Again, the NGS mapping step will perform alignment
of all NGS libraries.

The precise actions that are performed in the alignment are defined by the wrappers (e.g., the
``bwa`` or ``star``) wrappers.  Generally, this includes converting into BAM format, sorting
by coordinate, an indexing using a BAI file.  For short reads, this can include marking of
duplicates using Samblaster and depends on the actual configuration (see below for the default
configuration).

==========
Properties
==========

overall stability

    **stable**

applicable to

    germline and somatic read alignment

generally applicable to

    short and long read DNA and RNA sequencing

.. _ngs_mapping_step_input:

==========
Step Input
==========

For each library defined in all sample sheets, the instances of this step will search for the input
files according to the configuration.  The found read files will be linked into
``work/input_links/{library_name}`` (status quo, not a output path, thus path not guaranteed
to be stable between minor versions).

This is different to the other steps that use the output of previous steps for their input.

----------------------
Data Set Configuration
----------------------

Consider the following data set definition from the main configuration file.

.. code-block:: yaml

    data_sets:
      first_batch:
        file: 01_first_batch.json
        search_patterns:
          # Note that currently only "left" and "right" key known
          - {'left': '*/L???/*_R1.fastq.gz', 'right': '*/L???/*_R2.fastq.gz'}
        search_paths: ['../input/01_first_batch']

Here, the data set ``first_batch`` is defined.  The sample sheet file is named
``01_first_batch.json`` and looked for in the relative path to the configuration file.  The input
search will be start in the (one, but could be more than one) path ``../input/01_first_batch``
(relative to the directory containing the configuration file).  The sample sheet provides a
``folderName`` ``extraInfo`` entry for each NGS library.  This folder name is searched for (e.g.,
``P001-N1-DNA1-WES``).  Once such a folder is found, the patterns in the values of the dict
``search_patterns`` are used for locating the paths of the actual files.

Currently, the only supported keys in the ``search_patterns`` dict are ``"left"`` and ``"right""``
(the lattern can be omitted when only searching for single-end reads).

Consider the following example::

    ../input/
    `-- 01_first_batch
        |-- P001-N1-DNA1-WES1
        |   `-- 42KF5AAXX
        |       `-- L001
        |           |-- P001-N1-DNA1-WES1_R1.fastq.gz
        |           |-- P001-N1-DNA1-WES1_R1.fastq.gz.md5
        |           |-- P001-N1-DNA1-WES1_R2.fastq.gz
        |           `-- P001-N1-DNA1-WES1_R2.fastq.gz.md5
        [...]

Here, the folder ``01_first_batch`` will be searched for a directory named ``P001-N1-DNA1-WES``.
After finding, the relative paths ``42KF5AAXX/L001/P001-N1-DNA1-WES1_R1.fastq.gz`` and
``42KF5AAXX/L001/P001-N1-DNA1-WES1_R2.fastq.gz`` will be found and used for the left/right parts of
a paired read set.

Mixing Single-End and Paired-End Reads
======================================

By default, it is checked that for each ``search_pattern``, the same number of matching files
has to be found, otherwise directories are ignored.  The reason is to reduce the number of
possible errors when linking in files.  You can change this behaviour by specifying
``mixed_se_pe: True`` in the data set information.  Then, it will be allowed to have the matches
for the ``right`` entry to be empty.  However, you will need to consistently have either SE or
PE data for each library; it is allowed to mix SE and PE libraries within one project but not
to have PE and SE data for one library.

===========
Step Output
===========

For each NGS library with name ``library_name`` and each read mapper ``mapper`` that the library
has been aligned with, the pipeline step will create a directory
``output/{mapper}.{library_name}/out`` with symlinks of the following names to the resulting
sorted BAM files with corresponding BAI and MD5 files.

- ``{mapper}.{library_name}.bam``
- ``{mapper}.{library_name}.bam.bai``
- ``{mapper}.{library_name}.bam.md5``
- ``{mapper}.{library_name}.bam.bai.md5``

In addition, several tools are used to automatically generate reports based on the BAM and BAI
files.  See the Reports section below for more details

The BAM files are only postprocessed if configured so.

.. note::

    In contrast to other pipeline steps, the NGS mapping step will also generate the BAM files for
    the background data sets as there are currently problems with Snakemake sub workflows and
    input functions.

====================
Global Configuration
====================

- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_ngs_mapping.rst

======================
Available Read Mappers
======================

The following read mappers are available for the alignment of DNA-seq and RNA-seq reads.

- (short/Illumina) DNA
    - ``"bwa"``
    - ``"external"``

- (short/Illumina) RNA-seq
    - ``"star"``
    - ``"external"``

- (long/PacBio/Nanopore) DNA
    - ``"minimap2"``
    - ``"external"``

=======
Reports
=======

Currently, the following reports are generated based on the BAM and BAI file output by this step.

General Alignment Statistics (.txt)
  The tools ``samtools bamstats``, ``samtools flagstats`` and ``samtools idxstats`` are always
  called by default, and are linked out into the ``output/{mapper}.{library_name}/report/bam_qc``
  directory. The file names for these reports (and their MD5s) use the following naming convention:

  - ``{mapper}.{library_name}.bamstats.txt``
  - ``{mapper}.{library_name}.flagstats.txt``
  - ``{mapper}.{library_name}.idxstats.txt``
  - ``{mapper}.{library_name}.bamstats.txt.md5``
  - ``{mapper}.{library_name}.flagstats.txt.md5``
  - ``{mapper}.{library_name}.idxstats.txt.md5``

For example, it will look as follows for the example bam files shown above:

::

    output/
    +-- bwa.P001-N1-DNA1-WES1
    |   |-- out
    |   |   |-- bwa.P001-N1-DNA1-WES1.bam
    |   |   |-- bwa.P001-N1-DNA1-WES1.bam.bai
    |   |   |-- bwa.P001-N1-DNA1-WES1.bam.bai.md5
    |   |   `-- bwa.P001-N1-DNA1-WES1.bam.md5
    |   `-- report
    |       `-- bam_qc
    |           |-- bwa.P001-N1-DNA1-WES1.bam.bamstats.txt
    |           |-- bwa.P001-N1-DNA1-WES1.bam.bamstats.txt.md5
    |           |-- bwa.P001-N1-DNA1-WES1.bam.flagstats.txt
    |           |-- bwa.P001-N1-DNA1-WES1.bam.flagstats.txt.md5
    |           |-- bwa.P001-N1-DNA1-WES1.bam.idxstats.txt
    |           `-- bwa.P001-N1-DNA1-WES1.bam.idxstats.txt.md5
    [...]

Target Coverage Report (.txt)
  If ``ngs_mapping/path_target_regions`` is set to a BED file with the target regions
  (either capture regions of capture kits in the case of targeted sequencing or exons for WES/WGS
  sequencing) a target coverage report is generated and linked out into the
  ``output/{mapper}.{library_name}/report/cov_qc``
  directory. The file names for these reports (and their MD5s) use the following naming convention:

  - ``{mapper}.{library_name}.txt``
  - ``{mapper}.{library_name}.txt.md5``

  For example, it will look as follows for the example bam files shown above:

::

    output/
    +-- bwa.P001-N1-DNA1-WES1
    |   `-- report
    |       |-- bam_qc
    |       |   [...]
    |       `-- cov_qc
    |           |-- bwa.P001-N1-DNA1-WES1.txt
    |           `-- bwa.P001-N1-DNA1-WES1.txt.md5
    [...]


Genome-wide Coverage Count (.bed.gz)
  If ``ngs_mapping/compute_coverage_bed`` to be set to ``true`` a report is generated
  that gives the depth at each base of the genome. (note: currently this report only appears
  in ``work/`` and is not yet linked out into the ``output/`` directory).

  (TODO: add file name rules and example)

::

    work/
    +-- bwa.P001-N1-DNA1-WES1
    |   `-- report
    |       `-- bam_qc
    |           |-- bwa.P001-N1-DNA1-WES1.bam.bamstats.d
    |           |   |-- acgt-cycles.gp
    |           |   |-- acgt-cycles.png
    |           |   |-- coverage.gp
    |           |   |-- coverage.png
    |           |   |-- gc-content.gp
    |           |   |-- gc-content.png
    |           |   |-- gc-depth.gp
    |           |   |-- gc-depth.png
    |           |   |-- indel-cycles.gp
    |           |   |-- indel-cycles.png
    |           |   |-- indel-dist.gp
    |           |   |-- indel-dist.png
    |           |   |-- index.html
    |           |   |-- insert-size.gp
    |           |   |-- insert-size.png
    |           |   |-- quals2.gp
    |           |   |-- quals2.png
    |           |   |-- quals3.gp
    |           |   |-- quals3.png
    |           |   |-- quals.gp
    |           |   |-- quals-hm.gp
    |           |   |-- quals-hm.png
    |           |   `-- quals.png
    |           [...]
    [...]


"""

from itertools import chain
import os
import re
import sys

from biomedsheets.shortcuts import GenericSampleSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import InvalidConfiguration, UnsupportedActionException
from snappy_pipeline.utils import DictQuery, dictify, flatten, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    ResourceUsage,
    get_ngs_library_folder_name,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# TODO: Need something smarter still for @RG

#: Extensions of files to create as main payload
EXT_VALUES = (".bam", ".bam.bai", ".bam.md5", ".bam.bai.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("bam", "bai", "bam_md5", "bai_md5")

#: Available read mappers for (short/Illumina) DNA-seq data
READ_MAPPERS_DNA = ("bwa", "bwa_mem2")

#: Available read mappers for (short/Illumina) RNA-seq data
READ_MAPPERS_RNA = ("star",)

#: Available read mappers for (long/PacBio/Nanopoare) DNA-seq data
READ_MAPPERS_DNA_LONG = ("minimap2",)

#: Default configuration
DEFAULT_CONFIG = r"""
step_config:
  ngs_mapping:
    # Aligners to use for the different NGS library types
    tools:
      dna: []      # Required if DNA analysis; otherwise, leave empty. Example: 'bwa'.
      rna: []      # Required if RNA analysis; otherwise, leave empty. Example: 'star'.
      dna_long: [] # Required if long-read mapper used; otherwise, leave empty. Example: 'minimap2'.
    path_link_in: ""   # OPTIONAL Override data set configuration search paths for FASTQ files
    # Thresholds for targeted sequencing coverage QC.  Enabled by specifying
    # the path_arget_regions setting above
    target_coverage_report:
      # Mapping from enrichment kit to target region BED file, for either computing per--target
      # region coverage or selecting targeted exons.
      #
      # The following will match both the stock IDT library kit and the ones
      # with spike-ins seen fromr Yale genomics.  The path above would be
      # mapped to the name "default".
      # - name: IDT_xGen_V1_0
      #   pattern: "xGen Exome Research Panel V1\\.0*"
      #   path: "path/to/targets.bed"
      path_target_interval_list_mapping: []
      # Maximal/minimal/warning coverage
      max_coverage: 200
      min_cov_warning: 20  # >= 20x for WARNING
      min_cov_ok: 50  # >= 50x for OK
      detailed_reporting: false  # per-exon details (cannot go into multiqc)
    # Depth of coverage collection, mainly useful for genomes.
    bam_collect_doc:
      enabled: false
      window_length: 1000
    # Compute fingerprints with ngs-chew
    ngs_chew_fingerprint:
      enabled: true
    # Configuration for BWA
    bwa:
      path_index: REQUIRED # Required if listed in ngs_mapping.tools.dna; otherwise, can be removed.
      num_threads_align: 16
      num_threads_trimming: 8
      num_threads_bam_view: 4
      num_threads_bam_sort: 4
      memory_bam_sort: 4G
      trim_adapters: false
      mask_duplicates: true
      split_as_secondary: false  # -M flag
    # Configuration for BWA-MEM2
    bwa_mem2:
      path_index: REQUIRED # Required if listed in ngs_mapping.tools.dna; otherwise, can be removed.
      bwa_mode: auto  # in ['auto', 'bwa-aln', 'bwa-mem']
      num_threads_align: 16
      num_threads_trimming: 8
      num_threads_bam_view: 4
      num_threads_bam_sort: 4
      memory_bam_sort: 4G
      trim_adapters: false
      mask_duplicates: true
      split_as_secondary: true  # -M flag
    # Configuration for STAR
    star:
      path_index: REQUIRED # Required if listed in ngs_mapping.tools.rna; otherwise, can be removed.
      path_features: ""    # Required for computing gene counts
      num_threads_align: 16
      num_threads_trimming: 8
      num_threads_bam_view: 4
      num_threads_bam_sort: 4
      memory_bam_sort: 4G
      genome_load: NoSharedMemory
      raw_star_options: ''
      align_intron_max: 1000000                # ENCODE option
      align_intron_min: 20                     # ENCODE option
      align_mates_gap_max: 1000000             # ENCODE option
      align_sjdb_overhang_min: 1               # ENCODE option
      align_sj_overhang_min: 8                 # ENCODE option
      out_filter_mismatch_n_max: 999           # ENCODE option
      out_filter_mismatch_n_over_l_max: 0.04   # ENCODE option
      out_filter_multimap_n_max: 20            # ENCODE option
      out_filter_type: BySJout                 # ENCODE option
      out_filter_intron_motifs: None    # or for cufflinks: RemoveNoncanonical
      out_sam_strand_field: None        # or for cufflinks: intronMotif
      transcriptome: false              # true to output transcript coordinate bam for RSEM
      trim_adapters: false
      mask_duplicates: false
      include_unmapped: true
    # Configuration for Minimap2
    minimap2:
      mapping_threads: 16
"""


class MappingGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for mapping steps"""

    tool_category = None

    def skip_result_files_for_library(self, library_name: str) -> bool:
        """Override function such that the mapper is applicable to the library"""
        extra_infos = self.parent.ngs_library_to_extra_infos[library_name]
        extraction_type = extra_infos.get("extractionType", "DNA").lower()
        if extra_infos["seqPlatform"] in ("ONT", "PacBio"):
            suffix = "_long"
        else:
            suffix = ""
        library_tool_category = f"{extraction_type}{suffix}"
        if self.tool_category not in ("__any__", library_tool_category):
            return True
        else:
            return self.name not in self.config["tools"][library_tool_category]

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return a list of all paths with prefix ``output/` that are
        returned by ``self.get_output_files()`` for all actions in ``self.actions``.

        You can override the ``skip_result_files_for_library()`` method to skip result files for
        a library.
        """
        # Skip if step part has a tool category and it is not enabled
        if (
            self.tool_category != "__any__"
            and self.name not in self.config["tools"][self.tool_category]
        ):
            return

        for action in self.actions:
            # Get list of result path templates.
            output_files_tmp = self.get_output_files(action)
            if isinstance(output_files_tmp, dict):
                output_files = output_files_tmp.values()
            else:
                output_files = output_files_tmp
            result_paths_tpls = list(
                filter(
                    lambda p: p.startswith("output/"),
                    flatten(output_files),
                )
            )
            #: Generate all concrete output paths.
            for path_tpl in result_paths_tpls:
                for library_name in self.parent.ngs_library_to_extra_infos.keys():
                    if not self.skip_result_files_for_library(library_name):
                        yield from expand(path_tpl, mapper=[self.name], library_name=library_name)


class ReportGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for report generation steps"""

    def skip_result_files_for_library(self, library_name: str) -> bool:
        if not getattr(self, "tool_categories", None):
            return False

        extra_infos = self.parent.ngs_library_to_extra_infos[library_name]
        extraction_type = extra_infos.get("extractionType", "DNA").lower()
        if extra_infos["seqPlatform"] in ("ONT", "PacBio"):
            suffix = "_long"
        else:
            suffix = ""
        library_tool_category = f"{extraction_type}{suffix}"
        return library_tool_category not in self.tool_categories

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return a list of all paths with prefix ``output/` that are
        returned by ``self.get_output_files()`` for all actions in ``self.actions``.

        You can override the ``skip_result_files_for_library()`` method to skip result files for
        a library.
        """
        # Skip if step part has a tool category and it is not enabled
        for action in self.actions:
            # Get list of result path templates.
            output_files_tmp = self.get_output_files(action)
            if isinstance(output_files_tmp, dict):
                output_files = output_files_tmp.values()
            else:
                output_files = output_files_tmp
            result_paths_tpls = list(
                filter(
                    lambda p: p.startswith("output/"),
                    flatten(output_files),
                )
            )
            #: Generate all concrete output paths.
            for library_name in self.parent.ngs_library_to_extra_infos.keys():
                for sub_step in self.parent.sub_steps.values():
                    if (
                        isinstance(sub_step, ReadMappingStepPart)
                        and not sub_step.skip_result_files_for_library(library_name)
                        and not self.skip_result_files_for_library(library_name)
                    ):
                        for path_tpl in result_paths_tpls:
                            yield path_tpl.format(mapper=sub_step.name, library_name=library_name)
                if action == "collect":
                    break  # only once


class ReadMappingStepPart(MappingGetResultFilesMixin, BaseStepPart):
    """Base class for read mapping step parts"""

    #: Class available actions. Read mapping step parts only support action "run".
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{mapper}.{{library_name}}/out/{mapper}.{{library_name}}{ext}"
        self.extensions = EXT_VALUES
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config["path_link_in"],
        )

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            result = {
                "input": {
                    "reads_left": list(
                        sorted(self._collect_reads(wildcards, wildcards.library_name, ""))
                    )
                },
                "sample_name": wildcards.library_name,
                "platform": "ILLUMINA",
            }
            reads_right = list(
                sorted(self._collect_reads(wildcards, wildcards.library_name, "right-"))
            )
            if reads_right:
                result["input"]["reads_right"] = reads_right
            return result

        assert action == "run", "Unsupported actions"
        return args_function

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            return "work/input_links/{library_name}/.done".format(**wildcards)

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return concrete output files that all read mapping sub steps will write"""
        assert action in self.actions
        # Obtain and yield the paths in the ``work/`` directory
        paths_work = self._get_output_files_run_work()
        yield from paths_work.items()
        # Return list of paths to the links that will be created in ``output/``
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(paths_work.values(), self.get_log_file(action).values())
        ]

    @dictify
    def _get_output_files_run_work(self):
        """Return dict of output files to write to the ``work/`` directory."""
        for ext in self.extensions:
            yield ext[1:].replace(".", "_"), self.base_path_out.format(mapper=self.name, ext=ext)
        for ext in (".bamstats.txt", ".flagstats.txt", ".idxstats.txt"):
            path = (
                "work/{mapper}.{{library_name}}/report/bam_qc/" "{mapper}.{{library_name}}.bam{ext}"
            ).format(mapper=self.name, ext=ext)
            yield "report_" + ".".join(ext.split(".")[1:3]).replace(".", "_"), path
        for ext in (
            ".bamstats.txt.md5",
            ".flagstats.txt.md5",
            ".idxstats.txt.md5",
        ):
            path = (
                "work/{mapper}.{{library_name}}/report/bam_qc/" "{mapper}.{{library_name}}.bam{ext}"
            ).format(mapper=self.name, ext=ext)
            yield "report_" + ".".join(ext.split(".")[1:3]).replace(".", "_") + "_md5", path

    @dictify
    def get_log_file(self, action):
        """Return dict of log files in the "log" directory."""
        _ = action
        mapper = self.__class__.name
        prefix = f"work/{mapper}.{{library_name}}/log/{mapper}.{{library_name}}.mapping"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config["path_link_in"]:
            folder_name = library_name
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        seen = []
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            path = os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)
            if path in seen:
                print("WARNING: ignoring path seen before %s" % path, file=sys.stderr)
            else:
                seen.append(path)
                yield path


class BwaStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using BWA"""

    #: Step name
    name = "bwa"

    #: Tool category
    tool_category = "dna"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        mem_mb = int(4.5 * 1024 * self.config["bwa"]["num_threads_align"])
        return ResourceUsage(
            threads=self.config["bwa"]["num_threads_align"],
            time="3-00:00:00",  # 3 days
            memory=f"{mem_mb}M",
        )

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute BWA are present in the
        configuration. It further checks that the provided index has all the expected file
        extensions. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if self.__class__.name not in self.config["tools"]["dna"]:
            return  # BWA not run, don't check configuration  # pragma: no cover

        # Check required configuration settings present
        self.parent.ensure_w_config(
            config_keys=("step_config", "ngs_mapping", "bwa", "path_index"),
            msg="Path to BWA index is required",
        )

        # Check that the path to the BWA index is valid.
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            expected_path = self.config["bwa"]["path_index"] + ext
            if not os.path.exists(expected_path):  # pragma: no cover
                tpl = "Expected BWA input path {expected_path} does not exist!".format(
                    expected_path=expected_path
                )
                raise InvalidConfiguration(tpl)


class BwaMem2StepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using BWA-MEM 2"""

    name = "bwa_mem2"
    tool_category = "dna"

    def get_resource_usage(self, action: str) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        mem_mb = int(4.5 * 1024 * self.config["bwa_mem2"]["num_threads_align"])
        return ResourceUsage(
            threads=self.config["bwa_mem2"]["num_threads_align"],
            time="3-00:00:00",  # 3 days
            memory=f"{mem_mb}M",
        )

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute BWA-MEM2 are present in the
        configuration. It further checks that the provided index has all the expected file
        extensions. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if self.__class__.name not in self.config["tools"]["dna"]:
            return  # BWA-MEM2 not run, don't check configuration  # pragma: no cover

        # Check required configuration settings present
        self.parent.ensure_w_config(
            config_keys=("step_config", "ngs_mapping", "bwa_mem2", "path_index"),
            msg="Path to BWA-MEM2 index is required",
        )

        # Check that the path to the BWA-MEM2 index is valid.
        for ext in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
            expected_path = self.config["bwa_mem2"]["path_index"] + ext
            if not os.path.exists(expected_path):  # pragma: no cover
                raise InvalidConfiguration(
                    f"Expected BWA-MEM2 input path {expected_path} does not exist!"
                )


class StarStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using STAR"""

    #: Step name
    name = "star"

    #: Tool category
    tool_category = "rna"

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute BWA are present in the
        configuration. It further checks that the provided index has all the expected file
        extensions. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if self.__class__.name not in self.config["tools"]["rna"]:
            return  # STAR not run, don't check configuration  # pragma: no cover

        # Check required configuration settings present
        self.parent.ensure_w_config(
            config_keys=("step_config", "ngs_mapping", "star", "path_index"),
            msg="Path to STAR index is required",
        )
        self.parent.ensure_w_config(
            config_keys=("static_data_config", "reference"),
            msg="No reference genome FASTA file given",
        )

        # Check validity of the STAR index
        full_path = self.config["star"]["path_index"]
        # a lot of files should be in this dir, justtest these
        for indfile in ("Genome", "SA", "SAindex"):
            expected_path = os.path.join(full_path, indfile)
            if not os.path.exists(expected_path):  # pragma: no cover
                tpl = "Expected STAR index file {expected_path} does not exist!".format(
                    expected_path=expected_path
                )
                raise InvalidConfiguration(tpl)

    @dictify
    def _get_output_files_run_work(self):
        """Override base class' function to make Snakemake aware of extra files for STAR."""
        output_files = super()._get_output_files_run_work()
        if self.config[self.name]["path_features"]:
            output_files["gene_counts"] = self.base_path_out.format(
                mapper=self.name, ext=".GeneCounts.tab"
            )
            output_files["gene_counts_md5"] = output_files["gene_counts"] + ".md5"
        if self.config[self.name]["transcriptome"]:
            output_files["transcriptome"] = self.base_path_out.format(
                mapper=self.name, ext=".toTranscriptome.bam"
            )
            output_files["transcriptome_md5"] = output_files["transcriptome"] + ".md5"
        return output_files

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        mem_gb = int(3.5 * self.config["star"]["num_threads_align"])
        return ResourceUsage(
            threads=self.config["star"]["num_threads_align"],
            time="2-00:00:00",  # 2 days
            memory=f"{mem_gb}G",
        )


class Minimap2StepPart(ReadMappingStepPart):
    """Support for performing long-read alignment using minimap2"""

    #: Step name
    name = "minimap2"

    #: Tool category
    tool_category = "dna_long"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        mem_gb = int(3.5 * self.config["minimap2"]["mapping_threads"])
        return ResourceUsage(
            threads=self.config["minimap2"]["mapping_threads"],
            time="2-00:00:00",  # 2 days
            memory=f"{mem_gb}G",
        )

    def get_params(self, action):
        assert action == "run", "Parameters only available for action 'run'."
        return getattr(self, "_get_params_run")

    def _get_params_run(self, wildcards):
        return {"extra_infos": self.parent.ngs_library_to_extra_infos[wildcards.library_name]}


class ExternalStepPart(ReadMappingStepPart):
    """Support for linking in external BAM files"""

    #: Step name
    name = "external"

    #: Use wildcard for tool library
    tool_category = "__any__"

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute BWA are present in the
        configuration. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if "external" not in self.config["tools"]["dna"]:
            return  # External not run, don't check configuration  # pragma: no cover

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            return {
                "input": self._collect_bams(wildcards, wildcards.library_name),
                "sample_name": wildcards.library_name,
                "platform": "EXTERNAL",
            }

        assert action == "run", "Unsupported actions"
        return args_function

    @listify
    def _collect_bams(self, wildcards, library_name):
        """Yield the path to bam files"""
        _ = library_name
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        for _, path_infix, filename in self.path_gen.run(folder_name, ("bam",)):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        return ResourceUsage(
            threads=1,
            time="00:10:00",  # 10 minutes
            memory="1G",
        )


class TargetCoverageReportStepPart(ReportGetResultFilesMixin, BaseStepPart):
    """Build target coverage report"""

    #: Step name
    name = "target_coverage_report"

    #: Run for "dna" and "dna_long" only
    tool_categories = ("dna", "dna_long")

    #: Class available actions
    actions = ("run", "collect")

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return required input files"""
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards):
        mapper_lib = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", f"work/{mapper_lib}/out/{mapper_lib}.bam"
        yield "bai", f"work/{mapper_lib}/out/{mapper_lib}.bam.bai"

    @listify
    def _get_input_files_collect(self, wildcards):
        _ = wildcards
        for sheet in self.parent.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                extraction_type = ngs_library.test_sample.extra_infos.get("extractionType", "DNA")
                if ngs_library.extra_infos["seqPlatform"] in ("ONP", "PacBio"):
                    suffix = "_long"
                else:
                    suffix = ""
                for mapper in self.config["tools"][extraction_type.lower() + suffix]:
                    if (
                        self.parent.default_kit_configured
                        or ngs_library.name in self.parent.ngs_library_to_kit
                    ):
                        yield self._get_output_files_run_work()["txt"].format(
                            mapper=mapper, library_name=ngs_library.name
                        )

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        self._validate_action(action)
        # Obtain and yield the paths in the ``work/`` directory
        paths_work = getattr(self, f"_get_output_files_{action}_work")()
        yield from paths_work.items()
        # Return list of paths to the links that will be created in ``output/``
        paths_all = list(paths_work.values())
        paths_logs = self.get_log_file(action)
        if isinstance(paths_logs, dict):
            paths_all += list(paths_logs.values())
        elif isinstance(paths_logs, list):
            paths_all += paths_logs
        elif isinstance(paths_logs, str):
            paths_all.append(paths_logs)
        yield "output_links", [re.sub(r"^work/", "output/", work_path) for work_path in paths_all]

    @dictify
    def _get_output_files_run_work(self):
        yield "txt", "work/{mapper}.{library_name}/report/cov_qc/{mapper}.{library_name}.txt"
        yield "txt_md5", "work/{mapper}.{library_name}/report/cov_qc/{mapper}.{library_name}.txt.md5"

    @dictify
    def _get_output_files_collect_work(self):
        yield "txt", "work/target_cov_report/out/target_cov_report.txt"
        yield "txt_md5", "work/target_cov_report/out/target_cov_report.txt.md5"

    @dictify
    def get_log_file(self, action):
        self._validate_action(action)
        if action == "run":
            prefix = "work/{mapper}.{library_name}/log/{mapper}.{library_name}.target_cov_report"
            key_ext = (
                ("log", ".log"),
                ("conda_info", ".conda_info.txt"),
                ("conda_list", ".conda_list.txt"),
                ("wrapper", ".wrapper.py"),
                ("env_yaml", ".environment.yaml"),
            )
            for key, ext in key_ext:
                yield key, prefix + ext
                yield key + "_md5", prefix + ext + ".md5"
        else:
            yield "log", "work/target_cov_report/log/snakemake.target_coverage.log"

    def get_params(self, action):
        assert action == "run", "Parameters only available for action 'run'."
        return getattr(self, "_get_params_run")

    def _get_params_run(self, wildcards):
        # Find bed file associated with library kit
        library_name = wildcards.library_name
        path_targets_bed = ""
        kit_name = self.parent.ngs_library_to_kit.get(library_name, "__default__")
        for item in self.config["target_coverage_report"]["path_target_interval_list_mapping"]:
            if item["name"] == kit_name:
                path_targets_bed = item["path"]
                break

        return {
            "path_targets_bed": path_targets_bed,
            "max_coverage": self.config["target_coverage_report"]["max_coverage"],
            "min_cov_warning": self.config["target_coverage_report"]["min_cov_warning"],
            "min_cov_ok": self.config["target_coverage_report"]["min_cov_ok"],
            "detailed_reporting": self.config["target_coverage_report"]["detailed_reporting"],
        }

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="04:00:00",  # 4 hours
            memory="20G",
        )


class BamCollectDocStepPart(ReportGetResultFilesMixin, BaseStepPart):
    """Generate depth of coverage files."""

    #: Step name
    name = "bam_collect_doc"

    #: Run for "dna" and "dna_long" only
    tool_categories = ("dna", "dna_long")

    #: Class available actions
    actions = ("run",)

    def skip_result_files_for_library(self, library_name: str) -> bool:
        return not self.config["bam_collect_doc"]["enabled"]

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return required input files"""
        self._check_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def _check_action(self, action):
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)

    @dictify
    def _get_input_files_run(self):
        yield "bam", "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
        yield "bai", "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam.bai"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        self._check_action(action)
        # Obtain and yield the paths in the ``work/`` directory
        paths_work = self._get_output_files_run_work()
        yield from paths_work.items()
        # Return list of paths to the links that will be created in ``output/``
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(paths_work.values(), self.get_log_file(action).values())
        ]

    @dictify
    def _get_output_files_run_work(self):
        yield "vcf", "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz"
        yield "vcf_md5", "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.md5"
        yield "vcf_tbi", "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.tbi"
        yield "vcf_tbi_md5", "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.vcf.gz.tbi.md5"
        yield "bw", "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.bw"
        yield "bw_md5", "work/{mapper}.{library_name}/report/cov/{mapper}.{library_name}.cov.bw.md5"

    @dictify
    def get_log_file(self, action):
        _ = action
        prefix = "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_collect_doc"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._check_action(action)
        return ResourceUsage(
            threads=1,
            time="04:00:00",
            memory="2G",
        )


class NgsChewStepPart(ReportGetResultFilesMixin, BaseStepPart):
    """Analyze BAM File with ``ngs-chew``, e.g., ``fingerprint``"""

    #: Step name
    name = "ngs_chew"

    #: Class available actions
    actions = ("fingerprint",)

    def __init__(self, parent):
        super().__init__(parent)

    def skip_result_files_for_library(self, library_name: str) -> bool:
        return not self.config["ngs_chew_fingerprint"]["enabled"]

    def get_input_files(self, action):
        """Return required input files"""
        self._check_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def _check_action(self, action):
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)

    @dictify
    def _get_input_files_fingerprint(self):
        yield "bam", "work/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        self._check_action(action)
        # Obtain and yield the paths in the ``work/`` directory
        paths_work = self._get_output_files_fingerprint_work()
        yield from paths_work.items()
        # Return list of paths to the links that will be created in ``output/``
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(paths_work.values(), self.get_log_file(action).values())
        ]

    @dictify
    def _get_output_files_fingerprint_work(self):
        yield "npz", "work/{mapper}.{library_name}/report/fingerprint/{mapper}.{library_name}.npz"
        yield "npz_md5", "work/{mapper}.{library_name}/report/fingerprint/{mapper}.{library_name}.npz.md5"

    def get_log_file(self, action):
        self._check_action(action)
        return getattr(self, "_get_log_files_{action}".format(action=action))()

    @dictify
    def _get_log_files_fingerprint(self):
        prefix = "work/{mapper}.{library_name}/log/{mapper}.{library_name}.ngs_chew_fingerprint"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._check_action(action)
        return ResourceUsage(
            threads=1,
            time="04:00:00",
            memory="2G",
        )


class NgsMappingWorkflow(BaseStep):
    """Perform NGS Mapping"""

    #: Step name
    name = "ngs_mapping"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(workflow, config, config_lookup_paths, config_paths, workdir)
        self.register_sub_step_classes(
            (
                BwaStepPart,
                BwaMem2StepPart,
                ExternalStepPart,
                LinkInStep,
                Minimap2StepPart,
                StarStepPart,
                TargetCoverageReportStepPart,
                BamCollectDocStepPart,
                NgsChewStepPart,
            )
        )
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit, self.default_kit_configured = self._build_ngs_library_to_kit()
        # Create shortcut from library to all extra infos.
        self.ngs_library_to_extra_infos = self._build_ngs_library_to_extra_infos()
        # Validate project
        self.validate_project(config_dict=self.config, sample_sheets_list=self.shortcut_sheets)

    def _build_ngs_library_to_extra_infos(self):
        result = {}
        for donor in self._all_donors():
            for bio_sample in donor.bio_samples.values():
                for test_sample in bio_sample.test_samples.values():
                    for library in test_sample.ngs_libraries.values():
                        result.setdefault(library.name, {}).update(test_sample.extra_infos)
                        result.setdefault(library.name, {}).update(library.extra_infos)
        return result

    def _build_ngs_library_to_kit(self):
        cov_config = DictQuery(self.w_config).get("step_config/ngs_mapping/target_coverage_report")
        # Build mapping.
        default_kit_configured = False
        regexes = {}
        for item in cov_config["path_target_interval_list_mapping"]:
            if item["name"] == "__default__":
                default_kit_configured = True
            else:
                regexes[item["pattern"]] = item["name"]
        result = {}
        for donor in self._all_donors():
            for bio_sample in donor.bio_samples.values():
                for test_sample in bio_sample.test_samples.values():
                    for library in test_sample.ngs_libraries.values():
                        if library.extra_infos.get("libraryKit"):
                            library_kit = library.extra_infos.get("libraryKit")
                            for pattern, name in regexes.items():
                                if re.match(pattern, library_kit):
                                    result[library.name] = name
        return result, default_kit_configured

    @listify
    def _all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = filter(is_not_background, sheets)
        for sheet in sheets:
            for entity in sheet.bio_entities.values():
                yield entity

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample sheets.
        """
        for sub_step in self.sub_steps.values():
            if sub_step.name not in (LinkInStep.name,):
                yield from sub_step.get_result_files()

    def validate_project(self, config_dict, sample_sheets_list):
        """Validates project.

        Method compares sample information included in the sample sheet and the configuration. If
        sheet contains 'DNA' samples, a DNA mapper should be defined. Similarly, if it contains
        'RNA', a RNA mapper should be defined.

        :param config_dict: Dictionary with configurations as found in the project's yaml file.
        :type config_dict: dict

        :param sample_sheets_list: List with biomedical sample sheets.
        :type sample_sheets_list: list
        """
        # Initialise variables
        dna_bool_list = []
        rna_bool_list = []

        # Get tools dictionary
        tools_dict = config_dict["tools"]

        # Iterate over sheets
        for sheet in sample_sheets_list:
            dna_present, rna_present = self.extraction_type_check(sample_sheet=sheet)
            # Append to respective lists
            dna_bool_list.append(dna_present)
            rna_bool_list.append(rna_present)

        # Evaluate type of project
        dna_analysis = any(dna_bool_list)
        rna_analysis = any(rna_bool_list)

        # Validate DNA project
        dna_tool_list = tools_dict.get("dna", [])
        if dna_analysis and not dna_tool_list:
            raise InvalidConfiguration(
                "Sample sheet contains DNA but configuration has no DNA "
                "mapper defined in tool list."
            )
        # Validate RNA project
        rna_tool_list = tools_dict.get("rna", [])
        if rna_analysis and not rna_tool_list:
            raise InvalidConfiguration(
                "Sample sheet contains RNA but configuration has no RNA "
                "mapper defined in tool list."
            )

    @staticmethod
    def extraction_type_check(sample_sheet):
        """Retrieve extraction type from biomedsheet.

        Method crawls through all bio entities in the biomedsheet and checks if there are DNA
        and/or RNA extraction types. In both cases, the test will be consider True if at least one
        test sample contains the extraction type (i.e., DNA or RNA).

        :param sample_sheet: Sample sheet.
        :type sample_sheet: biomedsheets.models.Sheet

        :return: Returns tuple with boolean for DNA, RNA extraction types: (DNA extraction type
        present, RNA extraction type present).
        """
        # Initialise variables
        contains_rna_extraction = False
        contains_dna_extraction = False

        # Crawl over bio entities until test_sample
        for _, entity in sample_sheet.bio_entities.items():
            for _, bio_sample in entity.bio_samples.items():
                for _, test_sample in bio_sample.test_samples.items():
                    extraction_type = test_sample.extra_infos.get("extractionType")
                    if extraction_type.lower() == "dna":
                        contains_dna_extraction = True
                    elif extraction_type.lower() == "rna":
                        contains_rna_extraction = True

        # Return
        return contains_dna_extraction, contains_rna_extraction
