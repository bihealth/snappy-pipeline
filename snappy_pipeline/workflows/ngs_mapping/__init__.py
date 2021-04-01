# -*- coding: utf-8 -*-
"""Implementation of the ``ngs_mapping`` step

The ngs_mapping step allows for the alignment of NGS data with standard read mappers, such as BWA
for DNA data and STAR for RNA data.  Also, it provides functionality to compute post-alignment
statistics, such as the coverage of target (e.g., exome or panel) regions and genome-wide per base
pair coverage.

There is a distinction made between "normal" DNA reads (short reads from Illumina) and "long"
DNA reads, such as PacBio/Oxford Nanopore.  Again, the NGS mapping step will perform alignment
of all NGS libraries.

The precise actions that are performed in the alignment are defined by the wrappers (e.g., the
``bwa`` or ``star``) wrappers.  Generally, this includes converting into BAM format, sorting
by coordinate, an indexing using a BAI file.  For short reads, this can include marking of
duplicates using Samblaster and depends on the actual configuration (see below for the default
configuration).

==========
Stability
==========

The ngs_mapping pipeline step is considered stable for short Illumina reads (DNA and RNA).

The long read mapping steps may not be stable as they are currently considered experimental.

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
    - ``"ngmlr"``
    - ``"external"``

=======
Reports
=======

Currently, the following reports are generated based on the BAM and BAI file output by this step.

General Alignment Statistics (.txt)
  The tools ``samtools bamstats``, ``samtools flagstats`` and ``samtools idxstats`` are always
  called by default, and are linked out into the ``output/{mapper}.{library_name}/report/bam_qc``
  directory. The file names for these reports (and their MD5s) use the following naming convention:

  - ``{mapper}.{library_name}.bamstats.html``
  - ``{mapper}.{library_name}.bamstats.txt``
  - ``{mapper}.{library_name}.flagstats.txt``
  - ``{mapper}.{library_name}.idxstats.txt``
  - ``{mapper}.{library_name}.bamstats.html.md5``
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
    |           |-- bwa.P001-N1-DNA1-WES1.bam.bamstats.html
    |           |-- bwa.P001-N1-DNA1-WES1.bam.bamstats.html.md5
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

Bamstats Plots (.gp .png)
  Output of ``plot-bamstats`` is always generated by default.  If you specify the ``ref_gc_stats``
  setting then GC-dependent metrics will be computed as well. These plots are not linked out into
  the ``output/`` directory and are saved in the
  ``work/{mapper}.{library_name}/report/bam_qc/{mapper}.{library_name}.bam.bamstats.d``
  directory. The plots have default file names that do not refer to the input bam or upstream
  processing.

  For example, it will look as follows for the example bam files shown above:

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

import os
import re
import sys
import textwrap

from biomedsheets.shortcuts import GenericSampleSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.workflows.abstract import (
    STDERR_TO_LOG_FILE,
    BaseStepPart,
    BaseStep,
    LinkOutStepPart,
    LinkInStep,
    LinkInPathGenerator,
    get_ngs_library_folder_name,
)
from snappy_pipeline.utils import dictify, listify, DictQuery

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# TODO: Need something smarter still for @RG

#: Extensions of files to create as main payload
EXT_VALUES = (".bam", ".bam.bai", ".bam.md5", ".bam.bai.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("bam", "bai", "bam_md5", "bai_md5")

#: Value for GATK BAM postprocessing
POSTPROC_GATK_POST_BAM = "gatk_post_bam"

#: Available read mappers for (short/Illumina) DNA-seq data
READ_MAPPERS_DNA = ("bwa",)

#: Available read mappers for (short/Illumina) RNA-seq data
READ_MAPPERS_RNA = ("star",)

#: Available read mappers for (long/PacBio/Nanopoare) DNA-seq data
READ_MAPPERS_DNA_LONG = ("minialign", "ngmlr", "ngmlr_chained")

#: Default configuration
DEFAULT_CONFIG = r"""
step_config:
  ngs_mapping:
    # Aligners to use for the different NGS library types
    tools:
      dna: []      # Required if DNA analysis; otherwise, leave empty. Example: 'bwa'.
      rna: []      # Required if RNA analysis; otherwise, leave empty. Example: 'star'.
      dna_long: [] # Required if long-read mapper used; otherwise, leave empty. Example: 'ngmlr'.
    # Whether or not to compute coverage BED file
    compute_coverage_bed: false
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
    # Enable Picard HS metrics by setting both paths
    picard_hs_metrics:
      path_targets_interval_list: null
      path_baits_interval_list: null
    # Configuration for BWA
    bwa:
      path_index: REQUIRED # Required if listed in ngs_mapping.tools.dna; otherwise, can be removed.
      ref_gc_stats: null  # optional
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
      ref_gc_stats: null
      num_threads_align: 16
      num_threads_trimming: 8
      num_threads_bam_view: 4
      num_threads_bam_sort: 4
      memory_bam_sort: 4G
      trim_adapters: false
      mask_duplicates: false
      raw_star_options: ''
      align_intron_max: 1000000
      align_intron_min: 20
      align_mates_gap_max: 1000000
      align_sjdb_overhang_min: 1
      align_sj_overhang_min: 8
      genome_load: NoSharedMemory
      out_filter_intron_motifs: None # or for cufflinks: RemoveNoncanonical
      out_filter_mismatch_n_max: 999
      out_filter_mismatch_n_over_l_max: 0.04
      out_filter_multimap_n_max: 20
      out_filter_type: BySJout
      out_sam_strand_field: None # or for cufflinks: intronMotif
      include_unmapped: true
      quant_mode: ''
    # Configuration for Minialign
    minialign:
      # `path_index`: Required if listed in ngs_mapping.tools.dna_long; otherwise, can be removed.
      path_index: REQUIRED
      ref_gc_stats: null    # Optional
      mapping_threads: 16
      num_threads_bam_view: 4
    # Configuration for NGMLR
    ngmlr:
      # `path_index`: Required if listed in ngs_mapping.tools.dna_long; otherwise, can be removed.
      path_index: REQUIRED
      ref_gc_stats: null    # Optional
    # Select postprocessing method, only for DNA alignment
    postprocessing: null # optional, {'gatk_post_bam'}
    # Configuration for GATK BAM postprocessing
    gatk_post_bam:
      do_realignment: true
      do_recalibration: true
      realigned_infix: realigned
      recalibrated_infix: recalibrated
      path_known_sites_vcf: REQUIRED # is list
      time_multiplicator_wgs: 4
"""


class ReadMappingStepPart(BaseStepPart):
    """Base class for read mapping step parts"""

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{mapper}.{{library_name}}/out/{mapper}.{{library_name}}{ext}"
        self.extensions = EXT_VALUES
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir, self.parent.data_set_infos, self.parent.config_lookup_paths
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
        """Return output files that all read mapping sub steps must return
        (BAM + BAI file)
        """
        assert action == "run"
        for ext in self.extensions:
            yield ext[1:].replace(".", "_"), self.base_path_out.format(mapper=self.name, ext=ext)
        for ext in (".bamstats.html", ".bamstats.txt", ".flagstats.txt", ".idxstats.txt"):
            path = (
                "work/{mapper}.{{library_name}}/report/bam_qc/" "{mapper}.{{library_name}}.bam{ext}"
            ).format(mapper=self.name, ext=ext)
            yield "report_" + ".".join(ext.split(".")[1:3]).replace(".", "_"), path
        for ext in (
            ".bamstats.html.md5",
            ".bamstats.txt.md5",
            ".flagstats.txt.md5",
            ".idxstats.txt.md5",
        ):
            path = (
                "work/{mapper}.{{library_name}}/report/bam_qc/" "{mapper}.{{library_name}}.bam{ext}"
            ).format(mapper=self.name, ext=ext)
            yield "report_" + ".".join(ext.split(".")[1:3]).replace(".", "_") + "_md5", path

    @dictify
    def _get_log_file(self, _action):
        """Return dict of log files."""
        prefix = "work/{mapper}.{{library_name}}/log/{mapper}.{{library_name}}".format(
            mapper=self.__class__.name
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def _collect_reads(self, wildcards, _library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
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

    name = "bwa"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_bwa_run"] = {
            "mem": int(4.5 * 1024 * self.config["bwa"]["num_threads_align"]),
            "time": "140:00",
            "ntasks": self.config["bwa"]["num_threads_align"],
        }

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


class StarStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using STAR"""

    name = "star"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_star_run"] = {
            "mem": int(3.7 * 1024 * 16),
            "time": "40:00",
            "ntasks": 16,
        }

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


class Minimap2StepPart(ReadMappingStepPart):
    """Support for performing long-read alignment using minimap2"""

    name = "minimap2"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_minimap2_run"] = {
            "mem": int(3.7 * 1024 * 16),
            "time": "96:00",
            "ntasks": 16,
        }


class NgmlrStepPart(ReadMappingStepPart):
    """Support for performing PacBio alignment using NGMLR without chaining"""

    name = "ngmlr"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_ngmlr_run"] = {
            "mem": int(3.7 * 1024 * 16),
            "time": "96:00",
            "ntasks": 16,
        }

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute BWA are present in the
        configuration. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if not (set(self.config["tools"]["dna_long"]) & {"ngmlr", "ngmlr_chained"}):
            return  # NGLMR not run, don't check configuration  # pragma: no cover

        # Check required configuration settings present
        self.parent.ensure_w_config(
            config_keys=("step_config", "ngs_mapping", "ngmlr", "path_index"),
            msg="Path to NGMLR index is required",
        )


class ExternalStepPart(ReadMappingStepPart):
    """Support for linking in external BAM files"""

    name = "external"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_external_run"] = {"mem": 1024, "time": "04:00", "ntasks": 1}

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
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        for _, path_infix, filename in self.path_gen.run(folder_name, ("bam",)):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)


class GatkPostBamStepPart(BaseStepPart):
    """Support for supporting BAM postprocessing with GATK

    This uses the snappy-gatk_post_bam wrapper that performs read realignment and base
    recalibration. Note that in particular the base recalibration step takes quite a long time
    for large files.
    """

    name = "gatk_post_bam"

    def __init__(self, parent):
        super().__init__(parent)
        self.path_tpl = (
            "work/{{mapper}}.{{library_name}}/out/{{mapper}}.{{library_name}}{infix}{ext}"
        )

    def check_config(self):
        """Check parameters in configuration.

        Method checks that all parameters required to execute BWA are present in the
        configuration. If invalid configuration, it raises InvalidConfiguration exception.
        """
        # Check if tool is at all included in workflow
        if "gatk_post_bam" not in (self.config["postprocessing"] or []):  # pylint: disable=C0325
            return  # GATK BAM postprocessing not enabled, skip

        # Check required configuration settings present
        self.parent.ensure_w_config(
            config_keys=("step_config", "ngs_mapping", "gatk_post_bam", "paths_known_sites"),
            msg="Known sites list cannot be empty for GATK BAM postprocessing",
        )
        self.parent.ensure_w_config(
            config_keys=("static_data_config", "reference", "path"),
            msg="Path to reference FASTA required for GATK BAM postprocessing",
        )

    def get_input_files(self, action):
        """Return required input files"""
        assert action == "run", "Unsupported action"
        return self.path_tpl.format(infix="", ext=".bam")

    @dictify
    def get_output_files(self, action):
        """Return output files that are generated by snappy-gatk_post_bam"""
        assert action == "run", "Unsupported action"
        realigned_infix = "." + self.config["gatk_post_bam"]["realigned_infix"]
        recalibrated_infix = "." + self.config["gatk_post_bam"]["recalibrated_infix"]
        if self.config["gatk_post_bam"]["do_realignment"]:
            for ext_name, ext in zip(EXT_NAMES, EXT_VALUES):
                yield ext_name + "_realigned", self.path_tpl.format(infix=realigned_infix, ext=ext)
            recalibrated_infix = realigned_infix + recalibrated_infix
        if self.config["gatk_post_bam"]["do_recalibration"]:
            for ext_name, ext in zip(EXT_NAMES, EXT_VALUES):
                yield ext_name + "_recalibrated", self.path_tpl.format(
                    infix=recalibrated_infix, ext=ext
                )

    def get_log_file(self, action):
        return "work/{mapper}.{library_name}/log/snakemake.gatk_post_bam.log"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_gatk_post_bam_run"] = {
            "mem": int(3.7 * 1024 * 16),
            "time": "40:00",
            "ntasks": 16,
        }


class LinkOutBamStepPart(BaseStepPart):
    """Link out the read mapping results

    Depending on the configuration, the files are linked out after
    postprocessing
    """

    # TODO: do not attempt realignment of RNA-seq data

    name = "link_out_bam"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = (
            "work/{wildcards.mapper}.{wildcards.library_name}/out/"
            "{wildcards.mapper}.{wildcards.library_name}{{postproc}}{{ext}}"
        )
        self.base_pattern_out = (
            "output/{{mapper}}.{{library_name}}/out/{{mapper}}.{{library_name}}{ext}"
        )
        self.base_path_out = self.base_pattern_out.replace(",[^.]", "")
        self.extensions = EXT_VALUES

    def get_input_files(self, action):
        """Return required input files"""

        def input_function(wildcards):
            """Helper rapper function"""
            return expand(
                self.base_path_in.format(wildcards=wildcards),
                postproc=[self._get_postproc_token()],
                ext=self.extensions,
            )

        assert action == "run", "Unsupported action"
        return input_function

    def get_output_files(self, action):
        """Return output files that are generated by snappy-gatk_post_bam"""
        assert action == "run", "Unsupported action"
        return expand(self.base_pattern_out, ext=self.extensions)

    def get_shell_cmd(self, action, wildcards):
        """Return call for linking out postprocessed (or not) files"""
        assert action == "run", "Unsupported action"
        ins = expand(
            self.base_path_in.format(wildcards=wildcards),
            postproc=[self._get_postproc_token()],
            ext=self.extensions,
        )
        outs = [s.format(**wildcards) for s in expand(self.base_path_out, ext=self.extensions)]
        assert len(ins) == len(outs)
        return "\n".join(
            (
                "test -L {out} || ln -sr {in_} {out}".format(in_=in_, out=out)
                for in_, out in zip(ins, outs)
            )
        )

    def _get_postproc_token(self):
        """Return file name token for result files for the postprocessing

        The value depends on the configuration.  This way, the workflow
        controls whether to execute postprocessing or not (and which
        postprocessing if there was more than one option).
        """
        if self.config["postprocessing"] == "gatk_post_bam":
            do_realignment = self.config["gatk_post_bam"]["do_realignment"]
            do_recalibration = self.config["gatk_post_bam"]["do_recalibration"]
        else:
            do_realignment, do_recalibration = False, False
        realigned_infix = self.config["gatk_post_bam"]["realigned_infix"]
        recalibrated_infix = self.config["gatk_post_bam"]["recalibrated_infix"]
        return {
            (False, False): "",
            (False, True): "." + recalibrated_infix,
            (True, False): "." + realigned_infix,
            (True, True): "." + realigned_infix + "." + recalibrated_infix,
        }[(do_realignment, do_recalibration)]


class PicardHsMetricsStepPart(BaseStepPart):
    """Build target report from Picard HsMetrics"""

    name = "picard_hs_metrics"

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return required input files"""
        assert action == "run", "Unsupported action"
        return {
            "bam": "work/{mapper_lib}/out/{mapper_lib}.bam",
            "bai": "work/{mapper_lib}/out/{mapper_lib}.bam.bai",
        }

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        assert action == "run", "Unsupported action"
        yield "txt", "work/{mapper_lib}/report/picard_hs_metrics/{mapper_lib}.txt"
        yield "txt_md5", "work/{mapper_lib}/report/picard_hs_metrics/{mapper_lib}.txt.md5"

    def get_log_file(self, action):
        return "work/{mapper_lib}/log/snakemake.picard_hs_metrics.log"

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_picard_hs_metrics_run"] = {
            "mem": int(19 * 1024 * 2),
            "time": "04:00",
            "ntasks": 2,
        }


class TargetCoverageReportStepPart(BaseStepPart):
    """Build target coverage report"""

    name = "target_coverage_report"
    actions = ("run", "collect")

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return required input files"""
        assert action in self.actions, "Invalid action"
        return getattr(self, "_get_input_files_{action}".format(action=action))

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "bam", "work/{mapper_lib}/out/{mapper_lib}.bam".format(**wildcards)
        yield "bai", "work/{mapper_lib}/out/{mapper_lib}.bam.bai".format(**wildcards)

    @listify
    def _get_input_files_collect(self, wildcards):
        for mapper in self.config["tools"]["dna"]:
            for sheet in self.parent.shortcut_sheets:
                for library in sheet.all_ngs_libraries:
                    if library.name in self.parent.ngs_library_to_kit:
                        kv = {"mapper_lib": "{}.{}".format(mapper, library.name)}
                        yield self._get_output_files_run()["txt"].format(**kv)

    def get_output_files(self, action):
        """Return output files"""
        assert action in self.actions, "Invalid action"
        return getattr(self, "_get_output_files_{action}".format(action=action))()

    @dictify
    def _get_output_files_run(self):
        yield "txt", "work/{mapper_lib}/report/cov_qc/{mapper_lib}.txt"
        yield "txt_md5", "work/{mapper_lib}/report/cov_qc/{mapper_lib}.txt.md5"

    @dictify
    def _get_output_files_collect(self):
        yield "txt", "work/target_cov_report/out/target_cov_report.txt"
        yield "txt_md5", "work/target_cov_report/out/target_cov_report.txt.md5"

    def get_log_file(self, action):
        if action == "run":
            return "work/{mapper_lib}/log/snakemake.target_coverage.log"
        else:
            return "work/target_cov_report/log/snakemake.target_coverage.log"

    def update_cluster_config(self, cluster_config):
        for action in self.actions:
            cluster_config["ngs_mapping_target_coverage_report_{}".format(action)] = {
                "mem": int(19 * 1024 * 2),
                "time": "10:00",
                "ntasks": 16,
            }


class GenomeCoverageReportStepPart(BaseStepPart):
    """Build genome-wide per-base coverage report"""

    name = "genome_coverage_report"

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action):
        """Return required input files"""
        assert action == "run", "Unsupported action"
        return {
            "bam": "work/{mapper_lib}/out/{mapper_lib}.bam",
            "bai": "work/{mapper_lib}/out/{mapper_lib}.bam.bai",
        }

    def get_output_files(self, action):
        """Return output files"""
        assert action == "run", "Unsupported action"
        return {
            "bed": "work/{mapper_lib}/report/coverage/{mapper_lib}.bed.gz",
            "tbi": "work/{mapper_lib}/report/coverage/{mapper_lib}.bed.gz.tbi",
        }

    def get_log_file(self, action):
        return "work/{mapper_lib}/log/snakemake.genome_coverage.log"

    # TODO(holtgrewe): can this be removed?
    def get_shell_cmd(self, action, wildcards):
        """Return bash script to execute"""
        assert action == "run", "Unsupported action"
        return (
            STDERR_TO_LOG_FILE
            + textwrap.dedent(
                r"""
            module load SAMtools/1.2-foss-2015a
            module load HTSlib/1.2.1-foss-2015a

            samtools depth {input.bam} \
            | awk -F $'\t' 'BEGIN {{ OFS=FS; }} {{ print $1, $2-1, $2, $3; }}' \
            | awk -F $'\t' '
                    BEGIN {{ OFS=FS; pr = -1; pb = -1; pe = -1; pc = -1; }}
                    {{
                        if (pr == $1 && $2 == pe && $4 == pc) {{
                            pe = $3;
                        }} else {{
                            if (pr != -1)
                                print pr, pb, pe, pc;
                            pr = $1;
                            pb = $2;
                            pe = $3;
                            pc = $4;
                        }}
                    }}
                    END {{ if (pr != -1) {{
                        print pr, pb, pe, pc;
                    }} }}' \
            | bgzip -c \
            > {output.bed}

            tabix {output.bed}
            """
            ).lstrip()
        )

    def update_cluster_config(self, cluster_config):
        cluster_config["ngs_mapping_genome_coverage_report_run"] = {
            "mem": int(3.7 * 1024),
            "time": "04:00",
            "ntasks": 1,
        }


class NgsMappingWorkflow(BaseStep):
    """Perform NGS Mapping"""

    name = "ngs_mapping"
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
        )
        self.register_sub_step_classes(
            (
                BwaStepPart,
                ExternalStepPart,
                GatkPostBamStepPart,
                GenomeCoverageReportStepPart,
                LinkInStep,
                LinkOutBamStepPart,
                LinkOutStepPart,
                Minimap2StepPart,
                NgmlrStepPart,
                PicardHsMetricsStepPart,
                StarStepPart,
                TargetCoverageReportStepPart,
            )
        )
        self.sub_steps["link_out"].disable_patterns = expand("**/*{ext}", ext=EXT_VALUES)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        cov_config = DictQuery(self.w_config).get("step_config/ngs_mapping/target_coverage_report")
        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in cov_config["path_target_interval_list_mapping"]
            if item["name"] != "__default__"
        }
        result = {}
        for donor in self._all_donors():
            for bio_sample in donor.bio_samples.values():
                for test_sample in bio_sample.test_samples.values():
                    for library in test_sample.ngs_libraries.values():
                        if library.extra_infos.get("libraryKit"):
                            library_kit = library.extra_infos.get("libraryKit")
                            for pattern, name in regexes.items():
                                if re.match(pattern, library_kit):
                                    yield library.name, name
        return result

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

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        token = "{mapper}.{ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", token, "out", token + "{ext}"), ext=EXT_VALUES
        )
        yield from self._yield_result_files(
            os.path.join("output", token, "log", "{mapper}.{ngs_library.name}.{ext}"),
            ext=(
                "log",
                "conda_info.txt",
                "conda_list.txt",
                "log.md5",
                "conda_info.txt.md5",
                "conda_list.txt.md5",
            ),
        )
        yield from self._yield_result_files(
            os.path.join("output", token, "report", "bam_qc", token + ".bam.{report}.txt"),
            report=("bamstats", "flagstats", "idxstats"),
        )
        yield from self._yield_result_files(
            os.path.join("output", token, "report", "bam_qc", token + ".bam.{report}.txt.md5"),
            report=("bamstats", "flagstats", "idxstats"),
        )
        yield from self._yield_result_files(
            os.path.join("output", token, "report", "bam_qc", token + ".bam.bamstats.html")
        )
        yield from self._yield_result_files(
            os.path.join("output", token, "report", "bam_qc", token + ".bam.bamstats.html.md5")
        )

        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                if ngs_library.name in self.ngs_library_to_kit:
                    extraction_type = ngs_library.test_sample.extra_infos["extractionType"]
                    suffix = (
                        "_long"
                        if ngs_library.extra_infos["seqPlatform"] in ("PacBio", "ONP")
                        else ""
                    )
                    # Per-sample target coverage report.
                    yield from expand(
                        os.path.join("output", token, "report", "cov_qc", token + ".{ext}"),
                        mapper=self.config["tools"][extraction_type.lower() + suffix],
                        ngs_library=[ngs_library],
                        ext=["txt", "txt.md5"],
                    )
        yield "output/target_cov_report/out/target_cov_report.txt"
        yield "output/target_cov_report/out/target_cov_report.txt.md5"
        if (
            self.config["picard_hs_metrics"]["path_targets_interval_list"]
            and self.config["picard_hs_metrics"]["path_baits_interval_list"]
        ):
            yield from self._yield_result_files(
                os.path.join("output", token, "report", "picard_hs_metrics", token + ".txt")
            )
            yield from self._yield_result_files(
                os.path.join("output", token, "report", "picard_hs_metrics", token + ".txt.md5")
            )
        if self.config["compute_coverage_bed"]:
            yield from self._yield_result_files(
                os.path.join("output", token, "report", "coverage", token + "{ext}"),
                ext=(".bed.gz", ".bed.gz.tbi"),
            )
        else:
            print(
                "Genome-wide coverage BED generation disabled", file=sys.stderr
            )  # pragma: no cover

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                extraction_type = ngs_library.test_sample.extra_infos["extractionType"]
                if ngs_library.extra_infos["seqPlatform"] in ("ONP", "PacBio"):
                    suffix = "_long"
                else:
                    suffix = ""
                yield from expand(
                    tpl,
                    mapper=self.config["tools"][extraction_type.lower() + suffix],
                    ngs_library=[ngs_library],
                    **kwargs
                )

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
