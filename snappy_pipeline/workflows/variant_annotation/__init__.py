# -*- coding: utf-8 -*-
"""Implementation of the ``variant_annotation`` step

The ``variant_annotation`` step takes as the input the results of the ``variant_calling`` step
(called germline variants in vcf.gz format), a transcript database, and (optionally) a variety of
other snp database resources.  The results are variant files with additional annotations describing
their functional impact on transcripts and additional details from the provided snp resources.

==========
Stability
==========

This pipeline step currently only supports Jannovar annotations, and is considered stable.

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``variant_calling`` step.

===========
Step Output
===========

Variant annotation will be performed on all input vcfs, separately for each configured variant
annotator.  The name of the primary DNA NGS library of the index will be used as an identification
token in the output file.  For each read mapper, variant caller, variant annotator, and pedigree
the following files will be generated:


- ``{mapper}.{var_caller}.{var_annotator}.{lib_name}.vcf.gz``
- ``{mapper}.{var_caller}.{var_annotator}.{lib_name}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{var_annotator}.{lib_name}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{var_annotator}.{lib_name}.vcf.gz.tbi.md5``


For example, annotation of freebayes variant calls could look like this:

::

    output/
    +-- bwa.freebayes.jannovar_annotate_vcf.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.freebayes.jannovar_annotate_vcf.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.freebayes.jannovar_annotate_vcf.P001-N1-DNA1-WES1.vcf.gz.tbi
    |       |-- bwa.freebayes.jannovar_annotate_vcf.P001-N1-DNA1-WES1.vcf.gz.md5
    |       `-- bwa.freebayes.jannovar_annotate_vcf.P001-N1-DNA1-WES1.vcf.gz.tbi.md5
    [...]

====================
Global Configuration
====================

The ``static_data_config/<$resource>/path`` must be properly configured for all desired resources
and point to a file in .vcf.gz format.
Supported options include:

- ``static_data_config/dbsnp/path``
- ``static_data_config/exac/path``
- ``static_data_config/gnomad_exomes/path``
- ``static_data_config/gnomad_genomes/path``
- ``static_data_config/uk10k/path``
- ``static_data_config/clinvar/path``
- ``static_data_config/cosmic/path``

Also ``static_data_config/reference/path`` must be set appropriately to point to a reference
genome fasta file.

.. note::

    You must ensure agreement between all these global config ``static_data_config/*/path`` files,
    as well as the required parameter for ``step_config/variant_annotation/path_jannovar_ser``
    shown below. i.e. all these files correspond to the *same version of the reference genome* and
    use the *same chromosome names*.

=====================
Default Configuration
=====================

The default configuration is as follows. Note that the ``path_jannovar_ser`` parameter must
be set to point to the desired transcript annotations db as generated by ``jannovar download``.

.. include:: DEFAULT_CONFIG_variant_annotation.rst

============================
Available Variant Annotators
============================

The following variant annotator is currently available:

- ``"jannovar"`` : See the `software documentation <http://jannovar.readthedocs.io/en/master/>`_ \
for more details

=======
Reports
=======

Currently, no separate report files are generated.

The transcript-based annotations are added by jannovar to the INFO field of the vcf records as
follows:

    ANN
        Description: Functional annotations

        Value: Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|
        Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos / cDNA.length|CDS.pos / CDS.length|
        AA.pos / AA.length|Distance|ERRORS / WARNINGS / INFO

e.g. ::

    ANN=G|missense_variant|MODERATE|SELENON|57190|transcript|NM_020451.2|Coding|4/13|c.409A>G|
    p.(Thr137Ala)|464/18047|409/1773|137/591||;

If familial pedigree information is provided in the BioMedSheet, the following attributes will be
added to the INFO column of the vcf records as applicable:

    INHERITANCE
        Description: Compatible Mendelian inheritance modes.

        Possible values: AD, AR, XD, XR
    INHERITANCE_RECESSIVE_DETAIL
        Description: Extra annotation for recessive inheritance sub type.

        Possible values: AR_HOM_ALT, AR_COMP_HET, XR_HOM_ALT, XR_COMP_HET

e.g. ::

    INHERITANCE=AR;INHERITANCE_RECESSIVE_DETAIL=AR_HOM_ALT;

Annotations added based on snp databases will vary depending on the resource used. For COSMIC
annotations, the COSMIC database id will be added to the ID column and the following attributes
will be added to the INFO column, as appropriate.

    COSMIC_CNT
        Description: Number of samples in COSMIC having this mutation. (Integer)

    COSMIC_IDS
        Description: COSMIC IDs with matching alternative positions and alleles, for each
        alternative allele, separated by '|'.

    COSMIC_OVL_CNT
        Description: Number of samples in COSMIC having this mutation (requiring no genotype match,
        only position overlap).

    COSMIC_OVL_IDS
        Description: COSMIC IDs with overlapping alternative positions, not necessarily matching
        alleles, for each alternative allele, separated by '|'.

    COSMIC_OVL_SNP
        Description: Classified as SNP (polymorphism) in COSMIC. (Flag)

    COSMIC_SNP
        Description: Classified as SNP (polymorphism) in COSMIC. (Flag)

e.g. ::

    COSMIC_CNT=2;COSMIC_IDS=COSM1343868;COSMIC_OVL_CNT=2;COSMIC_OVL_IDS=COSM1343868;COSMIC_OVL_SNP;
    COSMIC_SNP;

Please see the header of the output vcfs for more details about the annotations provided when using
other databases.

==================
Parallel Execution
==================

Parallel execution of variant annotation is implemented similar as for ``variant_calling``.  See
the :ref:`variant_calling_parallel_execution` section of that pipeline step for more information.
"""

import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    WritePedigreeStepPart
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

# TODO: the number of restart times is high because tabix in HTSJDK/Jannovar is flaky...

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_annotation
step_config:
  variant_annotation:
    drmaa_snippet: ''         # value to pass in as additional DRMAA arguments
    window_length: 5000000    # split input into windows of this size, each triggers a job
    num_jobs: 100             # number of windows to process in parallel
    use_drmaa: true           # use DRMAA for parallel processing
    restart_times: 10         # number of times to re-launch jobs in case of failure
    max_jobs_per_second: 10   # throttling of job creation
    max_status_checks_per_second: 10   # throttling of status checks
    debug_trunc_tokens: 0     # truncation to first N tokens (0 for none)
    keep_tmpdir: never        # keep temporary directory, {always, never, onerror}
    job_mult_memory: 1        # memory multiplier
    job_mult_time: 1          # running time multiplier
    merge_mult_memory: 1      # memory multiplier for merging
    merge_mult_time: 1        # running time multiplier for merging
    ignore_chroms:            # patterns of chromosome names to ignore
    - NC_007605  # herpes virus
    - hs37d5     # GRCh37 decoy
    - chrEBV     # Eppstein-Barr Virus
    - '*_decoy'  # decoy contig
    - 'HLA-*'    # HLA genes
    use_advanced_ped_filters: false  # whether or not to use the advanced pedigree filters flag
    path_ngs_mapping: ../ngs_mapping
    path_variant_calling: ../variant_calling
    path_jannovar_ser: REQUIRED  # REQUIRED
    tools_ngs_mapping: []      # default to those configured for ngs_mapping
    tools_variant_calling: []  # default to those configured for variant_calling
    dbnsfp:  # configuration for default genome release, needs change if differing
      col_contig: 1
      col_pos: 2
      columns: []
    annotation_tracks_bed: []
    annotation_tracks_tsv: []
    annotation_tracks_vcf: []
"""


class JannovarAnnotateVcfStepPart(BaseStepPart):
    """Annotate VCF file using "Jannovar annotate-vcf" """

    name = "jannovar"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}/out/.done"
        )

    @dictify
    def get_input_files(self, action):
        """Return path to pedigree input file"""
        assert action == "annotate_vcf"
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        tpl = (
            "output/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}"
        )
        KEY_EXT = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["variant_calling"]
        for key, ext in KEY_EXT.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        assert action == "annotate_vcf"
        prefix = (
            "work/{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}"
        )
        KEY_EXT = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        for key, ext in KEY_EXT.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        assert action == "annotate_vcf"
        prefix = (
            "work/{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}/log/"
            "{mapper}.{var_caller}.jannovar_annotate_vcf.{index_ngs_library}"
        )

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["variant_annotation_jannovar_annotate_vcf"] = {
            "mem": 7 * 1024 * 2,
            "time": "100:00",
            "ntasks": 2,
        }


class VariantAnnotationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    name = "variant_annotation"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow,
            config,
            cluster_config,
            config_lookup_paths,
            config_paths,
            workdir,
            (VariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, JannovarAnnotateVcfStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        self.register_sub_workflow("variant_calling", self.config["path_variant_calling"])
        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_variant_calling"]:
            self.config["tools_variant_calling"] = self.w_config["step_config"]["variant_calling"][
                "tools"
            ]

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        token = "{mapper}.{caller}.jannovar_annotate_vcf.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", token, "out", token + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files(
            os.path.join("output", token, "log", token + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:  # pragma: no cover
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                elif not pedigree.index.dna_ngs_library:  # pragma: no cover
                    msg = "INFO: pedigree index without DNA NGS library (names: {})"
                    print(
                        msg.format(  # pragma: no cover
                            list(sorted(d.name for d in pedigree.donors))
                        ),
                        file=sys.stderr,
                    )
                    continue  # pragma: no cover
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_annotation", "path_variant_calling"),
            ("Path to variant calling not configured but required for variant annotation"),
        )
        self.ensure_w_config(
            ("step_config", "variant_annotation", "path_jannovar_ser"),
            ("Path to serialized Jannovar database"),
        )
