# -*- coding: utf-8 -*-
"""Implementation of the ``variant_calling`` step

The ``variant_calling` step takes the output of the ``ngs_mapping`` step and performs small
variant calling on the read alignments.  The output are variant calls in VCF (and optionally gVCF)
files and quality control statistics on these data.

==========
Properties
==========

overall stability

    **stable**

applicable to

    germline variant calling

generally applicable to

    short read variant calling

==========
Step Input
==========

BAM files from the ``ngs_mapping`` step.

===========
Step Output
===========

Creates one output directory for each read mapper (from ``ngs_mapping``), each variant caller, and
each pedigree from the germline sample sheet.

**Primary Output**

- ``output/{mapper}.{caller}.{index_library}/out/{mapper}.{caller}.{index_library}.vcf.gz``

**Additional Output**

The callers implementing a gVCF workflow (currently only ``gatk4_hc_gvcf``) also create one output
gVCF file for the pedigree.

- ``output/{mapper}.{caller}.{index_library}/out/{mapper}.{caller}.{index_library}.g.vcf.gz``

Further, each VCF and gVCF file gets an appropriate TBI index file ``{vcf_file}.tbi`` and each output
is gets an appropriate MD5 checksum file ``{file}.md5``.

====================
Global Configuration
====================

- If GATK HaplotypeCaller or GATK UnifiedGenotyper are activated then
  ``static_data_config/dbsnp/path`` must be properly configured
- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_calling.rst

===============
Variant Callers
===============

The following germline variant callers are currently available.

``gatk4_hc_gvcf``

    Variant calling with GATK v4 HaplotypeCaller using the gVCF workflow consisting of
    variant discovery with ``HaplotypeCaller``, merging of the gVCF files withing each
    pedigree with ``CombineGVCFs`` and genotyping with ``GenotypeGVCFs``.

    This is the mainly used variant caller and the only one enabled by default.

    The reason is this being the main advertised run mode by the GATK team and this workflow
    enables physical phasing information in the output VCF files.

``gatk4_hc_joint``

    Variant calling with the GATK v4 HaplotypeCaller using joint calling with direct VCF
    generation.

    This variant caller is provided as a fallback to explore problems with *de novo* variant
    calls that may have been introduced by the gVCF workflow.

    Disabled by default.

``gatk3_hc``

    Joint calling with GATK v3 HaplotypeCaller.

    This caller is provided for historical reasons as earlier versions of SNAPPY pipeline
    were based on this workflow.

    Disabled by default.

``gatk3_ug``

    Joint calling with GATK v3 UnifiedGenotyper.

    This caller is provided for historical reasons and to provide a vote in creating consensus
    sets of variant calls.

``bcftools_call``

    Variant calling with ``bcftools mpileup | bcftools call``.

    This caller is provided for establishing baseline variant calls in benchmark situations.
    BCFtools allows for fast and efficient variant calling at the cost of some sensitivity
    and specificity.

    Disabled by default.

``bcftools_call``

    Variant calling with Freebayes.

    This caller is provided as it is useful to genereate variant calls with the same artifacts
    as external pipelines.

    Disabled by default.

=======
Reports
=======

``jannovar_stats``

    Create statistics on variants using ``jannovar statsistics`` for each pedigree.

    ::

        report/jannovar_stats/{mapper}.{caller}.{index_library}.{donor_library}.txt

``bcftools_stats``

    Create statistics on variants using ``bcftools stats`` for each donor in each pedigree
    for each mapper and caller.

    ::

        report/bcftools_stats/{mapper}.{caller}.{index_library}.{donor_library}.txt

``baf_file_generation``

    Create one UCSC BigWig file for each individual in each pedigree for each mapper and caller
    with B-allele fraction.  These files can be used for to visually confirm structural variants
    or runs of homozygosity.

    ::

        report/baf/{mapper}.{caller}.{index_library}.{donor_library}.bw

``roh_calling``

    Perform run-of-homozygosity calling with ``bcftools roh``.

=========
Log Files
=========

For each variant caller and report generator, the following log files are created into the
``log`` directory.

``{file}.conda_info.txt``

    Output of ``conda info`` of the executing conda environment.

``{file}.conda_list.txt``

    Output of ``conda list`` of the executing conda environment with list of the full package
    list and exact versions.

``{file}.log``

    Log output of the execution.

``{file}.wrapper.py``

    The actual Snakemake wrapper file with all input / output / parameter values.

====================
Implementation Notes
====================

- All variant callers are parallelized using GNU parallel on genome-wide windows generated by
  GATK v4 ``PreprocessIntervals``.
- Each output file has an accompanying MD5 sum.

==============
Example Output
==============

Given a pedigree with index ``index`` and two more donors ``mother`` and ``father``, the following
files would be created into ``output/`` (each VCF file has a ``.tbi`` file and overall each file has
a ``.md5`` file).  In this case, the read mapper is ``bwa`` and the variant caller is ``gatk4_hc_gvcf``.


```
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.baf_file_generation_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.baf_file_generation_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.baf_file_generation_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.baf_file_generation_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.baf_file_generation_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.bcftools_stats_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.bcftools_stats_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.bcftools_stats_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.bcftools_stats_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.bcftools_stats_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.baf_file_generation_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.baf_file_generation_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.baf_file_generation_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.baf_file_generation_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.baf_file_generation_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.bcftools_stats_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.bcftools_stats_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.bcftools_stats_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.bcftools_stats_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.bcftools_stats_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.baf_file_generation_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.baf_file_generation_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.baf_file_generation_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.baf_file_generation_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.baf_file_generation_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.bcftools_stats_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.bcftools_stats_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.bcftools_stats_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.bcftools_stats_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.bcftools_stats_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.gatk4_hc_gvcf_genotype.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.gatk4_hc_gvcf_genotype.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.gatk4_hc_gvcf_genotype.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.gatk4_hc_gvcf_genotype.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.gatk4_hc_gvcf_genotype.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.jannovar_stats_run.conda_info.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.jannovar_stats_run.conda_list.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.jannovar_stats_run.environment.yaml
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.jannovar_stats_run.log
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/log/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.jannovar_stats_run.wrapper.py
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/out/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.g.vcf.gz
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/out/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.vcf.gz
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/baf/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.baf.bw
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/baf/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.baf.bw
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/baf/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.baf.bw
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/bcftools_stats/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.index-N1-DNA1-WES1.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/bcftools_stats/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.father-N1-DNA1-WES1.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/bcftools_stats/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.mother-N1-DNA1-WES1.txt
bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1/report/jannovar_stats/bwa.gatk4_hc_gvcf.index-N1-DNA1-WES1.txt
```
"""

from collections import OrderedDict
from itertools import chain
import re
import sys
import typing
import warnings

from biomedsheets.shortcuts import GermlineCaseSheet, Pedigree, is_not_background
from snakemake.io import Wildcards, expand

from snappy_pipeline.utils import dictify, flatten, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import (
    SnakemakeDict,
    SnakemakeDictItemsGenerator,
    SnakemakeListItemsGenerator,
)
from snappy_pipeline.workflows.abstract.exceptions import InvalidConfigurationException
from snappy_pipeline.workflows.abstract.warnings import InconsistentPedigreeWarning
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available germline variant callers
VARIANT_CALLERS = (
    "bcftools_call",
    "freebayes",
    "gatk3_hc",
    "gatk3_ug",
    "gatk4_hc_gvcf",
    "gatk4_hc_joint",
)

#: Default configuration for the variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_calling
step_config:
  variant_calling:
    # Common configuration
    path_ngs_mapping: ../ngs_mapping  # REQUIRED

    # Report generation
    baf_file_generation:
      enabled: true
      min_dp: 10  # minimal DP of variant, must be >=1
    bcftools_stats:
      enabled: true
    jannovar_stats:
      enabled: true
      path_ser: REQUIRED  # REQUIRED
    bcftools_roh:
      enabled: true
      path_targets: null  # REQUIRED; optional
      path_af_file: null  # REQUIRED
      ignore_homref: false
      skip_indels: false
      rec_rate: 1e-8

    # Variant calling tools and their configuration
    #
    # Common configuration
    tools: ['gatk4_hc_gvcf']  # REQUIRED
    ignore_chroms:
    - '^NC_007605$' # herpes virus
    - '^hs37d5$'    # GRCh37 decoy
    - '^chrEBV$'    # Eppstein-Barr Virus
    - '_decoy$'     # decoy contig
    - '^HLA-'       # HLA genes

    # Variant caller specific configuration
    bcftools_call:
      max_depth: 250
      max_indel_depth: 250
      window_length: 10000000
      num_threads: 16
    freebayes:
      use_standard_filters: true
      window_length: 10000000
      num_threads: 16
      min_alternate_fraction: 0.05  # FreeBayes default
      min_mapping_quality: 1        # FreeBayes default
      min_repeat_entropy: 1         # FreeBayes default
      haplotype_length: 3           # FreeBayes default
    gatk3_hc:
      num_threads: 16
      window_length: 10000000
      allow_seq_dict_incompatibility: false
    gatk3_ug:
      num_threads: 16
      window_length: 10000000
      allow_seq_dict_incompatibility: false
      downsample_to_coverage: 250
    gatk4_hc_joint:
      window_length: 10000000
      num_threads: 16
      allow_seq_dict_incompatibility: false
    gatk4_hc_gvcf:
      window_length: 10000000
      num_threads: 16
      allow_seq_dict_incompatibility: false
"""


class GetResultFilesMixin:
    """Mixin to provide ``get_result_files()`` function for variant calling and
    variant annotation steps.

    The function will use ``get_output_files()`` for all actions to obtain
    the files in the ``output/directory`` to expect from the given step.
    """

    @listify
    def get_result_files(self) -> SnakemakeListItemsGenerator:
        """Return concrete result file paths"""

        def strip_tpl(tpl):
            return tpl.replace(r",[^\.]+", "")

        index_dna_ngs_libraries = self._get_index_dna_ngs_libraries()
        for action in self.actions:
            output_files = self.get_output_files(action).values()
            result_paths_tpls = [
                strip_tpl(p) for p in flatten(output_files) if p.startswith("output/")
            ]
            for path_tpl in result_paths_tpls:
                for index_library_name, member_library_names in index_dna_ngs_libraries.items():
                    kwargs = {
                        "mapper": self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                    }
                    if "index_library_name" in path_tpl:
                        kwargs["index_library_name"] = [index_library_name]
                        kwargs["donor_library_name"] = member_library_names
                    else:
                        kwargs["library_name"] = [index_library_name]
                    for key, value in self.get_extra_kv_pairs().items():
                        if "{%s}" % key in path_tpl:
                            kwargs[key] = value
                    yield from expand(
                        path_tpl,
                        **kwargs,
                    )

    def get_extra_kv_pairs(self):
        return {"var_caller": self.parent.config["tools"]}

    @dictify
    def _get_index_dna_ngs_libraries(
        self,
    ) -> typing.Generator[typing.Tuple[str, typing.List[str]], None, None]:
        """Return ``dict`` that maps the index DNA library name to a list of all pedigree
        member's DNA library names.
        """
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if self._is_pedigree_good(pedigree):
                    index = pedigree.index.dna_ngs_library.name
                    donors = [
                        donor.dna_ngs_library.name
                        for donor in pedigree.donors
                        if donor.dna_ngs_library
                    ]
                    yield index, donors

    def _is_pedigree_good(self, pedigree: Pedigree) -> bool:
        """Check pedigrees for inconsistencies and issue warning for any.

        :return: ``True`` if there was no inconsistency reported
        """
        msg = None
        donor_names = list(sorted(d.name for d in pedigree.donors))
        if not pedigree.index:  # pragma: no cover
            msg = f"INFO: pedigree without index (name: {donor_names})"
        elif not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = f"INFO: pedigree index without DNA NGS library (names: {donor_names})"
        if msg:
            warnings.warn(InconsistentPedigreeWarning(msg))
        return not msg


class VariantCallingGetLogFileMixin:
    """Mixin to provide the generic ``get_log_file()`` function for variant calling"""

    @dictify
    def get_log_file(self, action) -> SnakemakeDictItemsGenerator:
        """Return dict of log files in the "log" directory."""
        _ = action
        token = f"{{mapper}}.{self.name}.{{library_name}}"
        prefix = f"work/{token}/log/{token}.{self.name}_{action}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, f"{prefix}{ext}"
            yield f"{key}_md5", f"{prefix}{ext}.md5"


class VariantCallingStepPart(GetResultFilesMixin, VariantCallingGetLogFileMixin, BaseStepPart):
    """Base class for germline variant calling step parts

    Variant calling is performed on a per-pedigree level.  The (single) index individual is used
    for naming the output file.
    """

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{index_library_name}}/out/"
            "{{mapper}}.{var_caller}.{{index_library_name}}{ext}"
        )
        self.base_path_tmp = self.base_path_out.replace("/out/", "/tmp/")
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action) -> SnakemakeDict:
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]

        if not pedigree.index or not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = "INFO: pedigree without index (names: {})"
            donor_names = list(sorted(d.name for d in pedigree.donors))
            print(msg.format(donor_names), file=sys.stderr)
            yield "bam", []
            yield "ped", None
        else:
            library_name = pedigree.index.dna_ngs_library.name
            yield "ped", f"work/write_pedigree.{library_name}/out/{library_name}.ped"

            bams = []
            for donor in pedigree.donors:
                if not donor.dna_ngs_library:
                    continue  # skip
                infix = f"{wildcards.mapper}.{donor.dna_ngs_library.name}"
                bams.append(ngs_mapping(f"output/{infix}/out/{infix}.bam"))
            yield "bam", bams

    def get_output_files(self, action) -> SnakemakeDict:
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        token = f"{{mapper}}.{self.name}.{{library_name}}"
        work_files = {
            "vcf": f"work/{token}/out/{token}.vcf.gz",
            "vcf_md5": f"work/{token}/out/{token}.vcf.gz.md5",
            "vcf_tbi": f"work/{token}/out/{token}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{token}/out/{token}.vcf.gz.tbi.md5",
        }
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]


class BcftoolsCallStepPart(VariantCallingStepPart):
    """Germline variant calling with bcftools"""

    #: Step name
    name = "bcftools_call"

    def get_resource_usage(self, action: str) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            time="2-00:00:00",
            memory=f"{int(3.75 * 1024 * 16)}M",
        )


class FreebayesStepPart(VariantCallingStepPart):
    """Germline variant calling with freebayes"""

    #: Step name
    name = "freebayes"

    def get_resource_usage(self, action):
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            time="2-00:00:00",  # 2 days
            memory=f"{int(3.75 * 1024 * 16)}M",
        )

    def get_params(self, action):
        """
        :param action: Action (i.e., step) in the workflow. Currently only available for 'run'.
        :type action: str
        :return: Returns get parameters function.
        :raises UnsupportedActionException: if action not 'run'.
        """
        self._validate_action(action)
        parameters_key_list = [
            "window_length",
            "min_alternate_fraction",
            "min_mapping_quality",
            "min_repeat_entropy",
            "haplotype_length",
        ]
        return {key: self.config["freebayes"][key] for key in parameters_key_list}


class GatkCallerStepPartBase(VariantCallingStepPart):
    """Base class for GATK v3/v4 variant callers"""

    def check_config(self):
        if self.__class__.name not in self.config["tools"]:
            return  # caller not enabled, skip  # pragma: no cover
        self.parent.ensure_w_config(
            ("static_data_config", "dbsnp", "path"),
            "dbSNP not configured but required for {}".format(self.__class__.name),
        )

    def get_resource_usage(self, action) -> ResourceUsage:
        self._validate_action(action)
        num_threads = self.config[self.name]["num_threads"]
        mem_per_thread = 5.5
        mem_total = int(mem_per_thread * num_threads + 0.5)
        return ResourceUsage(
            threads=num_threads,
            time="2-00:00:00",
            memory=f"{mem_total}G",
        )


class Gatk3HaplotypeCallerStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK v3 HaplotypeCaller"""

    #: Step name
    name = "gatk3_hc"


class Gatk3UnifiedGenotyperStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK v3 UnifiedGenotyper"""

    #: Step name
    name = "gatk3_ug"


class Gatk4HaplotypeCallerJointStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller doing joint calling per pedigree"""

    name = "gatk4_hc_joint"


class Gatk4HaplotypeCallerGvcfStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller and gVCF workflow"""

    name = "gatk4_hc_gvcf"

    actions = ("discover", "combine_gvcfs", "genotype")

    @dictify
    def _get_input_files_discover(self, wildcards):
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        bam_path = f"output/{infix}/out/{infix}.bam"
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        yield "bam", ngs_mapping(bam_path)

    @dictify
    def _get_input_files_combine_gvcfs(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]

        if not pedigree.index or not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = "INFO: pedigree without index (names: {})"
            donor_names = list(sorted(d.name for d in pedigree.donors))
            print(msg.format(donor_names), file=sys.stderr)
            yield "ped", None
            yield "gvcf", []
        else:
            library_name = pedigree.index.dna_ngs_library.name
            yield "ped", f"work/write_pedigree.{library_name}/out/{library_name}.ped"

            gvcfs = []
            for donor in pedigree.donors:
                if not donor.dna_ngs_library:
                    continue  # skip
                infix = f"{wildcards.mapper}.gatk4_hc_gvcf_discover.{donor.dna_ngs_library.name}"
                gvcfs.append(f"work/{infix}/out/{infix}.g.vcf.gz")
            yield "gvcf", gvcfs

    @dictify
    def _get_input_files_genotype(self, wildcards) -> SnakemakeDictItemsGenerator:
        infix = f"{wildcards.mapper}.gatk4_hc_gvcf_combine_gvcfs.{wildcards.library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_md5", f"work/{infix}/out/{infix}.g.vcf.gz.md5"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"
        yield "gvcf_tbi_md5", f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5"

    @dictify
    def _get_output_files_discover(self) -> SnakemakeDictItemsGenerator:
        infix = "{mapper}.gatk4_hc_gvcf_discover.{library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_md5", f"work/{infix}/out/{infix}.g.vcf.gz.md5"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"
        yield "gvcf_tbi_md5", f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5"
        yield "output_links", []

    @dictify
    def _get_output_files_combine_gvcfs(self) -> SnakemakeDictItemsGenerator:
        infix = "{mapper}.gatk4_hc_gvcf_combine_gvcfs.{library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_md5", f"work/{infix}/out/{infix}.g.vcf.gz.md5"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"
        yield "gvcf_tbi_md5", f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5"
        yield "output_links", []

    @dictify
    def _get_output_files_genotype(self) -> SnakemakeDictItemsGenerator:
        infix = "{mapper}.gatk4_hc_gvcf.{library_name}"
        result = {
            "gvcf": f"work/{infix}/out/{infix}.g.vcf.gz",
            "gvcf_md5": f"work/{infix}/out/{infix}.g.vcf.gz.md5",
            "gvcf_tbi": f"work/{infix}/out/{infix}.g.vcf.gz.tbi",
            "gvcf_tbi_md5": f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5",
            "vcf": f"work/{infix}/out/{infix}.vcf.gz",
            "vcf_md5": f"work/{infix}/out/{infix}.vcf.gz.md5",
            "vcf_tbi": f"work/{infix}/out/{infix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{infix}/out/{infix}.vcf.gz.tbi.md5",
        }
        yield from result.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(result.values(), self.get_log_file("genotype").values())
        ]


class ReportGetLogFileMixin:
    """Log file generation for reports"""

    #: Whether we generate per-donor files.
    report_per_donor = None

    @dictify
    def get_log_file(self, action: str) -> SnakemakeDictItemsGenerator:
        """Return dict of log files in the "log" directory."""
        self._validate_action(action)
        assert self.report_per_donor is not None
        token = "{mapper}.{var_caller}.{index_library_name}"
        prefix = f"work/{token}/log/{token}.{{donor_library_name}}.{self.name}_{action}"
        if not self.report_per_donor:
            prefix = prefix.replace("{donor_library_name}.", "")
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, f"{prefix}{ext}"
            yield f"{key}_md5", f"{prefix}{ext}.md5"


class BcftoolsStatsStepPart(GetResultFilesMixin, ReportGetLogFileMixin, BaseStepPart):
    """VCF statistics computation with ``bcftools stats``.

    Statistics are computed overall and per-sample
    """

    name = "bcftools_stats"

    actions = ("run",)

    report_per_donor = True

    def __init__(self, parent):
        super().__init__(parent)

    def get_input_files(self, action: str) -> SnakemakeDict:
        """Return required input files"""
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")()

    @dictify
    def _get_input_files_run(self) -> SnakemakeDictItemsGenerator:
        yield "vcf", (
            "work/{mapper}.{var_caller}.{index_library_name}/out/"
            "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
        )

    def get_output_files(self, action: str) -> SnakemakeDict:
        """Return step part output files"""
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        base_path = (
            "work/{mapper}.{var_caller}.{index_library_name}/report/bcftools_stats/"
            "{mapper}.{var_caller}.{index_library_name}.{donor_library_name}"
        )
        work_files = {key: f"{base_path}{ext}" for key, ext in ext_names.items()}
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]

    def get_resource_usage(self, action: str) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="02:00:00",
            memory="1024M",
        )


class BcftoolsRohStepPart(GetResultFilesMixin, ReportGetLogFileMixin, BaseStepPart):
    """ROH calling with `bcftools roh`."""

    name = "bcftools_roh"
    actions = ("run",)
    report_per_donor = False

    def get_input_files(self, action: str) -> SnakemakeDict:
        """Return required input files"""
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")()

    @dictify
    def _get_input_files_run(self) -> SnakemakeDictItemsGenerator:
        yield "vcf", (
            "output/{mapper}.{var_caller}.{index_library_name}/out/"
            "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
        )

    def get_output_files(self, action: str) -> SnakemakeDict:
        """Return step part output files"""
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        base_path = (
            "work/{mapper}.{var_caller}.{index_library_name}/report/roh/"
            "{mapper}.{var_caller}.{index_library_name}"
        )
        work_files = {key: f"{base_path}{ext}" for key, ext in ext_names.items()}
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]

    def get_resource_usage(self, action: str) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="01:00:00",
            memory="4G",
        )


class JannovarStatisticsStepPart(GetResultFilesMixin, ReportGetLogFileMixin, BaseStepPart):
    """Base class for VCF statistics computation with "jannovar statistics"

    Statistics are computed overall and per-sample
    """

    name = "jannovar_stats"
    actions = ("run",)
    report_per_donor = False

    @dictify
    def get_input_files(self, action) -> SnakemakeDictItemsGenerator:
        """Return path to input files"""
        self._validate_action(action)
        yield "vcf", (
            "work/{mapper}.{var_caller}.{index_library_name}/out/"
            "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
        )

    def get_output_files(self, action) -> SnakemakeDict:
        """Return step part output files"""
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        """Return output files that all germline variant calling sub steps must return (VCF +
        TBI file)
        """
        base_path = (
            "work/{mapper}.{var_caller}.{index_library_name}/report/jannovar_stats/"
            "{mapper}.{var_caller}.{index_library_name}"
        )
        ext_names = {"report": ".txt", "report_md5": ".txt.md5"}
        work_files = {key: f"{base_path}{ext}" for key, ext in ext_names.items()}
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]

    def get_resource_usage(self, action: str) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="04:00:00",
            memory=f"{int(3.75 * 1024 * 2)}M",
        )


class BafFileGenerationStepPart(GetResultFilesMixin, ReportGetLogFileMixin, BaseStepPart):
    """Class for computing B allele fraction files.

    One file is generated per sample in the output VCF files.
    """

    #: Step name
    name = "baf_file_generation"

    #: Class available actions
    actions = ("run",)

    #: Whether to report results per donor
    report_per_donor = True

    @dictify
    def get_input_files(self, action: str) -> SnakemakeDictItemsGenerator:
        self._validate_action(action)
        yield "vcf", (
            "work/{mapper}.{var_caller}.{index_library_name}/out/"
            "{mapper}.{var_caller}.{index_library_name}.vcf.gz"
        )

    @dictify
    def get_output_files(self, action: str) -> SnakemakeDictItemsGenerator:
        self._validate_action(action)
        base_path = (
            "{mapper}.{var_caller}.{index_library_name}/report/baf/"
            r"{mapper}.{var_caller}.{index_library_name}.{donor_library_name,[^\.]+}.baf"
        )
        ext_names = {"bw": ".bw", "bw_md5": ".bw.md5"}
        work_files = {}
        for key, ext in ext_names.items():
            work_files[key] = f"work/{base_path}{ext}"
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]

    def get_resource_usage(self, action: str) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="02:00:00",
            memory="1024M",
        )


class VariantCallingWorkflow(BaseStep):
    """Workflow implementation for germline variant calling"""

    name = "variant_calling"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls) -> str:
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
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                BcftoolsCallStepPart,
                Gatk3HaplotypeCallerStepPart,
                Gatk3UnifiedGenotyperStepPart,
                Gatk4HaplotypeCallerJointStepPart,
                Gatk4HaplotypeCallerGvcfStepPart,
                BcftoolsStatsStepPart,
                BcftoolsRohStepPart,
                JannovarStatisticsStepPart,
                BafFileGenerationStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self) -> SnakemakeListItemsGenerator:
        for tool in self.config["tools"]:
            yield from self.sub_steps[tool].get_result_files()
        for name in ("baf_file_generation", "bcftools_stats", "jannovar_stats", "bcftools_roh"):
            if self.w_config["step_config"]["variant_calling"][name]["enabled"]:
                yield from self.sub_steps[name].get_result_files()

    def check_config(self):
        # Checks for static data
        self.ensure_w_config(
            ("step_config", "variant_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for variant calling",
        )
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant calling",
        )
        # Check that only valid tools are selected
        selected = set(self.w_config["step_config"]["variant_calling"]["tools"])
        invalid = list(sorted(selected - set(VARIANT_CALLERS)))
        if invalid:
            raise InvalidConfigurationException(f"Invalid variant callers selected: {invalid}")
