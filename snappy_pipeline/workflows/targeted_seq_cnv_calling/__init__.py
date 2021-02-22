# -*- coding: utf-8 -*-
"""Implementation of the ``targeted_seq_cnv_calling`` step

This step allows for the detection of CNV events for germline samples from targeted sequenced
(e.g., exomes or large panels).  The wrapped tools start from the aligned reads (thus off
``ngs_mapping``) and generate CNV calls for germline variants.

The wrapped tools implement different strategies.  Some work "reference free" and just use the
germline BAM files for their input, others need the germlien BAM files and additionally a
set of BAM files generated by the same wet-lab and Bioinformatics protocols for their background.

.. note::

    Status: not implemented yet

==========
Step Input
==========

Germline CNV calling for targeted sequencing starts off the aligned reads, i.e.,
``ngs_mapping``.

===========
Step Output
===========

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_targeted_seq_cnv_calling.rst

=====================
Available CNV Callers
=====================

- ``xhmm``

"""

from collections import OrderedDict
import os
import re

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import glob_wildcards, expand, touch

from ..abstract import BaseStepPart, BaseStep, LinkOutStepPart
from ..ngs_mapping import NgsMappingWorkflow
from ...utils import listify, dictify, DictQuery

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload (VCF)
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Available WGS CNV callers
TARGETED_SEQ_CNV_CALLERS = ("xhmm", "gcnv")

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Default configuration for the targeted_seq_cnv_calling step
DEFAULT_CONFIG = r"""
# Default configuration targeted_seq_cnv_calling
step_config:
  targeted_seq_cnv_calling:
    # You can select select between the following modes of creating chunks.
    # The CNV calling will be performed on each chunk individually, all
    # individuals of a family will go into one chunk.
    #
    # single -- create one large chunk for all
    # evenly -- create evenly sized chunks of maximal size chunk_max_size
    # incremental -- iterate samples and create a new chunk after the maximal
    #       size has been reached; should be combined with background data
    #       to guarantee a minimal number of samples.
    chunk_mode: single
    # Maximal chunk size; background samples are not counted.
    chunk_max_size: 200
    # Sort samples by these attributes for chunk_mode=incremental to create
    # chunks.
    #
    # family_id -- the family ID from the sample sheet
    # index_name -- the name of the pedigree's index
    # batch_no -- the batch number from the sample sheet
    chunk_sort:
      - batch_no
      - family_id
      - index_name

    # Path to the ngs_mapping step.
    path_ngs_mapping: ../ngs_mapping
    path_variant_calling: ../variant_calling
    tools:
    - xhmm
    xhmm:
      path_target_interval_list: REQUIRED_OR_MAPPING # REQUIRED
      # The following allows to define one or more set of target intervals.
      path_target_interval_list_mapping: []
      # The following will match both the stock IDT library kit and the ones
      # with spike-ins seen fromr Yale genomics.  The path above would be
      # mapped to the name "default".
      # - name: IDT_xGen_V1_0
      #   pattern: "xGen Exome Research Panel V1\\.0*"
      #   path: "path/to/targets.bed"
    gcnv:
      path_target_interval_list: REQUIRED_OR_MAPPING # REQUIRED
      # The following allows to define one or more set of target intervals.
      path_target_interval_list_mapping: []
      # The following will match both the stock IDT library kit and the ones
      # with spike-ins seen fromr Yale genomics.  The path above would be
      # mapped to the name "default".
      # - name: IDT_xGen_V1_0
      #   pattern: "xGen Exome Research Panel V1\\.0*"
      #   path: "path/to/targets.bed"
      # Path to BED file with uniquely mappable regions.
      path_uniquely_mapable_bed: REQUIRED
"""


class XhmmStepPart(BaseStepPart):
    """Targeted seq. CNV calling with XHMM"""

    name = "xhmm"

    actions = (
        "coverage",
        "merge_cov",
        "ref_stats",
        "filter_center",
        "pca",
        "normalize",
        "zscore_center",
        "refilter",
        "discover",
        "genotype",
        "extract_ped",
    )

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to donor
        self.index_ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_donor.update(sheet.index_ngs_library_to_donor)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    def get_params(self, action):
        assert action == "coverage"

        def get_params(wildcards):
            return {"library_kit": self.ngs_library_to_kit[wildcards.library_name]}

        return get_params

    @dictify
    def _build_ngs_library_to_kit(self):
        xhmm_config = DictQuery(self.w_config).get("step_config/targeted_seq_cnv_calling/xhmm")
        if not xhmm_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.parent._all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"

        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in xhmm_config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.parent._all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    def get_input_files(self, action):
        """Return input function for XHMM rule"""
        assert action in self.actions
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def _get_input_files_coverage(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Yield input BAM and BAI file
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield key, ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    def _get_input_files_merge_cov(self, wildcards):
        token = "{mapper}.xhmm_coverage.{lib}"
        summaries = [
            "work/{token}/out/{token}.DATA.sample_interval_summary".format(token=token).format(
                lib=lib, **wildcards
            )
            for lib in sorted(self.index_ngs_library_to_donor)
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit
        ]
        return summaries

    @dictify
    def _get_input_files_filter_center(self, wildcards):
        yield "merge_cov", self._get_output_files_merge_cov()[0].format(**wildcards)
        yield "extreme_gc", self._get_output_files_ref_stats()["extreme_gc_targets"].format(
            **wildcards
        )

    def _get_input_files_pca(self, wildcards):
        return [self._get_output_files_filter_center()["centered"].format(**wildcards)]

    @dictify
    def _get_input_files_normalize(self, wildcards):
        yield "centered", self._get_output_files_filter_center()["centered"].format(**wildcards)
        yield "pca", self._get_output_files_pca()["pc"].format(**wildcards)

    def _get_input_files_zscore_center(self, wildcards):
        token = "{mapper}.xhmm_normalize.{library_kit}".format(**wildcards)
        return ["work/{token}/out/{token}".format(token=token)]

    @dictify
    def _get_input_files_refilter(self, wildcards):
        token = "{mapper}.xhmm_merge_cov.{library_kit}".format(**wildcards)
        yield "original", "work/{token}/out/{token}.RD.txt".format(token=token)
        for infix in ("filter_center", "zscore_center"):
            for kvs in (
                ("filtered_samples", ".filtered_samples.txt"),
                ("filtered_targets", ".filtered_targets.txt"),
            ):
                token = "{mapper}.xhmm_{infix}.{library_kit}".format(infix=infix, **wildcards)
                key = "{}_{}".format(kvs[0], infix)
                yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=kvs[1])

    @dictify
    def _get_input_files_discover(self, wildcards):
        token = "{mapper}.xhmm_zscore_center.{library_kit}".format(**wildcards)
        yield "center_zscore", "work/{token}/out/{token}".format(token=token)
        token = "{mapper}.xhmm_refilter.{library_kit}".format(**wildcards)
        yield "refilter_original", "work/{token}/out/{token}.RD.txt".format(token=token)

    @dictify
    def _get_input_files_genotype(self, wildcards):
        token = "{mapper}.xhmm_zscore_center.{library_kit}".format(**wildcards)
        yield "center_zscore", "work/{token}/out/{token}".format(token=token)
        token = "{mapper}.xhmm_refilter.{library_kit}".format(**wildcards)
        yield "refilter_original", "work/{token}/out/{token}.RD.txt".format(token=token)
        token = "{mapper}.xhmm_discover.{library_kit}".format(**wildcards)
        yield "discover_xcnv", "work/{token}/out/{token}.xcnv".format(token=token)

    @dictify
    def _get_input_files_extract_ped(self, wildcards):
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        token = "bwa.xhmm_filter_center.{library_kit}".format(library_kit=library_kit)
        yield (
            "filtered_samples",
            "work/{token}/out/{token}.filtered_samples.txt".format(token=token),
        )
        token = "{mapper}.xhmm_genotype.{library_kit}".format(library_kit=library_kit, **wildcards)
        for key, ext in (("vcf", ".vcf.gz"), ("tbi", ".vcf.gz.tbi")):
            yield key, "work/{token}/out/{token}{ext}".format(token=token, ext=ext)

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    def get_output_files(self, action):
        """Return output files that XHMM creates for the given action."""
        assert action in self.actions
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_coverage(self):
        exts = (
            "sample_interval_statistics",
            "sample_interval_summary",
            "sample_statistics",
            "sample_summary",
        )
        for ext in exts:
            token = "{mapper}.xhmm_coverage.{library_name}"
            yield ext, "work/{token}/out/{token}.DATA.{ext}".format(token=token, ext=ext)

    def _get_output_files_merge_cov(self):
        token = "{mapper}.xhmm_merge_cov.{library_kit}"
        return ["work/{token}/out/{token}.RD.txt".format(token=token)]

    @dictify
    def _get_output_files_ref_stats(self):
        token = "{mapper}.xhmm_ref_stats.{library_kit}"
        for infix in ("extreme_gc_targets",):
            yield infix, "work/{token}/out/{token}.{infix}.txt".format(token=token, infix=infix)

    @dictify
    def _get_output_files_filter_center(self):
        token = "{mapper}.xhmm_filter_center.{library_kit}"
        for infix in ("centered", "filtered_targets", "filtered_samples"):
            yield infix, "work/{token}/out/{token}.{infix}.txt".format(token=token, infix=infix)

    @dictify
    def _get_output_files_pca(self):
        token = "{mapper}.xhmm_pca.{library_kit}"
        kvs = (("pc_loading", ".PC_LOADINGS.txt"), ("pc_sd", ".PC_SD.txt"), ("pc", ".PC.txt"))
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    @dictify
    def _get_output_files_normalize(self):
        token = "{mapper}.xhmm_normalize.{library_kit}"
        kvs = (("normalized", ""), ("num_removed", ".num_removed_PC.txt"))
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    @dictify
    def _get_output_files_zscore_center(self):
        token = "{mapper}.xhmm_zscore_center.{library_kit}"
        kvs = (
            ("zscore_center", ""),
            ("filtered_samples", ".filtered_samples.txt"),
            ("filtered_targets", ".filtered_targets.txt"),
        )
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    def _get_output_files_refilter(self):
        token = "{mapper}.xhmm_refilter.{library_kit}"
        return ["work/{token}/out/{token}.RD.txt".format(token=token)]

    @dictify
    def _get_output_files_discover(self):
        token = "{mapper}.xhmm_discover.{library_kit}"
        kvs = (("xcnv", ".xcnv"), ("aux_xcnv", ".aux_xcnv"))
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    @dictify
    def _get_output_files_genotype(self):
        token = "{mapper}.xhmm_genotype.{library_kit}"
        kvs = (
            ("vcf", ".vcf.gz"),
            ("vcf_md5", ".vcf.gz.md5"),
            ("tbi", ".vcf.gz.tbi"),
            ("tbi_md5", ".vcf.gz.tbi.md5"),
        )
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    @dictify
    def _get_output_files_extract_ped(self):
        token = "{mapper}.xhmm.{library_name}"
        kvs = (
            ("vcf", ".vcf.gz"),
            ("vcf_md5", ".vcf.gz.md5"),
            ("tbi", ".vcf.gz.tbi"),
            ("tbi_md5", ".vcf.gz.tbi.md5"),
        )
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    def get_log_file(self, action):
        """Return path to log file"""
        if action == "coverage":
            return (
                "work/{{mapper}}.xhmm_{action}.{{library_name}}/log/"
                "snakemake.targeted_seq_cnv_calling.log"
            ).format(action=action)
        elif action == "extract_ped":
            return "work/{mapper}.xhmm.{library_name}/log/" "snakemake.targeted_seq_cnv_calling.log"
        else:
            return (
                "work/{{mapper}}.xhmm_{action}.{{library_kit}}/log/"
                "snakemake.targeted_seq_cnv_calling.log"
            ).format(action=action)

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration for XHMM CNV calling"""
        for action in self.actions:
            if action == "merge_cov":
                cluster_config["targeted_seq_cnv_calling_xhmm_{}".format(action)] = {
                    "mem": 12 * 1024,
                    "time": "24:00",
                    "ntasks": 1,
                }
            else:
                cluster_config["targeted_seq_cnv_calling_xhmm_{}".format(action)] = {
                    "mem": 12 * 1024,
                    "time": "08:00",
                    "ntasks": 1,
                }


class GcnvStepPart(BaseStepPart):
    """Targeted seq. CNV calling with GATK4 gCNV"""

    name = "gcnv"

    actions = (
        "preprocess_intervals",
        "annotate_gc",
        "filter_intervals",
        "scatter_intervals",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "merge_cohort_vcfs",
        "extract_ped",
    )

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to donor
        self.index_ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_donor.update(sheet.index_ngs_library_to_donor)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    def get_params(self, action):
        assert action == "coverage"

        def get_params(wildcards):
            return {"library_kit": self.ngs_library_to_kit[wildcards.library_name]}

        return get_params

    @dictify
    def _build_ngs_library_to_kit(self):
        gcnv_config = DictQuery(self.w_config).get("step_config/targeted_seq_cnv_calling/gcnv")
        if not gcnv_config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.parent._all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"

        # Build mapping.
        regexes = {
            item["pattern"]: item["name"]
            for item in gcnv_config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.parent._all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    def get_input_files(self, action):
        """Return input function for gCNV rule"""
        assert action in self.actions
        return getattr(self, "_get_input_files_{}".format(action))

    def _get_input_files_preprocess_intervals(self, wildcards):
        return {}

    @dictify
    def _get_input_files_annotate_gc(self, wildcards):
        token = "gcnv_preprocess_intervals.{wildcards.library_kit}".format(wildcards=wildcards)
        ext = "interval_list"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)

    @dictify
    def _get_input_files_filter_intervals(self, wildcards):
        yield from self._get_input_files_annotate_gc(wildcards).items()
        token = "gcnv_annotate_gc.{wildcards.library_kit}".format(wildcards=wildcards)
        ext = "tsv"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)
        key = "covs"
        covs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                token = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                covs.append("work/{token}/out/{token}.{ext}".format(token=token, ext="tsv"))
        yield key, covs

    @dictify
    def _get_input_files_scatter_intervals(self, wildcards):
        ext = "interval_list"
        token = "{mapper}.gcnv_filter_intervals.{library_kit}".format(**wildcards)
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards):
        ext = "interval_list"
        token = "{mapper}.gcnv_filter_intervals.{library_kit}"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                token = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                tsvs.append("work/{token}/out/{token}.{ext}".format(token=token, ext=ext))
        yield ext, tsvs

    @dictify
    def _get_input_files_coverage(self, wildcards):
        # Yield .interval list file.
        ext = "interval_list"
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        token = "gcnv_preprocess_intervals.{library_kit}".format(library_kit=library_kit)
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)
        # Yield input BAM and BAI file
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield key, ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_input_files_call_cnvs(self, wildcards):
        yield "interval_list_shard", "work/{token}/out/{token}/temp_{{shard}}/scattered.interval_list".format(
            token="{mapper}.gcnv_scatter_intervals.{library_kit}"
        )
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                token = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                tsvs.append("work/{token}/out/{token}.{ext}".format(token=token, ext=ext))
        yield ext, tsvs
        ext = "ploidy"
        token = "{mapper}.gcnv_contig_ploidy.{library_kit}".format(**wildcards)
        yield ext, "work/{token}/out/{token}/.done".format(token=token)
        key = "intervals"
        token = "gcnv_annotate_gc.{library_kit}"
        yield key, "work/{token}/out/{token}.{ext}".format(token=token, ext="tsv")

    @dictify
    def _get_input_files_post_germline_calls(self, wildcards, checkpoints):
        checkpoint = checkpoints.targeted_seq_cnv_calling_gcnv_scatter_intervals
        library_kit = self.ngs_library_to_kit.get(wildcards.library_name)
        scatter_out = checkpoint.get(library_kit=library_kit, **wildcards).output[0]
        shards = list(
            map(
                os.path.basename,
                glob_wildcards(os.path.join(scatter_out, "temp_{shard}/{file}")).shard,
            )
        )
        token = "{mapper}.gcnv_call_cnvs.{library_kit}".format(library_kit=library_kit, **wildcards)
        yield "calls", [
            "work/{token}.{shard}/out/{token}.{shard}/.done".format(token=token, shard=shard)
            for shard in shards
        ]
        ext = "ploidy"
        token = "{mapper}.gcnv_contig_ploidy.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        yield ext, "work/{token}/out/{token}/.done".format(token=token)

    @listify
    def _get_input_files_merge_cohort_vcfs(self, wildcards):
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                token = "{mapper}.gcnv_post_germline_calls.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                yield "work/{token}/out/{token}.vcf.gz".format(token=token)

    @dictify
    def _get_input_files_extract_ped(self, wildcards):
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        token = "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}".format(
            library_kit=library_kit, **wildcards
        )
        for key, ext in (("vcf", ".vcf.gz"), ("tbi", ".vcf.gz.tbi")):
            yield key, "work/{token}/out/{token}{ext}".format(token=token, ext=ext)

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    def get_output_files(self, action):
        """Return output files that gCNV creates for the given action."""
        assert action in self.actions
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_preprocess_intervals(self):
        ext = "interval_list"
        token = "gcnv_preprocess_intervals.{library_kit}"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)

    @dictify
    def _get_output_files_annotate_gc(self):
        ext = "tsv"
        token = "gcnv_annotate_gc.{library_kit}"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)

    @dictify
    def _get_output_files_filter_intervals(self):
        ext = "interval_list"
        token = "{mapper}.gcnv_filter_intervals.{library_kit}"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)

    def _get_output_files_scatter_intervals(self):
        return "work/{token}/out/{token}".format(
            token="{mapper}.gcnv_scatter_intervals.{library_kit}"
        )

    @dictify
    def _get_output_files_coverage(self):
        ext = "tsv"
        token = "{mapper}.gcnv_coverage.{library_name}"
        yield ext, "work/{token}/out/{token}.{ext}".format(token=token, ext=ext)

    @dictify
    def _get_output_files_contig_ploidy(self):
        ext = "done"
        token = "{mapper}.gcnv_contig_ploidy.{library_kit}"
        yield ext, touch("work/{token}/out/{token}/.{ext}".format(token=token, ext=ext))

    @dictify
    def _get_output_files_call_cnvs(self):
        ext = "done"
        token = "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}"
        yield ext, touch("work/{token}/out/{token}/.{ext}".format(token=token, ext=ext))

    @dictify
    def _get_output_files_post_germline_calls(self):
        token = "{mapper}.gcnv_post_germline_calls.{library_name}"
        pairs = {"ratio_tsv": ".ratio.tsv", "itv_vcf": ".interval.vcf.gz", "seg_vcf": ".vcf.gz"}
        for key, ext in pairs.items():
            yield key, touch("work/{token}/out/{token}{ext}".format(token=token, ext=ext))

    @dictify
    def _get_output_files_merge_cohort_vcfs(self):
        token = "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}"
        pairs = {
            "vcf": ".vcf.gz",
            "vcf_md5": ".vcf.gz.md5",
            "tbi": ".vcf.gz.tbi",
            "tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in pairs.items():
            yield key, "work/{token}/out/{token}{ext}".format(token=token, ext=ext)

    @dictify
    def _get_output_files_extract_ped(self):
        token = "{mapper}.gcnv.{library_name}"
        kvs = (
            ("vcf", ".vcf.gz"),
            ("vcf_md5", ".vcf.gz.md5"),
            ("tbi", ".vcf.gz.tbi"),
            ("tbi_md5", ".vcf.gz.tbi.md5"),
        )
        for key, suffix in kvs:
            yield key, "work/{token}/out/{token}{suffix}".format(token=token, suffix=suffix)

    def get_log_file(self, action):
        """Return path to log file"""
        if action in ("preprocess_intervals", "annotate_gc"):
            token = "gcnv_{action}.{{library_kit}}".format(action=action)
            return "work/{token}/log/{token}.log".format(token=token)
        elif action in (
            "filter_intervals",
            "contig_ploidy",
            "scatter_intervals",
            "merge_cohort_vcfs",
        ):
            token = "{{mapper}}.gcnv_{action}.{{library_kit}}".format(action=action)
            return "work/{token}/log/{token}.log".format(token=token)
        elif action == "call_cnvs":
            token = "{{mapper}}.gcnv_{action}.{{library_kit}}.{{shard}}".format(action=action)
            return "work/{token}/log/{token}.log".format(token=token)
        else:
            token = "{{mapper}}.gcnv_{action}.{{library_name}}".format(action=action)
            return "work/{token}/log/{token}.log".format(token=token)

    def update_cluster_config(self, cluster_config):
        """Update cluster configuration for gCNV CNV calling"""
        for action in self.actions:
            if action in ("call_cnvs", "post_germline_calls"):
                cluster_config["targeted_seq_cnv_calling_gcnv_{}".format(action)] = {
                    "mem": 12 * int(3.75 * 1024),
                    "time": "48:00",
                    "ntasks": 16,
                }
            else:
                cluster_config["targeted_seq_cnv_calling_gcnv_{}".format(action)] = {
                    "mem": 2 * int(3.75 * 1024),
                    "time": "04:00",
                    "ntasks": 1,
                }


class TargetedSeqCnvCallingWorkflow(BaseStep):
    """Perform germline targeted sequencing CNV calling"""

    name = "targeted_seq_cnv_calling"
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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((XhmmStepPart, GcnvStepPart, LinkOutStepPart))
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Build mapping from NGS DNA library to library kit.
        self.ngs_library_to_kit = self.sub_steps["xhmm"].ngs_library_to_kit

    @listify
    def _all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the germline targeted sequencing CNV calling workflow.

        If xhmm/path_target_interval_list_mapping is non-empty then we will use this mapping.  In
        this case, only the samples that have a ``libraryKit`` set with a matching entry in the
        mapping. Otherwise, we will create output files for the input primary DNA library of all
        donors.
        """
        # Get list of library kits and donors to use.
        library_kits, donors, kit_counts = self._pick_kits_and_donors()
        # Actually yield the result files.
        token = "{mapper}.{caller}.{index.dna_ngs_library.name}"
        callers = ("xhmm", "gcnv")
        vcf_tools = [t for t in self.config["tools"] if t in callers]
        yield from self._yield_result_files(
            os.path.join("output", token, "out", token + "{ext}"),
            donors,
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=vcf_tools,
            ext=EXT_VALUES,
        )
        if "xhmm" in self.config["tools"]:
            token = "{mapper}.xhmm_genotype.{library_kit}"
            min_kit_usages = 10
            chosen_kits = [kit for kit in library_kits if kit_counts.get(kit, 0) > min_kit_usages]
            chosen_donors = [
                donor
                for donor in donors
                if self.ngs_library_to_kit.get(donor.dna_ngs_library.name) in chosen_kits
            ]
            yield from expand(
                os.path.join("output", token, "out", token + "{ext}"),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=["xhmm"],
                library_kit=chosen_kits,
                ext=EXT_VALUES,
            )
        if "gcnv" in self.config["tools"]:
            token = "{mapper}.gcnv_merge_cohort_vcfs.{library_kit}"
            min_kit_usages = 10
            chosen_kits = [kit for kit in library_kits if kit_counts.get(kit, 0) > min_kit_usages]
            chosen_donors = [
                donor
                for donor in donors
                if self.ngs_library_to_kit.get(donor.dna_ngs_library.name) in chosen_kits
            ]
            yield from expand(
                os.path.join("output", token, "out", token + "{ext}"),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=["gcnv"],
                library_kit=chosen_kits,
                ext=EXT_VALUES,
            )

    def _pick_kits_and_donors(self):
        """Return ``(library_kits, donors)`` with the donors with a matching kit and the kits with a
        matching donor.
        """
        kit_counts = {name: 0 for name in self.ngs_library_to_kit.values()}
        for name in self.ngs_library_to_kit.values():
            kit_counts[name] += 1
        donors = [
            donor
            for donor in self._all_donors()
            if donor.dna_ngs_library and donor.dna_ngs_library.name in self.ngs_library_to_kit
        ]
        return list(sorted(set(self.ngs_library_to_kit.values()))), donors, kit_counts

    def _yield_result_files(self, tpl, donors, **kwargs):
        """Build output paths from path template and extension list.

        Will only yield the result files for pedigrees where the index is in ``donors``.
        """
        donor_names = {donor.name for donor in donors}
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if pedigree.index.name in donor_names:
                    yield from expand(tpl, index=[pedigree.index], **kwargs)

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        self.ensure_w_config(
            ("step_config", "targeted_seq_cnv_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for targeted seq. CNV calling",
        )
