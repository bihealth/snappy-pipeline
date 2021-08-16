# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_sv_calling`` step

The (germline) ``wgs_sv_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned NGS reads) and performs germline SV calling on them.  The result are called SVs in VCF
format.

In contrast to WGS CNV calling, WGS SV calling is able to identify more than just copy number
variation, e.g., inversions.  Large-range CNVs are often better detected by WGS CNV calling.

The WGS SV calling step is mostly followed by WGS SV filtration to reduce the FDR and remove
artifacts.

.. warning::

    Note that only one NGS library is currently supported per bio entity/donor, also if you have
    both Illumina and PacBio data for one patient, for example.  This means that you have to
    create a separate patient with different pk (and secondary id if you only use secondary id for
    file names) for the PacBio data set or create the patient in another project.

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step.

===========
Step Output
===========

For all pedigrees, variant calling will be performed on the primary DNA NGS libraries of all
members, separately for each configured read mapper and variant caller.  The name of the primary
DNA NGS library of the index will be used as an identification token in the output file.  For each
read mapper, variant caller, and pedigree, the following files will be generated:

- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.md5``
- ``{mapper}.{var_caller}.{lib_name}-{lib_pk}.vcf.gz.tbi.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.delly2.P001-N1-DNA1-WGS1-4
    |   `-- out
    |       |-- bwa.delly2.P001-N1-DNA1-WGS1-4.vcf.gz
    |       |-- bwa.delly2.P001-N1-DNA1-WGS1-4.vcf.gz.tbi
    |       |-- bwa.delly2.P001-N1-DNA1-WGS1-4.vcf.gz.md5
    |       `-- bwa.delly2.P001-N1-DNA1-WGS1-4.vcf.gz.tbi.md5
    [...]

Generally, these files will be unfiltered, i.e., contain low-quality variants.

--------------
Delly 2 Output
--------------

The Delly 2 workflow used in this pipeline step incorporates the variants of the whole cohort for
regenotyping.

.. note::

    The file will contain variant records for variants not present in the pedigree at hand. This
    will change in the future and variants not present in the pedigree will be removed.

====================
Global Configuration
====================

- ``static_data_config/reference/path`` must be set appropriately

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_wgs_sv_calling.rst

=============================
Available Germline SV Callers
=============================

The following germline SV callers are currently available

- ``"dna"`` (Illumina)
    - ``"delly2"``
    - ``"manta"``

- ``"dna_long"`` (PacBio)
    - ``"pb_honey_spots"``
    - ``"sniffles"``

=======
Reports
=======

Currently, no reports are generated.
"""

# TODO: remove variants not in pedigree after the final merge step
# TODO: assumption: same platform type
# TODO: only one primary NGS library!
# TODO: only WGS libraries!

from collections import OrderedDict
import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_wrappers.tools.genome_windows import yield_regions

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Available (short) DNA WGS SV callers
DNA_WGS_SV_CALLERS = ("delly2", "manta", "popdel")

#: Available (long) DNA WGS SV callers
LONG_DNA_WGS_SV_CALLERS = ("pb_honey_spots", "sniffles")

#: Default configuration for the wgs_sv_calling step
DEFAULT_CONFIG = r"""
# Default configuration
step_config:
  wgs_sv_calling:
    tools:
      dna:  # short
      - delly2
      dna_long:  # PacBio/Oxford Nanopore
      - sniffles
    path_ngs_mapping: ../ngs_mapping    # REQUIRED
    delly2:
      path_exclude_tsv: null  # optional
      max_threads: 16
      map_qual: 1
      geno_qual: 5
      qual_tra: 20
      mad_cutoff: 9
    manta:
      max_threads: 16
    popdel:
      window_size: 10000000
      max_sv_size: 20000  # == padding
      ignore_chroms:
      - NC_007605  # herpes virus
      - hs37d5     # GRCh37 decoy
      - chrEBV     # Eppstein-Barr Virus
      - '*_decoy'  # decoy contig
      - 'HLA-*'    # HLA genes
      - 'chrUn_*'  # unplaced contigs
    pb_honey_spots:
      num_threads: 16
    sniffles:
      num_threads: 16
"""


class Delly2StepPart(BaseStepPart):
    """WGS SV identification using Delly2

    Delly2 supports the calling based on whole cohorts.  The rough steps are as follows:

    - Perform variant calling on each sample individually ("delly2_call")
    - Merge called variants to get a cohort-wide site list ("delly2_merge_calls")
    - Perform genotyping of the variants in the cohort-wide site list in each sample
      ("delly2_genotype")
    - Merge cohort-wide site list ("delly2_merge_genotypes"); using bcftools
    - Reorder VCF and put pedigree in front; later on, non-pedigree variants should be removed.
    """

    name = "delly2"

    #: Actions in Delly 2 workflow
    actions = ("call", "merge_calls", "genotype", "merge_genotypes", "reorder_vcf")

    #: Directory infixes
    dir_infixes = {
        "call": "{mapper,[^\.]+}.delly2.call.{library_name,[^\.]+}",
        "merge_calls": "{mapper,[^\.]+}.delly2.merge_calls",
        "genotype": "{mapper,[^\.]+}.delly2.genotype.{library_name,[^\.]+}",
        "merge_genotypes": "{mapper,[^\.]+}.delly2.merge_genotypes",
        "reorder_vcf": r"{mapper,[^\.]+}.delly2.{index_ngs_library,[^\.]+}",
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{index_ngs_library}}/out/"
            "{{mapper}}.{var_caller}.{{index_ngs_library}}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)

    def get_library_extra_infos(self, wildcards):
        """Returns library extra infos for the given library name"""
        return self.library_name_to_library[wildcards.library_name].ngs_library.extra_infos

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        assert action in self.actions
        mapping = {
            "call": self._get_input_files_call,
            "merge_calls": self._get_input_files_merge_calls,
            "genotype": self._get_input_files_genotype,
            "merge_genotypes": self._get_input_files_merge_genotypes,
            "reorder_vcf": self._get_input_files_reorder_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_call(self, wildcards):
        """Return input files for "call" action"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @listify
    def _get_input_files_merge_calls(self, wildcards):
        """Return input files for "merge_calls" action"""
        infix = self.dir_infixes["call"]
        infix = infix.replace(r",[^\.]+", "")
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        for donor in self._donors_with_dna_ngs_library():
            yield tpl.format(library_name=donor.dna_ngs_library.name, **wildcards)

    @dictify
    def _get_input_files_genotype(self, wildcards):
        """Return input files for "genotype" action"""
        # Sites VCF file
        infix = self.dir_infixes["merge_calls"]
        infix = infix.replace(r",[^\.]+", "")
        yield "bcf", os.path.join("work", infix, "out", infix + ".bcf").format(**wildcards)
        # BAM files
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @listify
    def _get_input_files_merge_genotypes(self, wildcards):
        """Return input files for "merge_genotypes" action"""
        for donor in self._donors_with_dna_ngs_library():
            infix = self.dir_infixes["genotype"]
            infix = infix.replace(r",[^\.]+", "")
            tpl = os.path.join("work", infix, "out", infix + ".bcf")
            yield tpl.format(library_name=donor.dna_ngs_library.name, **wildcards)

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        """Return input files for "reorder_vcf" action"""
        infix = self.dir_infixes["merge_genotypes"]
        infix = infix.replace(r",[^\.]+", "")
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        yield "bcf", tpl.format(**wildcards)

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action in self.actions
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            infix = self.dir_infixes[action]
            infix2 = infix.replace(r",[^\.]+", "")
            if action != "reorder_vcf":  # generate bcf files internally
                name = name.replace("vcf", "bcf")
                ext = ext.replace("vcf.gz", "bcf")
                name = name.replace("tbi", "csi")
                ext = ext.replace(".tbi", ".csi")
            yield name, "work/" + infix + "/out/" + infix2 + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action in self.actions
        infix = self.dir_infixes[action]
        infix = infix.replace(r",[^\.]+", "")
        return "work/" + infix + "/log/snakemake.log"

    def update_cluster_config(self, cluster_config):
        for action in self.actions:
            if action in ("merge_genotypes", "merge_calls", "reorder_vcf"):  # cheap actions
                cluster_config["wgs_sv_calling_delly2_{}".format(action)] = {
                    "mem": 7 * 1024 * 2,
                    "time": "96:00",
                    "ntasks": 2,
                }
            else:
                cluster_config["wgs_sv_calling_delly2_{}".format(action)] = {
                    "mem": 20 * 1024 * 2,
                    "time": "168:00",
                    "ntasks": 2,
                }


class MantaStepPart(BaseStepPart):
    """WGS SV identification using Manta

    The work flow for Manta is very simple as it allows direct input of a pedigree's input.
    However, this has the drawback of not supporting any background information.
    """

    name = "manta"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.manta.{index_ngs_library}/out/{mapper}.manta.{index_ngs_library}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""

        @listify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to NGS mapping sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected pedigree.  The pedigree is selected
            # by the primary DNA NGS library of the index.
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            for donor in pedigree.donors:
                if donor.dna_ngs_library:
                    for ext in (".bam", ".bam.bai"):
                        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
                        yield ngs_mapping(
                            tpl.format(
                                library_name=donor.dna_ngs_library.name, ext=ext, **wildcards
                            )
                        )

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action == "run"
        infix = "{mapper}.manta.{index_ngs_library}"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action == "run"
        return "work/{mapper}.manta.{index_ngs_library}/log/snakemake.log"

    def update_cluster_config(self, cluster_config):
        cluster_config["wgs_sv_calling_manta_run"] = {
            "mem": int(3.75 * 1024 * 16),
            "time": "40:00",
            "ntasks": 16,
        }


class SvTkStepPart(BaseStepPart):
    """svtk implementation"""

    name = "svtk"

    #: Actions in svtk workflow
    actions = ("standardize",)

    #: Directory infixes
    dir_infixes = {
        "standardize": r"{mapper,[^\.]+}.{caller,[^\.]+}.svtk_standardize.{library_name,[^\.]+}",
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{caller}.svtk_{step}.{{index_ngs_library}}/out/"
            "{{mapper}}.{caller}.svtk_{step}.{{index_ngs_library}}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        assert action in self.actions
        mapping = {
            "standardize": self._get_input_files_standardize,
        }
        return mapping[action]

    @dictify
    def _get_input_files_standardize(self, wildcards):
        """Return input files for "standardize" action"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = (
            "output/{mapper}.{caller_infix}.{library_name}/out/"
            "{mapper}.{caller_infix}.{library_name}{ext}"
        )
        if wildcards.caller == "delly2":
            caller_infix = "delly2.call"
            name_ext = {"calls": ".bcf"}
        else:
            name_ext = {"calls": ".vcf.gz"}
        for name, ext in name_ext.items():
            yield name, tpl.format(ext=ext, caller_infix=caller_infix, **wildcards)

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action in self.actions
        infix = self.dir_infixes[action]
        infix2 = infix.replace(r",[^\.]+", "")
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            infix = self.dir_infixes[action].replace(r",[^\.]+", "")
            yield name, "work/" + infix + "/out/" + infix2 + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action in self.actions
        infix = self.dir_infixes[action]
        infix2 = infix.replace(r",[^\.]+", "")
        return "work/%s/log/%s.log" % (infix, infix2)

    def update_cluster_config(self, cluster_config):
        for action in self.actions:
            cluster_config["wgs_sv_svtk_%s" % action] = {
                "mem": 4 * 1024,
                "time": "4:00",
                "ntasks": 1,
            }


class PopDelStepPart(BaseStepPart):
    """WGS SV identification using PopDel.

    Implemented using chromosome-wise calling.
    """

    name = "popdel"

    #: Actions in PopDel workflow
    actions = ("profile", "call", "concat_calls", "reorder_vcf")

    #: Directory infixes
    dir_infixes = {
        "profile": "{mapper}.popdel.internal.profile.{index_ngs_library}",
        "call": "{mapper}.popdel.internal.call.{chrom}-{begin}-{end}",
        "concat_calls": "{mapper}.popdel.internal.concat_calls",
        "reorder_vcf": r"{mapper}.popdel.{index_ngs_library}",
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{var_caller}.{{index_ngs_library}}/out/"
            "{{mapper}}.{var_caller}.{{index_ngs_library}}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)

    def get_library_extra_infos(self, wildcards):
        """Returns library extra infos for the given library name"""
        return self.library_name_to_library[wildcards.library_name].ngs_library.extra_infos

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        assert action in self.actions
        mapping = {
            "profile": self._get_input_files_profile,
            "call": self._get_input_files_call,
            "concat_calls": self._get_input_files_concat_calls,
            "reorder_vcf": self._get_input_files_reorder_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_profile(self, wildcards):
        """Return input files for "call" action"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{index_ngs_library}/out/{mapper}.{index_ngs_library}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    def _get_input_files_call(self, wildcards):
        """Return input files for "call" action"""
        tpl = os.path.join(
            "work", self.dir_infixes["profile"], "out", self.dir_infixes["profile"] + ".profile"
        )
        result = {"profile": []}
        for donor in self._donors_with_dna_ngs_library():
            result["profile"].append(
                tpl.format(index_ngs_library=donor.dna_ngs_library.name, **wildcards)
            )
        return result

    def _get_input_files_concat_calls(self, wildcards):
        window_size = self.config["popdel"]["window_size"]
        padding = self.config["popdel"]["max_sv_size"]
        window_size
        """Return input files for "concat_calls" action"""
        tpl = os.path.join(
            "work", self.dir_infixes["call"], "out", self.dir_infixes["call"] + ".vcf.gz"
        )
        result = {"vcf": []}
        with open(self.get_fai_path(), "rt") as fai_file:
            for region in yield_regions(
                fai_file,
                window_size=window_size,
                padding=padding,
                ignore_chroms=self.get_ignore_chroms(),
            ):
                if region.begin == 0:
                    region.begin = 1
                result["vcf"].append(
                    tpl.format(chrom=region.chrom, begin=region.begin, end=region.end, **wildcards)
                )
        return result

    def get_fai_path(self):
        return self.w_config["static_data_config"]["reference"]["path"] + ".fai"

    def get_ignore_chroms(self):
        return self.config["popdel"]["ignore_chroms"]

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        """Return input files for "reorder_vcf" action"""
        infix = self.dir_infixes["concat_calls"].format(**wildcards)
        yield "vcf", "work/" + infix + "/out/" + infix + ".vcf.gz"

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action in self.actions
        infix = self.dir_infixes[action]
        if action == "profile":
            yield "profile", "work/" + infix + "/out/" + infix + ".profile"
            yield "profile_md5", "work/" + infix + "/out/" + infix + ".profile.md5"
        else:
            for name, ext in zip(EXT_NAMES, EXT_VALUES):
                yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action in self.actions
        infix = self.dir_infixes[action]
        return "work/" + infix + "/log/" + infix + ".log"

    def update_cluster_config(self, cluster_config):
        for action in self.actions:
            cluster_config["wgs_sv_calling_popdel_{}".format(action)] = {
                "mem": int(12 * 1024 * 2),
                "time": "96:00",
                "ntasks": 2,
            }


class PbHoneySpotsStepPart(BaseStepPart):
    """WGS SV identification using PB Honey Spots"""

    name = "pb_honey_spots"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.pb_honey_spots.{index_ngs_library}/out/{mapper}."
            "manta.{index_ngs_library}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""

        @listify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to NGS mapping sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected pedigree.  The pedigree is selected
            # by the primary DNA NGS library of the index.
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            for donor in pedigree.donors:
                for ext in (".bam", ".bam.bai"):
                    tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
                    yield ngs_mapping(
                        tpl.format(library_name=donor.dna_ngs_library.name, ext=ext, **wildcards)
                    )

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action == "run"
        infix = "{mapper}.pb_honey_spots.{index_ngs_library}"
        yield "bed", "work/" + infix + "/out/" + infix + ".bed"
        yield "bed_md5", "work/" + infix + "/out/" + infix + ".bed.md5"

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action == "run"
        return "work/{mapper}.pb_honey_spots.{index_ngs_library}/log/snakemake.log"

    def update_cluster_config(self, cluster_config):
        cluster_config["wgs_sv_calling_pb_honey_spots_run"] = {
            "mem": int(3.75 * 1024 * 16),
            "time": "40:00",
            "ntasks": 16,
        }


class SnifflesStepPart(BaseStepPart):
    """WGS SV identification using Sniffles"""

    name = "sniffles"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.sniffles.{index_ngs_library}/out/{mapper}."
            "sniffles.{index_ngs_library}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""

        @listify
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to NGS mapping sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected pedigree.  The pedigree is selected
            # by the primary DNA NGS library of the index.
            pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
            for donor in pedigree.donors:
                for ext in (".bam", ".bam.bai"):
                    tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
                    yield ngs_mapping(
                        tpl.format(library_name=donor.dna_ngs_library.name, ext=ext, **wildcards)
                    )

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        assert action == "run"
        infix = "{mapper}.sniffles.{index_ngs_library}"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        assert action == "run"
        return "work/{mapper}.sniffles.{index_ngs_library}/log/snakemake.log"

    def update_cluster_config(self, cluster_config):
        cluster_config["wgs_sv_calling_sniffles_run"] = {
            "mem": int(3.75 * 1024 * 16),
            "time": "40:00",
            "ntasks": 16,
        }


class WgsSvCallingWorkflow(BaseStep):
    """Perform (germline) WGS SV calling"""

    name = "wgs_sv_calling"
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
        self.register_sub_step_classes(
            (
                Delly2StepPart,
                MantaStepPart,
                SvTkStepPart,
                PopDelStepPart,
                PbHoneySpotsStepPart,
                SnifflesStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_pattern = "{mapper}.{caller}.{index_library.name}"
        # Illumina DNA WGS SV calling
        yield from self._yield_result_files_dna_short(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            caller=self.config["tools"]["dna"],
            ext=EXT_VALUES,
        )
        # Long Read DNA WGS SV calling
        bed_tools = set(["pb_honey_spots"])
        yield from self._yield_result_files_dna_long(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna_long"],
            caller=set(self.config["tools"]["dna_long"]) - bed_tools,
            ext=EXT_VALUES,
        )
        # Long Read WGS SV Calling (BED output)
        yield from self._yield_result_files_dna_long(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna_long"],
            caller=set(self.config["tools"]["dna_long"]) & bed_tools,
            ext=(".bed", ".bed.md5"),
        )

    def _yield_result_files_dna_short(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                ngs_library = pedigree.index.dna_ngs_library
                if not ngs_library:
                    msg = "WARNING: index of pedigree has no NGS library (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                if (
                    ngs_library.extra_infos["libraryType"] != "WGS"
                    or ngs_library.extra_infos["seqPlatform"] != "Illumina"
                ):
                    continue  # not WGS or no long read
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def _yield_result_files_dna_long(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                ngs_library = pedigree.index.dna_ngs_library
                if not ngs_library:
                    msg = "WARNING: index of pedigree has no NGS library (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                if (
                    ngs_library.extra_infos["libraryType"] != "WGS"
                    or ngs_library.extra_infos["seqPlatform"] != "PacBio"
                ):
                    continue  # not WGS or no long read
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "wgs_sv_calling", "path_ngs_mapping"),
            ("Path to NGS mapping not configured but required for variant calling"),
        )
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            ("Path to reference FASTA not configured but required for variant calling"),
        )
        # Check that only valid tools are selected
        selected = set(self.w_config["step_config"]["wgs_sv_calling"]["tools"]["dna"])
        invalid = selected - set(DNA_WGS_SV_CALLERS)
        if invalid:
            raise Exception(
                "Invalid short-read WGS SV caller selected: {}".format(list(sorted(invalid)))
            )
        selected = set(self.w_config["step_config"]["wgs_sv_calling"]["tools"]["dna_long"])
        invalid = selected - set(LONG_DNA_WGS_SV_CALLERS)
        if invalid:
            raise Exception(
                "Invalid long-read WGS SV caller selected: {}".format(list(sorted(invalid)))
            )
