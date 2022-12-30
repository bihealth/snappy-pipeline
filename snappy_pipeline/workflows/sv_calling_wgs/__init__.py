"""Implementation of the ``sv_calling_wgs`` step
"""

from collections import OrderedDict
import itertools
import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand, touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_wrappers.tools.genome_windows import yield_regions

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available (short) DNA WGS SV callers
DNA_WGS_SV_CALLERS = ("delly2", "manta", "popdel", "melt")

#: Available (long) DNA WGS SV callers
LONG_DNA_WGS_SV_CALLERS = ("pb_honey_spots", "sniffles", "sniffles2")

#: Default configuration for the sv_calling_wgs step
DEFAULT_CONFIG = r"""
# Default configuration
step_config:
  sv_calling_wgs:
    tools:
      dna: [delly2] # Required if short-read mapper used; otherwise, leave empty. Example: 'delly2'.
      dna_long: []  # Required if long-read mapper used (PacBio/Oxford Nanopore); otherwise, leave empty. Example: 'sniffles'.

    path_ngs_mapping: ../ngs_mapping    # REQUIRED

    # Short-read SV calling tool configuration
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
    gcnv:
      # Path to gCNV model - will execute analysis in CASE MODE.
      #
      # Example of precomputed model:
      # - library: "Agilent SureSelect Human All Exon V6"  # Library name
      #   contig_ploidy: /path/to/ploidy-model         # Output from `DetermineGermlineContigPloidy`
      #   model_pattern: /path/to/model_*              # Output from `GermlineCNVCaller`
      precomputed_model_paths: []  # REQUIRED

      # Path to BED file with uniquely mappable regions.
      path_uniquely_mapable_bed: null  # REQUIRED
    melt:
      me_refs_infix: 1KGP_Hg19
      me_types:
      - ALU
      - LINE1
      - SVA
      genes_file: add_bed_files/1KGP_Hg19/hg19.genes.bed  # adjust, e.g., Hg38/Hg38.genes.bed

    # Long-read SV calling tool configuration
    pb_honey_spots:
      num_threads: 16
    sniffles:
      num_threads: 16
    sniffles2:
      num_threads: 16
      tandem_repeats: /fast/groups/cubi/work/projects/biotools/sniffles2/trf/GRCh37/human_hs37d5.trf.bed  # REQUIRED

    # Common configuration
    ignore_chroms:
    - NC_007605  # herpes virus
    - hs37d5     # GRCh37 decoy
    - chrEBV     # Eppstein-Barr Virus
    - '*_decoy'  # decoy contig
    - 'HLA-*'    # HLA genes
    - 'chrUn_*'  # unplaced contigs
"""


class Delly2StepPart(BaseStepPart):
    """WGS SV identification using Delly2

    Delly2 supports the calling based on pedigrees.  The rough steps are as follows:

    - Perform variant calling on each sample individually ("delly2_call")
    - Merge called variants to get a cohort-wide site list ("delly2_merge_calls")
    - Perform genotyping of the variants in the cohort-wide site list in each sample
      ("delly2_genotype")
    - Merge cohort-wide site list ("delly2_merge_genotypes"); using bcftools
    - Reorder VCF and put pedigree in front; later on, non-pedigree variants should be removed.
    """

    #: Step name
    name = "delly2"

    #: Actions in Delly 2 workflow
    actions = ("call", "merge_calls", "genotype", "merge_genotypes")

    #: Directory infixes
    dir_infixes = {
        "call": r"{mapper,[^\.]+}.delly2.call.{library_name,[^\.]+}",
        "merge_calls": r"{mapper,[^\.]+}.delly2.merge_calls.{index_ngs_library,[^\.]+}",
        "genotype": r"{mapper,[^\.]+}.delly2.genotype.{library_name,[^\.]+}",
        "merge_genotypes": r"{mapper,[^\.]+}.delly2.{index_ngs_library,[^\.]+}",
    }

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "cheap_action": ResourceUsage(
            threads=2,
            time="4-00:00:00",  # 4 days
            memory=f"{7 * 1024 * 2}M",
        ),
        "default": ResourceUsage(
            threads=2,
            time="7-00:00:00",  # 7 days
            memory=f"{20 * 1024 * 2}M",
        ),
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
        # Build shortcut from donor library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)

    def get_library_extra_infos(self, wildcards):
        """Returns library extra infos for the given library name"""
        return self.library_name_to_library[wildcards.library_name].ngs_library.extra_infos

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "call": self._get_input_files_call,
            "merge_calls": self._get_input_files_merge_calls,
            "genotype": self._get_input_files_genotype,
            "merge_genotypes": self._get_input_files_merge_genotypes,
        }
        return mapping[action]

    @dictify
    def _get_input_files_call(self, wildcards):
        """Return input files for "call" action"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    def _get_input_files_merge_calls_or_genotype(self, wildcards, step):
        assert step in ("call", "genotype")
        # Create path template to per-sample call/genotype BCF
        infix = self.dir_infixes[step]
        infix = infix.replace(r",[^\.]+", "")
        tpl = os.path.join("work", infix, "out", infix + ".bcf")
        # Yield paths to pedigree's per-sample call BCF files
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                yield tpl.format(library_name=donor.dna_ngs_library.name, **wildcards)

    @listify
    def _get_input_files_merge_calls(self, wildcards):
        """Return input files for "merge_calls" action"""
        yield from self._get_input_files_merge_calls_or_genotype(wildcards, "call")

    @dictify
    def _get_input_files_genotype(self, wildcards):
        """Return input files for "genotype" action"""
        pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
        # Per-pedigree site BCF file
        infix = self.dir_infixes["merge_calls"]
        infix = infix.replace(r",[^\.]+", "")
        infix = infix.format(
            mapper=wildcards.mapper, index_ngs_library=pedigree.index.dna_ngs_library.name
        )
        yield "bcf", os.path.join("work", infix, "out", infix + ".bcf").format(**wildcards)
        # BAM files
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @listify
    def _get_input_files_merge_genotypes(self, wildcards):
        """Return input files for "merge_genotypes" action"""
        yield from self._get_input_files_merge_calls_or_genotype(wildcards, "genotype")

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
        # Validate action
        self._validate_action(action)
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            infix = self.dir_infixes[action]
            infix2 = infix.replace(r",[^\.]+", "")
            if action != "merge_genotypes":  # generate bcf files internally
                name = name.replace("vcf", "bcf")
                ext = ext.replace("vcf.gz", "bcf")
                name = name.replace("tbi", "csi")
                ext = ext.replace(".tbi", ".csi")
            yield name, "work/" + infix + "/out/" + infix2 + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        # Validate action
        self._validate_action(action)
        infix = self.dir_infixes[action]
        infix = infix.replace(r",[^\.]+", "")
        return "work/" + infix + "/log/snakemake.log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        if action in ("merge_genotypes", "merge_calls"):  # cheap actions
            return self.resource_usage_dict.get("cheap_action")
        else:
            return self.resource_usage_dict.get("default")


class MantaStepPart(BaseStepPart):
    """WGS SV identification using Manta

    The work flow for Manta is very simple as it allows direct input of a pedigree's input.
    However, this has the drawback of not supporting any background information.
    """

    #: Step name
    name = "manta"

    #: Class available actions
    actions = ("run",)

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
        # Validate action
        self._validate_action(action)

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

        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        # Validate action
        self._validate_action(action)
        infix = "{mapper}.manta.{index_ngs_library}"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        # Validate action
        self._validate_action(action)
        return "work/{mapper}.manta.{index_ngs_library}/log/snakemake.log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=16,
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{int(3.75 * 1024 * 16)}M",
        )


class GcnvWgsStepPart(RunGcnvStepPart):
    """WGS CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        # No mapping given as WGS, we will use the "default" one for all.
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library:
                yield donor.dna_ngs_library.name, "default"


class PopDelStepPart(BaseStepPart):
    """WGS SV identification using PopDel.

    Implemented using chromosome-wise calling.
    """

    #: Step name
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
        # Validate action
        self._validate_action(action)
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
        return self.config["ignore_chroms"]

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
        # Validate action
        self._validate_action(action)
        infix = self.dir_infixes[action]
        if action == "profile":
            yield "profile", "work/" + infix + "/out/" + infix + ".profile"
            yield "profile_md5", "work/" + infix + "/out/" + infix + ".profile.md5"
        else:
            for name, ext in zip(EXT_NAMES, EXT_VALUES):
                yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        # Validate action
        self._validate_action(action)
        infix = self.dir_infixes[action]
        return "work/" + infix + "/log/" + infix + ".log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="4-00:00:00",  # 4 days
            memory=f"{12 * 1024 * 2}M",
        )


class MeltStepPart(BaseStepPart):
    """MEI calling using Melt

    We perform cohort-wide genotyping of the mobile elements.  This requires a somewhat complicated
    workflow.  The methods input/output/etc. delegate to appropriate privated functions, thus
    this class is somewhat lengthy.  In the, however, it's not that complex.
    """

    #: Step name
    name = "melt"

    #: Class available actions
    actions = (
        "preprocess",
        "indiv_analysis",
        "group_analysis",
        "genotype",
        "make_vcf",
        "merge_vcf",
        "reorder_vcf",
    )

    def __init__(self, parent):
        super().__init__(parent)
        #: All individual's primary NGS libraries
        self.all_dna_ngs_libraries = []
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    self.all_dna_ngs_libraries.append(donor.dna_ngs_library.name)
        #: Linking NGS libraries to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        mapping = {
            "preprocess": self._get_input_files_preprocess,
            "indiv_analysis": self._get_input_files_indiv_analysis,
            "group_analysis": self._get_input_files_group_analysis,
            "genotype": self._get_input_files_genotype,
            "make_vcf": self._get_input_files_make_vcf,
            "merge_vcf": self._get_input_files_merge_vcf,
            "reorder_vcf": self._get_input_files_reorder_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_preprocess(self, wildcards):
        # Get shorcut to NGS mapping sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_input_files_indiv_analysis(self, wildcards):
        tpl = "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}{ext}"
        yield "orig_bam", tpl.format(ext=".bam", **wildcards)
        yield "orig_bai", tpl.format(ext=".bam.bai", **wildcards)
        yield "disc_bam", tpl.format(ext=".bam.disc", **wildcards)
        yield "disc_bai", tpl.format(ext=".bam.disc.bai", **wildcards)

    @listify
    def _get_input_files_group_analysis(self, wildcards):
        for library_name in self.all_dna_ngs_libraries:
            yield "work/{mapper}.melt.indiv_analysis.{me_type}/out/.done.{library_name}".format(
                library_name=library_name, **wildcards
            )

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield "done", "work/{mapper}.melt.group_analysis.{me_type}/out/.done".format(**wildcards)
        yield "bam", "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}.bam".format(
            **wildcards
        )

    @listify
    def _get_input_files_make_vcf(self, wildcards):
        yield "work/{mapper}.melt.group_analysis.{me_type}/out/.done".format(**wildcards)
        # Probably, best create a create-list intermediate target
        for library_name in self.all_dna_ngs_libraries:
            yield "work/{mapper}.melt.genotype.{me_type}/out/.done.{library_name}".format(
                library_name=library_name, **wildcards
            )

    @listify
    def _get_input_files_merge_vcf(self, wildcards):
        for me_type in self.config["melt"]["me_types"]:
            yield "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz".format(
                me_type=me_type, **wildcards
            )

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        yield "vcf", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz".format(
            **wildcards
        )
        yield "vcf_tbi", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi".format(
            **wildcards
        )

    def get_ped_members(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        mapping = {
            "preprocess": self._get_output_files_preprocess,
            "indiv_analysis": self._get_output_files_indiv_analysis,
            "group_analysis": self._get_output_files_group_analysis,
            "genotype": self._get_output_files_genotype,
            "make_vcf": self._get_output_files_make_vcf,
            "merge_vcf": self._get_output_files_merge_vcf,
            "reorder_vcf": self._get_output_files_reorder_vcf,
        }
        return mapping[action]()

    @dictify
    def _get_output_files_preprocess(self):
        # Note that mapper is not part of the output BAM file as MELT infers sample file from BAM
        # file name instead of using sample name from BAM header.
        tpl = "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}%s"
        yield "orig_bam", tpl % ".bam"
        yield "orig_bai", tpl % ".bam.bai"
        yield "disc_bam", tpl % ".bam.disc"
        yield "disc_bai", tpl % ".bam.disc.bai"
        yield "disc_fq", tpl % ".bam.fq"

    @dictify
    def _get_output_files_indiv_analysis(self):
        yield "done", touch("work/{mapper}.melt.indiv_analysis.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_output_files_group_analysis(self):
        yield "done", touch("work/{mapper}.melt.group_analysis.{me_type}/out/.done")

    @dictify
    def _get_output_files_genotype(self):
        yield "done", touch("work/{mapper}.melt.genotype.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_output_files_make_vcf(self):
        yield "list_txt", "work/{mapper}.melt.genotype.{me_type}/out/list.txt"
        yield "done", touch("work/{mapper}.melt.make_vcf.{me_type}/out/.done")
        yield "vcf", "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz"
        yield "vcf_tbi", "work/{mapper}.melt.merge_vcf.{me_type}/out/{me_type}.final_comp.vcf.gz.tbi"

    @dictify
    def _get_output_files_merge_vcf(self):
        yield "vcf", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz"
        yield "vcf_tbi", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi"
        yield "vcf_md5", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.md5"
        yield "vcf_tbi_md5", "work/{mapper}.melt.merge_vcf/out/{mapper}.melt.merge_vcf.vcf.gz.tbi.md5"

    @dictify
    def _get_output_files_reorder_vcf(self):
        tpl = "work/{mapper}.melt.{index_library_name}/out/{mapper}.melt.{index_library_name}.%s%s"
        key_ext = {"vcf": "vcf.gz", "vcf_tbi": "vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, tpl % (ext, "")
            yield key + "_md5", tpl % (ext, ".md5")

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        if action == "preprocess":
            return "work/{mapper}.melt.preprocess.{library_name}/log/snakemake.wgs_mei_calling.log"
        elif action in ("indiv_analysis", "genotype"):
            return (
                "work/{{mapper}}.melt.{action}.{{me_type}}/log/"
                "snakemake.wgs_mei_calling.{{library_name}}.log"
            ).format(action=action)
        elif action in ("indiv_analysis", "group_analysis", "genotype", "make_vcf"):
            return (
                "work/{{mapper}}.melt.{action}.{{me_type}}/log/" "snakemake.wgs_mei_calling.log"
            ).format(action=action)
        elif action == "merge_vcf":
            return "work/{mapper}.melt.merge_vcf/log/snakemake.wgs_mei_calling.log"
        elif action == "reorder_vcf":
            return (
                "work/{mapper}.melt.reorder_vcf.{index_library_name}/log/"
                "snakemake.wgs_mei_calling.log"
            )
        return None

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=6,
            time="5-07:00:00",  # ~5.3 days
            memory=f"{int(3.75 * 1024 * 6)}M",
        )


class PbHoneySpotsStepPart(BaseStepPart):
    """WGS SV identification using PB Honey Spots"""

    #: Step name
    name = "pb_honey_spots"

    #: Class available actions
    actions = ("run",)

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
        # Validate action
        self._validate_action(action)

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

        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        # Validate action
        self._validate_action(action)
        infix = "{mapper}.pb_honey_spots.{index_ngs_library}"
        yield "bed", "work/" + infix + "/out/" + infix + ".bed"
        yield "bed_md5", "work/" + infix + "/out/" + infix + ".bed.md5"

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        # Validate action
        self._validate_action(action)
        return "work/{mapper}.pb_honey_spots.{index_ngs_library}/log/snakemake.log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["pb_honey_spots"]["num_threads"],
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{int(3.75 * 1024 * 16)}M",
        )


class SnifflesStepPart(BaseStepPart):
    """WGS SV identification using Sniffles"""

    #: Step name
    name = "sniffles"

    #: Class available actions
    actions = ("run",)

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
        # Validate action
        self._validate_action(action)

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

        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        # Validate action
        self._validate_action(action)
        infix = "{mapper}.sniffles.{index_ngs_library}"
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, "work/" + infix + "/out/" + infix + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        # Validate action
        self._validate_action(action)
        return "work/{mapper}.sniffles.{index_ngs_library}/log/snakemake.log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["sniffles"]["num_threads"],
            time="1-16:00:00",  # 1 day and 16 hours
            memory=f"{int(3.75 * 1024 * 16)}M",
        )


class Sniffles2StepPart(BaseStepPart):
    """WGS SV identification using Sniffles 2"""

    #: Step name
    name = "sniffles2"

    #: Class available actions
    actions = ("bam_to_snf", "snf_to_vcf")

    #: Directory infixes
    dir_infixes = {
        "bam_to_snf": r"{mapper,[^\.]+}.sniffles2.bam_to_snf.{library_name,[^\.]+}",
        "snf_to_vcf": r"{mapper,[^\.]+}.sniffles2.{index_ngs_library,[^\.]+}",
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.sniffles2.{index_ngs_library}/out/"
            "{mapper}.sniffles2.{index_ngs_library}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from donor library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)

    def get_input_files(self, action):
        """Return appropriate input function for the given action"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "bam_to_snf": self._get_input_files_bam_to_snf,
            "snf_to_vcf": self._get_input_files_snf_to_vcf,
        }
        return mapping[action]

    @dictify
    def _get_input_files_bam_to_snf(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for name, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield name, ngs_mapping(tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_input_files_snf_to_vcf(self, wildcards):
        # Create path template to per-sample call/genotype BCF
        infix = self.dir_infixes["bam_to_snf"]
        infix = infix.replace(r",[^\.]+", "")
        tpl = os.path.join("work", infix, "out", infix + ".snf")
        # Yield paths to pedigree's per-sample call BCF files
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        yield "snf", [
            tpl.format(library_name=donor.dna_ngs_library.name, **wildcards)
            for donor in pedigree.donors
            if donor.dna_ngs_library
        ]

    @dictify
    def get_output_files(self, action):
        """Return output paths for the given action; include wildcards"""
        # Validate action
        self._validate_action(action)
        ext_names = EXT_NAMES
        ext_values = EXT_VALUES
        if action == "bam_to_snf":
            ext_names = list(itertools.chain(ext_names, ["snf", "snf_md5"]))
            ext_values = list(itertools.chain(ext_values, [".snf", ".snf.md5"]))
        for name, ext in zip(ext_names, ext_values):
            infix = self.dir_infixes[action]
            infix2 = infix.replace(r",[^\.]+", "")
            yield name, "work/" + infix + "/out/" + infix2 + ext

    def get_log_file(self, action):
        """Return log file path for the given action; includes wildcards"""
        # Validate action
        self._validate_action(action)
        infix = self.dir_infixes[action]
        infix = infix.replace(r",[^\.]+", "")
        return "work/" + infix + "/log/snakemake.log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config["sniffles2"]["num_threads"], time="0-02:00:00", memory="4G"
        )


class SvCallingWgsWorkflow(BaseStep):
    """Perform (germline) WGS SV calling"""

    name = "sv_calling_wgs"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
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
                Delly2StepPart,
                MantaStepPart,
                PopDelStepPart,
                GcnvWgsStepPart,
                PbHoneySpotsStepPart,
                SnifflesStepPart,
                Sniffles2StepPart,
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
        bed_tools = {"pb_honey_spots"}
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
                seq_platform = ngs_library.extra_infos["seqPlatform"].lower()
                library_type = ngs_library.extra_infos["libraryType"].lower()
                if not ngs_library:
                    msg = "WARNING: index of pedigree has no NGS library (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                if library_type != "wgs" or seq_platform != "illumina":
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
            ("step_config", "sv_calling_wgs", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for variant calling",
        )
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant calling",
        )
        # Check that only valid tools are selected
        selected = set(self.w_config["step_config"]["sv_calling_wgs"]["tools"]["dna"])
        invalid = selected - set(DNA_WGS_SV_CALLERS)
        if invalid:
            raise Exception(
                "Invalid short-read WGS SV caller selected: {}".format(list(sorted(invalid)))
            )
        selected = set(self.w_config["step_config"]["sv_calling_wgs"]["tools"]["dna_long"])
        invalid = selected - set(LONG_DNA_WGS_SV_CALLERS)
        if invalid:
            raise Exception(
                "Invalid long-read WGS SV caller selected: {}".format(list(sorted(invalid)))
            )
