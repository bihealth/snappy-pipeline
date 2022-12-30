"""Implementation of the ``sv_calling_wgs`` step
"""

from collections import OrderedDict

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import (
    ForwardResourceUsageMixin,
    ForwardSnakemakeFilesMixin,
)
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.sv_calling import (
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_wrappers.tools.genome_windows import yield_regions

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Available (short) DNA WGS SV callers
DNA_WGS_SV_CALLERS = ("delly2", "manta", "popdel", "melt", "gcnv")

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
      jar_file: REQUIRED
      genes_file: add_bed_files/1KGP_Hg19/hg19.genes.bed  # adjust, e.g., Hg38/Hg38.genes.bed

    # Long-read SV calling tool configuration
    sniffles2:
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


class GcnvWgsStepPart(RunGcnvStepPart):
    """WGS CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        # No mapping given for WGS; we will use the "wgs" one for all.
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library:
                yield donor.dna_ngs_library.name, "wgs"


class PopDelStepPart(
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    BaseStepPart,
):
    """WGS SV identification using PopDel.

    PopDel requires whole-cohort wide calling.  We implement cohort-wise calling to speed
    things up a bit more.
    """

    name = "popdel"
    actions = ("profile", "call", "concat_calls", "reorder_vcf")

    _resource_usage = ResourceUsage(
        threads=2,
        time="4-00:00:00",
        memory="24G",
    )
    resource_usage_dict = {
        "profile": _resource_usage,
        "call": _resource_usage,
        "concat_calls": _resource_usage,
        "reorder_vcf": _resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)
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

    @dictify
    def _get_input_files_profile(self, wildcards):
        """Return input files for "call" action"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.index_ngs_library}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")

    @dictify
    def _get_output_files_profile(self):
        infix = "{mapper}.popdel_profile.{index_library_name}"
        yield "profile", f"work/{infix}/out/{infix}.profile"
        yield "profile_md5", f"work/{infix}/out/{infix}.profile.md5"

    @dictify
    def _get_input_files_call(self, wildcards):
        paths = []
        for donor in self._donors_with_dna_ngs_library():
            infix = f"{wildcards.mapper}.popdel_profile.{donor.dna_ngs_library.name}"
            paths.append(f"work/{infix}/out/{infix}.profile")
        yield "profile", paths

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    @dictify
    def _get_output_files_call(self):
        infix = "{mapper}.popdel_call.{chrom}-{begin}-{end}"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"

    @dictify
    def _get_input_files_concat_calls(self, wildcards):
        window_size = self.config["popdel"]["window_size"]
        padding = self.config["popdel"]["max_sv_size"]
        vcfs = {}
        with open(self._get_fai_path(), "rt") as fai_file:
            for r in yield_regions(
                fai_file,
                window_size=window_size,
                padding=padding,
                ignore_chroms=self._get_ignore_chroms(),
            ):
                if r.begin == 0:
                    r.begin = 1
                infix = f"{wildcards.mapper}.popdel.call.{r.chrom}-{r.begin}-{r.end}"
                vcfs.append(f"work/{infix}/out/{infix}.vcf.gz")
        yield "vcf", vcfs

    def _get_fai_path(self):
        return self.w_config["static_data_config"]["reference"]["path"] + ".fai"

    def _get_ignore_chroms(self):
        return self.config["ignore_chroms"]

    @dictify
    def _get_output_files_concat_calls(self):
        infix = "{mapper}.popdel_concat_calls"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        infix = f"{wildcards.mapper}.popdel.internal.concat_calls"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"

    @dictify
    def _get_output_files_reorder_vcf(self):
        infix = "{mapper}.popdel.{index_ngs_library}"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"

    def get_ped_members(self, wildcards):
        """Used in Snakefile to rule ``sv_calling_wgs_popdel_reorder_vcf``"""
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )


class MeltStepPart(
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    BaseStepPart,
):
    """MEI calling using MELT

    We implement the workflow as per-pedigree calling.  Generally, this leads to consistent
    positions within each pedigree but not necessarily across the whole cohort.

    Note that MELT is not free software, so further setup is needed.
    """

    name = "melt"
    actions = (
        "preprocess",
        "indiv_analysis",
        "group_analysis",
        "genotype",
        "make_vcf",
        "merge_vcf",
        "reorder_vcf",
    )

    _resource_usage = ResourceUsage(
        threads=6,
        time="1-00:00:00",
        memory=f"{int(3.75 * 1024 * 6)}M",
    )
    resource_usage_dict = {
        "preprocess": _resource_usage,
        "indiv_analysis": _resource_usage,
        "group_analysis": _resource_usage,
        "genotype": _resource_usage,
        "make_vcf": _resource_usage,
        "merge_vcf": _resource_usage,
        "reorder_vcf": _resource_usage,
    }

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

    @dictify
    def _get_input_files_preprocess(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")

    @dictify
    def _get_output_files_preprocess(self):
        # Note that mapper is not part of the output BAM file as MELT infers sample file from BAM
        # file name instead of using sample name from BAM header.
        prefix = "work/{mapper}.melt.preprocess.{library_name}/out/{library_name}"
        yield "orig_bam", f"{prefix}.bam"
        yield "orig_bai", f"{prefix}.bam.bai"
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "disc_bai", f"{prefix}.bam.disc.bai"
        yield "disc_fq", f"{prefix}.bam.fq"

    @dictify
    def _get_input_files_indiv_analysis(self, wildcards):
        infix = f"{wildcards.mapper}.melt.preprocess.{wildcards.library_name}"
        yield "orig_bam", f"work/{infix}/out/{wildcards.library_name}.bam"
        yield "disc_bam", f"work/{infix}/out/{wildcards.library_name}.bam.disc"

    @dictify
    def _get_output_files_indiv_analysis(self):
        infix = "{mapper}.melt.indiv_analysis.{me_type}"
        yield "done", touch(f"work/{infix}/out/.done.{{library_name}}")

    @listify
    def _get_input_files_group_analysis(self, wildcards):
        for library_name in self.all_dna_ngs_libraries:
            infix = f"{wildcards.mapper}.melt.indiv_analysis.{wildcards.me_type}"
            yield f"work/{infix}/out/.done.{library_name}"

    @dictify
    def _get_output_files_group_analysis(self):
        yield "done", touch("work/{mapper}.melt.group_analysis.{me_type}/out/.done")

    @dictify
    def _get_input_files_genotype(self, wildcards):
        infix_done = f"{wildcards.mapper}.melt.group_analysis.{wildcards.me_type}"
        yield "done", f"work/{infix_done}/out/.done".format(**wildcards)
        infix_bam = f"{wildcards.mapper}.melt.preprocess.{wildcards.library_name}"
        yield "bam", f"work/{infix_bam}/out/{wildcards.library_name}.bam"

    @dictify
    def _get_output_files_genotype(self):
        yield "done", touch("work/{mapper}.melt.genotype.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_input_files_make_vcf(self, wildcards):
        infix = f"{wildcards.mapper}.melt.group_analysis.{wildcards.me_type}"
        yield "group_analysis", f"work/{infix}/out/.done"
        paths = []
        for library_name in self.all_dna_ngs_libraries:
            infix = f"{wildcards.mapper}.melt.genotype.{wildcards.me_type}"
            yield f"work/{infix}/out/.done.{library_name}"
        yield "genotype", paths

    @dictify
    def _get_output_files_make_vcf(self):
        infix = "{mapper}.melt.genotype.{me_type}"
        yield "list_txt", f"work/{infix}/out/list.txt"
        yield "done", touch(f"work/{infix}/out/.done")
        yield "vcf", f"work/{infix}.final_comp.vcf.gz"
        yield "vcf_tbi", f"work/{infix}.final_comp.vcf.gz.tbi"

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        vcfs = []
        for me_type in self.config["melt"]["me_types"]:
            infix = f"{wildcards.mapper}.melt.merge_vcf.{wildcards.me_type}"
            vcfs.append(f"work/{infix}/out/{me_type}.final_comp.vcf.gz")
        yield "vcf", vcfs

    @dictify
    def _get_output_files_merge_vcf(self):
        infix = "{mapper}.melt.merge_vcf"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        infix = f"{wildcards.mapper}.melt.merge_vcf"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"

    @dictify
    def _get_output_files_reorder_vcf(self):
        infix = "{mapper}.melt.{index_library_name}"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"


class Sniffles2StepPart(BaseStepPart):
    """WGS SV identification using Sniffles 2"""

    name = "sniffles2"
    actions = ("bam_to_snf", "snf_to_vcf")

    _resource_usage = ResourceUsage(threads=2, time="0-02:00:00", memory="4G")
    resource_usage_dict = {
        "bam_to_snf": _resource_usage,
        "snf_to_vcf": _resource_usage,
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

    @dictify
    def _get_input_files_bam_to_snf(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")

    @dictify
    def _get_output_files_bam_to_snf(self):
        infix = "{mapper}.sniffles2_bam_to_snf.{library_name}"
        yield "snf", f"work/{infix}/out/{infix}.snf"

    @dictify
    def _get_input_files_snf_to_vcf(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        snfs = []
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                infix = f"{wildcards.mapper}.sniffles2_bam_to_snf.{donor.dna_ngs_library.name}.snf"
                snfs.append(f"work/{infix}/out/{infix}.snf")
        yield "snf", snfs

    @dictify
    def _get_output_files_snf_to_vcf(self):
        infix = "{mapper}.sniffles2.{index_ngs_library}"
        yield "snf", f"work/{infix}/out/{infix}.snf"


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
                MeltStepPart,
                Sniffles2StepPart,
                WritePedigreeStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            filter(is_not_background, sheets)
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample sheets.
        """
        for sub_step in self.sub_steps.values():
            yield from sub_step.get_result_files()

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
