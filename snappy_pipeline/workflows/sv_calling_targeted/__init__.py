"""SV calling for targeted sequencing

Based on the output of ``ngs_mapping``, call structural variants from depth of coverage,
read pair, and split read signal.
"""

import re

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import SvCallingTargeted as SvCallingTargetedConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload (VCF)
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available SV callers
SV_CALLERS = ("gcnv", "delly2", "manta", "melt")

#: Minimum number of samples per kit to apply gCNV calling criteria to be analyzed
GCNV_MIN_KIT_SAMPLES = 10

#: Default configuration for the sv_calling_targeted step
DEFAULT_CONFIG = SvCallingTargetedConfigModel.default_config_yaml_string()


class GcnvTargetedStepPart(RunGcnvStepPart):
    """Targeted seq. CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self.parent.ngs_library_to_kit


class SvCallingTargetedWorkflow(BaseStep):
    """Perform germline targeted sequencing CNV calling"""

    #: Workflow name
    name = "sv_calling_targeted"

    sheet_shortcut_class = GermlineCaseSheet

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SvCallingTargetedConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Build mapping from NGS library name to kit
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                GcnvTargetedStepPart,
                Delly2StepPart,
                MantaStepPart,
                MeltStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)
        # Build dictionary with sample count per library kit
        _, _, self.library_kit_counts_dict = self.pick_kits_and_donors()

    @dictify
    def _build_ngs_library_to_kit(self):
        config = self.w_config.step_config["sv_calling_targeted"].gcnv
        if not config.path_target_interval_list_mapping:
            # No mapping given, we will use the "default" one for all.
            for donor in self.all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"
        # Build mapping
        regexes = {item.pattern: item.name for item in config.path_target_interval_list_mapping}
        result = {}
        for donor in self.all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    @classmethod
    def default_config_yaml(cls):
        """Default configuration.

        :return: Returns default config YAML, to be overwritten by project-specific one.
        """
        return DEFAULT_CONFIG

    def get_library_count(self, library_kit):
        """Get library count.

        :param library_kit: Library kit name.
        :type library_kit: str

        :return: Returns number of samples with inputted library kit. If library name not defined,
        it returns zero.
        """
        return self.library_kit_counts_dict.get(library_kit, 0)

    @listify
    def all_donors(self, include_background=True):
        """Get all donors.

        :param include_background: Boolean flag to defined if background should be included or not.
        Default: True, i.e., background will be included.

        :return: Returns list of all donors in sample sheet.
        """
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
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

    def pick_kits_and_donors(self):
        """Return ``(library_kits, donors)`` with the donors with a matching kit and the kits with a
        matching donor.
        """
        kit_counts = {name: 0 for name in self.ngs_library_to_kit.values()}
        for name in self.ngs_library_to_kit.values():
            kit_counts[name] += 1
        donors = [
            donor
            for donor in self.all_donors()
            if donor.dna_ngs_library and donor.dna_ngs_library.name in self.ngs_library_to_kit
        ]
        return list(sorted(set(self.ngs_library_to_kit.values()))), donors, kit_counts
