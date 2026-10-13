"""Code shared by steps in ``sv_calling_{targeted,wgs}``"""

import re

from snakemake.io import expand

from snappy_pipeline.utils import dictify, flatten, listify


#: Output keys and paths of the VCF that the SV callers publish
VCF_OUTPUTS = {
    "vcf": "work/{library_name}/out/{library_name}.vcf.gz",
    "vcf_md5": "work/{library_name}/out/{library_name}.vcf.gz.md5",
    "vcf_tbi": "work/{library_name}/out/{library_name}.vcf.gz.tbi",
    "vcf_tbi_md5": "work/{library_name}/out/{library_name}.vcf.gz.tbi.md5",
}


class SvCallingGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for SV calling steps"""

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return the paths of ``self.get_output_files()`` of the last
        action, with the prefix ``work/`` replaced by ``output/``.
        """
        if self.name != self.config.tool:
            return  # tool not enabled, no result files

        ngs_mapping_config = self.parent.get_task_config("alignments")
        ngs_mapping_tool = ngs_mapping_config.tool
        is_dna = (
            ngs_mapping_tool.is_dna()
            if hasattr(ngs_mapping_tool, "is_dna")
            else str(ngs_mapping_tool) in {"bwa", "bwa_mem2", "minimap2", "mbcs"}
        )
        if is_dna:
            # Get list of result path templates.
            output_files_tmp = self.get_output_files(self.actions[-1])
            if isinstance(output_files_tmp, dict):
                output_files = output_files_tmp.values()
            else:
                output_files = output_files_tmp
            result_paths_tpls = [
                re.sub(r"^work/", "output/", p)
                for p in flatten(output_files)
                if p.startswith("work/")
            ]
            #: Generate all concrete output paths.
            for path_tpl in result_paths_tpls:
                for library_name in self.index_ngs_library_to_pedigree.keys():
                    if cfg := getattr(self.config, self.name, None):
                        if library_name not in cfg.skip_libraries:
                            yield from expand(path_tpl, library_name=library_name)


class SvCallingGetLogFileMixin:
    """Mixin that provides the log file of SV calling steps"""

    @dictify
    def _get_log_file(self, action: str):
        """Return dict of log files in the "log" directory"""
        if hasattr(self, f"_get_log_file_infix_{action}"):
            infix = getattr(self, f"_get_log_file_infix_{action}")()
        else:
            infix = "{library_name}"
        yield "log", f"work/{infix}/log/{infix}.{action}.log"
