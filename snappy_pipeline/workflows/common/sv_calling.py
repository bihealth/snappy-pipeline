"""Code shared by steps in ``sv_calling_{targeted,wgs}``"""

from snakemake.io import expand

from snappy_pipeline.utils import DictQuery, dictify, flatten, listify


class SvCallingGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for SV calling steps"""

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return a list of all paths with prefix ``output/` that are
        returned by ``self.get_output_files()`` for all actions in ``self.actions``.
        """
        ngs_mapping_config = DictQuery(self.w_config).get("step_config/ngs_mapping")
        for mapper in ngs_mapping_config["tools"]["dna"]:
            # Get list of result path templates.
            output_files_tmp = self.get_output_files(self.actions[-1])
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
                for library_name in self.index_ngs_library_to_pedigree.keys():
                    yield from expand(path_tpl, mapper=[mapper], library_name=library_name)


class SvCallingGetLogFileMixin:
    """Mixin that provides ``get_log_files()`` for SV calling steps"""

    @dictify
    def get_log_file(self, action):
        """Return dict of log files in the "log" directory"""
        _ = action
        if action != "merge_genotypes":
            infix = f"delly2_{action}"
        else:
            infix = "delly2"
        prefix = f"work/{{mapper}}.{infix}.{{library_name}}/log/{{mapper}}.{infix}.{{library_name}}.sv_calling"
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
