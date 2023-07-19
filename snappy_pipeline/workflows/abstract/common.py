"""Commonly used code and types"""

from itertools import chain
import re
import typing

from snakemake.io import Wildcards

from snappy_wrappers.resource_usage import ResourceUsage

#: Type for a Snaekmake key/path dict.
SnakemakeDict = typing.Dict[str, typing.Union[typing.List[str], str]]

#: Type for a Snakemake input function.
SnakemakeInputFunc = typing.Callable[
    [Wildcards],
    SnakemakeDict,
]

#: Type for generating pairs with key/path(s) mapping for our workflows.
SnakemakeDictItemsGenerator = typing.Generator[
    typing.Tuple[str, typing.Union[typing.List[str], str]],
    None,
    None,
]

#: Type for generating path(s) mapping for our workflows.
SnakemakeListItemsGenerator = typing.Generator[
    str,
    None,
    None,
]


class ForwardSnakemakeFilesMixin:
    """Mixin for forwarding calls to ``get_{input,output,log}_file(s)``"""

    def get_input_files(self, action: str) -> SnakemakeInputFunc:
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def get_output_files(self, action: str) -> SnakemakeDict:
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    def get_log_file(self, action: str) -> SnakemakeDict:
        self._validate_action(action)
        return getattr(self, f"_get_log_file_{action}")()


class ForwardResourceUsageMixin:
    """Mixin for forwarding calls to ``get_resource_usage`` to ``resource_usage_dict``."""

    #: Resource usage definitions
    resource_usage_dict: typing.Optional[typing.Dict[str, ResourceUsage]] = None

    def get_resource_usage(self, action: str) -> ResourceUsage:
        self._validate_action(action)
        assert self.resource_usage_dict is not None, "resource_usage_dict not set!"
        assert action in self.resource_usage_dict, f"No resource usage entry for {action}"
        return self.resource_usage_dict[action]


def augment_work_dir_with_output_links(
    work_dir_dict: SnakemakeDict, log_files: typing.Optional[typing.List[str]] = None
) -> SnakemakeDict:
    """Augment a dictionary with key/value pairs to work directory with ``"output_links"`` key.

    Optionally, the output files will be augmented from the paths in ``log_files``.
    """
    result = dict(work_dir_dict)
    result["output_links"] = [
        re.sub(r"^work/", "output/", work_path)
        for work_path in chain(work_dir_dict.values(), log_files or [])
    ]
    return result
