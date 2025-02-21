import enum
import re

from typing import Annotated, Any, Union

from pydantic import Field, AliasGenerator

from snappy_pipeline.models import SnappyModel

# TODO: These regex need to be checked & solidified!!
VCF_TAG_PATTERN = re.compile(
    r"""
^(
(?P<info>INFO/([A-Za-z_][0-9A-Za-z_.]*|1000G))            # INFO allows 1000G as tag (see https://samtools.github.io/hts-specs/VCFv4.5.pdf, section 1.6.1.8)
|
(?P<fmt_or_id>((FORMAT|FMT)/)?([A-Za-z_][0-9A-Za-z_.]*))  # FORMAT doesn't (see section 1.6.2), & IDs without INFO or FORMAT are allowed by bcftools
)$
""",
    re.VERBOSE,
)
ANNOTATION_VCF_TAG_PATTERN = re.compile(
    r"""
^
(?P<prefix>([\^=~\+\.-]|\.\+)?)                                                                 # Prefix: ^, =~, +, -, .+
(?P<target>(INFO/([A-Za-z_][0-9A-Za-z_.]*|1000G)|((FORMAT|FMT)/)?([A-Za-z_][0-9A-Za-z_.]*)))    # VCF_TAG_PATTERN
(:=(?P<src>(INFO/([A-Za-z_][0-9A-Za-z_.]*|1000G)|((FORMAT|FMT)/)?([A-Za-z_][0-9A-Za-z_.]*))))?  # Source VCF_TAG PATTERN (when source & dest are different tags)
$""",
    re.VERBOSE,
)
# TODO: create the pattern from BcftoolsMarkSitePresentAbsent enum
PRESENT_ABSENT_PATTERN = re.compile("^([+-])(.+)$")


def BcftoolsModelAliasGenerator(option: str, exceptions: dict[str, str] = {}) -> str:
    """
    The alias generator should replace underscores with dashes for an option,
    unless it belongs to the list of exceptions.
    In that case, the alias is generated from the exceptions dict.

    Currently, there is no replacement at all, because the automatic alias generation
    seems to interfere with the re-validation occurring when a sub-workflow is registered.
    """
    if option in exceptions:
        return exceptions[option]
    return option.replace("_", "_")


class BcftoolsModel(SnappyModel):
    model_config = {
        "alias_generator": AliasGenerator(
            serialization_alias=lambda x: BcftoolsModelAliasGenerator(x)
        )
    }

    extra_args: Annotated[
        dict[str, Any],
        Field(
            examples=[
                {"'max-depth'": 5000},
                {"columns": ["FORMAT/DP", "-INFO/SCBZ"]},
                {"'illumina1.3+'": True},
            ],
            alias="extra_args",
        ),
    ] = {}
    """
    Placeholder for arguments available to the tool, but not present in the model.
    The arguments must be input as 'option: value'. Values should not contain the apostrophe (') character.

    The user is responsible to check that the arguments make sense
    (for example that the maximum allowed value is greater that the minimum allowed value).

    To turn off a flag set in the model, the user should add the flag to the extra arguments and set it to False.
    """


# class BcftoolsBamFlag(enum.StrEnum):
#     PAIRED = enum.auto()
#     PROPER_PAIR = enum.auto()
#     UNMAP = enum.auto()
#     MUNMAP = enum.auto()
#     REVERSE = enum.auto()
#     MREVERSE = enum.auto()
#     READ1 = enum.auto()
#     READ2 = enum.auto()
#     SECONDARY = enum.auto()
#     QCFAIL = enum.auto()
#     DUP = enum.auto()
#     SUPPLEMENTARY = enum.auto()
#
#
# def check_enum_name(x: Any) -> list[int]:
#     if isinstance(x, BcftoolsBamFlags):
#         return [x.value]
#     if isinstance(x, str):
#         if x not in BcftoolsBamFlags.__members__:
#             raise ValueError(f"Illegal flag '{x}'")
#         return [BcftoolsBamFlags[x].value]
#     if isinstance(x, int):
#         all_flags = 0
#         l = []
#         for v in BcftoolsBamFlags:
#             all_flags |= v.value
#             if x & v != 0:
#                 l.append(v.value)
#         if x & all_flags != 0:
#             raise ValueError(f"Unexpected flags in '{x}'")
#         return list(set(l))
#     else:
#         raise ValueError(f"Illegal object '{x}'")
#
#
# def convert_enum_names(x: list[Any] | BcftoolsBamFlags | str | int) -> list[int]:
#     if isinstance(x, list):
#         l = []
#         for v in x:
#             l += check_enum_name(v)
#         l = list(set(l))
#     else:
#         l = check_enum_name(x)
#     return l
#
#
# class BcftoolsBamModel(BcftoolsModel):
#     skip_all_set: list[BcftoolsBamFlag | int] | int = []
#     """Skip reads with all of the bits set"""
#     skip_any_set: list[BcftoolsBamFlag | int] | int = [
#         BcftoolsBamFlag.UNMAP,
#         BcftoolsBamFlag.SECONDARY,
#         BcftoolsBamFlag.QCFAIL,
#         BcftoolsBamFlag.DUP
#     ]
#     """Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP]"""
#     skip_all_unset: list[BcftoolsBamFlag | int] | int = []
#     """Skip reads with all of the bits unset"""
#     skip_any_unset: list[BcftoolsBamFlag | int] | int = []
#     ""Skip reads with any of the bits unset"""
#
#     validate_skip_all_set = field_validator("skip_all_set", mode="before")(convert_enum_names)
#     validate_skip_any_set = field_validator("skip_any_set", mode="before")(convert_enum_names)
#     validate_skip_all_unset = field_validator("skip_all_unset", mode="before")(convert_enum_names)
#     validate_skip_any_unset = field_validator("skip_any_unset", mode="before")(convert_enum_names)


class BcftoolsBamFlag(enum.Flag):
    NONE = 0
    PAIRED = 1
    PROPER_PAIR = 2
    UNMAP = 4
    MUNMAP = 8
    REVERSE = 16
    MREVERSE = 32
    READ1 = 64
    READ2 = 128
    SECONDARY = 256
    QCFAIL = 512
    DUP = 1024
    SUPPLEMENTARY = 2048


BcftoolsBamFlagMultipleTypes = Union[str | int | list[str | int] | BcftoolsBamFlag]


def transform_to_flag(
    value: BcftoolsBamFlagMultipleTypes = BcftoolsBamFlag.NONE,
) -> BcftoolsBamFlag:
    flags = BcftoolsBamFlag.NONE
    if isinstance(value, list):
        for v in value:
            flags |= transform_to_flag(v)
    else:
        if isinstance(value, BcftoolsBamFlag):
            flag = value
        elif isinstance(value, str):
            try:
                value = int(value)
                flag = BcftoolsBamFlag(value)
            except ValueError:
                flag = BcftoolsBamFlag(getattr(BcftoolsBamFlag, value))
        elif isinstance(value, int):
            flag = BcftoolsBamFlag(value)
        else:
            raise ValueError(f"Illegal flag '{value}'")
        flags |= flag
    return flags


def ensure_unique(x: list[BcftoolsBamFlag], set_type: str = "?"):
    if len(x) > len(set(x)):
        raise ValueError("Duplicated flags '{}' in set '{}'".format(", ".join(x), set_type))


def ensure_no_common_elements(
    x: list[BcftoolsBamFlag], y: list[BcftoolsBamFlag], set_type: str = "?"
):
    intersection = set(x) & set(y)
    if len(intersection) > 0:
        raise ValueError(
            "Flag(s) '{}' are common to '{}' sets".format(", ".join(intersection), set_type)
        )


def ensure_valid_bit_sets(model: BcftoolsModel):
    for a in ("skip_all_set", "skip_any_set", "skip_all_unset", "skip_any_unset"):
        try:
            getattr(model, a)
        except AttributeError:
            raise ValueError(f"Model doesn't contain option '{a}'")

    ensure_unique(model.skip_all_set, set_type="all-set")
    ensure_unique(model.skip_any_set, set_type="any-set")
    ensure_unique(model.skip_all_unset, set_type="all-unset")
    ensure_unique(model.skip_any_unset, set_type="any-unset")

    ensure_no_common_elements(model.skip_all_set, model.skip_any_set, set_type="set")
    ensure_no_common_elements(model.skip_all_unset, model.skip_any_unset, set_type="unset")

    ensure_no_common_elements(
        set(model.skip_all_set + model.skip_any_set),
        set(model.skip_all_unset + model.skip_any_unset),
        set_type="complete",
    )

    return model


class BcftoolsCaller(enum.StrEnum):
    CONSENSUS = "consensus"
    MULTIALLELIC = "multiallelic"


class BcftoolsPloidy(enum.StrEnum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    X = "X"
    Y = "Y"
    HAPLOID = "1"
    DIPLOID = "2"


class BcftoolsMarkSitePresentAbsent(enum.StrEnum):
    PRESENT = "+"
    ABSENT = "-"


def ensure_valid_mark_sites_tag(option_str: str):
    if option_str is not None:
        m = PRESENT_ABSENT_PATTERN.match(option_str)
        if not m:
            raise ValueError(
                f"Mark sites tag '{option_str} must be prefixed with "
                f"'{BcftoolsMarkSitePresentAbsent.PRESENT}' (mark present sites) or "
                f"'{BcftoolsMarkSitePresentAbsent.ABSENT}' (mark absent sites)"
            )
        if not VCF_TAG_PATTERN.match(m.group(2)):
            raise ValueError(f"Illegal Mark sites tag '{m.group(2)}'")
    return option_str


def ensure_valid_tags(tags: list[str]):
    for tag in tags:
        if not ANNOTATION_VCF_TAG_PATTERN.match(tag):
            raise ValueError(f"Illegal annotation tag pattern '{tag}'")
    return tags
