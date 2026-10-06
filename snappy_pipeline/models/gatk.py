from typing import Annotated

from pydantic import AfterValidator, Field

from snappy_pipeline.models import SnappyModel


def argument(args: list[str]) -> list[str]:
    def _is_valid_argument(arg: str) -> bool:
        return arg.startswith("--")

    if invalid_args := list(filter(lambda x: not _is_valid_argument(x), args)):
        raise ValueError(f"invalid arguments: {invalid_args}")
    return args


class GATK(SnappyModel):
    extra_arguments: Annotated[
        list[str],
        AfterValidator(argument),
        Field(
            examples=[
                "--max-depth-per-sample 1000",
                "--maximum-population-allele-frequency 0.2",
            ]
        ),
    ] = []
    """
    List additional GetPileupSummaries arguments.
    Each additional argument must be of the form:
    "--<argument name> <argument value>"
    """

    java_options: Annotated[str, Field(examples=["-Xms4000m -Xmx8000m"])] = ""
    """Options for the java run-time compiler"""
