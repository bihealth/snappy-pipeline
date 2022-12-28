"""Commonly used code and types"""

import typing

#: Type for a Snaekmake key/path dict.
SnakemakeDict = typing.Dict[str, typing.Union[typing.List[str], str]]

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
