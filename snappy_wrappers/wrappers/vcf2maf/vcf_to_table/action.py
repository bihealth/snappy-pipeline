from enum import Enum


class Action(Enum):
    ERROR = 0
    STOP = 1
    SKIP = 2
    IGNORE = 3
    DEFAULT = 4
    OK = 5
