class IllegalValue(Exception):
    """Illegal value found during parsing"""

    pass


class MissingValue(Exception):
    """Missing value found during parsing"""

    pass


class ActionStop(Exception):
    """Interrupting parsing upon missing value or parsing error"""

    pass


class FunctionError(Exception):
    """Error encountered during function"""

    pass
