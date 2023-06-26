import importlib.util
import inspect
import logging
import os
import re
import typing

# from action import Action
from common_functions import CommonFunctions
import exceptions


class Functions:
    mapping_files = {}
    TAB = re.compile("\t")

    def __init__(self, config: typing.Dict[str, typing.Any]):
        self.compiled = Functions.importCode(config["functions"]) if "functions" in config else {}
        self.compiled["identity"] = Functions.identity
        self.compiled["map"] = Functions.mapper

        self.compiled["minimize_mutation"] = CommonFunctions.minimize_mutation
        self.compiled["variant_classification"] = CommonFunctions.variant_classification
        self.compiled["variant_type"] = CommonFunctions.variant_type
        self.compiled["parse_protein_mutation"] = CommonFunctions.parse_protein_mutation
        self.compiled["strip_sequence_version"] = CommonFunctions.strip_sequence_version

        for col_name, col_def in config["output"].items():
            if "function" in col_def and not col_def["function"] in self.compiled:
                raise exceptions.MissingValue(
                    'Function "{}" requested to process column "{}" not found'.format(
                        col_def["function"], col_name
                    )
                )

    @staticmethod
    def importCode(functions, module_name="extra_functions"):
        """Read & compile processing functions from python files

        :param functions: list of python files to compile
        :return: a dict of compiled functions as values, with the function names as keys
        """
        to_compile = []
        for filename in functions:
            with open(filename, "rt") as f:
                code = f.read()
            to_compile.append(code)
            logging.info('Loaded additional function(s) from file "{}"'.format(filename))

        spec = importlib.util.spec_from_loader(module_name, loader=None)
        module = importlib.util.module_from_spec(spec)
        exec("\n\n".join(to_compile), module.__dict__)
        globals()[module_name] = module

        hooks = dict(inspect.getmembers(module, inspect.isfunction))
        logging.debug("Functions: {}".format(hooks.keys()))
        return hooks

    @staticmethod
    def identity(value, args: typing.Dict[str, typing.Any] = None):
        if value is None:
            return None
        assert isinstance(value, list)
        assert len(value) == 1
        result = value[0]

        return result

    @staticmethod
    def mapper(value, args: typing.Dict[str, typing.Any] = None):
        if args is None or not isinstance(args, dict) or "filename" not in args:
            raise exceptions.MissingValue("Missing mapping filename argument")
        filename = os.path.realpath(os.path.join(os.path.dirname(__file__), args["filename"]))

        if value is None:
            return None
        assert isinstance(value, list)
        assert len(value) == 1
        value = value[0]
        if value is None:
            return None

        if filename not in Functions.mapping_files:
            m = {}
            with open(filename, "rt") as f:
                for line in f:
                    try:
                        k, v = Functions.TAB.split(line.strip(), 1)
                        m[k] = v
                    except Exception:
                        pass
            Functions.mapping_files[filename] = m
        m = Functions.mapping_files[filename]

        return [m.get(v, None) for v in value]

    @staticmethod
    def max_length(previous: int = 0, current: int = 0) -> int:
        if current == 0:
            return previous
        if previous < 2:
            return current
        if current == 1:
            return previous
        if current != previous:
            raise exceptions.IllegalValue(
                "Record expansion impossible: previous length = {}, current length = {}".format(
                    previous, current
                )
            )
        return current

    def run(
        self,
        values,
        fct: str = "identity",
        args: typing.Dict[str, typing.Any] = None,
        default: str = None,
    ):
        if values is None:
            return None
        assert isinstance(values, list)
        alts = {}
        for value in values:
            if value is None:
                continue
            assert isinstance(value, dict)
            for k, v in value.items():
                if k not in alts:
                    alts[k] = 0
                if v is None:
                    continue
                assert isinstance(v, list) and len(v) > 0
                alts[k] = Functions.max_length(alts[k], len(v))

        f = self.compiled[fct]

        results = {}
        for alt in alts:
            x = []
            for value in values:
                assert value is None or value[alt] is None or isinstance(value[alt], list)
                x.append(value[alt])
            try:
                results[alt] = f(x, args)
            except Exception as e:
                logging.error('Function "{}" has returned error {}'.format(fct, e))
                # raise exceptions.FunctionError(
                #     'Error in function "{}" for value {} and arguments {}'.format(fct, x, args)
                # )
                results[alt] = None
            if results[alt] is None:
                results[alt] = default
            if results[alt] is not None and not isinstance(results[alt], list):
                results[alt] = [results[alt]]

        return results
