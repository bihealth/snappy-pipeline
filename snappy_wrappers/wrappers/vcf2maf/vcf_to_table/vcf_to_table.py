import argparse
import logging
import sys
import typing

import exceptions
import vcfpy
from action import Action
from functions import Functions
from parser import VcfParser
from ruamel.yaml import YAML


def _get_command_line_parser():
    parser = argparse.ArgumentParser(description="VCF to MAF conversion")

    parser.add_argument("--config", help="configuration to parse the vcf file")

    parser.add_argument(
        "--select-alt",
        help="select one allele for output (generally set to 0), otherwise all alleles are output",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--unique",
        help="uniquify output (useful for repeated missing annotations)",
        action="store_true",
    )

    parser.add_argument("--vcf-tumor-id", help="Tumor sample id as in the vcf column name")
    parser.add_argument("--vcf-normal-id", help="Normal sample id as in the vcf column name")
    parser.add_argument("--tumor-id", help="Tumor sample id output in the maf file")
    parser.add_argument("--normal-id", help="Normal sample id output in the maf file")

    parser.add_argument(
        "--NCBI_Build", help="Genome version (GRCh37 & GRCh38 are accepted by cBioPortal)"
    )
    parser.add_argument("--Center", help="Name of the center for the MAF")

    parser.add_argument("--title", action="store_true", help="Write column names on output")

    parser.add_argument("--debug", action="store_true", help="Turn on debugging logs")

    parser.add_argument("vcf", help="Input vcf filename", nargs="?", default=None, type=str)
    parser.add_argument(
        "maf", help="Output maf file", nargs="?", type=argparse.FileType("w"), default=sys.stdout
    )

    return parser


def _check_config(config: typing.Dict[str, typing.Any]) -> typing.Dict[str, typing.Any]:
    assert isinstance(config, dict), "Illegal (or empty) configuration file"
    assert "output" in config, '"output" missing from configuration file'
    assert isinstance(config["output"], dict), "Illegal (or empty) output definition"
    for k, v in config["output"].items():
        assert isinstance(v, dict), 'Illegal (or empty) definition of output column "{}"'.format(k)
        assert "input" in v or not isinstance(v["input"], list), (
            'Missing inputs for output column "{}"'.format(k)
        )
        for i in range(len(v["input"])):
            x = v["input"][i]
            assert isinstance(x, dict), (
                "Illegal (or empty) element {} of inputs definition for output column {}".format(
                    i + 1, k
                )
            )
            assert "set" in x and "column" in x, (
                '"set" and/or "column" definition missing from element {} of inputs definition for output column {}'.format(
                    i + 1, k
                )
            )
            assert x["set"] in (
                "annotation",
                "constant",
                "fixed",
                "format",
                "info",
                "variant",
            ), (
                'Illegal set "{}" for input column "{}" of inputs definition for output column {}, must be one of "annotation", "constant", "fixed", "format", "info" or "variant"'.format(
                    x["set"], x["column"], k
                )
            )
            if x["set"] == "format":
                assert "sample" in x, (
                    'Sample definition is missing from set "{}" for input column "{}" of inputs definition for output column {}, must be present for "format" sets'.format(
                        x["set"], x["column"], k
                    )
                )
            if "on_missing" in x:
                assert x["on_missing"] in (
                    "stop",
                    "skip",
                    "default",
                    "ignore",
                ), (
                    'Illegal value of requested action for missing input "{}" for output column "{}", must be either "stop", "skip", "default" or "ignore"'.format(
                        x["column"], k
                    )
                )
            else:
                config["output"][k]["input"][i]["on_missing"] = "ignore"
            if x["on_missing"] == "default":
                assert "default" in x and x["default"] is not None, (
                    'Missing default value for input "{}" for output column "{}"'.format(
                        x["column"], k
                    )
                )
            else:
                config["output"][k]["input"][i]["default"] = None
        if "on_missing" in v:
            assert v["on_missing"] in (
                "stop",
                "skip",
                "default",
                "ignore",
            ), (
                'Illegal value of requested action for output column "{}", must be either "stop", "skip", "default" or "ignore"'.format(
                    k
                )
            )
        else:
            config["output"][k]["on_missing"] = "ignore"
        if v["on_missing"] == "default":
            assert "default" in v and v["default"] is not None, (
                'Missing default value for output column "{}"'.format(k)
            )
        else:
            config["output"][k]["default"] = None
    return config


def _add_constants(parser: VcfParser, args):
    constants = ("NCBI_Build", "Center")
    for constant in constants:
        parser.addConstant(constant, vars(args).get(constant, None))
    parser.addConstant("tumor_id", args.tumor_id if args.tumor_id else args.vcf_tumor_id)
    parser.addConstant("normal_id", args.normal_id if args.normal_id else args.vcf_normal_id)


def _get_column_values(record, parser, description, functions):
    arguments = []
    for i in range(len(description["input"])):
        inputs = description["input"][i]
        value = parser.getInputValue(
            record=record,
            theSet=inputs["set"],
            theColumn=inputs["column"],
            theSample=inputs.get("sample", None),
            default=inputs.get("default", None),
        )

        if description.get("action", "ignore") == "stop":
            if value is None:
                raise exceptions.ActionStop(
                    'Missing mandatory input "{}" of set "{}"'.format(
                        inputs["column"], inputs["set"]
                    )
                )
            for alt in record.ALT:
                if alt.value not in value or value[alt.value] is None:
                    raise exceptions.ActionStop(
                        'Missing mandatory input "{}" of set "{}" for alt variant "{}"'.format(
                            inputs["column"], inputs["set"], alt.value
                        )
                    )
                for i in range(len(value[alt.value])):
                    if value[alt.value][i] is None:
                        raise exceptions.ActionStop(
                            'Missing mandatory element {} from input "{}" of set "{}" for alt variant "{}"'.format(
                                i + 1, inputs["column"], inputs["set"], alt.value
                            )
                        )

        arguments.append(value)
    return functions.run(
        arguments,
        description.get("function", "identity"),
        description.get("args", None),
        description.get("default", None),
    )


def _check_result(
    result: typing.Dict[str, typing.List[typing.Any]],
    action: str,
    default: str,
    alts: typing.List[str],
) -> typing.Dict[str, typing.List[typing.Any]]:
    assert action != "default" or default is not None

    if result is None:
        if action == "stop":
            raise exceptions.ActionStop("Missing result")
        elif action == "default":
            result = {alt: [default] for alt in alts}
        else:
            pass
    else:
        assert isinstance(result, dict)
        for alt in alts:
            if alt not in result or result[alt] is None:
                if action == "stop":
                    raise exceptions.ActionStop('Missing result for allele "{}"'.format(alt))
                elif action == "default":
                    result[alt] = [default]
                else:
                    pass
            else:
                assert isinstance(result[alt], list)
                for i in range(len(result[alt])):
                    if result[alt][i] is None:
                        if action == "stop":
                            raise exceptions.ActionStop(
                                'Missing element {} from result for allele "{}"'.format(i + 1, alt)
                            )
                        elif action == "default":
                            result[alt][i] = default
                        else:
                            pass
    return result


def _update_one_status(status, result, action):
    if status == Action.SKIP:
        return status
    assert isinstance(status, list)

    if result is None:
        if action == "stop":
            raise exceptions.ActionStopException("Missing value")
        if action == "skip":
            return Action.SKIP
        return status
    assert isinstance(result, list)

    if len(status) == 1 and len(result) > 1:
        status = [status[0] for i in range(len(result))]
    if len(status) > 1 and len(result) == 1:
        result = [result[0] for i in range(len(status))]
    if len(result) != len(status):
        raise exceptions.IllegalValue(
            "Expansion impossible, previous length = {}, current length = {}".format(
                len(status), len(result)
            )
        )

    for i in range(len(result)):
        if result[i] is not None:
            continue
        if action == "stop":
            raise exceptions.ActionStopException("Missing value")
        if action == "skip":
            status[i] = Action.SKIP
    return status


def _update_status(status, result, description):
    if status == Action.SKIP:
        return status

    if result is None:
        if description.get("on_missing", "ignore") == "skip":
            return Action.SKIP
        else:
            return status
    assert isinstance(result, dict)

    for alt, s in status.items():
        status[alt] = _update_one_status(
            s, result.get(alt, None), description.get("on_missing", "ignore")
        )
    return status


def _add_to_output(results, status, alts, record):
    if status != Action.SKIP:
        assert isinstance(status, dict)
        for alt in alts:
            if status[alt] == Action.SKIP:
                logging.info("Skipped record {}".format(vcf_record_to_string(record)))
                continue
            for i in range(len(status[alt])):
                if status[alt][i] == Action.SKIP:
                    logging.info("Skipped record {}".format(vcf_record_to_string(record)))
                    continue
                row = []
                for col, result in results.items():
                    r = None
                    if alt in result and result[alt]:
                        assert isinstance(result[alt], list)
                        if len(result[alt]) == 1:
                            r = result[alt][0]
                        else:
                            r = result[alt][i]
                    row.append(r)
                yield [str(x) if x is not None else "" for x in row]
    else:
        logging.info("Skipped record {}".format(vcf_record_to_string(record)))


def vcf_record_to_string(record):
    fields = [
        record.CHROM,
        str(record.POS),
        record.REF,
        ",".join([alt.value for alt in record.ALT]),
    ]
    info = list()
    for k, v in record.INFO.items():
        if isinstance(v, list) and not isinstance(v, str):
            v = ",".join([str(el) for el in v])
        info.append("{}={}".format(k, str(v)))
    fields.append(";".join(info))
    for call in record.calls:
        fmt = list()
        for v in call.data.values():
            if isinstance(v, list) and not isinstance(v, str):
                v = ",".join([str(el) for el in v])
            fmt.append(str(v))
        fields.append(";".join(fmt))
    return "\t".join(fields)


def main():
    parser = _get_command_line_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        level=logging.DEBUG if args.debug else logging.WARNING,
    )

    if not args.vcf_tumor_id and not args.vcf_normal_id:
        logging.error("Missing tumor and normal vcf ids")
        return -1
    samples = {"tumor": args.vcf_tumor_id, "normal": args.vcf_normal_id}

    with open(args.config, "r") as f:
        logging.info("Loading configuration from {}".format(args.config))
        config = YAML().load(f)
    try:
        config = _check_config(config)
    except AssertionError as e:
        logging.error(e)
        return -1

    functions = Functions(config)

    if args.vcf is None:
        logging.info("Reading VCF file from stdin")
        reader = vcfpy.Reader.from_stream(sys.stdin)
    else:
        logging.info("Reading VCF file {}".format(args.vcf))
        reader = vcfpy.Reader.from_path(args.vcf)
    vcf = VcfParser(config, reader.header, samples)

    _add_constants(vcf, args)

    if args.title:
        print("\t".join(config["output"].keys()), file=args.maf)

    unique = set()
    try:
        iRecord = 0
        for record in reader:
            iRecord += 1

            if args.select_alt is not None:
                if args.select_alt >= len(record.ALT):
                    logging.info(
                        "Record {} ignored, number of alt alleles ({}) smaller than requested ({})".format(
                            iRecord, len(record.ALT), args.select_alt + 1
                        )
                    )
                    continue
                record = vcf.subset(record, args.select_alt)
            vcf.annotation.set_minimizer(record)
            alts = [alt.value for alt in record.ALT]

            status = {alt: [Action.OK] for alt in alts}
            results = {}
            for column, description in config["output"].items():
                result = _get_column_values(record, vcf, description, functions)

                results[column] = _check_result(
                    result,
                    description.get("on_missing", "ignore"),
                    description.get("default", None),
                    alts,
                )

                status = _update_status(status, results[column], description)
                if status == Action.SKIP:
                    logging.info("Skipped record {}".format(vcf_record_to_string(record)))
                    break
            for row in _add_to_output(results, status, alts, record):
                s = "\t".join(row)
                if args.unique and s in unique:
                    continue
                print(s, file=args.maf)
                unique.add(s)
    except exceptions.FunctionError as e:
        logging.error(
            'Error during compution of values for column request "{}" of record {}:{}'.format(
                column, record.CHROM, record.POS
            )
        )
        logging.error("Message: {}".format(e))
        return -1
    except exceptions.ActionStop as e:
        logging.error(
            'Missing mandatory value for column request "{}" of record {}:{}'.format(
                column, record.CHROM, record.POS
            )
        )
        logging.error("Message: {}".format(e))
        return -1

    logging.info("Conversion complete, processed {} VCF records".format(iRecord))
    return 0


if __name__ == "__main__":
    sys.exit(main())
