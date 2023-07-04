import re


def minimize_mutation(x, args):
    assert isinstance(x, list) and len(x) == 3
    assert isinstance(x[0], list) and len(x[0]) == 1
    assert isinstance(x[1], list) and len(x[1]) == 1
    assert isinstance(x[2], list) and len(x[2]) == 1
    pos = int(x[0][0])
    ref = x[1][0]
    alt = x[2][0]

    assert args is not None and isinstance(args, dict) and "return_value" in args.keys()

    # print("DEBUG- on input: pos = {}, ref = {}, alt = {}".format(pos, ref, alt))

    i = 0
    # Trim common left part
    while i < len(ref) and i < len(alt):
        if ref[i] != alt[i]:
            break
        pos += 1
        i += 1
    if i == len(ref) and i == len(alt):
        raise Exception("Reference and alternate sequences are identical")
    ref = ref[i:]
    alt = alt[i:]

    # print("DEBUG- on after left trim: pos = {}, ref = {}, alt = {}".format(pos, ref, alt))

    # Trim common right part
    i = 0
    while i < len(ref) and i < len(alt):
        if ref[-i - 1] != alt[-i - 1]:
            if i > 0:
                ref = ref[:-i]
                alt = alt[:-i]
            break
        i += 1

    end = pos + len(ref) - 1

    # print("DEBUG- on output: pos = {}, ref = {}, alt = {}".format(pos, ref, alt))

    if ref == "":
        ref = "-"
        end = pos
        pos -= 1
    if alt == "":
        alt = "-"

    return_value = args["return_value"]
    if return_value == "pos":
        return [pos]
    elif return_value == "ref":
        return [ref]
    elif return_value == "alt":
        return [alt]
    elif return_value == "end":
        return [end]
    else:
        raise Exception('Unexpected return value request "{}"'.format(return_value))


def calc_end_pos(x, args):
    assert isinstance(x, list) and len(x) == 3
    assert isinstance(x[0], list) and len(x[0]) == 1
    assert isinstance(x[1], list) and len(x[1]) == 1
    assert isinstance(x[2], list) and len(x[2]) == 1
    pos = int(x[0][0])
    ref = x[1][0]
    alt = x[2][0]
    d = len(alt) - len(ref)
    if d <= 0:
        end = pos + len(ref) - 1
    else:
        end = pos + 1
    return [end]


def variant_type(x, args=None):
    assert isinstance(x, list) and len(x) == 2
    assert isinstance(x[0], list) and len(x[0]) == 1
    assert isinstance(x[1], list) and len(x[1]) == 1
    ref = x[0][0].replace("-", "")
    alt = x[1][0].replace("-", "")
    if len(ref) < len(alt):
        variant_type = "INS"
    elif len(ref) > len(alt):
        variant_type = "DEL"
    elif len(ref) == 0:
        raise Exception("Empty reference & alternative allele sequences")
    elif len(ref) == 1:
        variant_type = "SNP"
    elif len(ref) == 2:
        variant_type = "DNP"
    elif len(ref) == 3:
        variant_type = "TNP"
    else:
        variant_type = "ONP"
    return [variant_type]


strip_sequence_version_pattern = re.compile("\.[0-9]+$")  # noqa: W605


def strip_sequence_version(x, args):
    if x is None:
        return None
    assert isinstance(x, list) and len(x) == 1
    if x[0] is None:
        return None
    no_version = list()
    for el in x[0]:
        if el is None:
            no_version.append(None)
        else:
            no_version.append(strip_sequence_version_pattern.sub("", el))
    return no_version
