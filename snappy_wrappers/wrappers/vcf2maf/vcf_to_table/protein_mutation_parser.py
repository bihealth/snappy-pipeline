import functools
import re


class ProteinMutationFormatException(Exception):
    pass


aa_codes_short = (
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "B",
    "Z",
    "X",
)

aa_codes_long = (
    "Ala",
    "Cys",
    "Asp",
    "Glu",
    "Phe",
    "Gly",
    "His",
    "Ile",
    "Lys",
    "Leu",
    "Met",
    "Asn",
    "Pro",
    "Gln",
    "Arg",
    "Ser",
    "Thr",
    "Val",
    "Trp",
    "Tyr",
    "Asx",
    "Glx",
    "Xaa",
)

long_to_short = {
    "Ala": "A",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y",
    "Asx": "B",
    "Glx": "Z",
    "Xaa": "X",
    "Ter": "*",
    "A": "A",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "Y": "Y",
    "B": "B",
    "Z": "Z",
    "X": "X",
    "*": "*",
    "=": "=",
    "?": "?",
}

short_to_long = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
    "B": "Asx",
    "Z": "Glx",
    "X": "Xaa",
    "*": "Ter",
    "Ala": "Ala",
    "Cys": "Cys",
    "Asp": "Asp",
    "Glu": "Glu",
    "Phe": "Phe",
    "Gly": "Gly",
    "His": "His",
    "Ile": "Ile",
    "Lys": "Lys",
    "Leu": "Leu",
    "Met": "Met",
    "Asn": "Asn",
    "Pro": "Pro",
    "Gln": "Gln",
    "Arg": "Arg",
    "Ser": "Ser",
    "Thr": "Thr",
    "Val": "Val",
    "Trp": "Trp",
    "Tyr": "Tyr",
    "Asx": "Asx",
    "Glx": "Glx",
    "Xaa": "Xaa",
    "Ter": "Ter",
    "=": "=",
    "?": "?",
}


@functools.lru_cache
def _build_protein_pattern():
    prefix = "^(([A-z0-9_\.\(\)-]+):)?p\.\(?"  # noqa: W605
    postfix = "\)?$"  # noqa: W605

    aa = "(" + "|".join(aa_codes_short) + "|" + "|".join(aa_codes_long) + ")"
    aaTer = (
        "("
        + "|".join(aa_codes_short)
        + "|"
        + "|".join(aa_codes_long)
        + "|\*|Ter"  # noqa: W605
        + ")"
    )
    aaAll = (
        "("
        + "|".join(aa_codes_short)
        + "|"
        + "|".join(aa_codes_long)
        + "|\*|Ter|=|\?"  # noqa: W605
        + ")"
    )
    nb = "([0-9]+)"
    nb_unknown = "([0-9]+|\?)"  # noqa: W605

    interval = aa + nb + "_" + aa + nb
    one_or_interval = aa + nb + "(_" + aa + nb + ")?"

    substitution = aaTer + nb + aaAll

    deletion = one_or_interval + "del"
    duplication = one_or_interval + "dup"

    insertion = interval + "ins" + "(" + aa + "*)" + aaTer

    delins = one_or_interval + "delins" + "(" + aa + "*)" + aaTer

    frameshift = aaTer + nb + aa + "fs" + "(Ter|\*)" + nb_unknown  # noqa: W605

    extensionN = "(Met|M)1" + "ext(-[0-9]+|\?)"  # noqa: W605
    extensionC = "(Ter|\*)" + nb + aa + "ext" + "(Ter|\*)" + nb_unknown  # noqa: W605

    pattern = (
        prefix
        + "("
        + "|".join(
            [
                "(0|=|\?)",  # noqa: W605
                substitution,  # 4: ref, 5: position, 6: alt
                duplication,  # 7: start, 8: start pos, 10: end, 11: end pos
                deletion,  # 12: start, 13: start pos, 15: end, 16: end pos
                insertion,  # 17: start, 18: start pos, 19: end, 20: end pos (=start pos + 1), 21+23: insertion
                delins,  # 24: start, 25: start pos, 27: end, 28: end pos, 29+31: insertion
                frameshift,  # 32: start, 34: start pos, 34: end, 36: ter pos
                extensionN,  # 38: pos (negative number)
                extensionC,  # 40: pos, 41: alt, 43: ter pos
            ]
        )
        + ")"
        + postfix
    )

    return re.compile(pattern)


@functools.lru_cache
def _build_silent_dinucleotide():
    prefix = "^(([A-z0-9_\.\(\)-]+):)?p\.\(?"  # noqa: W605
    postfix = "\)?$"  # noqa: W605

    aa = "(" + "|".join(aa_codes_short) + "|" + "|".join(aa_codes_long) + ")"
    aaTer = (
        "("
        + "|".join(aa_codes_short)
        + "|"
        + "|".join(aa_codes_long)
        + "|\*|Ter"  # noqa: W605
        + ")"
    )
    nb = "([0-9]+)"

    return re.compile(prefix + aa + aaTer + "*" + nb + "=" + postfix)


def _aa(aa, want_short=True, want_long=False):
    if want_short and aa in long_to_short.keys():
        return long_to_short[aa]
    if want_long and aa in short_to_long.keys():
        return short_to_long[aa]
    return aa


def _is_long(aa):
    if aa not in long_to_short.keys():
        raise ProteinMutationFormatException("Unexpected AA code {}".format(aa))
    return len(aa) == 3


def _many_aas(many_aas, want_short=True, want_long=False, is_long=False):
    if not many_aas:
        return ""
    if not want_short and not want_long:
        return many_aas
    n = 3 if is_long else 1
    if want_short:
        mapping = long_to_short
    elif want_long:
        mapping = short_to_long
    else:
        mapping = None
    # print("DEBUG- _many_aas source = {} by {}".format(many_aas, n))
    i = 0
    converted = ""
    while i < len(many_aas):
        converted += mapping[many_aas[i : i + n]]
        i += n
    return converted


def _interval(groups, i0, want_short=True, want_long=False):
    first_aa = _aa(groups[i0 + 0], want_short=want_short, want_long=want_long)
    first_n = groups[i0 + 1]
    last_aa = _aa(groups[i0 + 2], want_short=want_short, want_long=want_long)
    last_n = groups[i0 + 3]
    return first_aa + first_n + "_" + last_aa + last_n


def _one_or_interval(groups, i0, want_short=True, want_long=False):
    one_aa = _aa(groups[i0 + 0], want_short=want_short, want_long=want_long)
    one_n = groups[i0 + 1]
    rslt = one_aa + one_n
    if groups[i0 + 2] is not None:
        interval_aa = _aa(groups[i0 + 3], want_short=want_short, want_long=want_long)
        interval_n = groups[i0 + 4]
        rslt += "_" + interval_aa + interval_n
    return rslt


def _parse_mutation_groups(groups, want_short=True, want_long=False):
    iGroup = 3
    if groups[iGroup] and groups[iGroup] in ("0", "=", "?"):
        return groups[iGroup]
    iGroup += 1

    if groups[iGroup] is not None:
        # Substitutions
        # print("DEBUG- substitution: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(4, 7)])))
        return (
            _aa(groups[iGroup], want_short=want_short, want_long=want_long)
            + groups[iGroup + 1]
            + _aa(groups[iGroup + 2], want_short=want_short, want_long=want_long)
        )
    iGroup += 3

    if groups[iGroup] is not None:
        # Duplications
        # print("DEBUG- duplication: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(8, 12)])))
        return _one_or_interval(groups, iGroup, want_short=want_short, want_long=want_long) + "dup"
    iGroup += 5

    if groups[iGroup] is not None:
        # Deletions
        # print("DEBUG- deletion: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(12, 17)])))
        return _one_or_interval(groups, iGroup, want_short=want_short, want_long=want_long) + "del"
    iGroup += 5

    if groups[iGroup] is not None:
        # Insertions
        # print("DEBUG- insertion: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(17, 24)])))
        is_long = _is_long(groups[iGroup])
        return (
            _interval(groups, iGroup, want_short=want_short, want_long=want_long)
            + "ins"
            + _many_aas(
                groups[iGroup + 4], want_short=want_short, want_long=want_long, is_long=is_long
            )
            + _aa(groups[iGroup + 6], want_short=want_short, want_long=want_long)
        )
    iGroup += 7

    if groups[iGroup] is not None:
        # Deletion-insertions
        # print("DEBUG- delins: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(24, 32)])))
        is_long = _is_long(groups[iGroup])
        return (
            _one_or_interval(groups, iGroup, want_short=want_short, want_long=want_long)
            + "delins"
            + _many_aas(
                groups[iGroup + 5], want_short=want_short, want_long=want_long, is_long=is_long
            )
            + _aa(groups[iGroup + 7], want_short=want_short, want_long=want_long)
        )
    iGroup += 8

    if groups[iGroup] is not None:
        # Frameshift
        # print("DEBUG- frameshift: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(32, 37)])))
        return (
            _aa(groups[iGroup], want_short=want_short, want_long=want_long)
            + groups[iGroup + 1]
            + _aa(groups[iGroup + 2], want_short=want_short, want_long=want_long)
            + "fs"
            + _aa(groups[iGroup + 3], want_short=want_short, want_long=want_long)
            + groups[iGroup + 4]
        )
    iGroup += 5

    if groups[iGroup] is not None:
        # Extension N terminus
        # print("DEBUG- extension N: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(37, 39)])))
        return (
            _aa(groups[iGroup], want_short=want_short, want_long=want_long)
            + "1ext"
            + groups[iGroup + 1]
        )
    iGroup += 2

    if groups[iGroup] is not None:
        # Extension C terminus
        # print("DEBUG- extension C: groups = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(39, 44)])))
        return (
            _aa(groups[iGroup], want_short=want_short, want_long=want_long)
            + groups[iGroup + 1]
            + _aa(groups[iGroup + 2], want_short=want_short, want_long=want_long)
            + "ext"
            + _aa(groups[iGroup + 3], want_short=want_short, want_long=want_long)
            + groups[iGroup + 4]
        )
    iGroup += 5

    return ""


def _parse_one_mutation(mut, want_short=True, want_long=False, add_sequence=False):
    if mut is None:
        return ""

    # print("DEBUG- received mutation {}".format(mut))
    m = _build_protein_pattern().match(mut)
    if not m:
        raise ProteinMutationFormatException(
            'Protein mutation "{}" doesn\'t follow the HGVS nomenclature'.format(mut)
        )
    groups = m.groups()
    # print("DEBUG- mutation {} is valid, found groups are:".format(mut))
    # for i in range(len(groups)):
    #     if groups[i] is not None:
    #         print("    {}: {}".format(i, groups[i]))

    rslt = "p." + _parse_mutation_groups(groups, want_short=want_short, want_long=want_long)

    if add_sequence:
        rslt = groups[1] + ":" + rslt

    return rslt


def parse_protein_mutation(x, args):
    if x is None:
        return None
    assert isinstance(x, list) and len(x) == 1
    if x[0] is None:
        return None
    x = x[0]

    add_sequence = False
    want_short = True
    want_long = False
    illegal = True
    if args is not None and isinstance(args, dict):
        add_sequence = args.get("add_sequence", add_sequence)
        want_short = args.get("want_short", want_short)
        want_long = args.get("want_long", want_long)
        illegal = args.get("illegal", illegal)
    if want_short and want_long:
        raise Exception("Cannot force both short & long amino-acide representations")

    rslts = []
    for mut in x:
        try:
            rslt = _parse_one_mutation(
                mut, want_short=want_short, want_long=want_long, add_sequence=add_sequence
            )
            rslts.append(rslt)
        except ProteinMutationFormatException as e:
            if not illegal:
                m = _build_silent_dinucleotide().match(mut)
                if not m:
                    rslts.append("")
                else:
                    raise Exception(e)
            else:
                raise Exception(e)

    return rslts


def _test_re(x):
    m = _build_protein_pattern().match(x)
    if not m:
        print('"{}" doesn\'t match'.format(x))
    else:
        print('Groups of "{}":'.format(x))
        gs = m.groups()
        for i in range(len(gs)):
            if gs[i] is not None:
                print("    {}: {}".format(i, gs[i]))


# with open("t", "rt") as f:
#     for line in f:
#         test_re(line.strip())
