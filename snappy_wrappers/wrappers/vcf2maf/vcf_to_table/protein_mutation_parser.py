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
                "0",  # No group, but seq in 2
                substitution,  # 3: ref, 4: position, 5: alt
                duplication,  # 6: start, 7: start pos, 9: end, 10: end pos
                deletion,  # 11: start, 12: start pos, 14: end, 15: end pos
                insertion,  # 16: start, 17: start pos, 18: end, 19: end pos (=start pos + 1), 20+22: insertion
                delins,  # 23: start, 24: start pos, 26: end, 27: end pos, 28+30: insertion
                frameshift,  # 31: start, 33: start pos, 33: end, 35: ter pos
                extensionN,  # 37: pos (negative number)
                extensionC,  # 39: pos, 40: alt, 42: ter pos
            ]
        )
        + ")"
        + postfix
    )

    return re.compile(pattern)


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
    if groups[2] == "0":
        rslt = "0"

    elif groups[3] is not None:
        # Substitutions
        rslt = (
            _aa(groups[3], want_short=want_short, want_long=want_long)
            + groups[4]
            + _aa(groups[5], want_short=want_short, want_long=want_long)
        )
        # print("DEBUG- substitution: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(3, 6)]), rslt))

    elif groups[6] is not None:
        # Duplications
        rslt = _one_or_interval(groups, 6, want_short=want_short, want_long=want_long) + "dup"
        # print("DEBUG- duplication: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(6, 11)]), rslt))

    elif groups[11] is not None:
        # Deletions
        rslt = _one_or_interval(groups, 11, want_short=want_short, want_long=want_long) + "del"
        # print("DEBUG- deletion: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(11, 16)]), rslt))

    elif groups[16] is not None:
        # Insertions
        is_long = _is_long(groups[16])
        rslt = (
            _interval(groups, 16, want_short=want_short, want_long=want_long)
            + "ins"
            + _many_aas(groups[20], want_short=want_short, want_long=want_long, is_long=is_long)
            + _aa(groups[22], want_short=want_short, want_long=want_long)
        )
        # print("DEBUG- insertion: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(16, 23)]), rslt))

    elif groups[23] is not None:
        # Deletion-insertions
        is_long = _is_long(groups[23])
        rslt = (
            _one_or_interval(groups, 23, want_short=want_short, want_long=want_long)
            + "delins"
            + _many_aas(groups[28], want_short=want_short, want_long=want_long, is_long=is_long)
            + _aa(groups[30], want_short=want_short, want_long=want_long)
        )
        # print("DEBUG- delins: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(23, 31)]), rslt))

    elif groups[31] is not None:
        # Frameshift
        rslt = (
            _aa(groups[31], want_short=want_short, want_long=want_long)
            + groups[32]
            + _aa(groups[33], want_short=want_short, want_long=want_long)
            + "fs"
            + _aa(groups[34], want_short=want_short, want_long=want_long)
            + groups[35]
        )
        # print("DEBUG- frameshift: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(31, 36)]), rslt))

    elif groups[36] is not None:
        # Extension N terminus
        rslt = _aa(groups[36], want_short=want_short, want_long=want_long) + "1ext" + groups[37]
        # print("DEBUG- extension N: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(36, 38)]), rslt))

    elif groups[38] is not None:
        # Extension C terminus
        rslt = (
            _aa(groups[38], want_short=want_short, want_long=want_long)
            + groups[39]
            + _aa(groups[40], want_short=want_short, want_long=want_long)
            + "ext"
            + _aa(groups[41], want_short=want_short, want_long=want_long)
            + groups[42]
        )
        # print("DEBUG- extension C: groups = {}, rslt = {}".format(", ".join(["{}={}".format(i, groups[i]) for i in range(38, 43)]), rslt))

    return rslt


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
    if args is not None and isinstance(args, dict):
        add_sequence = args.get("add_sequence", add_sequence)
        want_short = args.get("want_short", want_short)
        want_long = args.get("want_long", want_long)
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
