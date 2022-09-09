# -*- coding: utf-8 -*-
"""Quickly draw a Venn Diagram from a shared counts file

A "shared counts" file is a text file, each line describing the shared by some separator (default
is TAB).  The last entry is a number given the number of counts between sets.  Each line starts
with a list of set names, separated shared elements.  A number alone on a line is the number of
elements in none of the given sets.

Quickvenn can plot up to 5 sets.

Example (using space as the separator):
::

    one 3
    two 4
    one two 4
    one two three 5
    6

"""
import argparse
import itertools
import math
import sys

import matplotlib  # isort:skip

# Force using Agg backend, has to go before import of pyplot
matplotlib.use("Agg")  # noqa
# pylint: disable=wrong-import-position
import matplotlib.patches as mp  # noqa # isort:skip
import matplotlib.pyplot as plt  # noqa # isort:skip

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Transparency alpha value to use
ALPHA = 0.5
#: Names of the colours to use in the diagrams
COLORS = ["blue", "red", "yellow", "green", "violet"]


class VennDiagramDescription:
    """Describe a Venn diagram in terms of shared elements"""

    def __init__(self, names=None, shared_counts=None):
        self.names = list(names or [])
        self.shared_counts = dict(shared_counts or {})

    def __str__(self):
        tpl = "VennDiagramDescription(names={}, shared_counts={})"
        return tpl.format(self.names, self.shared_counts)

    def __repr__(self):
        return str(self)


def powerset(iterable):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))


def load_description(f, separator):
    """Load set overlap description from ``f`` using ``separator``"""
    names = []
    shared_counts = []
    for line in f:
        if line.startswith("#"):
            continue
        arr = line.strip().split(separator)
        keys, num = arr[:-1], arr[-1]
        shared_counts.append((tuple(sorted(keys)), num))
        names += keys
    names = tuple(sorted(set(names)))
    return VennDiagramDescription(names, shared_counts)


def plot_venn_diagram_one(desc, args):
    """Create a Venn diagram for VennDiagramDescription having one set"""
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    e = mp.Ellipse(xy=(50, 50), width=80, height=80)
    e.set_clip_box(ax.bbox)
    e.set_alpha(ALPHA)
    e.set_facecolor(COLORS[0])
    ax.add_artist(e)

    key = desc.names[0]
    ax.annotate(
        desc.shared_counts.get(key, args.default_value),
        fontsize=16,
        xy=(50, 50),
        horizontalalignment="center",
        verticalalignment="center",
    )

    key = tuple()
    if desc.shared_counts.get(key):
        ax.annotate(
            desc.shared_counts[key],
            fontsize=16,
            xy=(10, 85),
            horizontalalignment="center",
            verticalalignment="center",
        )

    ax.annotate(
        desc.names[0],
        fontsize=16,
        xy=(50, 95),
        horizontalalignment="center",
        verticalalignment="center",
    )

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    return fig


def plot_venn_diagram_two(desc, args):
    """Create a Venn diagram for VennDiagramDescription having two sets"""
    ES = [
        mp.Ellipse(xy=(50, 50), width=80, height=80, angle=0),
        mp.Ellipse(xy=(90, 50), width=80, height=80, angle=0),
    ]
    VAL_COORDS = [(10, 85), (35, 50), (110, 50), (70, 50)]
    LABEL_COORDS = [(45, 95), (90, 95)]

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    for i, e in enumerate(ES):
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_alpha(ALPHA)
        e.set_facecolor(COLORS[i])

    for i, key in enumerate(powerset(desc.names)):
        ax.annotate(
            desc.shared_counts.get(key, args.default_value),
            fontsize=16,
            xy=VAL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    for i, c in enumerate(LABEL_COORDS):
        ax.annotate(
            desc.names[i],
            fontsize=16,
            xy=LABEL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    ax.set_xlim(10, 140)
    ax.set_ylim(0, 100)
    ax.axis("off")

    return fig


def plot_venn_diagram_three(desc, args):
    """Create a Venn diagram for VennDiagramDescription having three sets"""
    ES = [
        mp.Ellipse(xy=(50, 50), width=80, height=80, angle=0),
        mp.Ellipse(xy=(90, 50), width=80, height=80, angle=0),
        mp.Ellipse(xy=(75, 90), width=80, height=80, angle=0),
    ]
    VAL_COORDS = [(10, 85), (35, 45), (110, 45), (70, 100), (70, 40), (50, 75), (100, 75), (70, 65)]
    LABEL_COORDS = [(45, 5), (90, 5), (75, 135)]

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    for i, e in enumerate(ES):
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_alpha(ALPHA)
        e.set_facecolor(COLORS[i])

    for i, key in enumerate(powerset(desc.names)):
        ax.annotate(
            desc.shared_counts.get(key, args.default_value),
            fontsize=16,
            xy=VAL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    for i, c in enumerate(LABEL_COORDS):
        ax.annotate(
            desc.names[i],
            fontsize=16,
            xy=LABEL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    ax.set_xlim(10, 140)
    ax.set_ylim(0, 150)
    ax.axis("off")

    return fig


def plot_venn_diagram_four(desc, args):
    """Create a Venn diagram for VennDiagramDescription having four sets"""
    ES = [
        mp.Ellipse(xy=(50, 50), width=52, height=100, angle=45),  # A
        mp.Ellipse(xy=(73, 62), width=52, height=100, angle=45),  # B
        mp.Ellipse(xy=(73, 62), width=52, height=100, angle=-45),  # C
        mp.Ellipse(xy=(95, 50), width=52, height=100, angle=-45),  # D
    ]
    VAL_COORDS = [
        (20, 20),  # <null>
        (25, 65),  # A
        (50, 90),  # B
        (90, 90),  # C
        (120, 65),  # D
        (42.5, 78),  # AB
        (45, 39),  # AC
        (72.5, 20),  # AD,
        (72.5, 80),  # BC
        (102.5, 39),  # BD
        (103.5, 78),  # CD
        (57, 63),  # ABC
        (83.5, 30),  # ABD
        (62, 30),  # ACD
        (90, 63),  # BCD
        (71.5, 45),  # ABCD
    ]
    LABEL_COORDS = [(20, 100), (40, 110), (100, 110), (120, 100)]

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    for i, e in enumerate(ES):
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_alpha(ALPHA)
        e.set_facecolor(COLORS[i])

    for i, key in enumerate(powerset(desc.names)):
        ax.annotate(
            desc.shared_counts.get(key, args.default_value),
            fontsize=16,
            xy=VAL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    for i, c in enumerate(LABEL_COORDS):
        ax.annotate(
            desc.names[i],
            fontsize=16,
            xy=LABEL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    ax.set_xlim(10, 140)
    ax.set_ylim(10, 120)
    ax.axis("off")

    return fig


def plot_venn_diagram_five(desc, args):
    """Create a Venn diagram for VennDiagramDescription having five sets"""
    VAL_COORDS = [
        (-35, 40),  # <null>
        (4, 43),  # A
        (-40, 20),  # B
        (-28, -35),  # C
        (22, -40),  # D
        (40, 10),  # E
        (-12, 33),  # AB
        (13, 31),  # AC
        (10, -32),  # AD
        (29, 20),  # AE
        (-35, 0),  # BC
        (-28, 20),  # BD
        (35, -5),  # BE
        (-11, -34),  # CD
        (-30, -21),  # CE
        (30, -20),  # DE
        (5, 26),  # ABC
        (70, 70),  # ABD
        (28, 5),  # ABE
        (-3, -30),  # ACD
        (20, 23),  # ACE
        (14, -24),  # ADE
        (-24, 12),  # BCD
        (-32, -5),  # BCE
        (30, -14),  # BDE
        (-20, -19),  # CDE
        (-10, 20),  # ABCD
        (18, 15),  # ABCE
        (20, -8),  # ABDE
        (-5, -24),  # ACDE
        (-22, 3),  # BCDE
        (0, 0),  # ABCDE
    ]
    LABEL_COORDS = [(5, 60), (-55, 35), (-40, -50), (40, -50), (50, 30)]

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    for i in range(5):
        x = 8 * math.sin((i + 2) / 5 * 2 * math.pi)
        y = 8 * math.cos((i + 2) / 5 * 2 * math.pi)
        el = mp.Ellipse(xy=(x, -y), width=50, height=90, angle=i * 360 / 5)
        el.set_clip_box(ax.bbox)
        el.set_alpha(ALPHA)
        el.set_facecolor(COLORS[i])
        ax.add_artist(el)

    for i, key in enumerate(powerset(desc.names)):
        ax.annotate(
            desc.shared_counts.get(key, args.default_value),
            fontsize=16,
            xy=VAL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    for i, c in enumerate(LABEL_COORDS):
        ax.annotate(
            desc.names[i],
            fontsize=16,
            xy=LABEL_COORDS[i],
            horizontalalignment="center",
            verticalalignment="center",
        )

    ax.set_xlim(-55, 55)
    ax.set_ylim(-55, 65)
    ax.axis("off")

    return fig


def plot_venn_diagram(desc, args):
    """Perform plotting of the venn diagram plotting"""
    plotters = {
        1: plot_venn_diagram_one,
        2: plot_venn_diagram_two,
        3: plot_venn_diagram_three,
        4: plot_venn_diagram_four,
        5: plot_venn_diagram_five,
    }
    if len(desc.names) not in plotters:
        raise RuntimeError("I cannot create Venn diagram for {} sets".format(len(desc.names)))
    return plotters[len(desc.names)](desc, args)


def run(args):
    """Program entry point after parsing command line arguments"""
    desc = load_description(args.input_shared_counts, args.separator)
    fig = plot_venn_diagram(desc, args)
    fig.savefig(args.output_image, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    """Program entry point"""
    parser = argparse.ArgumentParser(description="Draw a quick Venn diagram")

    parser.add_argument(
        "--input-shared-counts",
        type=argparse.FileType("rt"),
        required=True,
        help="Set description file",
    )
    parser.add_argument(
        "--separator",
        type=str,
        default="\t",
        help="Separator to use in --input-shared-counts, defaults to TAB",
    )
    parser.add_argument(
        "--default-value", type=str, default="", help="value to use for empty entries"
    )
    parser.add_argument("--output-image", type=str, required=True, help="Path to output image file")

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
