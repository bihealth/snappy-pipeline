# -*- coding: utf-8 -*-
"""Wrapper for merging multiple tables in R on shared columns"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

import pandas as pd
import re
import sys

pattern = re.compile("(.+)-[DR]NA[0-9]+-(WES|WGS|mRNA_seq)[0-9]+$")


def readTable(fn, dataint=False):
    df = pd.read_csv(fn, sep="\t", header=0)
    df = df.dropna()
    df = df.astype({"Entrez_Gene_Id": "str"})
    df = df.groupby(["Hugo_Symbol", "Entrez_Gene_Id"]).median()
    df = df.reset_index()
    # WARNING- pandas rounds towards 0, not to nearest int
    if dataint:
        df = df.astype({df.columns[2]: "int32"})
    return df


joined = None
for fn in snakemake.input:
    df = readTable(fn, dataint=(snakemake.params.datatype == "int"))
    if joined is None:
        joined = df
    else:
        joined = pd.merge(joined, df, on=["Hugo_Symbol", "Entrez_Gene_Id"], how="outer")

joined.columns = [pattern.sub("\\1", x) for x in joined.columns]

if snakemake.params.datatype == "int":
    joined.to_csv(str(snakemake.output), sep="\t", na_rep="NA", index=False, float_format="%.0f")
else:
    joined.to_csv(str(snakemake.output), sep="\t", na_rep="NA", index=False)
