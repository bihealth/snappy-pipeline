from functools import cache

import pandas as pd


@cache
def samplesheet(pipeline: str) -> pd.DataFrame:
    return pd.read_csv(
        config["pipeline-configuration"][pipeline]["samplesheet"],
        sep="\t",
        comment="#",
        skiprows=12,
    )


def samples(pipeline: str) -> list[str]:
    return list(samplesheet(pipeline)["sampleName"].unique())


def folder_names(pipeline: str) -> list[str]:
    return list(samplesheet(pipeline)["folderName"].unique())


def reference_path() -> str:
    if chrom := config["reference"].get("chromosome", None):
        return f"resources/refs/{chrom}.fa.gz"
    else:
        return "resources/refs/genome.fa.gz"


def reference_faidx_region_string(wildcards) -> str:
    chrom = config["reference"].get("chromosome", None)
    region = config["reference"].get("region", None)
    match chrom, region:
        case None, None:
            return f""
        case chrom, None:
            return f"{chrom}"
        case chrom, region:
            return f"{chrom}:{region}"
