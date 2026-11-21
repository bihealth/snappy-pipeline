"""Code shared between workflows"""


def vcf_outputs(prefix: str) -> dict[str, str]:
    """The output keys of a bgzipped VCF result file: the file, its tabix index and their md5 files.

    ``prefix`` is the path without the ``.vcf.gz`` suffix, such as ``work/{library_name}/out/{library_name}``.
    """
    return {
        "vcf": f"{prefix}.vcf.gz",
        "vcf_tbi": f"{prefix}.vcf.gz.tbi",
        "vcf_md5": f"{prefix}.vcf.gz.md5",
        "vcf_tbi_md5": f"{prefix}.vcf.gz.tbi.md5",
    }
