import enum

from snappy_pipeline.models import SnappyModel


class VepTxFlag(enum.StrEnum):
    gencode_basic = "gencode_basic"
    refseq = "refseq"
    merged = "merged"


class VepPickOrder(enum.StrEnum):
    biotype = "biotype"
    """Biotype of transcript (protein_coding preferred)"""
    mane = "mane"
    """MANE transcript status"""
    appris = "appris"
    """APPRIS isoform annotation"""
    tsl = "tsl"
    """Transcript support level"""
    ccds = "ccds"
    """CCDS status of transcript"""
    canonical = "canonical"
    """Canonical status of transcript"""
    rank = "rank"
    """Consequence rank"""
    length = "length"
    """Translated, transcript or feature length (longer preferred)"""
    mane_select = "mane_select"
    """MANE Select status (available from version 103)"""
    mane_plus_clinical = "mane_plus_clinical"
    """MANE Plus Clinical transcript status (available from version 103)"""
    ensembl = "ensembl"
    """Undocumented (not available in 102)"""
    refseq = "refseq"
    """Undocumented (not available in 102)"""


class VepOutputOptions(enum.StrEnum):
     everything = "everything"
     """sift, polyphen, ccds, uniprot, hgvs, symbol, numbers, domains"""
     sift = "sift b"
     polyphen = "polyphen b"
     ccds = "ccds"
     uniprot = "uniprot"
     hgvs = "hgvs"
     symbol = "symbol"
     numbers = "numbers"
     domains = "domains"
     regulatory = "regulatory"
     canonical = "canonical"
     protein = "protein"
     biotype = "biotype"
     tsl = "tsl"
     appris = "appris"
     gene_phenotype = "gene_phenotype"
     af = "af"
     af_1kg = "af_1kg"
     af_esp = "af_esp"
     max_af = "max_af"
     pubmed = "pubmed"
     mane = "mane"
     variant_class = "variant_class"
     var_synonyms = "var_synonyms"
     """Removed since version 106"""
     af_gnomad = "af_gnomad"
     """
     Superseded by 'af_gnomade' & 'af_gnomadg' for version 107.
     'af_gnomad' has the same function as 'af_gnomade'
     """
     af_gnomade = "af_gnomade"
     af_gnomadg = "af_gnomadg"
     mirna = "mirna"
     """Available from version 109"""


class Vep(SnappyModel):
    cache_dir: str = ""
    """Defaults to $HOME/.vep Not a good idea on the cluster, because the database is very large."""

    species: str = "homo_sapiens"

    assembly: str = "GRCh38"

    cache_version: str = "102"
    """WARNING- this must match the wrapper's vep version!"""

    tx_flag: VepTxFlag = VepTxFlag.gencode_basic
    """The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged"."""

    pick_order: list[VepPickOrder] = [
        VepPickOrder.biotype,
        VepPickOrder.mane,
        VepPickOrder.mane_plus_clinical,
        VepPickOrder.appris,
        VepPickOrder.tsl,
        VepPickOrder.ccds,
        VepPickOrder.canonical,
        VepPickOrder.rank,
        VepPickOrder.length,
    ]
    """
    Ranking of transcripts returned by VEP.
    Important when only one transcript is selected for annotation.
    The default order is different from the default order proposed by VEP.
    Here, variants in protein coding transcripts will be annotated before MANE transcripts.
    """

    num_threads: int = 8
    buffer_size: int = 1000
    output_options: list[VepOutputOptions] = [VepOutputOptions.everything]
