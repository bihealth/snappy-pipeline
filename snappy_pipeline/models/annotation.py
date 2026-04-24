import enum
from typing import Literal

from pydantic import Field

from snappy_pipeline.models import SnappyModel


class VepTxFlag(enum.StrEnum):
    gencode_basic = "gencode_basic"
    refseq = "refseq"
    merged = "merged"


class Vep(SnappyModel):
    cache_dir: str = ""
    """Defaults to $HOME/.vep Not a good idea on the cluster"""

    species: str = "homo_sapiens"

    assembly: str = "GRCh38"

    cache_version: str = "102"
    """WARNING- this must match the wrapper's vep version!"""

    tx_flag: VepTxFlag = VepTxFlag.gencode_basic
    """The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged"."""

    pick_order: list[str] = [
        "biotype",
        "mane_select",
        "mane_plus_clinical",
        "appris",
        "tsl",
        "ccds",
        "canonical",
        "rank",
        "length",
    ]
    num_threads: int = 8
    buffer_size: int = 1000
    output_options: list[str] = ["everything"]


class Mehari(SnappyModel):
    threads: int = 1
    """
    Number of threads to use for annotation.
    A sweet spot regarding trade-off between I/O and CPU-bound tasks is around 5.
    """

    assembly: str = "GRCh38"
    """Assembly to use"""

    reference: str
    """Reference genome FASTA file (with accompanying index)"""

    transcripts: list[str] = Field(default_factory=list, min_length=1)
    """Transcript database(s) containing the transcript information."""

    frequencies: str | None = None
    """Frequency database."""

    clinvar: str | None = None
    """ClinVar database."""

    report_most_severe_consequence_by: Literal["gene", "transcript", "allele"] | None = None

    pick_transcript: (
        Literal[
            "mane-select",
            "mane-select-backport",
            "mane-plus-clinical",
            "mane-plus-clinical-backport",
            "length",
            "ensembl-canonical",
            "ensembl-canonical-backport",
            "ref-seq-select",
            "ref-seq-select-backport",
            "gencode-primary",
            "gencode-primary-backport",
            "basic",
            "basic-backport",
        ]
        | None
    ) = None
    """
    Which kind of transcript to pick / restrict to. Default is not to pick at all.
    Depending on `--pick-transcript-mode`, if multiple transcripts match the selection, either the first one is kept or all are kept.
    """

    pick_transcript_mode: Literal["first", "all"] | None = None
    """
    Determines how to handle multiple transcripts. Default is to keep all.
    When transcript picking is enabled via `--pick-transcript`, either keep the first one found or keep all that match.
    """

    keep_intergenic: bool | None = None
    """
    Whether to keep intergenic variants
    """

    discard_utr_splice_variants: bool | None = None
    """
    Whether to report splice variants in UTRs
    """

    in_memory_reference: bool | None = None
    """Read the reference genome into memory"""

    report_cdna_sequence: Literal["none", "reference", "alternative", "both"] | None = None
    """Whether to report cDNA sequence"""

    report_protein_sequence: Literal["none", "reference", "alternative", "both"] | None = None
    """Whether to report protein sequence"""

    enable_compound_variants: bool | None = None
    """
    Experimental: Enable variant grouping to evaluate the compound effect of multiple variants on the same transcript.
    When disabled, Mehari evaluates each variant independently
    """

    phasing_strategy: Literal["strict", "relaxed", "ignore"] | None = None
    """
    Experimental: The strategy used to evaluate grouped variants for compound effects

    Possible values:
    - strict:  Variants are only grouped if explicitly phased ('|') and sharing a Phase Set (PS). Unphased variants are evaluated independently
    - relaxed: Respects explicit phasing, but treats homozygous variants as universally phased Unphased heterozygous variants remain independent
    - ignore:  Completely ignores phasing metadata and _assumes_ all variants on the transcript are on the same haplotype
    """
