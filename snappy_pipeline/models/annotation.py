import enum
from typing import Literal

from pydantic import model_validator

from snappy_pipeline.models import ConfigPath, SnappyModel, ToolModel


class VepTxFlag(enum.StrEnum):
    gencode_basic = "gencode_basic"
    refseq = "refseq"
    merged = "merged"


class VepPlugin(SnappyModel):
    name: str
    """The plugin's name, as for VEP's ``--plugin``"""

    path: ConfigPath
    """The plugin's ``.pm`` file, e.g. ``Frameshift.pm`` and ``Wildtype.pm`` of pVACtools"""


class Vep(ToolModel):
    species: str = "homo_sapiens"

    assembly: str = "GRCh38"

    cache_version: str = "115"
    """WARNING- this must match the wrapper's vep version!"""

    tx_flag: VepTxFlag = VepTxFlag.gencode_basic
    """The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged"."""

    pick: bool = True
    """Annotate one consequence per variant, chosen by ``pick_order`` (VEP's ``--pick``);
    ``false`` annotates all transcripts, for example for neoepitope prediction"""

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
    """VEP flags (without ``--``); ``no_stats`` switches off the summary report"""
    plugins: list[VepPlugin] = []
    """VEP plugins to run"""

    @model_validator(mode="after")
    def stats_only_through_output_options(self):
        """The summary is an output of the rule, which only ``output_options`` switches off"""
        if any(a.split()[0] in ("--no_stats", "--stats_file") for a in self.extra_arguments if a):
            raise ValueError("use output_options: [no_stats] instead of --no_stats or --stats_file")
        return self


#: The tags of a transcript that mehari picks or restricts to
MehariTranscriptTag = Literal[
    "mane-select",
    "mane-select-backport",
    "mane-plus-clinical",
    "mane-plus-clinical-backport",
    "refseq-select",
    "refseq-select-backport",
    "ensembl-canonical",
    "ensembl-canonical-backport",
    "gencode-primary",
    "gencode-primary-backport",
    "basic",
    "basic-backport",
]

#: A criterion of mehari's transcript ranking: a tag, or a property of the transcript
MehariPickCriterion = Literal[
    "mane-select",
    "mane-select-backport",
    "mane-plus-clinical",
    "mane-plus-clinical-backport",
    "refseq-select",
    "refseq-select-backport",
    "ensembl-canonical",
    "ensembl-canonical-backport",
    "gencode-primary",
    "gencode-primary-backport",
    "basic",
    "basic-backport",
    "coding",
    "unflagged",
    "severity",
    "length",
]


class Mehari(ToolModel):
    """``mehari annotate seqvars`` (mehari 0.51); every field is a command line option, unset
    fields leave mehari's default"""

    threads: int = 1
    """Number of threads to use for annotation; mehari calls 16 a good value"""

    assembly: str = "GRCh38"
    """Assembly to use"""

    pick: Literal["gene", "transcript", "allele"] | None = None
    """Report only one annotation per gene, transcript or allele, chosen by ``pick_order``
    (``--pick``); default: report every transcript"""

    pick_order: list[MehariPickCriterion] = []
    """How to rank transcripts for ``pick``, ``restrict_fallback`` and the row of a gene that a
    change contains whole (``--pick-order``); each criterion only breaks the ties left by the ones
    before it. Empty for mehari's default: mane-select, mane-plus-clinical, refseq-select,
    ensembl-canonical, gencode-primary, basic, coding, unflagged, severity, length"""

    restrict_transcripts: list[MehariTranscriptTag] = []
    """Report only transcripts that have one of these tags (``--restrict-transcripts``); a gene
    without such a transcript gets no annotation, unless ``restrict_fallback`` is set"""

    restrict_fallback: bool | None = None
    """Report the top-ranked transcript of a gene that has none of the ``restrict_transcripts``
    tags, marked with the message ``INFO_RESTRICT_FALLBACK``"""

    flag_pick: bool | None = None
    """Report all annotations and mark the picked ones in the ``pick`` column, instead of
    dropping the others"""

    keep_intergenic: bool | None = None
    """Whether to keep intergenic variants"""

    discard_utr_splice_variants: bool | None = None
    """Drop the splice terms of variants in UTRs"""

    hgvs_parentheses: Literal["predicted", "never"] | None = None
    """Parentheses around predicted changes, e.g. ``p.(Gly2Arg)``, as HGVS recommends:
    ``predicted`` (every ``p.`` except ``p.?`` and ``p.Met1?``, mehari's default) or ``never``"""

    in_memory_reference: bool | None = None
    """Read the reference genome into memory"""

    plugins: list[
        Literal["tx-sequence-ref", "tx-sequence-alt", "aa-sequence-ref", "aa-sequence-alt", "nmd"]
    ] = []
    """Plugins that add ``ANN`` columns (``--plugin``): the reference or alternative transcript
    and amino acid sequences, and nonsense-mediated decay"""

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

    compound_window_size: int | None = None
    """Experimental: Maximum distance in bp between neighbouring variants of a compound group;
    unlimited if not set, a good value for peptide analyses is 33"""

    @model_validator(mode="after")
    def check_option_dependencies(self):
        """The dependencies that mehari's command line enforces, checked when the project loads."""
        if self.restrict_fallback and not self.restrict_transcripts:
            raise ValueError("restrict_fallback needs restrict_transcripts")
        if self.flag_pick and not self.pick:
            raise ValueError("flag_pick needs pick")
        return self
