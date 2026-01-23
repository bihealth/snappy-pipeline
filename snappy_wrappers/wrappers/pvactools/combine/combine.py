import argparse
import enum
import logging
import os
import re
import sys

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import vcfpy

VERSION = "0.1.0"

# vep description format for INFO/CSQ. The annotation names are in "titles", separated by ANNOTATION_SEPARATOR
ANNOTATION_DESCRIPTION_REGEX = "^Consequence annotations from Ensembl VEP. Format: (?P<titles>.+)$"
# Separator character use to split vep annotations and annotation names
ANNOTATION_SEPARATOR = r"\|"
# vep annotation name for gene ENSEMBL ids
ANNOTATION_GENE_ID = "Gene"
# vep annotation name for transcript ENSEMBL ids
ANNOTATION_TRANSCRIPT_ID = "Feature"

# ENSEMBL ids for human gene, transcripts & protein features
ENSEMBL_ID_PATTERN = re.compile(
    r"^(?P<id>ENS(?P<type>[GTP])[0-9]{11})((?P<version>\.[0-9]+)(?P<extra>.+)?)?$"
)

DEFAULT_VCF_IDS = {
    "annotation": "CSQ",  # Variant annotation produced by vep (INFO/CSQ)
    "rna_allele_depth": "AD",  # Allelic depth produced by bcftools (FORMAT/AD)
    "out_tpm_gene_id": "GX",  # Gene expression required by pVACseq (FORMAT/GX)
    "out_tpm_transcript_id": "TX",  # Transcript expression required by pVACseq (FORMAT/TX)
    "out_allele_depth_id": "RAD",  # RNA allelic depth required by pVACseq (FORMAT/RAD)
    "out_depth_id": "RDP",  # RNA total depth required by pVACseq (FORMAT/RDP)
}


class MissingDataException(Exception):
    pass


class CantOverwriteDataException(Exception):
    pass


class FeatureType(enum.StrEnum):
    """Character used in ENSEMBL ids to differentiate feature types (human only)"""

    Gene = "G"
    Transcript = "T"
    Protein = "P"
    Any = ""


@dataclass
class TPMFileFormat:
    """Description of generic TPM file format files"""

    feature_column: str | int = "Name"  # Feature column name (salmon)
    tpm_column: str | int = "TPM"  # TPM column name (salmon)
    sep: re.Pattern = re.compile(r"\t")  # Column separator character (usually either tab or comma)
    has_column_names: bool = True  # Has the file a header row with column names?
    ensembl_feature: bool = True  # Are the features following ENSEMBL id format?
    ignore_version: bool = False  # Is the feature version (*.[0-9]+) to be considered?


class TPM:
    """
    Class to read & store expression data (TPMs)

    The class is meant & store either gene or transcript-based TPMs,
    and retrieve the TPM value given a gene or transcript ID.
    If the feature ID is ENSEMBL, then feature id removal is under user control.
    This is necessary to ensure mapping between the feature ids provided in
    annotations (typically by vep) and those provided by the expression quantification
    software (typically salmon). By default, the former ids have no version, while
    the latter do have them. In that case, the format definition should set "ignore_version"
    to "True".
    """

    def __init__(self, feature_type: FeatureType):
        self.feature_type = feature_type

        self.filename: Optional[str | Path] = None
        self.tpm: dict[str, float] = {}

    def _process_feature(self, feature_id: str, fmt: TPMFileFormat) -> str:
        if fmt.ensembl_feature:
            m = ENSEMBL_ID_PATTERN.match(feature_id)
            assert m, f"Unexpected feature id '{feature_id}', not ENSEMBL format"
            assert self.feature_type == FeatureType.Any or self.feature_type == m.group("type"), (
                "Wrong feature type in {feature_id}, expected type {feature_type}"
            )
            if fmt.ignore_version:
                feature_id = m.group("id")
        return feature_id

    def read_file(self, filename: str | Path, fmt: TPMFileFormat) -> dict[str, float]:
        """Flexible function to extract TPM (any float value, really) for each feature"""
        self.filename = filename
        self.tpms = {}
        with open(self.filename, "rt") as f:
            # Read optional title line, and set the feature & value columns as column indices rather than names
            if fmt.has_column_names:
                titles = fmt.sep.split(f.readline().strip())
                if isinstance(fmt.feature_column, str):
                    try:
                        feature_column = titles.index(fmt.feature_column)
                    except ValueError:
                        raise MissingDataException(
                            f"Column {fmt.feature_column} missing from file {self.filename}"
                        )
                if isinstance(fmt.tpm_column, str):
                    try:
                        tpm_column = titles.index(fmt.tpm_column)
                    except ValueError:
                        raise MissingDataException(
                            f"Column {fmt.tpm_column} missing from file {self.filename}"
                        )
            else:
                feature_column = fmt.feature_column
                tpm_column = fmt.tpm_column

            # Loop over remainder of rows
            for line in f:
                tokens = fmt.sep.split(line.strip())
                try:
                    feature_id = tokens[feature_column]
                    value = tokens[tpm_column]
                except IndexError:
                    raise MissingDataException(
                        f"Can't extract feature or TPM from line {line.strip()} of file {self.filename}"
                    )

                # Post-process feature id (currently only for ENSEMBL)
                feature_id = self._process_feature(feature_id, fmt)

                assert feature_id not in self.tpms.keys(), (
                    f"Duplicate feature '{feature_id}' in file '{self.filename}'"
                )
                try:
                    self.tpms[feature_id] = float(value)
                except ValueError:
                    raise MissingDataException(
                        f"TPM value {value} not numeric for feature '{feature_id}'"
                    )

        logging.info(
            f"{len(self.tpms)} features read from file {self.filename}, summing to {sum(self.tpms.values())}"
        )
        return self.tpms

    def get_tpm_value(self, feature: str, rec_pos: str = "") -> str:
        """Extract TPM value from TPM dict for the query feature"""
        assert self.filename, "Can't extract TPM value before reading data"
        if feature:
            try:
                tpm = "{feature}|{value}".format(feature=feature, value=str(self.tpms.get(feature)))
            except KeyError:
                logging.warning(f"Unknown feature id {feature} in record {rec_pos}")
                tpm = "."
        else:
            logging.info(f"Empty feature in {rec_pos}")
            tpm = "."
        return tpm


def aggregate_counts(
    records: list[vcfpy.Record],
    sample: str,
    ref: str,
    alt: str,
    allele_depth: str = DEFAULT_VCF_IDS["rna_allele_depth"],
) -> list[int]:
    """
    Aggregates allelic counts over a list of records.

    The records all overlap one single position, which reference allele is in ref,
    and the target alt allele in alt.
    The counts are aggregated for one sample only, and the allelic depth FORMAT ID is in
    allele_depth.
    On output, the counts for the reference is in 0, for the target alt in 1 & for other alts in 2.
    """
    ed = [0, 0, 0]
    first = True
    for record in records:
        assert record.REF == ref, (
            f"Mismatch in reference ({ref} != {record.REF}) at {record.CHROM}:{record.POS}"
        )
        ad = record.call_for_sample.get(sample, None)
        if ad and allele_depth in ad.data:
            ad = ad.data[allele_depth]
            assert len(ad) == len(record.ALT) + 1, (
                f"Mismatch in {allele_depth} length at {record.CHROM}:{record.POS}"
            )
            if first:
                ed[0] = ad[0]
                first = False
            if len(record.ALT) > 0:
                for i in range(len(record.ALT)):
                    if record.ALT[i].serialize() == alt:
                        assert ed[1] == 0, (
                            f"Alt allele {alt} appears multiple times in {record.CHROM}:{record.POS}"
                        )
                        ed[1] = ad[i + 1]
                    else:
                        ed[2] += ad[i + 1]
    return ed


def get_annotation_titles(
    dna: vcfpy.Reader,
    annotation_id: str = DEFAULT_VCF_IDS["annotation"],
    annotation_regex: re.Pattern = re.compile(ANNOTATION_DESCRIPTION_REGEX),
    annotation_separator: re.Pattern = re.compile(ANNOTATION_SEPARATOR),
) -> list[str]:
    """Returns the annotation labels of the annotation INFO description (CSQ for vep, ANN for jannovar)"""
    if annotation_id not in dna.header.info_ids():
        raise MissingDataException(
            f"No annotation {annotation_id} in vcf {dna.path}, the TPM values cannot be added"
        )
    annotation_titles = dna.header.get_info_field_info(annotation_id).description
    m = annotation_regex.match(annotation_titles)
    if not m:
        raise MissingDataException(
            f"Annotation description '{annotation_titles}' cannot be parsed using regular expression {annotation_regex.pattern}"
        )
    annotation_titles = annotation_separator.split(m.group("titles"))
    return annotation_titles


def open_vcf_with_requests(
    vcf_fn: str | Path, infos: list[str] = [], fmts: list[str] = [], samples: list[str] = []
):
    """Open a vcf file with tabix, ensuring presence of required INFO, FORMAT & samples"""
    tabix_path = None
    for ext in (".csi", ".tbi"):
        if os.path.exists(vcf_fn + ext):
            tabix_path = vcf_fn + ext
            break
    if not tabix_path:
        raise MissingDataException(
            f"Mandatory index of vcf file {vcf_fn} is missing, please run tabix on the file"
        )

    vcf = vcfpy.Reader.from_path(vcf_fn, tabix_path=tabix_path)

    missing = set(infos) - set(vcf.header.info_ids())
    if len(missing) > 0:
        raise MissingDataException(f"Missing requested INFO ID(s) {missing} in vcf file {vcf_fn}")
    missing = set(fmts) - set(vcf.header.format_ids())
    if len(missing) > 0:
        raise MissingDataException(f"Missing requested FORMAT ID(s) {missing} in vcf file {vcf_fn}")
    missing = set(samples) - set(vcf.header.samples.names)
    if len(missing) > 0:
        raise MissingDataException(f"Missing requested samples {missing} in vcf file {vcf_fn}")

    return vcf


def add_header_lines(
    dna: vcfpy.Reader,
    gene_tpms: Optional[TPM],
    tx_tpms: Optional[TPM],
    rna_fn: Optional[str | Path],
    add_symbolic: bool = True,
    out_tpm_gene_id: str = DEFAULT_VCF_IDS["out_tpm_gene_id"],
    out_tpm_tx_id: str = DEFAULT_VCF_IDS["out_tpm_transcript_id"],
    out_depth_id: str = DEFAULT_VCF_IDS["out_depth_id"],
    out_allele_depth_id: str = DEFAULT_VCF_IDS["out_allele_depth_id"],
):
    """Add required lines to vcf header (program, and if available FORMAT for pileup & expression)"""
    dna.header.add_line(vcfpy.HeaderLine("combineVersion", VERSION))

    dna.header.add_line(
        vcfpy.HeaderLine(
            "combineCommand",
            "{gene_tpms} {tx_tpms} {pileups} {dna_fn}; Date={today}".format(
                gene_tpms=f"--gene-tpms {gene_tpms.filename}" if gene_tpms else "",
                tx_tpms=f"--transcript-tpms {tx_tpms.filename}" if tx_tpms else "",
                pileups=f"--pileup {rna_fn}" if rna_fn else "",
                dna_fn=dna.path if dna.path else "-",
                today=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            ),
        )
    )

    if gene_tpms:
        if out_tpm_gene_id in dna.header.format_ids():
            logging.info(f"FORMAT ID {out_tpm_gene_id} already present in vcf")
        else:
            dna.header.add_format_line(
                {
                    "ID": out_tpm_gene_id,
                    "Number": ".",
                    "Type": "String",
                    "Description": "Gene-based TPM",
                }
            )

    if tx_tpms:
        if out_tpm_tx_id in dna.header.format_ids():
            logging.info(f"FORMAT ID {out_tpm_tx_id} already present in vcf")
        else:
            dna.header.add_format_line(
                {
                    "ID": out_tpm_tx_id,
                    "Number": ".",
                    "Type": "String",
                    "Description": "Transcript-based TPM",
                }
            )

    if rna_fn:
        if out_depth_id in dna.header.format_ids():
            logging.info(f"FORMAT ID {out_depth_id} already present in vcf")
        else:
            dna.header.add_format_line(
                {
                    "ID": out_depth_id,
                    "Number": 1,
                    "Type": "Integer",
                    "Description": "Read depth in expression",
                }
            )
        if out_allele_depth_id in dna.header.format_ids():
            logging.info(f"FORMAT ID {out_allele_depth_id} already present in vcf")
        else:
            dna.header.add_format_line(
                {
                    "ID": out_allele_depth_id,
                    "Number": 3 if add_symbolic else "R",
                    "Type": "Integer",
                    "Description": "Number of reads supporting the reference, the alternative allele & other alleles",
                }
            )


def annotation_to_dict(
    info: str,
    annotations: list[str],
    annotation_separator: re.Pattern = re.compile(ANNOTATION_SEPARATOR),
) -> dict[str, str]:
    """Split annotations from INFO & annotation name & value dict"""
    ann = annotation_separator.split(info)
    assert len(annotations) == len(ann), f"Annotation length mismatch for annotation '{info}'"
    return dict(zip(annotations, ann))


def run(
    sample_id: str,
    dna_fn: Optional[str | Path],
    combined_fn: Optional[str | Path],
    gene_tpms: Optional[TPM],
    transcript_tpms: Optional[TPM],
    rna_fn: Optional[str | Path],
    vcf_ids: dict[str, str] = DEFAULT_VCF_IDS,
    annotation_names: dict[str, str] = {"gene": "Gene", "transcript": "Feature"},
    annotation_regex: re.Pattern = re.compile(ANNOTATION_DESCRIPTION_REGEX),
    annotation_separator: re.Pattern = re.compile(ANNOTATION_SEPARATOR),
    add_symbolic: bool = False,
) -> int:
    """
    Merge a vcf containing variants with values from RNA (TPM & pileups)

    - The input vcf must be normalised, and pileup values are only added for SNVs.
    - The input vcf must be annotated so that the TPM value can be extracted for the correct feature.
    - Multiple samples are allowed in the vcfs (somatic variants in dna_fn & rna pileups in rna_fn),
      but the merging is done for one sample only (sample_id), which must obviously be present in
      the somatic variant vcf file, and in the pileup file when requested.
    - All vcf & annotation identifiers are under user control, with resonable defaults.

    Args:
        sample_id -- sample id, must be present in vcf files
        dna_fn -- somatic variant vcf file path, or use stdin when missing
        combined_fn -- output vcf file path, or use stdout when missing
        gene_tpms -- optional TPM object containing gene ids to tpm expression values mappings
        transcript_tpms -- optional TPM object containing transcript ids to tpm expression values mappings
        rna_fn -- optional RNA allelic depth vcf file path
        vcf_ids -- INFO & FORMAT vcf ids for all the fields (with reasonable defaults)
            - annotation (default "CSQ"): somatic variant vcf INFO field containing variant annotations
            - rna_allele_depth (default "AD"): rna pileup vcf FORMAT field containing
              allelic depths in RNA data at somatic variant loci.
            - out_tpm_gene_id (default "GX"): output vcf FORMAT field containing TPM value for annotated gene
            - out_tpm_transcript_id (default "TX"): output vcf FORMAT field containing TPM value for annotated transcript
            - out_allele_depth_id (default "RAD"): rna allelic depths in output vcf FORMAT field
            - out_depth_id (default "RDP"): total RNA depth in output vcf FORMAT field
        annotation_names -- gene id & transcript id names in annotation
        annotation_regexp -- regular expression used to parse annotation names from INFO description
        annotation_separator -- regular expression used to split annotation nales and annotation values
            in INFO descrition and in INFO fields
        add_symbolic -- output allelic depth of SNV variants other than the reference allele or
            the alt allele described in the somatic variant file

    Returns:
        0 for success & 1 for failure
    """
    if dna_fn and dna_fn != "-":
        dna = vcfpy.Reader.from_path(dna_fn)
    else:
        dna = vcfpy.Reader.from_stream(sys.stdin)

    if sample_id not in dna.header.samples.names:
        logging.error(f"Sample {sample_id} not in vcf file {dna.path}")
        return -1

    # Ensure that the CSQ annotation contains Gene & Feature
    if gene_tpms or transcript_tpms:
        annotation_titles = get_annotation_titles(
            dna, vcf_ids["annotation"], annotation_regex, annotation_separator
        )
        if gene_tpms and annotation_names["gene"] not in annotation_titles:
            logging.error(
                f"Missing {annotation_names['gene']} annotation in file {dna.path} INFO/{vcf_ids['annotation']}"
            )
            return 1
        if transcript_tpms and annotation_names["transcript"] not in annotation_titles:
            logging.error(
                f"Missing {annotation_names['transcript']} annotation in file {dna.path} INFO/{vcf_ids['annotation']}"
            )
            return 1

    if rna_fn:
        try:
            rna = open_vcf_with_requests(
                rna_fn, fmts=[vcf_ids["rna_allele_depth"]], samples=[sample_id]
            )
        except MissingDataException as e:
            logging.error(e)
            return 1

    add_header_lines(
        dna,
        gene_tpms=gene_tpms,
        tx_tpms=transcript_tpms,
        rna_fn=rna_fn if rna_fn else None,
        add_symbolic=add_symbolic,
        out_tpm_gene_id=vcf_ids["out_tpm_gene_id"],
        out_tpm_tx_id=vcf_ids["out_tpm_transcript_id"],
        out_allele_depth_id=vcf_ids["out_allele_depth_id"],
    )

    if combined_fn:
        writer = vcfpy.Writer.from_path(combined_fn, dna.header)
    else:
        writer = vcfpy.Writer.from_stream(sys.stdout, dna.header)

    iRecord = 0
    nTPM = 0
    nPileup = 0
    for record in dna:
        iRecord += 1

        if len(record.ALT) != 1:
            logging.error(
                f"Multiple alternative alleles in record {iRecord}. Not allowed, the vcf must be normalised first"
            )
            return 1

        ref = record.REF
        alt = record.ALT[0].serialize()
        rec_pos = f"{record.CHROM}:{record.POS}_{ref}/{alt}"
        data = record.call_for_sample[sample_id].data

        if gene_tpms or transcript_tpms:
            csqs = [
                annotation_to_dict(a, annotation_titles, annotation_separator)
                for a in record.INFO[vcf_ids["annotation"]]
            ]
            if gene_tpms:
                record.add_format(vcf_ids["out_tpm_gene_id"], [])
                data[vcf_ids["out_tpm_gene_id"]] = [
                    gene_tpms.get_tpm_value(csq[annotation_names["gene"]], rec_pos) for csq in csqs
                ]
            if transcript_tpms:
                record.add_format(vcf_ids["out_tpm_transcript_id"], [])
                data[vcf_ids["out_tpm_transcript_id"]] = [
                    transcript_tpms.get_tpm_value(csq[annotation_names["transcript"]], rec_pos)
                    for csq in csqs
                ]
            nTPM += len(csqs)

        if rna_fn and record.is_snv():
            try:
                expr = list(
                    filter(
                        lambda r: r.is_snv(),
                        rna.fetch(record.CHROM, record.POS - 1, record.POS),
                    )
                )
                eds = aggregate_counts(
                    expr, sample_id, ref, alt, allele_depth=vcf_ids["rna_allele_depth"]
                )
                record.add_format(vcf_ids["out_depth_id"])
                record.add_format(vcf_ids["out_allele_depth_id"], [])
                data[vcf_ids["out_depth_id"]] = sum(eds)
                if not add_symbolic:
                    eds.pop()
                data[vcf_ids["out_allele_depth_id"]] = eds
                nPileup += 1

            except ValueError:
                # No expression at this position -> rna.fetch raises a ValueError
                pass

        writer.write_record(record)

    writer.close()
    if rna_fn:
        rna.close()
    dna.close()

    logging.info(f"{iRecord} variants have been processed")
    if gene_tpms or transcript_tpms:
        logging.info(f"TPM values have been added for {nTPM} features")
    if rna_fn:
        logging.info(f"RNA allelic depth has been added to {nPileup} variants")

    return 0


def _get_expression_format(args: argparse.Namespace, feature_type: FeatureType) -> TPMFileFormat:
    """Parse command line argument to create a TPM file format definition"""
    if args.format:
        if args.header_line or args.column_separator or args.feature_column or args.tpm_column:
            raise argparse.ArgumentError(
                "--format option incompatible with other expression file format arguments"
            )
        match (args.format, feature_type):
            case ("salmon", _):
                feature_column = "Name"
                tpm_column = "TPM"
                sep = r"\t"
                has_column_name = True
            case ("kallisto", "gene"):
                feature_column = "gene"
                tpm_column = "abundance"
                sep = r"\t"
                has_column_name = True
            case ("kallisto", "transcript"):
                feature_column = "target_id"
                tpm_column = "tpm"
                sep = r"\t"
                has_column_name = True
            case ("cufflinks", _):
                feature_column = "tracking_id"
                tpm_column = "FPKM"
                sep = r"\t"
                has_column_name = True
            case ("stringtie", "gene"):
                feature_column = "Gene ID"
                tpm_column = "TPM"
                sep = r"\t"
                has_column_name = True
            case ("stringtie", "transcript"):
                feature_column = "transcript_id"
                tpm_column = "TPM"
                sep = r"\t"
                has_column_name = True
            case ("unimplemented", _):
                raise argparse.ArgumentError("--format 'unimplemented' is unimplemented")
            case _:
                raise argparse.ArgumentError(f"--format '{args.format}' is unimplemented")
    else:
        if not args.header_line:
            try:
                feature_column = int(args.feature_column)
                tpm_column = int(args.tpm_column)
            except ValueError:
                raise argparse.ArgumentError(
                    "In absence of header line, both feature & TPM columns must be positive integer"
                )
            if feature_column < 0 or tpm_column < 0:
                raise argparse.ArgumentError(
                    "In absence of header line, both feature & TPM columns must be positive integer"
                )
        sep = args.column_separator
        has_column_name = args.header_line

    if not args.ensembl_id and args.use_ensembl_version:
        raise argparse.ArgumentError("--ensemb-id must be enabled to use --use-ensembl-version")

    return TPMFileFormat(
        feature_column=feature_column,
        tpm_column=tpm_column,
        sep=re.compile(sep),
        has_column_names=has_column_name,
        ensembl_feature=args.ensembl_id,
        ignore_version=not args.use_ensembl_version,
    )


def main() -> int:
    """
    Script to combine somatic variant vcf with expression data (TPMs) & RNA allelic depth.

    - The addition of expression and allelic depth is optional (but one of them should be present).
    - The expression values should be stored in a tab- or comma-delimited file, containing
      feature ids in one column, and TPM values in another. Column selection is under user control.
      When the file has a header row, then the user can select the columns by their name, otherwise
      column selection is done by numbers (starting from 0).
    - The expression file format definition is hard-coded in the script for common expression
      quantification software (currently only salmon).
    - Expression quantification can be provided for genes or transcripts or both.
    - When requested, allelic depth is attached to FORMAT records only for somatic SNVs.
    - The addition of expression and/or allelic depth is done for a single sample.
      This sample must be persent is the somatic variant and (when requested) allelic depth vcfs.
    - The input somatic variant must be annotated if expression values are requested.
    - The input somatic variant must be normalised (no multiple alt alleles), but
      multiple annotations are allowed. The TPM values are output as "<feature id>|<TPM value>".
    """
    parser = argparse.ArgumentParser(
        prog="combine",
        description="combine a vcf containing (somatic) variants with pileups of mRNA expression at the same loci, and TPM value from salmon",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Increase verbosity level")
    parser.add_argument("--version", action="version", version="%(prog)s " + VERSION)

    parser.add_argument("-s", "--sample", help="Sample ID")
    parser.add_argument(
        "-g",
        "--gene-tpms",
        help="Gene expression file (*.gene.sf), must contain the ENSEMBL feature ID in column 'Name', and the TPM value in column 'TPM'",
    )
    parser.add_argument(
        "-t",
        "--transcript-tpms",
        help="Transcript expression file (*.transcript.sf), must contain the ENSEMBL feature ID in column 'Name', and the TPM value in column 'TPM'",
    )
    parser.add_argument(
        "-p",
        "--pileup",
        help="Vcf file produced by bcftools mpileup at the position of somatic variants on the RNA expression data mapping by STAR",
    )
    parser.add_argument(
        "-a",
        "--add_symbolic",
        action="store_true",
        help="Add symbolic allele (neither ref nor alt) to pileup results",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Normalised vcf file with somatic variants including tpm & pileup values if available",
    )

    for k, v in DEFAULT_VCF_IDS.items():
        parser.add_argument(
            "--{key}".format(key=k.replace("_", "-")),
            default=v,
            help=f"VCF ID (can be INFO, FORMAT or annotation), default {v}",
        )

    parser.add_argument("--ensembl-id", action="store_true", help="Feature IDs are ENSEMBL IDs")
    parser.add_argument(
        "--use-ensembl-version",
        action="store_true",
        help="Use ENSEMBL feature version when mapping features between expression & vcf annotations",
    )

    group = parser.add_argument_group("Expression", "Definition of the expression file format")
    group.add_argument(
        "-f",
        "--format",
        choices=("salmon", "kallisto", "cufflinks", "stringtie", "unimplemented"),
        help="Expression file format (from the program that created it)",
    )
    group.add_argument(
        "--header-line", action="store_true", help="File contains a header line with column titles"
    )
    group.add_argument(
        "--feature-column", help="Name or index of the column containing feature ids (0-based)"
    )
    group.add_argument(
        "--tpm-column", help="Name or index of the column containing TPM values (0-based)"
    )
    group.add_argument("--column-separator", help="Column separator")

    parser.add_argument(
        "--annotation-description-regex",
        default=ANNOTATION_DESCRIPTION_REGEX,
        help="Regular expression used to parse annotation description. Field names must be in group named 'titles'",
    )
    parser.add_argument(
        "--annotation-separator",
        default=ANNOTATION_SEPARATOR,
        help="Character separating annotations",
    )
    parser.add_argument(
        "--annotation-gene-id",
        default=ANNOTATION_GENE_ID,
        help="Name of the gene ID in the annotation description",
    )
    parser.add_argument(
        "--annotation-transcript-id",
        default=ANNOTATION_TRANSCRIPT_ID,
        help="Name of the transcript ID in the annotation description",
    )

    parser.add_argument(
        "variants",
        nargs="?",
        default=None,
        help="Normalised vcf file with somatic variants",
    )

    args = parser.parse_args()
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(msg)s",
        level=logging.DEBUG if args.verbose else logging.WARNING,
    )

    if not (args.gene_tpms or args.transcript_tpms or args.pileup):
        logging.error("Nothing to combine DNA data with")
        return -1

    if args.gene_tpms:
        expression_format = _get_expression_format(args, FeatureType.Gene)
        gene_tpms = TPM(feature_type=FeatureType.Gene)
        gene_tpms.read_file(args.gene_tpms, fmt=expression_format)
    else:
        gene_tpms = None

    if args.transcript_tpms:
        expression_format = _get_expression_format(args, FeatureType.Transcript)
        transcript_tpms = TPM(feature_type=FeatureType.Transcript)
        transcript_tpms.read_file(args.transcript_tpms, fmt=expression_format)
    else:
        transcript_tpms = None

    return run(
        args.sample,
        args.variants,
        args.output,
        gene_tpms,
        transcript_tpms,
        args.pileup,
        vcf_ids={k: getattr(args, k) for k in DEFAULT_VCF_IDS.keys()},
        annotation_names={
            "gene": args.annotation_gene_id,
            "transcript": args.annotation_transcript_id,
        },
        annotation_regex=re.compile(args.annotation_description_regex),
        annotation_separator=re.compile(args.annotation_separator),
        add_symbolic=args.add_symbolic,
    )


if __name__ == "__main__":
    sys.exit(main())
