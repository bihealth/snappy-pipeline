import argparse
import logging
import re
import sys

import pandas as pd
import pyranges as pr
from vcfpy import OrderedDict, Reader, Writer


########################## Define tool parameters ##########################
def define_parser():
    parser = argparse.ArgumentParser(
        "comb_rna.py",
        description="A tool that will add the data from several expression tools' output files"
        + "and allelic depths from mRNA sequencing data for snvs"
        + "to the VCF INFO column. Supported tools are StringTie, Kallisto, Star, snappy_custom"
        + "and Cufflinks. There also is a ``custom`` option to annotate with data "
        + "from any tab-delimited file.",
    )

    parser.add_argument("--input-vcf", help="A VEP-annotated VCF file")
    parser.add_argument("--expression-file", help="A expression file")
    parser.add_argument("--genecode", help="A genecode file for calculate TPM from star gene count")
    parser.add_argument(
        "--rna-vcf",
        help="A VCF file of RNA-seq data",
    )
    parser.add_argument(
        "--format",
        choices=["kallisto", "stringtie", "cufflinks", "star", "snappy_custom", "custom"],
        help="The file format of the expression file to process. "
        + "Use `custom` to process file formats not explicitly supported. "
        + "The `custom` option requires the use of the --id-column and --expression-column arguments.",
    )
    parser.add_argument(
        "--mode",
        choices=["gene", "transcript"],
        help="The type of expression data in the expression_file",
    )
    parser.add_argument(
        "-i",
        "--id-column",
        help="The column header in the expression_file for the column containing gene/transcript ids. Required when using the `custom` format.",
    )
    parser.add_argument(
        "-e",
        "--expression-column",
        help="The column header in the expression_file for the column containing expression data. Required when using the `custom` and `star` format.",
    )
    parser.add_argument(
        "-s",
        "--sample-name",
        help="If the input_vcf contains multiple samples, the name of the sample to annotate.",
    )
    parser.add_argument(
        "-o",
        "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
        + "written next to the input VCF file with a .tx.vcf or .gx.vcf file ending.",
    )
    parser.add_argument(
        "--ignore-ensembl-id-version",
        help='Assumes that the final period and number denotes the Ensembl ID version and ignores it (i.e. for "ENST00001234.3" - ignores the ".3").',
        action="store_true",
    )

    return parser


########################## Arguments check ##########################
def args_check(args):
    if args.format == "custom":
        if args.id_column is None:
            raise Exception(
                "--id-column is not set. This is required when using the `custom` format."
            )
        if args.expression_column is None:
            raise Exception(
                "--expression-column is not set. This is required when using the `custom` format."
            )
    if args.format == "star":
        if args.genecode is None:
            raise Exception("--genecode is not set. This is required when using the `star` format")


########################## Get expected column names ##########################
def get_expected_column_names(format, mode, id_column, exp_column):
    match (format, mode):
        case ("star", _):
            try:
                exp = int(exp_column)
                if 1 < exp and exp <= 4:
                    return [0, exp - 1]
                else:
                    raise Exception("Expression column is invalid")
            except ValueError:
                raise ValueError("Expression column in star should be a number")
        case ("snappy_custom", _):
            return [0, 1]
        case ("kallisto", "gene"):
            return ["gene", "abundance"]
        case ("kallisto", "transcript"):
            return ["target_id", "tpm"]
        case ("cufflinks", _):
            return ["tracking_id", "FPKM"]
        case ("stringtie", "gene"):
            return ["Gene ID", "TPM"]
        case ("stringtie", "transcript"):
            return ["transcript_id", "TPM"]
        case ("custom", _):
            return [id_column, exp_column]
        case _:
            raise Exception("We don't support {format} with mode as {mode} yet.")


########################## Genecode for star parser ##########################
def count_to_TPM(genecode_path, exp_df, no_version=False):
    gc = pr.read_gtf(genecode_path, as_df=True)
    if no_version:
        gc.loc[:, "gene_id"] = gc.loc[:, "gene_id"].apply(lambda x: re.sub(r"\.[0-9]+$", "", x))
    gc = gc[gc["gene_id"].isin(exp_df.iloc[:, 0])]
    gc = gc[(gc.Feature == "gene")]
    exon = gc[["Chromosome", "Start", "End", "gene_id", "gene_name"]]
    # Convert columns to proper types.
    exon.loc[:, "Start"] = exon.loc[:, "Start"].astype(int)
    exon.loc[:, "End"] = exon.loc[:, "End"].astype(int)
    # Sort in place.
    exon = exon.sort_values(by=["Chromosome", "Start", "End"])
    # Group the rows by the Ensembl gene identifier.
    groups = exon.groupby("gene_id")
    lengths = groups.apply(count_bp)
    # Create a new DataFrame with gene lengths and EnsemblID.
    gene_id = lengths.index
    if len(gene_id) != len(exp_df):
        raise Exception("Duplicated gene id in expression file")
    ldf = pd.DataFrame({"length": lengths.values, "gene_id": gene_id}).sort_values(by="gene_id")
    # Calculate TPM
    tpm = calculate_TPM(exp_df, ldf)
    tpm_dict = {k: v for k, v in zip(exp_df.iloc[:, 0], tpm)}
    return tpm_dict


def calculate_TPM(exp_df, ldf):
    "TPM = (reads_transcript / transcript_length_kb) / scaled_total_mapped_reads"
    transcript_length_kb = ldf["length"] / 1000
    RPK = exp_df.iloc[:, 1] / transcript_length_kb
    scale_factor_permil = sum(RPK) / (10**6)
    return RPK / scale_factor_permil


def count_bp(df):
    """Given a DataFrame with the exon coordinates from Gencode for a single
    gene, return the total number of coding bases in that gene.
    Example:
        >>> import numpy as np
        >>> n = 3
        >>> r = lambda x: np.random.sample(x) * 10
        >>> d = pd.DataFrame([np.sort([a,b]) for a,b in zip(r(n), r(n))], columns=['start','end']).astype(int)
        >>> d
           start  end
        0      6    9
        1      3    4
        2      4    9
        >>> count_bp(d)
        7
    Here is a visual representation of the 3 exons and the way they are added:
          123456789  Length
        0      ----       4
        1   --            2
        2    ------       6
            =======       7
    """
    start = df.Start.min()
    end = df.End.max()
    bp = [False] * (end - start + 1)
    for i in range(df.shape[0]):
        s = df.iloc[i]["Start"] - start
        e = df.iloc[i]["End"] - start + 1
        bp[s:e] = [True] * (e - s)
    return sum(bp)


########################## Expression file parser ##########################
def colnames_check(names, exp_names):
    return set(names).issubset(set(exp_names))


def expression_file_parser(
    path, format, mode, id_column, exp_column, no_version=False, genecode=""
):
    col_names = get_expected_column_names(format, mode, id_column, exp_column)
    if isinstance(col_names[0], int):
        exp_df = pd.read_csv(path, skiprows=4, sep="\t", header=None)
        exp_df = exp_df.iloc[:, col_names]
        exp_df = exp_df.sort_values(by=0)
        if no_version:
            exp_df.iloc[:, 0] = exp_df.iloc[:, 0].apply(lambda x: re.sub(r"\.[0-9]+$", "", x))
        tpm_dict = count_to_TPM(genecode, exp_df, no_version)
    else:
        if format == "stringtie" and mode == "transcript":
            exp_df = pr.read_gtf(path, as_df=True)
            exp_df = exp_df[exp_df["feature"] == "transcript"]
        else:
            exp_df = pd.read_csv(path, sep="\t", header=0)

        if colnames_check(col_names, exp_df.columns):
            exp_df = exp_df.loc[:, col_names]
        else:
            raise Exception(
                "Gene id column or expression column is not presense in expression file"
            )

        if no_version:
            exp_df.iloc[:, 0] = exp_df.iloc[:, 0].apply(lambda x: re.sub(r"\.[0-9]+$", "", x))
        tpm_dict = {k: v for k, v in zip(exp_df.iloc[:, 0], exp_df.iloc[:, 1])}
    return tpm_dict


########################## Dna Vcf parser ##########################
def dna_vcf_parser(path, sample_name):
    vcf_reader = Reader.from_path(path)
    is_multi_sample = len(vcf_reader.header.samples.names) > 1
    if is_multi_sample and sample_name is None:
        vcf_reader.close()
        raise Exception(
            "ERROR: VCF {} contains more than one sample. Please use the -s option to specify which sample to annotate.".format(
                path
            )
        )
    elif is_multi_sample and sample_name not in vcf_reader.header.samples.names:
        vcf_reader.close()
        raise Exception(
            "ERROR: VCF {} does not contain a sample column for sample {}.".format(
                path, sample_name
            )
        )
    if "CSQ" not in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception(
            "ERROR: VCF {} is not VEP-annotated. Please annotate the VCF with VEP before running this tool.".format(
                path
            )
        )
    return vcf_reader, is_multi_sample


########################## Rna Vcf parser ##########################
def rna_vcf_parser(path):
    rna_vcf_reader = Reader.from_path(path)
    return rna_vcf_reader


########################## Add expression values functions ##########################
def to_array(dictionary):
    array = []
    for key, value in dictionary.items():
        array.append("{}|{}".format(key, value))
    return sorted(array)


def add_Expression_to_vcf(
    entry,
    is_multi_sample,
    sample_name,
    exp_dict,
    genes,
    tag,
    no_version,
    missing_expressions_count,
    entry_count,
):
    expressions = {}
    for gene in genes:
        entry_count += 1
        if no_version:
            gene = re.sub(r"\.[0-9]+$", "", gene)
        if gene in exp_dict:
            expressions[gene] = exp_dict[gene]
        else:
            missing_expressions_count += 1
    if is_multi_sample:
        entry.FORMAT += [tag]
        entry.call_for_sample[sample_name].data[tag] = to_array(expressions)
    else:
        entry.add_format(tag, to_array(expressions))
    return (entry, missing_expressions_count, entry_count)


def add_AD_to_vcf(entry, is_multi_sample, sample_name, rna_vcf):
    RAD_temp = []
    ROT = []
    try:
        rna_records = rna_vcf.fetch(entry.CHROM, entry.affected_start, entry.affected_end)
        rna_record = next(rna_records)
    except (StopIteration, ValueError):
        # if there is no snp. RAD and RoT will be set to 0
        if is_multi_sample:
            entry.FORMAT += ["RAD"]
            entry.call_for_sample[sample_name].data["RAD"] = 0
            entry.FORMAT += ["ROT"]
            entry.call_for_sample[sample_name].data["ROT"] = 0
        else:
            entry.add_format("RAD", 0)
            entry.add_format("ROT", 0)
    else:
        rna_alt = [alt.value for alt in rna_record.ALT]
        rna_AD = [c.data.get("AD") for c in rna_record.calls][0]
        ROT_temp = sum(rna_AD) - rna_AD[0]
        if rna_alt:
            rna_AD = [c.data.get("AD") for c in rna_record.calls][0]
            for dna_alt in entry.ALT:
                if dna_alt.value in rna_alt:
                    index = rna_record.ALT.index(dna_alt)
                    RAD_temp.append(rna_AD[0])  # AD of REF
                    RAD_temp.append(rna_AD[index + 1])  # AD of ALT
                    ROT_temp = ROT_temp - rna_AD[index + 1]
        ROT.append(ROT_temp)
        if is_multi_sample:
            entry.FORMAT += ["RAD"]
            entry.call_for_sample[sample_name].data["RAD"] = RAD_temp
            entry.FORMAT += ["ROT"]
            entry.call_for_sample[sample_name].data["ROT"] = ROT_temp
        else:
            entry.add_format("RAD", RAD_temp)
            entry.add_format("ROT", ROT_temp)
    return entry


########################## Add fields to header ##########################
def create_vcf_writer(path, vcf_reader, mode, rna_vcf, output_vcf):
    (head, _, tail) = path.rpartition(".vcf")
    new_header = vcf_reader.header.copy()
    if mode == "gene":
        new_header.add_format_line(
            OrderedDict(
                [
                    ("ID", "GX"),
                    ("Number", "."),
                    ("Type", "String"),
                    ("Description", "Gene Expressions"),
                ]
            )
        )
        output_file = ("").join([head, ".gx.vcf", tail])
    elif mode == "transcript":
        new_header.add_format_line(
            OrderedDict(
                [
                    ("ID", "TX"),
                    ("Number", "."),
                    ("Type", "String"),
                    ("Description", "Transcript Expressions"),
                ]
            )
        )
        output_file = ("").join([head, ".tx.vcf", tail])
    if rna_vcf:
        new_header.add_format_line(
            OrderedDict(
                [
                    ("ID", "RAD"),
                    ("Number", "A"),
                    ("Type", "Integer"),
                    (
                        "Description",
                        "Allelic depths from mRNA sequencing data of the same sample",
                    ),
                ]
            )
        )
        new_header.add_format_line(
            OrderedDict(
                [
                    ("ID", "ROT"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    (
                        "Description",
                        "Sum of allelic depths of other alternative allele from mRNA sequencing data of the same sample",
                    ),
                ]
            )
        )
        (head, _, tail) = output_file.rpartition(".vcf")
        output_file = ("").join([head, ".pileup.vcf", tail])
    if output_vcf:
        output_file = output_vcf
    return Writer.from_path(output_file, new_header)


########################## Generate new vcf file ##########################
def adding_extra_information(args):
    dna_vcf, is_multi_sample = dna_vcf_parser(args.input_vcf, args.sample_name)
    if args.rna_vcf:
        rna_vcf = rna_vcf_parser(args.rna_vcf)
    # CONTIGS = [x.id for x in dna_vcf.header.get_lines("contig")]
    format_pattern = re.compile("Format: (.*)")
    csq_format = (
        format_pattern.search(dna_vcf.header.get_info_field_info("CSQ").description)
        .group(1)
        .split("|")
    )
    if args.ignore_ensembl_id_version:
        no_version = True

    vcf_writer = create_vcf_writer(
        args.input_vcf, dna_vcf, args.mode, args.rna_vcf, args.output_vcf
    )
    exp_dict = expression_file_parser(
        args.expression_file,
        args.format,
        args.mode,
        args.id_column,
        args.expression_column,
        no_version,
        args.genecode,
    )

    missing_expressions_count = 0
    entry_count = 0

    for entry in dna_vcf:
        # Add expression data
        transcript_ids = set()
        genes = set()
        if "CSQ" not in entry.INFO:
            logging.warning(
                "Variant is missing VEP annotation. INFO column doesn't contain CSQ field for variant {}".format(
                    entry
                )
            )
            vcf_writer.write_record(entry)
            continue
        for transcript in entry.INFO["CSQ"]:
            for key, value in zip(csq_format, transcript.split("|")):
                if key == "Feature" and value != "" and not value.startswith("ENSR"):
                    transcript_ids.add(value)
                if key == "Gene" and value != "":
                    genes.add(value)

        if args.mode == "gene":
            genes = list(genes)
            if len(genes) > 0:
                (entry, missing_expressions_count, entry_count) = add_Expression_to_vcf(
                    entry,
                    is_multi_sample,
                    args.sample_name,
                    exp_dict,
                    genes,
                    "GX",
                    no_version,
                    missing_expressions_count,
                    entry_count,
                )
        elif args.mode == "transcript":
            transcript_ids = list(transcript_ids)
            if len(transcript_ids) > 0:
                (entry, missing_expressions_count, entry_count) = add_Expression_to_vcf(
                    entry,
                    is_multi_sample,
                    args.sample_name,
                    exp_dict,
                    transcript_ids,
                    "TX",
                    no_version,
                    missing_expressions_count,
                    entry_count,
                )
        # Add RAD and ROT
        if rna_vcf is not None:
            if entry.is_snv():
                entry = add_AD_to_vcf(entry, is_multi_sample, args.sample_name, rna_vcf)
        vcf_writer.write_record(entry)

    dna_vcf.close()
    vcf_writer.close()

    if missing_expressions_count > 0:
        logging.warning(
            "{} of {} {}s did not have an expression entry for their {} id.".format(
                missing_expressions_count, entry_count, args.mode, args.mode
            )
        )


########################## MAIN ##########################
def main(args_input=sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)
    adding_extra_information(args)


if __name__ == "__main__":
    main()
