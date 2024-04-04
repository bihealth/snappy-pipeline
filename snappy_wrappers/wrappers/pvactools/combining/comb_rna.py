import argparse
import sys
import vcfpy
import pyranges as pr, pandas as pd
import numpy as np
import re
from collections import OrderedDict
import logging

########################## Reading vcf files ##########################
def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if args.rna_vcf :
        rna_vcf_reader = vcfpy.Reader.from_path(args.rna_vcf)
    else : rna_vcf_reader=None
    #expression part
    is_multi_sample = len(vcf_reader.header.samples.names) > 1
    if is_multi_sample and args.sample_name is None:
        vcf_reader.close()
        raise Exception("ERROR: VCF {} contains more than one sample. Please use the -s option to specify which sample to annotate.".format(args.input_vcf))
    elif is_multi_sample and args.sample_name not in vcf_reader.header.samples.names:
        vcf_reader.close()
        raise Exception("ERROR: VCF {} does not contain a sample column for sample {}.".format(args.input_vcf, args.sample_name))
    if 'CSQ' not in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is not VEP-annotated. Please annotate the VCF with VEP before running this tool.".format(args.input_vcf))
    if args.mode == 'gene' and 'GX' in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is already gene expression annotated. GX format header already exists.".format(args.input_vcf))
    elif args.mode == 'transcript' and 'TX' in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is already transcript expression annotated. TX format header already exists.".format(args.input_vcf))
    #snv part

    return vcf_reader, rna_vcf_reader, is_multi_sample
########################## Adding extra fields into Format ##########################
def create_vcf_writer(args, vcf_reader):
    (head, sep, tail) = args.input_vcf.rpartition('.vcf')
    new_header = vcf_reader.header.copy()
    if args.mode == 'gene':
        new_header.add_format_line(vcfpy.OrderedDict([('ID', 'GX'), ('Number', '.'), ('Type', 'String'), ('Description', 'Gene Expressions')]))
        output_file = ('').join([head, '.gx.vcf', tail])
    elif args.mode == 'transcript':
        new_header.add_format_line(vcfpy.OrderedDict([('ID', 'TX'), ('Number', '.'), ('Type', 'String'), ('Description', 'Transcript Expressions')]))
        output_file = ('').join([head, '.tx.vcf', tail])
    if (args.rna_vcf) :
        new_header.add_format_line(vcfpy.OrderedDict([
        ('ID','RAD'),('Number', '.'), ('Type', 'String'),('Description', 'Allelic depths from mRNA sequencing data of the same sample')
        ]))
        new_header.add_format_line(vcfpy.OrderedDict([
        ('ID','ROT'),('Number', '.'), ('Type', 'String'),('Description', 'Sum of allelic depths of other alternative allele from mRNA sequencing data of the same sample')
        ]))
        (head,sep,tail) = output_file.rpartition('.vcf')
        output_file = ('').join([head, '.pileup.vcf', tail])
    if args.output_vcf:
        output_file = args.output_vcf
    return vcfpy.Writer.from_path(output_file, new_header)

########################## Auto assign gene id column based on different tools ##########################
def resolve_id_column(args):
    if args.format == 'cufflinks':
        return 'tracking_id'
    elif args.format == 'kallisto':
        if args.mode == 'gene':
            return 'gene'
        elif args.mode == 'transcript':
            return 'target_id'
    elif args.format == 'stringtie':
        if args.mode == 'gene':
            return 'Gene ID'
    elif args.format == 'star':
        return "gene_id"                  # NEED change
    elif args.format == 'custom':
        if args.id_column is None:
            raise Exception("ERROR: `--id-column` option is required when using the `custom` format")
        else:
            return args.id_column
########################## Auto assign expression column based on different tools ##########################
def resolve_expression_column(args):
    if args.format == 'cufflinks':
        return 'FPKM'
    elif args.format == 'kallisto':
        if args.mode == 'gene':
            return 'abundance'
        elif args.mode == 'transcript':
            return 'tpm'
    elif args.format == 'stringtie':
        return 'TPM'
    elif args.format == 'star':
        if args.expression_column is None:
            raise Exception("ERROR: `--expression-column` option is required when using the `star` format")
        else:
            return 'TPM_' + args.expression_column
    elif (args.format == 'custom'):
        if args.expression_column is None:
            raise Exception("ERROR: `--expression-column` option is required when using the `custom` format")
        else:
            return args.expression_column

def to_array(dictionary):
    array = []
    for key, value in dictionary.items():
        array.append("{}|{}".format(key, value))
    return sorted(array)

def resolve_stringtie_id_column(args, headers):
    if 'gene_id' in headers:
        return 'gene_id'
    else:
        return 'transcript_id'

#Star TPM calculation
#If there is genecount file, GENECODE file is required to calculate TPM
##########################Star TPM Calculation##########################
def genecount_reader(genecount_path):
    col_name=["gene_id","unstranded","rf","fr"]
    star_genecount = pd.read_csv(genecount_path,
                            skiprows=4,
                            sep="\t",
                            names=col_name,
                            header=None,
                            ).sort_values(by=["gene_id"])
    #star_genecount["Ensembl_gene_identifier"] = star_genecount.Ensembl_gene_identifier.map(lambda x: x.split(".")[0])
    star_genecount.unstranded=star_genecount.unstranded.astype(int,copy=None)
    star_genecount.rf=star_genecount.rf.astype(int,copy=None)
    star_genecount.fr=star_genecount.fr.astype(int,copy=None)
    #star_genecount.set_index("gene_id")
    return star_genecount

def genecode_reader(genecode_path,star_genecount,expression_column):
    gc = pr.read_gtf(genecode_path,as_df=True)
    gc = gc[gc["gene_id"].isin(star_genecount["gene_id"])]
    gc = gc[(gc.Feature=="gene")]
    exon = gc[["Chromosome","Start","End","gene_id","gene_name"]]
    # Convert columns to proper types.
    exon.Start = exon.Start.astype(int)
    exon.End = exon.End.astype(int)
    # Sort in place.
    exon.sort_values(by=['Chromosome','Start','End'], inplace=True)
    # Group the rows by the Ensembl gene identifier (with version numbers.)
    groups = exon.groupby('gene_id')

    lengths = groups.apply(count_bp)
    # Create a new DataFrame with gene lengths and EnsemblID.
    ensembl_no_version = lengths.index
    ldf = pd.DataFrame({'length': lengths.values,
                        'gene_id': ensembl_no_version})
    star_tpm = pd.merge(ldf,star_genecount, how='inner',on='gene_id')
    star_tpm["TPM_" + expression_column]= calculate_TPM(star_tpm,expression_column)
    return star_tpm

def calculate_TPM(df,strand):
    effLen = (df["length"])/1000
    rate = df[strand]/effLen
    denom = sum(rate)/(10**6)
    return rate/denom

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
        s = df.iloc[i]['Start'] - start
        e = df.iloc[i]['End'] - start + 1
        bp[s:e] = [True] * (e - s)
    return sum(bp)

##########################Passing expression file##########################
def parse_expression_file(args, vcf_reader, vcf_writer):
    if args.format == 'stringtie' and args.mode == 'transcript':
        df_all = pr.read_gtf(args.expression_file)
        df = df_all[df_all["feature"] == "transcript"]
        id_column = resolve_stringtie_id_column(args, df.columns.values)
    elif args.format == "star" :
        id_column = resolve_id_column(args)
        star_genecount = genecount_reader(args.expression_file)
        df = genecode_reader(args.genecode,star_genecount,args.expression_column)
    else:
        id_column = resolve_id_column(args)
        df = pd.read_csv(args.expression_file, sep='\t')
    if args.ignore_ensembl_id_version:
        df['transcript_without_version'] = df[id_column].apply(lambda x: re.sub(r'\.[0-9]+$', '', x))
    expression_column = resolve_expression_column(args)
    if expression_column not in df.columns.values:
        vcf_reader.close()
        vcf_writer.close()
        raise Exception("ERROR: expression_column header {} does not exist in expression_file {}".format(expression_column, args.expression_file))
    if id_column not in df.columns.values:
        vcf_reader.close()
        vcf_writer.close()
        raise Exception("ERROR: id_column header {} does not exist in expression_file {}".format(id_column, args.expression_file))
    return df, id_column, expression_column

##########################Adding RNA AD ##########################
def add_AD_rna_file(entry,is_multi_sample, sample_name,rna_vcf):
    RAD_temp=[]
    ROT = []
    for rna_record in rna_vcf.fetch(entry.CHROM,entry.affected_start,entry.affected_end):
        rna_alt= [alt.value for alt in rna_record.ALT]
        rna_AD = [c.data.get('AD') for c in rna_record.calls][0]
        ROT_temp=sum(rna_AD) - rna_AD[0]
        if not rna_alt:
            continue
        else:
            rna_AD = [c.data.get('AD') for c in rna_record.calls][0]
            for dna_alt in entry.ALT:
                if dna_alt.value in rna_alt:
                    index=rna_record.ALT.index(dna_alt)
                    RAD_temp.append(rna_AD[0])   #AD of REF
                    RAD_temp.append(rna_AD[index+1]) #AD of ALT
                    ROT_temp = ROT_temp - rna_AD[index+1]
        ROT.append(ROT_temp)
        if is_multi_sample:
            entry.FORMAT += ['RAD']
            entry.call_for_sample[sample_name].data['RAD'] = RAD_temp
            entry.FORMAT += ['ROT']
            entry.call_for_sample[sample_name].data['ROT'] = ROT
        else:
            entry.add_format('RAD', RAD_temp)
            entry.add_format('ROT', ROT)
    return entry

##########################Adding Expression file ##########################
def add_expressions(entry, is_multi_sample, sample_name, df, items, tag, id_column, expression_column, ignore_ensembl_id_version, missing_expressions_count, entry_count):
    expressions = {}
    for item in items:
        entry_count += 1
        if ignore_ensembl_id_version:
            subset = df.loc[df['transcript_without_version'] == re.sub(r'\.[0-9]+$', '', item)]
        else:
            subset = df.loc[df[id_column] == item]
        if len(subset) > 0:
            expressions[item] = subset[expression_column].sum()
        else:
            missing_expressions_count += 1
    if is_multi_sample:
        entry.FORMAT += [tag]
        entry.call_for_sample[sample_name].data[tag] = to_array(expressions)
    else:
        entry.add_format(tag, to_array(expressions))
    return (entry, missing_expressions_count, entry_count)

########################## Define tool parameters ##########################
def define_parser():
    parser = argparse.ArgumentParser(
        "comb_rna.py",
        description = "A tool that will add the data from several expression tools' output files" +
                      "and allelic depths from mRNA sequencing data for snvs" +
                      "to the VCF INFO column. Supported tools are StringTie, Kallisto, Star" +
                      "and Cufflinks. There also is a ``custom`` option to annotate with data " +
                      "from any tab-delimited file."
    )

    parser.add_argument(
        "--input_vcf",
        help="A VEP-annotated VCF file"
    )
    parser.add_argument(
        "--expression_file",
        help="A TSV file containing expression estimates"
    )
    parser.add_argument(
        "--genecode",
        help="A genecode file for calculate TPM from star gene count"
    )
    parser.add_argument(
        "--rna-vcf",
        help="A VCF file of RNA-seq data",
    )
    parser.add_argument(
        "--format",
        choices=['kallisto','stringtie','cufflinks','star','custom'],
        help="The file format of the expression file to process. "
            +"Use `custom` to process file formats not explicitly supported. "
            +"The `custom` option requires the use of the --id-column and --expression-column arguments."
    )
    parser.add_argument(
        "--mode",
        choices=['gene', 'transcript'],
        help="The type of expression data in the expression_file"
    )
    parser.add_argument(
        "-i", "--id-column",
        help="The column header in the expression_file for the column containing gene/transcript ids. Required when using the `custom` format."
    )
    parser.add_argument(
        "-e", "--expression-column",
        help="The column header in the expression_file for the column containing expression data. Required when using the `custom` and `star` format."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="If the input_vcf contains multiple samples, the name of the sample to annotate."
    )
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .tx.vcf or .gx.vcf file ending."
    )
    parser.add_argument(
        "--ignore-ensembl-id-version",
        help='Assumes that the final period and number denotes the Ensembl ID version and ignores it (i.e. for "ENST00001234.3" - ignores the ".3").',
        action="store_true"
    )

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)
    CONTIGS = ['chr{}'.format(i) for i in range(1, 23)] + ['chrX', 'chrY']
    if args.format == 'custom':
        if args.id_column is None:
            raise Exception("--id-column is not set. This is required when using the `custom` format.")
        if args.expression_column is None:
            raise Exception("--expression-column is not set. This is required when using the `custom` format.")
    if args.format == 'star':
        if args.genecode is None:
            raise Exception("--genecode is not set. This is required when using the `star` format")
        exp_col = ["unstranded","rf","fr"]
        if args.expression_column not in exp_col:
            raise Exception("""
                            --expression-column is not found. Please choose :
                            unstranded - Assumes a unstranded library.
                            rf - Assumes a stranded library fr-firststrand.
                            fr - Assumes a stranded library fr-secondstrand.
                            """)
    (vcf_reader,rna_vcf_reader, is_multi_sample) = create_vcf_reader(args)
    format_pattern = re.compile('Format: (.*)')
    csq_format = format_pattern.search(vcf_reader.header.get_info_field_info('CSQ').description).group(1).split('|')

    vcf_writer = create_vcf_writer(args, vcf_reader)
    (df, id_column, expression_column) = parse_expression_file(args, vcf_reader, vcf_writer)
    missing_expressions_count = 0
    entry_count = 0
    for entry in vcf_reader:
        #Add expression data
        transcript_ids = set()
        genes = set()
        if 'CSQ' not in entry.INFO:
            logging.warning("Variant is missing VEP annotation. INFO column doesn't contain CSQ field for variant {}".format(entry))
            vcf_writer.write_record(entry)
            continue
        for transcript in entry.INFO['CSQ']:
            for key, value in zip(csq_format, transcript.split('|')):
                if key == 'Feature' and value != '' and not value.startswith('ENSR'):
                    transcript_ids.add(value)
                if key == 'Gene' and value != '':
                    genes.add(value)

        if args.mode == 'gene':
            genes = list(genes)
            if len(genes) > 0:
                (entry, missing_expressions_count, entry_count) = add_expressions(entry, is_multi_sample, args.sample_name, df, genes, 'GX', id_column, expression_column, args.ignore_ensembl_id_version, missing_expressions_count, entry_count)
        elif args.mode == 'transcript':
            transcript_ids = list(transcript_ids)
            if len(transcript_ids) > 0:
                (entry, missing_expressions_count, entry_count) = add_expressions(entry, is_multi_sample, args.sample_name, df, transcript_ids, 'TX', id_column, expression_column, args.ignore_ensembl_id_version, missing_expressions_count, entry_count)
        #Add RAD and ROT
        if  rna_vcf_reader is not None:
            if (entry.CHROM in CONTIGS) and (entry.is_snv()):
                entry = add_AD_rna_file(entry, is_multi_sample, args.sample_name,rna_vcf_reader)
        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

    if missing_expressions_count > 0:
        logging.warning("{} of {} {}s did not have an expression entry for their {} id.".format(missing_expressions_count, entry_count, args.mode, args.mode))

if __name__ == '__main__':
    main()
