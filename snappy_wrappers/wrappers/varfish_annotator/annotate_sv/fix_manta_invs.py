"""Helper script to fix MANTA inversions."""

import argparse
import enum
import subprocess
import sys

import vcfpy


def looks_like_manta(header):
    """Checks VCF header whether the file looks like coming from MANTA."""
    for line in header.lines:
        if line.key == "source" and line.value.startswith("GenerateSVCandidates"):
            return True
    return False


@enum.unique
class InversionType(enum.Enum):
    """Inversion type."""

    INV3 = "INV3"
    INV5 = "INV5"
    NONE = "NONE"


def analyze_inv(record):
    def get_mate_info(alt, c):
        arr = alt.split(c)
        mate_chrom, mate_pos = arr[1].split(":")
        return mate_chrom, int(mate_pos)

    alt = record.ALT[0].serialize()
    if alt.startswith("["):
        mate_chrom, mate_pos = get_mate_info(alt, "[")
        if mate_chrom == record.CHROM:
            return mate_chrom, mate_pos, InversionType.INV5
    elif alt.endswith("]"):
        mate_chrom, mate_pos = get_mate_info(alt, "]")
        if mate_chrom == record.CHROM:
            return mate_chrom, mate_pos, InversionType.INV3
    return None, None, InversionType.NONE


def collect_inv_bnds(reader):
    result = {}
    for record in reader:
        _, _, inv_type = analyze_inv(record)
        if inv_type != InversionType.NONE:
            if record.ID[0] in result:
                result[record.ID[0]] = record
            else:
                result[record.INFO["MATEID"][0]] = None
    return result


def samtools_faidx(reference_fasta, chrom, start_pos, end_pos):
    region = f"{chrom}:{start_pos}-{end_pos}"
    samtools_out = subprocess.check_output(
        ["samtools", "faidx", reference_fasta, region], text=True
    )
    ref_seq = []
    for seq in samtools_out.split("\n"):
        if not seq.startswith(">"):
            ref_seq.append(seq)
    return "".join(ref_seq).upper()


def convert_manta_inversions(args, reader, writer, inv_bnds):
    for record in reader:
        if record.ID[0] in inv_bnds:
            continue  # skip mate record
        _, mate_pos, inv_type = analyze_inv(record)
        if inv_type != InversionType.NONE:
            if inv_type == InversionType.INV5:
                # Adjust POS for 5' inversion
                record.POS -= 1
                mate_pos -= 1
                record.REF = samtools_faidx(
                    args.reference_fasta, record.CHROM, record.POS, record.POS
                )

            id_suffix = record.ID[0].split("MantaBND")[1]
            idx = id_suffix.rfind(":")
            record.ID = ["MantaINV%s" % id_suffix[:idx]]

            record.ALT = [vcfpy.SymbolicAllele("INV")]

            record.INFO["END"] = mate_pos
            record.INFO["SVTYPE"] = "INV"
            record.INFO["SVLEN"] = [mate_pos - record.POS]
            if record.INFO.get("IMPRECISE"):
                mate_id = record.INFO["MATEID"][0]
                record.INFO["CIEND"] = inv_bnds[mate_id].INFO["CIPOS"]
            elif record.INFO.get("HOMLEN"):
                record.INFO["CIEND"] = [-record.INFO["HOMLEN"][0], 0]

            if inv_type == InversionType.INV5 and "HOMSEQ" in record.INFO:
                cipos = record.INFO["CIPOS"]
                hom_seq_start = record.POS + cipos[0] + 1
                hom_seq_end = record.POS + cipos[1]
                record.INFO["HOMSEQ"] = [
                    samtools_faidx(args.reference_fasta, record.CHROM, hom_seq_start, hom_seq_end)
                ]

            # remove BND-specific info
            for key in ("MATEID", "BND_DEPTH", "MATE_BND_DEPTH"):
                if key in record.INFO:
                    record.INFO.pop(key)

            # update INFO/EVENT
            if "EVENT" in record.INFO:
                eid_suffix = record.INFO["EVENT"].split("MantaBND")[1]
                idx = eid_suffix.rfind(":")
                record.INFO["EVENT"] = "MantaINV%s" % eid_suffix[:idx]

            # set INV3/INV5 tags
            if inv_type == InversionType.INV3:
                record.INFO["INV3"] = True
            elif inv_type == InversionType.INV5:
                record.INFO["INV5"] = True

            record.INFO = dict(sorted(record.INFO.items()))

        writer.write_record(record)


def run_manta(args, reader):
    """Run MANTA conversion."""
    # First, collect all inversion BND records.
    with vcfpy.Reader.from_path(args.input_vcf) as reader2:
        inv_bnds = collect_inv_bnds(reader2)
    # Go through file a second time and convert inversions.
    header = reader.header.copy()
    header.add_info_line(
        {
            "ID": "INV3",
            "Number": "0",
            "Type": "Flag",
            "Description": "Inversion breakends open 3' of reported location",
        }
    )
    header.add_info_line(
        {
            "ID": "INV5",
            "Number": "0",
            "Type": "Flag",
            "Description": "Inversion breakends open 5' of reported location",
        }
    )
    header.add_line(
        vcfpy.AltAlleleHeaderLine.from_mapping(
            {
                "ID": "INV",
                "Description": "Inversion",
            }
        )
    )
    with vcfpy.Writer.from_path(args.output_vcf, header) as writer:
        convert_manta_inversions(args, reader, writer, inv_bnds)


def run(args):
    """Main entry point after parsing command line.

    Non-MANTA VCF files are handled directly by copying.  MANTA files are processed with
    run_manta().
    """
    with vcfpy.Reader.from_path(args.input_vcf) as reader:
        if not looks_like_manta(reader.header):  # simply copy
            with vcfpy.Writer.from_path(args.output_vcf) as writer:
                for record in reader:
                    writer.write_record(record)
        else:
            run_manta(args, reader)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA file.")
    parser.add_argument("--input-vcf", required=True, help="Input VCF file.")
    parser.add_argument("--output-vcf", required=True, help="Output VCF file.")

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
