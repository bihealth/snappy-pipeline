#!/bin/bash
# More safety, by turning some bugs into errors.
set -o errexit -o pipefail -o noclobber -o nounset

# =================================================================================================
#
# Mimicks C implementation of alleleCounter (https://github.com/cancerit/alleleCount)
# which is required by ascat (https://github.com/VanLoo-lab/ascat), but not in bioconda nor conda-forge.
#
# The replacement isn't perfect, and we are aware of two types of differences:
# 1. The loci not covered by any read at all are absent from the replacement script
# 2. Some variants have lower coverage in the original code. After extensive checks with samtools,
#    we haven't been able to understand why this coverage should be lower. The output from the script
#    seems in line with samtools mpileup.
#
# =================================================================================================

version=0.1.0

# Default values for optional arguments
default_min_base_qual=20    # Value taken from https://github.com/cancerit/alleleCount/blob/dev/c/src/bam_access.c#L39
default_min_map_qual=35     # Value taken from https://github.com/cancerit/alleleCount/blob/dev/c/src/bam_access.c#L40
default_required_flag=3     # Value taken from https://github.com/cancerit/alleleCount/blob/dev/c/src/bam_access.c#L41
default_filtered_flag=3852  # Value taken from https://github.com/cancerit/alleleCount/blob/dev/c/src/bam_access.c#L42

# Default max depth for bcftools (not original alleleCounter parameter)
max_depth=8000

# =================================================================================================
#
# Argument checking (mirrors alleleCounter)
#
# =================================================================================================

# Mandatory arguments
loci_file=""
hts_file=""
output_file=""

# Optional arguments
ref_file=""
min_base_qual=$default_min_base_qual
min_map_qual=$default_min_map_qual
contig=""
required_flag=$default_required_flag
filtered_flag=$default_filtered_flag

# Display help function
print_help() {
  echo "Usage: alleleCounter -l loci_file.txt -b sample.bam -o output.txt [-m int] [-r ref.fa.fai]"
  echo
  echo " -l  --loci-file [file]           Path to loci file."
  echo " -b  --hts-file [file]            Path to sample HTS file."
  echo " -o  --output-file [file]         Path write output file."
  echo
  echo "Optional"
  echo " -r  --ref-file [file]           Path to reference fasta index file."
  echo "                                 NB. If cram format is supplied via -b and the reference listed in the cram header"
  echo "                                     can't be found alleleCounter may fail to work correctly."
  echo " -m  --min-base-qual [int]       Minimum base quality [Default: ${default_min_base_qual}]."
  echo " -q  --min-map-qual [int]        Minimum mapping quality [Default: ${default_min_map_qual}]."
  echo " -c  --contig [string]           Limit calling to named contig."
  echo " -d  --dense-snps                Unimplemented, ignored."
  echo " -x  --is-10x                    Unimplemented, leads to error."
  echo " -f  --required-flag [int]       Flag value of reads to retain in allele counting default: [${default_required_flag}]."
  echo "                                 N.B. if the proper-pair flag is are selected, alleleCounter will assume paired-end"
  echo "                                      and filter out any proper-pair flagged reads not in F/R orientation."
  echo " -F  --filtered-flag [int]       Flag value of reads to exclude in allele counting default: [${default_filtered_flag}]."
  echo " -v  --version                   Display version number."
  echo " -h  --help                      Display this usage information."
  echo
}

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -l|--loci-file)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing loci file"
        print_help
        exit 1
      fi
      loci_file="$1"
      shift # past value
      ;;
    -b|--hts-file)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing loci file"
        print_help
        exit 1
      fi
      hts_file="$1"
      shift # past value
      ;;
    -o|--output-file)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing output file"
        print_help
        exit 1
      fi
      output_file="$1"
      shift # past value
      ;;
    -r|--ref-file)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing genome reference file"
        print_help
        exit 1
      fi
      ref_file="$1"
      shift # past value
      ;;
    -m|--min-base-qual)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing minimum base quality"
        print_help
        exit 1
      fi
      min_base_qual=$1
      shift # past value
      ;;
    -q|--min-map-qual)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing minimum mapping quality"
        print_help
        exit 1
      fi
      min_map_qual=$1
      shift # past value
      ;;
    -c|--contig)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing contig name"
        print_help
        exit 1
      fi
      contig="$1"
      shift # past value
      ;;
    -f|--required-flag)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing required flag"
        print_help
        exit 1
      fi
      required_flag=$1
      shift # past value
      ;;
    -F|--filtered-flag)
      shift # past argument
      if [[ $# -eq 0 ]]
      then
        echo "Error: Missing filtered flag"
        print_help
        exit 1
      fi
      filtered_flag=$1
      shift # past value
      ;;
    -d|--dense-snps)
      # Unimplemented
      echo "Warning: Option '-d|--dense-snps' is not implemented"
      shift # past argument
      ;;
    -x|--is-10x)
      # Unimplemented
      echo "Error: Option '-x|--is-10x' is not implemented"
      print_help
      exit 1
      ;;
    -v|--version)
      echo $version
      exit 0
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    -*|--*)
      echo "Error: Unknown option $1"
      print_help
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

if [[ ${#POSITIONAL_ARGS[@]} -ne 0 ]]
then
  echo "Error: Unknown additional arguments ${POSITIONAL_ARGS[@]}"
  print_help
  exit 1
fi

if [[ -z "$loci_file" ]]
then
  echo "Error: Missing mandatory input 'loci-file' (location of SNPs)"
  print_help
  exit 1
fi

if [[ -z "$hts_file" ]]
then
  echo "Error: Missing mandatory input 'hts-file' (mapped reads)"
  print_help
  exit 1
fi

if [[ -z "$output_file" ]]
then
  echo "Error: Missing mandatory input 'output-file'"
  print_help
  exit 1
fi

# =================================================================================================
#
# Start of the actions mimicking alleleCounter
#
# 1. Pileup at SNPs loci (excluding indels)
# 2. Transform vcf to alleleCounter's "proprietary" format (see below)
# 3. Copies pileup vcf file for debugging purposes
#
# =================================================================================================

if [[ -z "$ref_file" ]]
then
  echo "Error: Missing mandatory input 'genome reference file' (not required by original alleleCounter)"
  print_help
  exit 1
fi

d=$(mktemp -d -p $TMPDIR)

# -------------------------------------------------------------------------------------------------
# Pileup at ${loci_file} locations for SNPs only, with limits on read & base mapping quality
#
bcftools mpileup \
    --no-BAQ \
    --max-depth ${max_depth} \
    --skip-indels \
    --min-MQ ${min_map_qual} \
    --min-BQ ${min_base_qual} \
    --skip-any-set ${filtered_flag} --skip-any-unset ${required_flag} \
    --fasta-ref "${ref_file}" \
    --regions-file "${loci_file}" \
    "${hts_file}" \
| bcftools annotate \
    --remove 'INFO,^FORMAT/AD' \
| bcftools view \
    --output-type z --output $d/sites.vcf.gz --write-index=tbi

# -------------------------------------------------------------------------------------------------
# Reformat vcf file to mimic files produced by `alleleCounter`:
# Tab-separated file with one row per SNP. The columns are:
# 1. Chromosome name
# 2. Position
# 3. Number of reads supporting an "A" at the position
# 4. Number of reads supporting an "C" at the position
# 5. Number of reads supporting an "G" at the position
# 6. Number of reads supporting an "T" at the position
# 7. Total number of reads at the position
#
rm -f ${output_file}

bcftools query --format "%CHROM\t%POS\t%REF\t%ALT\t[%AD]" $d/sites.vcf.gz \
| awk -F '\t' '
BEGIN {
    print "#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth"
}
{
    if ($3!~/^[ACGT]$/ || $4!~/^([ACGT],)*([ACGT]|<\*>)$/ || $3==$4) {
        printf "Error in %s\n", $0
        exit 1
    }

    a=0; c=0; g=0; t=0;

    patsplit($4, alts, /([ACGT]|<\*>)/);
    patsplit($5, counts, /[0-9]+/)

    if (length(alts)+1 != length(counts)) {
        printf "Error in %s\n", $0
        exit 1
    }

    if ($3=="A") { a += counts[1] };
    if ($3=="C") { c += counts[1] };
    if ($3=="G") { g += counts[1] };
    if ($3=="T") { t += counts[1] };
    
    s = counts[1];
    for (i=1; i<=length(alts); ++i) {
        if (alts[i]=="A") { a += counts[i+1] };
        if (alts[i]=="C") { c += counts[i+1] };
        if (alts[i]=="G") { g += counts[i+1] };
        if (alts[i]=="T") { t += counts[i+1] };
        s += counts[i+1];
    }
    
    printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", $1, $2, a, c, g, t, s)
}' \
> $output_file

# -------------------------------------------------------------------------------------------------
# Save vcf file if possible (not from original alleleCounter, to make debugging easier)
#
if [[ "$output_file" =~ \.txt$ ]]
then
    vcf=$(echo "$output_file" | sed -e "s/\.txt$/\.vcf/")
    if [[ ! -e $vcf.gz ]] && [[ ! -e $vcf.gz.tbi ]]
    then
        mv $d/sites.vcf.gz $vcf.gz
        mv $d/sites.vcf.gz.tbi $vcf.gz.tbi
    fi
fi

