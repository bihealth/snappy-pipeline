# Helper bash functions.

# Function get_index() --------------------------------------------------------
#
# numeric_index get_index(vcf_path, sample_name)
#
# Get (and print to stdout) the numeric index of $sample_name in the VCF file
# at $vcf_path.

get_index()
{
    [[ "$#" -ne 2 ]] && return

    vcf=$1
    name=$2

    pat=$(
        bcftools view --header-only $vcf \
        | tail -n 1 \
        | cut -f 10- \
        | tr '\t' '|')

    set +o pipefail
    bcftools view --header-only $vcf \
    | grep '^#CHROM' \
    | tr '\t' '\n' \
    | cat -n \
    | grep "\s$name$" \
    | egrep -w "$pat" \
    | awk '{ print $1 - 10; }'
}

# Function samples_vcf_ped() --------------------------------------------------
#
# samples samples_vcf_ped(vcf_path, ped_path)
#
# Get sample names (line by line) that are both in the PED and the VCF file.

samples_vcf_ped()
{
    [[ "$#" -ne 2 ]] && return

    vcf=$1
    ped=$2

    pat=$(
        bcftools view --header-only $vcf \
        | tail -n 1 \
        | cut -f 10- \
        | tr '\t' '|')

    cut -f 2 $ped \
    | egrep -w "$pat"
}
