# -*- coding: utf-8 -*-

import subprocess
import textwrap


def test_add_ped_header():
    ret = subprocess.check_output(
        textwrap.dedent(
            r"""
        source snappy_wrappers/tools/add_ped_header.sh

        set -euxo pipefail
        in_vcf=$(mktemp --suffix=.vcf)
        in_ped=$(mktemp --suffix=.ped)
        out_vcf=$(mktemp --suffix=.vcf.gz)

        cat > $in_vcf <<EOF
        ##fileformat=VCFv4.3
        ##INFO=<ID=TEST1,Number=1,Type=Integer,Description="Test 1">
        ##INFO=<ID=TEST2,Number=1,Type=Integer,Description="Test 2">
        ##INFO=<ID=TEST3,Number=1,Type=Integer,Description="Test 3">
        ##INFO=<ID=TEST4,Number=1,Type=Integer,Description="Test 4">
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
        1	1000	rs1005	C	T	.	.	TEST1=1
        1	2000	rs1006	C	T	.	.	TEST2=2
        1	3000	rs1007	C	T	.	.	TEST3=3
        1	4000	rs1008	C	T	.	.	TEST4=4
        EOF
        bgzip $in_vcf

        cat > $in_ped <<EOF
        # comment
        FAM	II-1	I-1	I-2	1	2
        FAM I-1     0	0	1	1
        FAM I-2     0	0	2	1
        EOF

        snappy-add_ped_header \
            $in_ped \
            $in_vcf.gz \
            $out_vcf

        zcat $out_vcf
        """
        ).lstrip(),
        shell=True,  # nosec
        executable="bash",
    )

    ret = ret.decode("utf-8").split("\n")
    expected = [
        (8, "##META=<ID=Sex,Type=String,Number=1,Values=[Unknown, Male, Female]>"),
        (9, "##META=<ID=Disease,Type=String,Number=1,Values=[Unknown, Unaffected, Affected]>"),
        (10, "##SAMPLE=<ID=II-1,Sex=Male,Disease=Affected>"),
        (11, "##SAMPLE=<ID=I-1,Sex=Male,Disease=Unaffected>"),
        (12, "##SAMPLE=<ID=I-2,Sex=Female,Disease=Unaffected>"),
        (13, "##PEDIGREE=<ID=II-1,Family=FAM,Father=I-1,Mother=I-2>"),
        (14, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"),
    ]

    for e in expected:
        assert ret[e[0]] == e[1]
