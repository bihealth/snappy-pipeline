rule count_alignments:
    input:
        bam="{file}.bam",
    output:
        counts="{file}.alignment_counts.txt",
    threads: workflow.cores
    log:
        "logs/count_alignments/{file}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """samtools view -@{threads} -c {input.bam} > {output.counts} 2> {log}"""
