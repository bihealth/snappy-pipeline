rule get_single_chromosome:
    output:
        reference_path(),
    params:
        species=config["reference"].get("species", "homo_sapiens"),
        datatype=config["reference"].get("datatype", "dna"),
        build=config["reference"].get("build", "GRCh37"),
        release=config["reference"].get("release", 111),
        chromosome=[config["reference"].get("chromosome", "12")],
    log:
        "logs/get_genome.log",
    # cache: "omit-software"
    wrapper:
        "v3.5.3/bio/reference/ensembl-sequence"


rule bgzip:
    input:
        "{prefix}",
    output:
        "{prefix}.gz",
    params:
        extra="",  # optional
    threads: 1
    log:
        "logs/bgzip/{prefix}.log",
    wrapper:
        "v3.5.3/bio/bgzip"


rule samtools_index:
    input:
        "{sample}",
    output:
        "{sample}.fai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v3.5.3/bio/samtools/faidx"


# restrict the reference to a small region on chr12 containing KRAS
rule reference_subregion:
    input:
        reference_path(),
    output:
        "resources/refs/subregion.fa",
    log:
        "logs/restrict_reference.log",
    # cache: "omit-software"
    params:
        region=reference_faidx_region_string,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} {params.region} > {output} 2> {log}"


rule contig_list:
    input:
        "{reference}.fa",
    output:
        "{reference}.fa.genome",
    conda:
        "../envs/bioawk.yaml"
    log:
        "logs/{reference}.contig_list.log",
    shell:
        "bioawk -c fastx '{{print $name;}}' {input} > {output} 2> {log}"


rule gatk_dict:
    input:
        "{reference}.fa",
    output:
        "{reference}.dict",
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/{reference}.gatk_dict.log",
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"


rule bwa_index:
    input:
        "{genome}.fa",
    output:
        idx=multiext("{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{genome}.log",
    params:
        extra="",
    # cache: "omit-software"
    wrapper:
        "v3.5.3/bio/bwa/index"
