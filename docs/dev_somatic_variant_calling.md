(dev_somatic_variant_calling)=

# Somatic Variant Calling Dissection

```{note}
Before reading this chapter, you should have read {ref}`dev_intro` and {ref}`dev_conventions`.
```

There is no `somatic_variant_calling` step any more.
Somatic small variant calling with mutect2 is a tool of the `variant_calling` step ({ref}`step_variant_calling`).
The code is in `snappy_pipeline/workflows/variant_calling/`.
This chapter shows what is specific to a tool with many actions and with several outputs.

## Task configuration

```yaml
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: resources/refs/genome.fa}

  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]
      scan:
        search_paths: [raw]
        search_patterns:
          - left: '(?P<readgroup>.+)\.R1\.fastq\.gz'
            right: '(?P<readgroup>.+)\.R2\.fastq\.gz'

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: resources/refs/genome}

  - name: mapping
    step: ngs_mapping
    config:
      depends_on:
        reads: reads
        reference: genome
        index: bwa_index
      tool: bwa
      bwa: {}

  - name: somatic_calls
    step: variant_calling
    config:
      depends_on:
        alignments: mapping
        reference: genome
      tool: mutect2
      mutect2:
        contamination: {}
```

An optional `depends_on.panel_of_normals` names a `panel_of_normals` task with `tool: mutect2`.
The `run` action then reads the panel VCF and its index through `vcf_inputs("panel_of_normals")`, which takes the key `vcf` of that task from `get_upstream_paths()` and adds the index as a second rule input.

## Step parts and actions

`VariantCallingWorkflow` registers `WritePedigreeStepPart`, the step part of the selected tool and the report step parts.
`Mutect2StepPart` implements the scatter-gather pattern with these actions:

| Action          | Job                                                                   |
| --------------- | --------------------------------------------------------------------- |
| `scatter`       | splits the genome into regions                                        |
| `run`           | calls variants in one region (`{scatteritem}`)                        |
| `gather`        | joins the calls and the statistics, learns the read orientation model |
| `pileup_normal` | pileup of the normal library, for the contamination estimate          |
| `pileup_tumor`  | pileup of the tumor library, for the contamination estimate           |
| `contamination` | estimates the contamination                                           |
| `filter`        | filters the calls                                                     |

```{graphviz}
:alt: The mutect2 actions: scatter, run per region, gather and filter; the two pileups feed contamination, which feeds filter

digraph mutect2 {
    rankdir=LR; nodesep=0.12; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781", fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    scatter [label=<<B>scatter</B>>];
    run [label=<<B>run</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">one job per region ({scatteritem})</FONT><BR ALIGN="LEFT"/>>];
    gather [label=<<B>gather</B>>];
    pileup_normal [label=<<B>pileup_normal</B>>];
    pileup_tumor [label=<<B>pileup_tumor</B>>];
    contamination [label=<<B>contamination</B>>];
    filter [label=<<B>filter</B>>];
    scatter -> run -> gather -> filter;
    pileup_normal -> contamination;
    pileup_tumor -> contamination;
    contamination -> filter;
}
```

The fragments of `scatter` and `run` are listed in `temp_outputs`.
They are `temp()` and deleted after `gather`.
The outputs of `gather`, the pileups and the contamination tables are listed in `intermediate_outputs`.
They are deleted only if the task's `delete_intermediates` names `intermediates`.

The Snakefile declares `scattergather: mutect2=...` with `mutect2.num_jobs` as the number of regions.
The step part inherits `RegionScatter`, which reaches Snakemake's `scatter` and `gather` functions through `self.parent.workflow.globals`.

## Libraries of a pair

The wildcard of every rule is the tumor `library_name`.
The workflow's `default_relationships` adds the column `matched_normal_lib` to the library table: the normal DNA library of the same donor.
`Mutect2StepPart.get_normal_lib_name` reads that column through `BaseStep.matched_normal`.
Both BAM files come from `get_upstream_paths("alignments", library_name=...)`.

## Several outputs

`VariantCallingWorkflow.task_produces` returns the signature of the task's tool.
For mutect2, it returns two signatures: one for all calls, and the same with the tag `filtered` for the PASS calls.

The `filter` action writes the PASS calls to `{library_name}.vcf.gz` and all calls to `{library_name}.full.vcf.gz`.
`get_output_paths` returns the file that matches the signature it receives.
A consumer that requires `-filtered` gets the `.full.vcf.gz` file.
All other consumers get the PASS calls (see {ref}`dev_contracts`).
