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

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: resources/refs/genome}

  - name: mapping
    step: ngs_mapping
    config:
      depends_on:
        reads: data_sets
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
The `run` action then reads the panel VCF and its index from the keys `vcf` and `vcf_tbi` of that task, through `get_upstream_paths()`.

## Step parts and actions

`VariantCallingWorkflow` registers `WritePedigreeStepPart`, the step part of the selected tool and the report step parts.
`Mutect2StepPart` implements the scatter-gather pattern with these actions:

| Action          | Job                                                      |
| --------------- | -------------------------------------------------------- |
| `scatter`       | splits the genome into regions                           |
| `run`           | calls variants in one region (`{scatteritem}`)           |
| `gather`        | joins the calls, the statistics and the read orientation |
| `pileup_normal` | pileup of the normal library, for the contamination      |
| `pileup_tumor`  | pileup of the tumor library, for the contamination       |
| `contamination` | estimates the contamination                              |
| `filter`        | filters the calls                                        |

The fragments of `scatter` and `run` are listed in `temp_outputs`.
They are `temp()` and deleted after `gather`.
The outputs of `gather`, the pileups and the contamination tables are listed in `intermediate_outputs`.
They are deleted only if the task sets `delete_intermediates`.
The Snakefile declares `scattergather: mutect2=...` with `mutect2.num_jobs` as the number of regions.
The step part reaches the scatter and gather functions through `self.parent.workflow.globals`.

## Libraries of a pair

The wildcard of every rule is the tumor `library_name`.
`default_relationships` of the workflow adds the column `matched_normal_lib`, found through the donor and the role `normal`.
`SomaticVariantCallingStepPart.get_normal_lib_name` reads it.
Both BAM files come from `get_upstream_paths("alignments", library_name=...)`.

## Several outputs

`VariantCallingWorkflow.task_produces` returns one signature for the germline tools.
For mutect2, it returns two: all calls, and the same signature with the tag `filtered`.
The `filter` action writes the PASS calls to `{library_name}.vcf.gz` and all calls to `{library_name}.full.vcf.gz`.
`get_output_paths` returns the file that matches the signature it receives.
A consumer that requires `-filtered` gets the `.full.vcf.gz` file.
All other consumers get the PASS calls (see {ref}`dev_contracts`).
