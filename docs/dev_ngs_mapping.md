(dev_ngs_mapping)=

# NGS Mapping Dissection

```{note}
Before reading this chapter, you should

- know the pipeline from the user's perspective (start at {ref}`usage`),

- have read {ref}`dev_intro` and {ref}`dev_conventions`.

After reading this chapter, you should

- know how a step reads FASTQ files,

- know how a step reads a reference and an index from `external_data` tasks,

- know how a step declares what it produces.
```

`ngs_mapping` ({ref}`step_ngs_mapping`) is one of the few steps that read raw FASTQ files.
Most steps work on its alignments or on data further downstream.
`umi_mapping` ({ref}`step_umi_mapping`) reads FASTQ files the same way and produces alignments under the same contract, so this chapter applies to it as well.
The code is in `snappy_pipeline/workflows/ngs_mapping/`.

## Tasks and configuration

A mapping task needs a FASTQ source, a reference and an index of the selected tool.
The reference and the index come from `external_data` tasks.
Their `produces` entries carry the tags that the `depends_on` fields of `ngs_mapping` require (`NgsMappingDependsOn` in `model.py`):

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
```

`tool` selects the mapper.
`BaseStep.__init__` registers the step part of that tool from `tool_parts` (`BwaStepPart`, `BwaMem2StepPart`, `BwaMem3StepPart`, `Minimap2StepPart` or `StarStepPart`).
`NgsMappingWorkflow.__init__` adds the strandedness part for STAR, else the `bam_collect_doc` part when it is enabled, and the BAM QC summary part.
The Snakefile defines the rules of the selected tool only.

## Reading FASTQ files

`depends_on.reads` names the FASTQ source.
It is `data_sets`, an `external_data` task or a task that writes FASTQ files such as `adapter_trimming`.
`BaseStep.read_groups(library_name)` finds the files of a library and returns them as read groups, ordered by the path of the left file.
`BaseStep.reads_input(library_name)` returns the rule inputs that provide them.

`ReadMappingStepPart` hands the reads to the mapping rule with `reads_input_files` and `reads_params` (`snappy_pipeline/workflows/common/reads.py`):

- When the files exist before the run, the rule inputs are the FASTQ files themselves (`reads_left` and `reads_right`).
- When another task writes them, the rule input is that task's `out/` directory.
  The wrappers then read the file lists from `params.args.input`.

There are no link-in jobs.

## Index and reference

A step part reads the index and the reference through `get_upstream_paths`, not through the config.
`NgsMappingWorkflow.get_index_path(tool)` asks `depends_on.index` for an index with the tag of the tool.
`get_index_files(tool)` expands it to the files of the index, which are inputs of the mapping rule:

```python
def _get_input_files_run(self, wildcards):
    inputs = super()._get_input_files_run(wildcards)
    return inputs | {"index": self.parent.get_index_files("bwa")}
```

## Outputs and the contract

`NgsMappingWorkflow.task_produces` declares `alignments [dna]` for the DNA mappers and `alignments [rna]` for STAR.
`load_project()` checks these signatures for every task that names a mapping task in `depends_on.alignments`.

`NgsMappingWorkflow.get_output_paths` returns the `work/` paths of the BAM file, its index and the QC reports.
Consumers get them from `get_upstream_paths("alignments", library_name=...)`.
The mapping rule writes its BAM file to `work/{library_name}/out/{library_name}.bam`, without the tool name in the path.
