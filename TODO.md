# TODO

## Ground Rules For This Refactor

- [ ] Treat this as a clean break: no deprecations, no compatibility layer, no migration shim.
- [ ] Prefer task-level composition over internal workflow multiplexing.
- [ ] Keep only truly variable wildcards in rule/output patterns.

## 1) Converge Somatic/Germline Variant Workflows

### 1.1 Unified annotation workflow

- [ ] Build one annotation workflow for both contexts (`somatic` and `germline`) with signature-driven behavior.
- [ ] Make both annotation tools available in all contexts: `vep` and `mehari`.
- [ ] Replace context-specific branching with shared step parts where possible (I/O, logging, output contracts).
- [ ] Keep context-specific sheet interpretation behind adapters, not in path templates.
- [ ] Remove internal state flags where task composition/tags can model the state instead.

Current findings:
- `variant_annotation` and `somatic_variant_annotation` already expose compatible `get_output_paths()` keys (`vcf`, `vcf_tbi`).
- Main divergence is in sheet model handling and tool/config model shape.

### 1.2 Unified filtration approach

- [ ] Enforce one filter per task (single-tool/single-filter workflow task design).
- [ ] Replace filter chains inside one task with explicit chained tasks in `tasks:` + `depends_on`.
- [ ] Define a shared filtered-variants contract for downstream consumers.
- [ ] Remove wildcard-heavy internal chain naming in favor of task names for composition.

Current findings:
- Germline filtration is currently a fixed multi-stage chain.
- Somatic filtration is currently an internal configurable chain.
- Both need to move to task-level composition with one filter per task.

### 1.3 Variant calling convergence (stretch)

- [ ] After annotation/filtration convergence, evaluate unifying `variant_calling` and `somatic_variant_calling` under one signature-driven model.
- [ ] Preserve caller-specific options while standardizing contracts (`consumes`/`produces`, expected paths).

## 2) Implement Polars-Based Sample Data Model (Biomedsheets/SODAR)

- [ ] Use `polars` as the internal dataframe backend (single backend; convert to pandas only at external boundaries if needed).
- [ ] Add a normalized, workflow-agnostic table model in abstract workflow code.
- [ ] Ensure the schema can represent at least:
  - [ ] patient/donor identifiers
  - [ ] sample identifiers
  - [ ] relationships (pedigree and tumor/normal pairing)
  - [ ] library identifiers
  - [ ] library type / assay type (e.g., DNA vs RNA)
  - [ ] unit-level data (lane, read group, repetition/run)
  - [ ] assets (`fq1`, `fq2`, `bam`, other file types)
  - [ ] biological/analysis tags (somatic vs germline)
  - [ ] optional metadata such as purity
- [ ] Normalize ingestion from both sources into the same schema:
  - [ ] biomedsheets TSV/JSON
  - [ ] SODAR-derived content (`snappy pull-sheet` path)
- [ ] Add tests validating parity between legacy `shortcut_sheets` traversal and dataframe-driven selectors.

Current findings:
- `DataSetInfo._load_sheet()` currently returns biomedsheets object graphs used directly by workflows.
- There is no central dataframe abstraction yet.

## 3) Add Frozen-Upstream Execution Mode

- [ ] Add first-class frozen-upstream mode in `snappy run` (orchestrator-level feature, not only user-passed Snakemake args).
- [ ] In frozen mode, enforce:
  - [ ] if upstream outputs exist: consume as immutable inputs
  - [ ] if upstream outputs are missing: fail fast with actionable error
  - [ ] never schedule upstream regeneration
- [ ] Evaluate implementation options and pick one:
  - [ ] Orchestrator-level module registration control (do not register upstream producer modules in frozen mode)
  - [ ] Snakemake flag strategy (`--omit-from`/`--until`/other DAG-limiting options) wrapped by CLI
  - [ ] Hybrid approach (CLI preset + orchestrator guardrails)
- [ ] Add integration tests for the scratch-space lifecycle case (upstream raw data gone, downstream task rerun still possible with frozen outputs).

Current findings:
- `snappy run --task <name>` currently changes target selection only; it does not itself freeze upstream execution.

