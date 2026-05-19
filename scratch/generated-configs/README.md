# Generated Task Configs

This directory is for auto-generated task-based `config.yaml` files.

## Generate

```zsh
cd /home/till/projects/snappy-pipeline
python scripts/generate_task_configs.py
```

By default this writes:

- `scratch/generated-configs/all-workflows/config.yaml`

The generator takes `static_data_config` and `data_sets` from:

- `.tests/test-workflow/pipelines/snappy-cancer_wes/.snappy_pipeline/config.yaml`

## Dry Run

```zsh
cd /home/till/projects/snappy-pipeline
snappy-snake -d scratch/generated-configs/all-workflows -n -- --cores 1
```

