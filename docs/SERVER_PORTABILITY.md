# Server Portability

`SDSpipeline/` is meant to move across servers without changing tracked code for every storage layout.

## Rules

- Keep reusable workflow logic in git.
- Keep per-server absolute paths in `config/paths.env`.
- Keep large inputs, outputs, caches, logs, and one-off repair scripts in `SDSworkspace/`.
- Do not commit server-specific mount points into `scripts/`, `templates/`, or `benchmark/` helpers.

## Expected local setup

1. Copy `config/paths.env.example` to `config/paths.env`.
2. Fill in the local roots for `SDSlog`, `SDSpipeline`, `SDSworkspace`, runtime envs, and external tools.
3. Run `scripts/admin/bootstrap_workspace.sh` if the workspace tree does not exist yet.
4. Put server-local helper code under `SDSworkspace/oneoff/`.

## Path policy

- Default script outputs should resolve inside `SDSworkspace/results/` or `SDSworkspace/runs/`.
- Default script inputs should resolve through `common_env.sh` helpers such as `find_population_vcf`, `find_population_sample_list`, and `find_default_g_file`.
- Legacy paths may remain as fallback lookup only during migration. They are not the canonical contract.

## Provenance policy

Every real run should record:

- which server it ran on
- which local paths were used
- which files were unavailable and replaced
- which one-off scripts were introduced locally

That record belongs under `SDSworkspace/provenance/`, with the high-level conclusion summarized in `SDSlog/`.
