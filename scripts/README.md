# Canonical Scripts

This directory is the first canonical snapshot of reusable SDS pipeline scripts.

It is intentionally flat for compatibility with the current relative-path assumptions among scripts. Logical ownership is tracked in:

- `../manifests/canonical_scripts.tsv`

Rules:

- reusable edits happen here first
- one-off scripts stay in `SDSworkspace/oneoff/`
- legacy copies under `../sds/scripts` should be treated as transitional once this repo becomes active
- every script should have a detailed document for manually command implementation, especially the required arguments and default values

Demography support now also lives here:

- `run_smcpp_benchmark.py`
- `run_smcpp_benchmark_lsf.sh`
- `prepare_smcpp_runtime.sh`
- `run_relate_download_refs.sh`
