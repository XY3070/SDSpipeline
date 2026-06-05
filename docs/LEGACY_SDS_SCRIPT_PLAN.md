# SDS Script Arrangement And Implementation Plan

## Scope

This document is the implementation contract for SDS script cleanup in this workspace.

- No execution changes are performed yet.
- Future coding should follow this document unless we revise the document first.
- The canonical script directory is `scripts/`.
- `old_scripts/` is reference material only and should not be treated as the source of truth.

## User Constraints To Preserve

1. The SDS runtime environment is `/data/home/grp-wangyf/intern/miniforge3/envs/sds`.
2. Every chromosome must be divisible into at least two chunks: `p` and `q`.
3. `--test` must test a `2,000,000` bp window starting at the first observed variant position on one chromosome.

## Current Workspace State

### Current canonical scripts

- `scripts/run_sds_input.sh`
- `scripts/generate_b_file.sh`
- `scripts/generate_s_file.sh`
- `scripts/generate_t_file.sh`
- `scripts/generate_o_file.sh`
- `scripts/submit_sds_array.sh`
- `scripts/compute_SDS.py`

### Old reference scripts

- `old_scripts/mk_sds_input.sh`
- `old_scripts/mk_sds_input_login.sh`
- `old_scripts/submit_chunk.sh`
- `old_scripts/run_chunk_sds.sh`

### Current data layout detected

- Population sample lists:
  - `data/NCN.txt`
  - `data/SCN.txt`
- VCFs:
  - `data/vcf/NCN/UKBQC_NCN_chr*.vcf.gz`
  - `data/vcf/SCN/UKBQC_SCN_chr*.vcf.gz`
- Current processed output root:
  - `data/processed/sds_input/<POP>/`

### Important current gap

The new `scripts/` directory currently covers SDS input generation, but it does not yet include a new canonical compute-stage runner equivalent to the old chunked compute workflow in `old_scripts/submit_chunk.sh` and `old_scripts/run_chunk_sds.sh`.

That means the document must plan both:

- input generation arrangement
- compute-stage arrangement

## Target Script Arrangement

The target arrangement should keep `scripts/` as the only active script directory.

### Stage A: Input generation

- `scripts/run_sds_input.sh`
  - top-level input-generation controller for one `POP` + one `CHR`
  - owns `--test`, skip flags, temp directory creation, and p/q orchestration
- `scripts/generate_b_file.sh`
  - writes chromosome arm boundaries from centromere coordinates
- `scripts/generate_s_file.sh`
  - generates singleton positions for one chunk
- `scripts/generate_t_file.sh`
  - generates test SNP genotype table for one chunk
- `scripts/generate_o_file.sh`
  - generates observability vector aligned to sample order

### Stage B: Compute

Planned canonical compute-stage arrangement:

- `scripts/run_sds_compute.sh`
  - top-level compute controller for one `POP` + one `CHR`
  - prepares chunk list, runs `compute_SDS.py`, merges results
- `scripts/run_sds_compute_chunk.sh`
  - worker wrapper for a single compute chunk
- `scripts/compute_SDS.py`
  - unchanged core numerical engine unless a bug fix is needed
- `scripts/submit_sds_array.sh`
  - scheduler entrypoint
  - should become a thin wrapper that can submit either input generation, compute, or both in sequence

## Design Rules

### Rule 1: `scripts/` is the only maintained implementation

- New behavior must be implemented in `scripts/`.
- `old_scripts/` remains only for reference and comparison.
- Any logic copied from `old_scripts/` must be normalized to the new path layout and naming convention.

### Rule 2: environment activation must use the exact env path

Do not rely on environment names like `SDS`.

Planned requirement:

- activate `/data/home/grp-wangyf/intern/miniforge3/envs/sds` by path
- all scheduler and shell entrypoints must use the same environment rule

Planned preference:

- use explicit activation via Miniforge's `conda.sh`, then `conda activate /data/home/grp-wangyf/intern/miniforge3/envs/sds`
- avoid `mamba activate SDS`
- avoid hidden dependence on the user's shell startup files

## Chunking Policy

### Minimum required chunking

Every chromosome must always have at least these two logical chunks:

- `p`
- `q`

This is mandatory even if one arm is much larger than the other.

### Input-generation chunking

For input generation:

- `b_file` defines the p/q arm boundaries
- `run_sds_input.sh` launches `generate_s_file.sh` and `generate_t_file.sh` separately for `p` and `q`
- p and q are processed in parallel inside one chromosome job when resources allow
- outputs from `p` and `q` are merged into chromosome-level final files

### Compute-stage chunking

For compute:

- the minimum compute split must also respect `p` and `q`
- the plan is to avoid a design where a chromosome is treated as one undivided compute unit
- if later performance requires more than two chunks, extra subchunks may be created inside `p` or `q`
- but the top-level visible chunk model remains `chrN_p` and `chrN_q`

### Naming convention

Planned chunk naming:

- arm-level:
  - `chr1_p`
  - `chr1_q`
- optional subchunk-level, only if needed later:
  - `chr1_p_0001`
  - `chr1_q_0001`

## `--test` Behavior Contract

### Required semantics

`--test` must mean:

- run only on one chosen chromosome
- find the first observed variant on that chromosome in the input VCF
- process only the `2,000,000` bp window starting from that first observed variant position
- do not run the full chromosome
- do not silently expand beyond that variant-anchored 2 Mb window

### Planned detailed behavior

When `--test` is enabled:

- query the selected VCF and detect the first variant position on `chr<CHR>`
- define:
  - `TEST_START = first_variant_pos`
  - `TEST_END = min(chr_length, TEST_START + 2000000 - 1)`
- test region is `chr<CHR>:TEST_START-TEST_END`
- this test region is treated as a single limited-range smoke test
- q-arm processing is disabled for that test run because the test window is a single contiguous window, not a true p/q split
- output files should still follow normal naming, but their contents represent only the test interval
- if the chromosome has no variants at all, the script should fail clearly with an explicit error

### About `chr22`

Planning decision:

- `--test` should not hard-code chromosome 22 inside the script
- the user must still choose the chromosome explicitly
- the chosen chromosome is acceptable as long as the script can detect a first variant and anchor the 2 Mb window there

Interpretation of the previous `chr22` issue:

- the old absolute `1-2000000` rule can yield empty test windows on these datasets
- the new variant-anchored rule avoids that failure mode

So the document treats `--test` strictly as a smoke-test mode, but one that is expected to contain real variants.

## Output Arrangement

### Input outputs

Per population under `data/processed/sds_input/<POP>/`:

- `chrN_b_file.txt`
- `chrN_s_file.txt`
- `chrN_t_file.txt`
- `chrN_o_file.txt`
- `logs/`
- `tmp/`

### Planned compute outputs

Per population under `data/processed/sds_output/<POP>/`:

- `chrN_p.sds.tsv`
- `chrN_q.sds.tsv`
- `chrN.sds.tsv`
- `chunks/chrN/` if temporary chunk files are needed
- `logs/`

Planned merge rule:

- merge p/q compute results into one chromosome result
- preserve one header in the final merged file

## Scheduler Plan

### Current scheduler script

`scripts/submit_sds_array.sh` currently:

- assumes `POP="NCN"`
- submits array `1-22`
- runs only `scripts/run_sds_input.sh`
- requests `2` CPUs, matching p/q parallelism

### Planned scheduler behavior

The scheduler script should become parameterizable and path-stable.

Planned responsibilities:

- accept `POP`
- accept stage selection:
  - input only
  - compute only
  - full pipeline
- activate the exact SDS environment path
- write logs into a repo-local or population-local log directory

### Resource plan

Input-generation stage:

- request at least `2` CPUs per chromosome job
- one CPU for `p`, one CPU for `q`

Compute stage:

- request at least `2` CPUs per chromosome job even for minimum p/q execution
- if later subchunking is introduced, CPU request may scale above `2`
- initial implementation should stay simple and make `2` CPUs the guaranteed baseline

## Detailed Implementation Schedule

### Phase 1: freeze interfaces

- normalize all entrypoints around `POP` and `CHR`
- define exact CLI options for input and compute controllers
- define environment bootstrap shared by scheduler and worker scripts

### Phase 2: finish input-stage contract

- review `run_sds_input.sh` against this document
- ensure `--test` means first 2 Mb only
- ensure p/q split is always explicit in normal mode
- ensure merged output order is deterministic

### Phase 3: create new compute-stage scripts

- add `scripts/run_sds_compute.sh`
- add `scripts/run_sds_compute_chunk.sh`
- reuse `scripts/compute_SDS.py`
- port necessary logic from `old_scripts/submit_chunk.sh` and `old_scripts/run_chunk_sds.sh`
- remove dependence on `uv run` if the environment already provides Python and required packages directly

### Phase 4: scheduler unification

- update `scripts/submit_sds_array.sh`
- make it environment-aware
- make it stage-selectable
- keep per-chromosome array submission

### Phase 5: verification

- smoke test with one chromosome and `--test`
- verify generated `b/s/t/o` files exist and align
- verify compute stage runs on p/q chunks
- verify merged chromosome result is produced

## Planned CLI Contract

### `scripts/run_sds_input.sh`

Retain or refine:

- `--pop`
- `--chr`
- `--skip-s`
- `--skip-t`
- `--skip-o`
- `--skip-b`
- `--test`
- `--force`

Behavioral contract:

- normal mode: process full chromosome through p/q arms
- test mode: find the first variant on the selected chromosome, then process only the following 2 Mb window

### `scripts/run_sds_compute.sh`

Planned options:

- `--pop`
- `--chr`
- `--test`
- `--force`
- `--keep-temp`
- `--chunk-lines <N>` only if line-based subchunking is needed

Behavioral contract:

- normal mode: compute on both p and q derived units
- test mode: compute only from the input files generated for the variant-anchored 2 Mb test interval

## Open Items To Resolve During Coding

1. Where `g_file.txt` should live in the new canonical layout.
   Current new `scripts/` flow does not define it, but `compute_SDS.py` still requires it.
2. Whether compute-stage splitting should be strictly arm-based first, or arm-based plus optional line splits from the start.
3. Final naming convention for chromosome SDS outputs.
4. Whether scheduler logs should stay under script-local `logs/` or move under processed output directories.

## Non-Goals For The First Coding Pass

- No redesign of the SDS numerical model in `compute_SDS.py`
- No expansion beyond chromosomes `1-22`
- No hidden fallback to legacy `old_scripts/`
- No change to the test-region size from `2,000,000` bp

## Immediate Coding Direction

When coding starts, the work should proceed in this order:

1. Make environment handling explicit and path-based.
2. Lock `run_sds_input.sh` behavior to this document's `--test` and p/q rules, including first-variant detection.
3. Add the missing new compute-stage scripts under `scripts/`.
4. Update the scheduler to drive the new arrangement.
5. Run smoke validation only after the above is complete.
