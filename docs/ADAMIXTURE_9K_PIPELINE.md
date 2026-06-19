# ADAMIXTURE on the 9k Chinese cohort

Reusable pipeline for population stratification on the frozen 9k Chinese
cohort, based on [ADAMIXTURE](https://github.com/AI-sandbox/ADAMIXTURE)
(Saurina-i-Ricos et al., bioRxiv 2026, DOI 10.64898/2026.02.13.700171).

Full proposal & method notes live in
[`SDSlog/logs/demography/2026-06-19_admixture_proposal_for_9k.md`](../../../SDSlog/logs/demography/2026-06-19_admixture_proposal_for_9k.md).

## Files

| File | Role |
| --- | --- |
| `scripts/setup_admixture_env.sh` | Bootstrap uv-managed Python env |
| `scripts/run_admixture_9k.sh`    | End-to-end driver (6 idempotent steps) |
| `config/paths.env`               | Server-local path overrides |

## One-time setup

```bash
./scripts/setup_admixture_env.sh
```

This creates a uv-managed venv at
`$SDS_ADMIXTURE_ENV_PREFIX` and installs ADAMIXTURE + scientific stack.
No conda/mamba required.

## Run

```bash
# Plan A default: 10 ready autosomes (chr4/5/7/14/17-22).
./scripts/run_admixture_9k.sh all

# To expand to more autosomes once chr1/chr6 are QC-ready:
CHRS_READY=(1 4 5 6 7 14 17 18 19 20 21 22) ./scripts/run_admixture_9k.sh all

# Step-by-step (each is idempotent via .done sentinels):
./scripts/run_admixture_9k.sh step1    # per-chr harmonize
./scripts/run_admixture_9k.sh step2    # concat
./scripts/run_admixture_9k.sh step3    # pgen + LD prune
./scripts/run_admixture_9k.sh step4    # build labels.txt
./scripts/run_admixture_9k.sh step5    # unsupervised K-sweep with CV
./scripts/run_admixture_9k.sh step6 6  # supervised at K=6
```

## Outputs

All under `$SDS_ADMIXTURE_RUN_ROOT`:

| Path | Content |
| --- | --- |
| `work/chrN.biallelic.vcf.gz` | per-chr intermediate |
| `work/9k_merged.vcf.gz` | merged autosome VCF |
| `work/9k.prune.in` | LD-prune SNP list |
| `input/9k_ldpruned.{pgen,psam,pvar}` | final input matrix |
| `input/9k_ldpruned.labels.txt` | Region labels (aligned to psam) |
| `unsupervised/9k.K*.Q` | per-K ancestry proportions (N×K) |
| `unsupervised/9k.K*.P` | per-K allele-frequency matrix |
| `unsupervised/cv_error.tsv` | CV curve for K selection |
| `unsupervised/*.png` | stacked-bar plots |
| `logs/admixture_unsupervised.log` | full run log |

## Portability

All server-local paths live in `config/paths.env`. To move to another server:

1. Copy `config/paths.env.example` → `config/paths.env` and edit.
2. Run `scripts/setup_admixture_env.sh`.
3. Run `scripts/run_admixture_9k.sh all`.
