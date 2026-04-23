# Reproducibility Guide

Code for the paper *General Frameworks for Conditional Two-Sample Testing*
(Lee, Cha, Kim 2024). All scripts assume the **package root** as the working
directory.

## Quick Start

```bash
Rscript setup.R           # one-time package install (CRAN + GitHub)
Rscript reproduce_all.R   # run every experiment registered in reproduce_all.R
```

## Directory Layout

```
.
├── R/                # library functions (sourced via R/load_all.R)
├── src/              # Rcpp accelerator (mmd_kernels.cpp)
├── data/             # superconductivity.csv (UCI mirror)
├── scripts/
│   ├── main/         # main-paper experiments
│   └── appendix/     # appendix experiments
├── reproduce_all.R   # master runner
└── setup.R           # dependency installer
```

## Individual Experiments

Each script is a standalone runner that writes one or more CSVs under
`results/<section>/`.

### Main paper

| Script | Output |
|---|---|
| `scripts/main/scenarios_1_3.R` | Scenarios 1–3 (S1U/B, S2U/B, S3U/B) — CIT + DRT + CK |
| `scripts/main/real_data.R` | Real data (Diamonds, Superconductivity) — CV DRT |

### Appendix

| Script | Output |
|---|---|
| `scripts/appendix/scenarios_1_3.R` | Scenarios 1–3 — non-CV DRT (appendix tables) |
| `scripts/appendix/scenarios_4_5.R` | Scenarios 4–5 — CIT + DRT + CK |
| `scripts/appendix/scenarios_4_5_single.R` | Scenarios 4–5 — single-split DRT |
| `scripts/appendix/real_data.R` | Real data — appendix tables (10 tests) |
| `scripts/appendix/algorithm1.R` | Algorithm 1 with/without sensitivity |
| `scripts/appendix/oracle_cit.R` | Oracle CIT vs Algorithm 1 |
| `scripts/appendix/split_ratio.R` | Split-ratio sensitivity |
| `scripts/appendix/epsilon.R` | Epsilon sensitivity |
| `scripts/appendix/bandwidth.R` | Bandwidth selection for MMD |
| `scripts/appendix/kstar.R` | Refined k\* bound comparison |
| `scripts/appendix/computation_cost.R` | Computational cost benchmarks |
| `scripts/appendix/density_ratio.R` | Density ratio estimation MSE |

## Running a Single Experiment

```bash
Rscript scripts/main/scenarios_1_3.R
```

Each runner sources `R/load_all.R` first, which:
- loads CRAN + GitHub package dependencies,
- compiles and links `src/mmd_kernels.cpp` (Rcpp), and
- registers a parallel cluster if `CDTST_CORES > 1` is set.

Common environment overrides supported by the runners:

| Variable | Effect |
|---|---|
| `CDTST_CORES` | Worker count (default 1; set 4–8 for typical machines) |
| `CDTST_N_REP` | Override per-cell repetition count (default 500) |
| `CDTST_SCENARIOS` | Comma list, e.g. `"S1U,S2U"` |
| `CDTST_N_VALUES` | Comma list, e.g. `"500,1000,2000"` |

## Compute Profile

- R ≥ 4.2, macOS or Linux
- ≥ 16 GB RAM recommended
- Full reproduction at `n_rep = 500` is on the order of **60–120 hours** on one
  core. Drop `CDTST_N_REP` for fast iteration.

## Output

- CSVs land under `results/main/...` and `results/appendix/...`.
- Both `results/` and `figures/` are git-ignored.

## Abbreviations

Used throughout file names, headers, and CSV columns.

| Abbrev. | Meaning |
|---|---|
| **C2ST** | Conditional two-sample testing |
| **CIT** | Conditional independence test (e.g. GCM, PCM, RCIT, WGSC) |
| **DRT** | Density-ratio-based test (e.g. LinearMMD, BlockMMD, CP, DCP) |
| **CK** | Conditional-kernel test (CMMD, CGED) |
| **CV** | Cross-validated variant of a DRT (e.g. `CV_LinearMMD`) |
| **MMD** | Maximum Mean Discrepancy (linear / block / quadratic) |
| **CP**, **DCP** | Classifier-based two-sample test family |
| **CMMD**, **CGED** | Conditional MMD / Conditional Generalised Energy Distance |
| **RCIT**, **GCM**, **PCM**, **WGSC** | Conditional independence test variants |
| **LL**, **KLR** | Density-ratio estimators: linear-logistic / kernel-logistic regression |
| **S1U/S1B/.../S5B** | Simulation scenarios — digit = scenario index, U = unbounded density ratio, B = bounded |
| **eps** (ε) | Algorithm 1 shrinkage parameter; default `1/sqrt(log(n))` |
| **k\*** | Subsample size returned by Algorithm 1 |
| **RR** | Empirical rejection rate (size or power) |
| **n_rep** | Monte-Carlo repetition count per cell (default 500) |
| **alpha** (α) | Test level (default 0.05) |
