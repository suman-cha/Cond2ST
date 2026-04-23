# Cond2ST

Reproducibility code for *General frameworks for conditional two-sample
testing* (Lee, Cha, Kim, 2024).
[arXiv:2410.16636](https://arxiv.org/abs/2410.16636).

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![R \>= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://cran.r-project.org/)
[![Cite](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](./CITATION.cff)

Given two samples \((X_j, Y_j) \sim P_j\), \(j = 1, 2\), test whether the
conditional distributions \(Y \mid X\) coincide. Two general frameworks,
both proved to control type-I error and to be consistent:

- **CIT route** — Algorithm 1 reduces the problem to a conditional
  independence test, so any CIT (GCM, PCM, RCIT, WGSC, ...) drops in.
- **Density-ratio route** — estimate \(p_1(x)/p_2(x)\) once, plug it into a
  weighted MMD or classifier-based two-sample statistic.

## Citation

```bibtex
@article{lee2024general,
  title   = {General frameworks for conditional two-sample testing},
  author  = {Lee, Seongchan and Cha, Suman and Kim, Ilmun},
  journal = {arXiv preprint arXiv:2410.16636},
  year    = {2024}
}
```

In R: `citation("Cond2ST")`. For the software citation, GitHub's *Cite
this repository* widget reads [`CITATION.cff`](./CITATION.cff).

## Installation

```sh
git clone https://github.com/suman-cha/Cond2ST.git
cd Cond2ST
Rscript setup.R
```

Or as an R package:

```r
remotes::install_github("suman-cha/Cond2ST")
```

`setup.R` installs CRAN dependencies plus five GitHub-only packages
(`RCIT`, `CondIndTests`, `KDist`, `GCM`, `vimp`). A C++17 toolchain is
required to compile `KDist` and the bundled Rcpp MMD kernels.

## Data

| Dataset | Source | License | Loaded via |
|---|---|---|---|
| `superconductivity.csv` | UCI ML Repository, dataset 464 ([Hamidieh, 2018](https://archive.ics.uci.edu/dataset/464/superconductivty+data)) | CC BY 4.0 | `data/superconductivity.csv` (21 263 × 82, shipped) |
| `diamonds` | `ggplot2::diamonds` (Wickham, 2016) | GPL-2 (via ggplot2) | Loaded at runtime; not stored in repo |

See [`data/README.md`](./data/README.md) for schema and preprocessing
notes. Both datasets are used as-is; the runners normalize features at
use time via `apply(X, 2, normalize)`.

## Repository layout

```
R/                # library functions, sourced via R/load_all.R
src/              # Rcpp accelerator for the MMD kernels
data/             # superconductivity.csv (UCI mirror)
scripts/main/     # main-paper experiments
scripts/appendix/ # appendix experiments
inst/             # CITATION (citation()), sessionInfo.txt (frozen env)
reproduce_all.R   # master runner
setup.R           # one-time dependency install
renv.lock         # frozen package versions (renv::restore())
```

## Reproduce

### Smoke test (~2 min)

```sh
CDTST_N_REP=20 CDTST_SCENARIOS=S1U Rscript scripts/main/scenarios_1_3.R
```

Outputs land in `results/main/`.

### Full reproduction (60–120 h, single core)

```sh
Rscript reproduce_all.R
```

`CDTST_CORES=N` (or `options(cdtst.cores = N)`) shards the Monte-Carlo
sweep across `N` workers; the serial and parallel paths are
bit-identical (per-rep RNG state is reset before any output-affecting
draw). Common environment overrides are documented in
[`scripts/README.md`](./scripts/README.md).

### Paper-to-script map

| Artifact | Script |
|---|---|
| Figs 1–3 | `scripts/main/scenarios_1_3.R` |
| Figs 4, 6 | `scripts/main/real_data.R` |
| Fig 5 | `scripts/appendix/density_ratio.R` |
| Figs 7–8 | `scripts/appendix/scenarios_4_5.R` |
| Fig 9 | `scripts/appendix/bandwidth.R` |
| Fig 10 | `scripts/appendix/computation_cost.R` |
| Table 2 | `scripts/appendix/algorithm1.R` |
| Table 3 | `scripts/appendix/epsilon.R` |
| Table 4 | `scripts/appendix/kstar.R` |
| App. A.10 | `scripts/appendix/oracle_cit.R` |
| App. tables (CIT split, non-CV DRT, 10-test real data) | `scripts/appendix/{split_ratio,scenarios_1_3,real_data,scenarios_4_5_single}.R` |

The full table with output filenames lives in
[`scripts/README.md`](./scripts/README.md).

## Computational requirements

- **R** ≥ 4.0, **C++17** toolchain (for `KDist` and the Rcpp MMD kernels).
- **RAM** ≥ 16 GB recommended; per-worker footprint stays under 4 GB.
- **No GPU.**
- **Wall time** at `n_rep = 500`: full sweep ≈ 60–120 h on one core.
  Embarrassingly parallel across reps (`CDTST_CORES=8` brings this to
  ~10–15 h on an 8-core workstation).
- Tested on macOS 14 (arm64) and Linux (x86_64).

## Reproducibility notes

- Every Monte-Carlo runner resets the per-rep RNG with
  `set.seed(seed_base + sim)` before generating data; every `*_test()`
  function starts with `if (!is.null(seed)) set.seed(seed)` before any
  output-affecting RNG draw. Reseeding makes the parallel path
  bit-identical to the serial path.
- Algorithm 1 default shrinkage: `eps = 1/sqrt(log(n))`.
- Frozen package versions are pinned in [`renv.lock`](./renv.lock);
  restore with `renv::restore()`. The exact session captured at paper
  submission is in [`inst/sessionInfo.txt`](./inst/sessionInfo.txt).

## License

[MIT](./LICENSE).

## Acknowledgements

National Research Foundation of Korea (2022R1A4A1033384) and the Korean
government MSIT (RS-2023-00211073).

Issues and questions: <https://github.com/suman-cha/Cond2ST/issues>.
