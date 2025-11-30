# C2ST vs Oracle CIT Comparison

## Purpose

This folder contains simulation studies addressing the reviewer comment:

> "The core idea for using CIT based methods for the conditional two-sample problem is to use subsampling where on average O(√(n log(1/ε))) many samples are discarded, which can affect the finite sample performance. [...] one can compare the power of the proposed test and the corresponding CIT in an oracle setting where we generate Z ∈ {1,2} and then generate (X,Y)|Z ∼ P_{XY}^{(Z)}."

## Experimental Design

### Two Settings Compared

1. **C2ST Setting (with Algorithm 1)**
   - Generate independent samples: $(X_i^{(1)}, Y_i^{(1)}) \sim P_{XY}^{(1)}$ and $(X_i^{(2)}, Y_i^{(2)}) \sim P_{XY}^{(2)}$
   - Apply Algorithm 1 (subsampling) to construct CIT data
   - Discards $O(\sqrt{n \log(1/\varepsilon)})$ samples on average
   - Apply CIT test on subsampled data

2. **Oracle CIT Setting (direct i.i.d. samples)**
   - Generate $Z_i \sim \text{Bernoulli}(n_1/(n_1+n_2))$
   - Generate $(X_i, Y_i) | Z_i \sim P_{XY}^{(Z_i)}$
   - Apply CIT test directly (no subsampling)
   - Represents the "ideal" CIT setting

### Data Generating Process

- **Group 1**: $X^{(1)} \sim N(0, I_p)$, $Y^{(1)} | X^{(1)} \sim N(X^{(1)\top}\beta, 1)$
- **Group 2**: $X^{(2)} \sim N(\mu_2, I_p)$, $Y^{(2)} | X^{(2)} \sim N(\delta + X^{(2)\top}\beta, 1)$

Where:
- $\delta = 0$ for null hypothesis (same conditional distributions)
- $\delta > 0$ for alternative hypothesis (different conditional distributions)

### Parameters Varied

| Parameter | Values |
|-----------|--------|
| Sample size ($n$) | 200, 500, 1000, 2000 |
| Effect size ($\delta$) | 0 (H0), 0.25, 0.5, 0.75 |
| Epsilon types ($\varepsilon$) | 1/n, 1/√log(n), 1/log(n), 1/√n |
| CIT methods | **GCM, PCM, RCIT, WGSC** |
| Regression methods | **linear, ranger (RF), xgboost** |
| Simulations per config | **500** |

### CIT Methods

| Method | Description | Regression Required |
|--------|-------------|---------------------|
| GCM | Generalized Covariance Measure | Yes (linear, RF, XGBoost) |
| PCM | Projected Covariance Measure | Yes (linear, RF, XGBoost) |
| RCIT | Randomized Conditional Independence Test | No (uses internal RFF) |
| WGSC | Weighted General Significance Criterion | Yes (linear, RF, XGBoost) |

## Files

### `c2st_vs_oracle_cit_comparison.R`

Main simulation script that:
1. Generates data in both settings
2. Applies all four CIT tests with three regression methods
3. Compares Type I error and Power
4. Quantifies samples discarded by Algorithm 1
5. Analyzes different epsilon types from paper appendix
6. Produces summary tables and figures

## Running the Experiment

```r
# From the repository root directory
source("experiments/appendix/CIT_oracle/c2st_vs_oracle_cit_comparison.R")
```

**Note**: This is a computationally intensive experiment. Consider running on a cluster or reducing `n_sims` for quick tests.

## Output

### Results Files (in `results/appendix/`)
- `c2st_vs_oracle_cit_detailed.csv`: Full simulation results
- `c2st_vs_oracle_cit_summary.csv`: Summary statistics

### Figures (in `figures/appendix/`)
- `c2st_vs_oracle_type1.pdf`: Type I error comparison by CIT method
- `c2st_vs_oracle_power.pdf`: Power curves comparison
- `c2st_vs_oracle_power_loss.pdf`: Power difference (Oracle - C2ST)
- `c2st_effective_sample_ratio.pdf`: Effective sample ratio (ñ/n)
- `c2st_samples_discarded.pdf`: Percentage of samples discarded
- `c2st_relative_efficiency.pdf`: Relative efficiency (C2ST Power / Oracle Power)
- `c2st_power_loss_by_reg_method.pdf`: Power loss comparison by regression method

## Key Findings

1. **Type I Error**: Both methods maintain correct Type I error control at α = 0.05 for all CIT methods

2. **Power Loss**: Oracle CIT consistently has higher power due to using all samples
   - Power gap is largest for small n and aggressive epsilon (1/√log(n))
   - Gap decreases as n increases
   - GCM and PCM show similar power characteristics
   - RCIT is competitive without requiring regression method specification
   - WGSC tends to be more conservative

3. **Epsilon Type Comparison**:
   - 1/n: Most conservative, fewest samples discarded
   - 1/√n: Moderate sample loss
   - 1/log(n): **Recommended default** - good balance
   - 1/√log(n): Most aggressive, highest sample loss

4. **Regression Method Impact**:
   - Linear: Fast, works well when model is correctly specified
   - Random Forest: Flexible, good default choice
   - XGBoost: High capacity, may overfit in small samples

5. **Practical Implication**: The power loss is modest, especially for large samples and conservative epsilon choices

## Theoretical Background

Algorithm 1 uses subsampling to ensure that the merged dataset has the same distribution as i.i.d. samples from $P_{XYZ}$. The expected number of discarded samples is:

$$\mathbb{E}[n - \tilde{n}] = O\left(\sqrt{n \log(1/\varepsilon)}\right)$$

This loss becomes relatively smaller as n increases (percentage → 0), but the absolute loss grows with √n.

## References

- Paper Section 3: Approach via Conditional Independence Testing
- Algorithm 1: Converting a Conditional Independence Test into a Conditional Two-Sample Test
- Theorem 2: Theoretical guarantees of Algorithm 1
- Paper Appendix: Various epsilon comparison experiments

