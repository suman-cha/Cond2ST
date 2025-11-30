# CIT Estimation Performance Simulations

## Purpose

This folder contains simulation studies addressing the reviewer comment:

> "Additionally, Example 1 assume knowledge of conditional expectations and Example 2 is an extreme case instability under estimation. I think a more detailed discussion is required on how estimation affects instability and as a consequence the alternative approach."

## Files

### 1. `cit_estimation_stability.R`

**Objective**: Quantify how regression estimation errors affect CIT test stability.

**Key Analyses**:
- Compare oracle GCM (known conditional expectations) vs estimated GCM
- Measure regression errors (MSE, RMSE, R²) for f(X) = E[Y|X] and g(X) = E[Z|X]
- Correlate regression quality with test statistic stability
- Compare multiple regression methods: linear, random forest, XGBoost

**Output**:
- `results/appendix/cit_estimation_stability_detailed.csv`: Full simulation results
- `results/appendix/cit_estimation_stability_summary.csv`: Summary statistics
- `figures/appendix/cit_type1_comparison.pdf`: Type I error comparison
- `figures/appendix/cit_error_vs_stability.pdf`: Regression error vs stability plot
- `figures/appendix/cit_stat_distribution.pdf`: Test statistic distributions
- `figures/appendix/cit_regression_r2.pdf`: Regression quality by method

### 2. `stable_vs_unstable_comparison.R`

**Objective**: Directly compare stable (Example 1) vs unstable (Example 2) estimation strategies.

**Methods Compared**:
1. **Oracle**: True conditional expectations (baseline)
2. **Stable (Sample Split)**: Well-specified model with sample splitting
3. **Unstable (Count-Dependent)**: Pathological estimator sensitive to sample composition
4. **Misspecified (Overfit RF)**: Realistic instability via overfitting

**Output**:
- `results/appendix/stable_vs_unstable_detailed.csv`: Full simulation results
- `results/appendix/stable_vs_unstable_summary.csv`: Summary statistics
- `figures/appendix/stable_vs_unstable_type1.pdf`: Type I error comparison
- `figures/appendix/stable_vs_unstable_distributions.pdf`: Test statistic distributions
- `figures/appendix/stable_vs_unstable_sd.pdf`: Stability comparison

## Running the Simulations

```r
# From the repository root directory
source("experiments/appendix/CIT_estimation_performance/cit_estimation_stability.R")
source("experiments/appendix/CIT_estimation_performance/stable_vs_unstable_comparison.R")
```

## Key Findings

1. **Stable estimation** (with sample splitting and correct model specification) yields test statistics that are asymptotically equivalent to the oracle, confirming Example 1.

2. **Unstable estimation** (count-dependent or pathologically overfit estimators) leads to:
   - Inflated Type I error rates
   - Higher variance in test statistics
   - Poor approximation to the limiting N(0,1) distribution

3. **Regression error is correlated with test instability**: Higher RMSE in estimating conditional expectations leads to larger deviations from oracle test statistics.

4. **Practical recommendations**:
   - Use sample splitting to separate estimation from testing
   - Employ regularization to prevent overfitting
   - Avoid estimators that depend on exact sample counts

## References

- Paper Section 3: Approach via Conditional Independence Testing
- Example 1 (Stable case): GCM with known f and g
- Example 2 (Unstable case): Pathological estimators dependent on sample composition

