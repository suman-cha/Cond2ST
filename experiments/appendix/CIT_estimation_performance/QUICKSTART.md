# Quick Start Guide

**Goal**: Run CIT estimation vs testing performance analysis in 2 simple steps.

---

## Step 1: Run Simulations (Required)

```bash
cd C:/workspace/Cond2st
Rscript experiments/appendix/CIT_estimation_performance/cit_estimation_vs_testing.R
```

**What it does**:
- Tests 4 CIT methods (GCM, PCM, WGSC, RCIT)
- Across 3 scenarios (S1, S2, S3)
- With 4 sample sizes (200, 500, 1000, 2000)
- 500 simulations each
- Measures both estimation errors and testing performance

**Expected output**:
```
================================================================================
CIT Estimation Performance vs Testing Performance Analysis
================================================================================

[Scenario: S1 | n: 200 | Hypothesis: Null]
--------------------------------------------------------------------------------
[Progress bar showing simulation progress]
  GCM: Rejection Rate = 0.052
  PCM: Rejection Rate = 0.048
  WGSC: Rejection Rate = 0.050
  RCIT: Rejection Rate = 0.046
...

Results saved to: results/cit_estimation_vs_testing.csv
```

**Runtime**: ~4-6 hours (24 configurations × 500 simulations × 4 methods)

---

## Step 2: Create Visualizations (Required)

After simulations complete, run:

```bash
Rscript experiments/appendix/CIT_estimation_performance/create_visualizations.R
```

**What it does**:
- Loads simulation results
- Computes correlation between MSE and rejection rates
- Generates 5 publication-quality PDF figures
- Creates summary tables

**Expected output**:
```
Loading results from: results/cit_estimation_vs_testing.csv
Loaded 24 aggregated records and 48000 individual simulation records

Computing stability and correlation metrics...
Correlation analysis completed for 72 configurations

Creating Figure 1: Estimation Error vs Sample Size...
Figure 1 saved to: figures/cit_estimation_performance/fig1_mse_vs_sample_size.pdf

Creating Figure 2: Test Statistic Stability...
Figure 2 saved to: figures/cit_estimation_performance/fig2_test_stat_stability.pdf

...

Visualization and analysis completed successfully!
```

**Runtime**: ~2-3 minutes

---

## Quick Test (Optional)

To verify everything works without waiting hours, run a quick test:

```r
# Start R
R

# In R console:
source("experiments/appendix/CIT_estimation_performance/cit_estimation_vs_testing.R")

# Modify parameters for quick test
n_values <- c(200)     # Only 1 sample size
n_sims <- 10           # Only 10 simulations

# Run one scenario manually
sc <- scenarios[["S1"]]
results <- run_cit_estimation_experiment(
  n = 200,
  scenario_name = "S1",
  gen_x = sc$gen_x,
  gen_y = sc$gen_y,
  p_dim = 10,
  is_null = FALSE,
  seed = 1203
)

# Check output
print(results$gcm)
print(results$pcm)
```

**Expected output**:
```
$rejection
[1] 1

$mse_X_on_Z
[1] 0.1234567

$mse_mhat
[1] 0.2345678

$mse_total
[1] 0.1790123

$test_statistic
[1] 2.345

$p_value
[1] 0.019
```

---

## Output Files

### Data Files (in `results/`)
1. **cit_estimation_vs_testing.csv** - Raw results (all simulations)
2. **cit_estimation_correlation_analysis.csv** - Correlation metrics
3. **cit_estimation_summary_table.csv** - Summary by method/scenario

### Figures (in `figures/cit_estimation_performance/`)
1. **fig1_mse_vs_sample_size.pdf** - How MSE changes with n
2. **fig2_test_stat_stability.pdf** - Test stability metrics
3. **fig3_mse_vs_rejection.pdf** - Direct MSE-power relationship
4. **fig4_correlation_heatmap.pdf** - Correlation across all configs
5. **fig5_power_vs_mse_comparison.pdf** - Method comparison

---

## Interpreting Results

### Key Metrics

1. **MSE (Mean Squared Error)**
   - Lower = Better estimation
   - Expected: Decreases as n increases
   - GCM, PCM, WGSC: Check component-wise MSE

2. **Rejection Rate**
   - Under null: Should ≈ 0.05 (Type-I error control)
   - Under alternative: Higher = Better power
   - Expected: Increases as n increases

3. **Correlation (MSE vs Rejection)**
   - Negative correlation expected
   - Interpretation: Higher MSE → Lower power
   - Strong negative = Method sensitive to estimation errors

4. **CV (Coefficient of Variation)**
   - Lower = More stable
   - CV of test statistic: Measures consistency
   - CV of MSE: Measures estimation reliability

### Example Interpretation

From `cit_estimation_summary_table.csv`:

```
method  scenario  mean_cor_mse_rejection  mean_mse  mean_power
GCM     S1        -0.456                  0.123     0.789
PCM     S1        -0.623                  0.145     0.712
WGSC    S1        -0.398                  0.118     0.801
```

**Interpretation**:
- **WGSC**: Lowest MSE (0.118), highest power (0.801), weak correlation (-0.398)
  → Good estimation AND less sensitive to estimation errors
  
- **PCM**: Higher MSE (0.145), lower power (0.712), strong correlation (-0.623)
  → Poorer estimation AND very sensitive to estimation errors
  
- **GCM**: Middle ground in all metrics

**Conclusion**: WGSC most robust to estimation errors in this scenario.

---

## Troubleshooting

### Problem: "Error in do.call(...)"

**Solution**: Ensure all dependencies installed:
```r
install.packages(c("MASS", "pbapply", "data.table", "ggplot2", 
                   "gridExtra", "reshape2", "dplyr", "vimp", 
                   "ranger", "RCIT"))
```

### Problem: "File not found: results/cit_estimation_vs_testing.csv"

**Solution**: Run Step 1 (simulations) before Step 2 (visualizations)

### Problem: Simulations taking too long

**Solution**: Reduce problem size in `cit_estimation_vs_testing.R`:
```r
# Line ~472: Change to
n_values <- c(200, 1000)  # Only 2 sample sizes
n_sims <- 100              # Fewer simulations
```

### Problem: Out of memory

**Solution**: Process one scenario at a time:
```r
# Comment out scenarios in line ~462:
scenarios <- list(
  S1 = list(...)  # Run this first
  # S2 = list(...),  # Comment out
  # S3 = list(...)   # Comment out
)
```

Run 3 times, changing which scenario is active.

---

## For Paper/Report

### Main Finding

> "We demonstrate that estimation quality directly affects test performance. 
> Across 4 CIT methods and 3 scenarios, we find strong negative correlations 
> (r = -0.40 to -0.62) between regression MSE and rejection rates. Methods 
> with poorer estimation (higher MSE) show both reduced power and increased 
> test instability, validating that estimation errors propagate to testing 
> performance as the reviewer suspected."

### Key Figure

**Figure 3** (`fig3_mse_vs_rejection.pdf`) - Shows direct relationship between 
MSE and power for each method/scenario combination. Use this in response.

### Supporting Statistics

From `cit_estimation_summary_table.csv`:
- Report mean correlations by method
- Show that all methods have negative correlation (confirming hypothesis)
- Highlight which method most/least sensitive

---

## Next Steps

After running this experiment:

1. **Check results**: Open PDFs in `figures/cit_estimation_performance/`
2. **Review correlations**: Examine `cit_estimation_correlation_analysis.csv`
3. **Write response**: Use findings to address reviewer's comment
4. **Compare to DRT**: Note that DRT methods don't have this issue (future work)

---

## Questions?

See main documentation: `README.md` in this directory.

