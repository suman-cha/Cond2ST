# reproduce_all.R — Run ALL experiments as subprocesses.
#
# Usage (from package root as working directory):
#   Rscript reproduce_all.R
#
# Prerequisites:
#   Rscript setup.R

experiments <- list(
  # Main paper
  list(script = "scripts/main/scenarios_1_3.R",
       desc   = "Scenarios 1-3 (S1U/B, S2U/B, S3U/B) with CIT + DRT + CK"),
  list(script = "scripts/main/real_data.R",
       desc   = "Real data (Diamonds + Superconductivity)"),
  # Appendix
  list(script = "scripts/appendix/scenarios_4_5.R",
       desc   = "Scenarios 4-5 (S4U/B, S5U/B)"),
  list(script = "scripts/appendix/scenarios_4_5_single.R",
       desc   = "Scenarios 4-5 single-split DRT (S4U/B, S5U/B)"),
  list(script = "scripts/appendix/scenarios_1_3.R",
       desc   = "Scenarios 1-3 appendix tables (non-CV DRT)"),
  list(script = "scripts/appendix/real_data.R",
       desc   = "Real data appendix tables (10 tests)"),
  list(script = "scripts/appendix/oracle_cit.R",
       desc   = "Oracle CIT vs Algorithm 1"),
  list(script = "scripts/appendix/split_ratio.R",
       desc   = "Split ratio sensitivity"),
  list(script = "scripts/appendix/algorithm1.R",
       desc   = "Algorithm 1 with/without sensitivity"),
  list(script = "scripts/appendix/epsilon.R",
       desc   = "Epsilon sensitivity"),
  list(script = "scripts/appendix/bandwidth.R",
       desc   = "Bandwidth selection for MMD"),
  list(script = "scripts/appendix/kstar.R",
       desc   = "Refined k* bound comparison"),
  list(script = "scripts/appendix/computation_cost.R",
       desc   = "Computational cost benchmarks"),
  list(script = "scripts/appendix/density_ratio.R",
       desc   = "Density ratio estimation MSE")
)

n_ok   <- 0L
n_fail <- 0L
failures <- character()

for (exp in experiments) {
  message("\n", strrep("=", 70))
  message("  [", exp$desc, "]")
  message("  ", exp$script)
  message(strrep("=", 70))

  if (!file.exists(exp$script)) {
    message("*** SKIPPED: file not found")
    n_fail <- n_fail + 1L
    failures <- c(failures, exp$script)
    next
  }

  ok <- tryCatch({
    source(exp$script, local = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) {
    message("*** FAILED: ", conditionMessage(e))
    FALSE
  })
  if (ok) {
    n_ok <- n_ok + 1L
    message("  OK")
  } else {
    n_fail <- n_fail + 1L
    failures <- c(failures, exp$script)
  }
}

message("\n", strrep("=", 70))
if (n_fail == 0L) {
  message("  ALL ", n_ok, " EXPERIMENTS COMPLETE")
} else {
  message("  ", n_ok, "/", n_ok + n_fail, " succeeded, ", n_fail, " failed:")
  for (f in failures) message("    - ", f)
}
message(strrep("=", 70))
message("Results: results/")
