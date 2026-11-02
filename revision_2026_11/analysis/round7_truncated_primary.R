# Round 7 (PRIMARY ANALYSIS): cumulative stabilised weights truncated at the strategy-specific 99th percentile, as prespecified in the registered protocol (OSF ZHP7C)
# cp_core::build() truncates when options(cp_trunc = TRUE); all other specifications are those of rounds 3-6.
#   1) main inference (round3_final.R; primary/full cohort; 1,000 bootstrap + 200 jackknife -> BCa; secondary outcomes; Bayesian) -> out/round7_report.txt
#   2) all sensitivity analyses (round4_sensitivity.R) -> out/round7s_report.txt
#   3) Figure 3 rows (adds the untruncated row using the round-3 BCa replicates) -> tables/fig3_forest_rows.csv
#   4) crude HR and Schoenfeld test -> out/round7_extras.txt
# Usage: Rscript round7_truncated_primary.R   (resumable; bootstrap files are named round7*)
setwd(Sys.getenv("CRRT_SCRIPT_DIR", getwd()))
options(cp_trunc = TRUE)
r7_step <- function(f, subs) { s <- paste(readLines(f, warn = FALSE), collapse = "\n")
  for (k in names(subs)) { stopifnot(grepl(k, s, fixed = TRUE)); s <- gsub(k, subs[[k]], s, fixed = TRUE) }
  eval(parse(text = s), envir = globalenv()); options(cp_trunc = TRUE) }   # must run in the global environment (cp_core objects such as L0 live there)
R7_T0 <- Sys.time()
# 1) main inference
r7_step("round3_final.R", list("round3_" = "round7_", "Round3 run" = "Round7 (99% truncation) run"))
# 2) sensitivity analyses (the "Weights truncated at 99th percentile" row is redundant here and not used)
r7_step("round4_sensitivity.R", list("round4_" = "round7s_", "Round4 run" = "Round7s (99% truncation) run"))
# 3) Figure 3 rows
OUTCSV0 <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "tables/fig3_forest_rows.csv")
if (!file.exists(sub(".csv", "_untruncated_primary.csv", OUTCSV0, fixed = TRUE))) file.copy(OUTCSV0, sub(".csv", "_untruncated_primary.csv", OUTCSV0, fixed = TRUE))
r7_step("fig3_rows_2026-11.R", list(
  "round3_boot_primary.rds" = "round7_boot_primary.rds", "round3_jk_primary.rds" = "round7_jk_primary.rds",
  "round3_boot_full.rds" = "round7_boot_full.rds", "round3_jk_full.rds" = "round7_jk_full.rds",
  "round4_boot_rrt_sofa.rds" = "round7s_boot_rrt_sofa.rds", "round4_jk_rrt_sofa.rds" = "round7s_jk_rrt_sofa.rds",
  "rob(G3, \"Weights truncated at 99th percentile\", P0, TV0, trunc = TRUE)" =
    "options(cp_trunc = FALSE); bootrow(G3, \"No weight truncation\", P0, TV0, \"round3_boot_primary.rds\", \"round3_jk_primary.rds\"); options(cp_trunc = TRUE)"))
# 4) crude HR and Schoenfeld test
sink("../out/round7_extras.txt", split = TRUE); r7_step("round6_results_extras.R", list("crude HR" = "crude HR")); sink()
cat("\nRound 7 complete in", round(as.numeric(difftime(Sys.time(), R7_T0, units = "mins"))), "minutes\n")
