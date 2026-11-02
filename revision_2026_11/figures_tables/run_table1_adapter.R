# Adapter: runs the original Table 1 generator (upstream/table1/table1_rich_sofa8_14_v2.R) unchanged except for the input cohort
# Replaced: input baseline.csv (corrected cohort: strict SOFA, baseline values at or before t0, early = [0, 24), two inconsistent records removed);
#           database cohort table -> VALUES subquery of the corrected cohort; count assertions and column headers. All other logic unchanged.
# Usage: Rscript run_table1_adapter.R
suppressPackageStartupMessages({library(dplyr); library(readr)})
ROOT <- Sys.getenv("CRRT_REVISION_DIR")
PRIV <- file.path(Sys.getenv("CRRT_REVISION_DIR"), "private_data"); OUTD <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "tables"); dir.create(OUTD, showWarnings = FALSE)
orig <- readLines(file.path(Sys.getenv("CRRT_REPO_DIR"), "upstream/table1/table1_rich_sofa8_14_v2.R"), warn = FALSE)
coh <- readRDS(file.path(PRIV, "cohort_round1.rds"))
for (pop in c("primary", "full")) {
  p <- if (pop == "primary") filter(coh, eq_strict) else coh
  b <- p %>% transmute(stay_id, hadm_id, subject_id, t0, sofa_t0 = sofa_t0_strict, early_crrt = early, died_28d = died28,
                       gender, age, apsiii, cci, map_t0, creatinine_t0, creatinine_baseline, lactate_t0, uo_rate_t0)
  d <- file.path(PRIV, paste0("table1_", pop)); dir.create(d, showWarnings = FALSE); write_csv(b, file.path(d, "baseline.csv"))
  vals <- paste(sprintf("(%d,%d,'%s'::timestamp,%s,%d)", b$stay_id, b$hadm_id, format(b$t0, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                        ifelse(b$early_crrt, "TRUE", "FALSE"), as.integer(b$sofa_t0)), collapse = ",")
  sub <- paste0("(SELECT * FROM (VALUES ", vals, ") AS vv(stay_id,hadm_id,t0,early_crrt,sofa_t0))")
  n <- nrow(b); ne <- sum(b$early_crrt); nd <- n - ne
  code <- orig
  code <- gsub('file.path(data_dir, "baseline.csv")', paste0('"', file.path(d, "baseline.csv"), '"'), code, fixed = TRUE)
  code <- gsub("JOIN data_extract_crrt_strictk3_v2.baseline b", paste0("JOIN ", sub, " b"), code, fixed = TRUE)
  code <- gsub("FROM data_extract_crrt_strictk3_v2.baseline", paste0("FROM ", sub, " AS cs"), code, fixed = TRUE)
  code <- gsub("WHERE sofa_t0 BETWEEN 8 AND 14", "", code, fixed = TRUE)
  code <- gsub("filter(sofa_t0 >= 8, sofa_t0 <= 14) %>%", "", code, fixed = TRUE)
  code <- gsub("stopifnot(nrow(df) == 1502, sum(df$early_crrt) == 146, sum(!df$early_crrt) == 1356)",
               sprintf("stopifnot(nrow(df) == %d, sum(df$early_crrt) == %d, sum(!df$early_crrt) == %d)", n, ne, nd), code, fixed = TRUE)
  code <- gsub("N=1502", paste0("N=", n), code, fixed = TRUE)
  code <- gsub("n=146)", paste0("n=", ne, ")"), code, fixed = TRUE)
  code <- gsub("n=1356)", paste0("n=", nd, ")"), code, fixed = TRUE)
  stopifnot(!any(grepl("data_extract_crrt_strictk3_v2", code)))
  tmp <- file.path(PRIV, paste0("table1_", pop)); Sys.setenv(CRRT_MANUSCRIPT_ROOT = tmp, CRRT_ANALYSIS_ROOT = tmp)
  dir.create(file.path(tmp, "table1_output", "results"), recursive = TRUE, showWarnings = FALSE)
  cat("\n######## Table 1 generator:", pop, " (N=", n, ", early=", ne, ") ########\n", sep = "")
  sink(file.path(OUTD, paste0("table1_", pop, "_run.log")), split = TRUE)
  eval(parse(text = code), envir = new.env())
  sink()
  r <- file.path(tmp, "table1_output", "results")
  file.copy(file.path(r, "table1_rich_sofa8_14_v2.csv"), file.path(OUTD, paste0("table1_", pop, ".csv")), overwrite = TRUE)
  file.copy(file.path(r, "table1_rich_sofa8_14_v2_bestbeffort_audit.csv"), file.path(OUTD, paste0("table1_", pop, "_besteffort.csv")), overwrite = TRUE)
}
cat("\nOutput: tables/ (aggregate tables only)\n")
