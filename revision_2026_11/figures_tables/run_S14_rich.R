# Supplemental Table S14 (full baseline characteristics): read-only recomputation with the original rich Table 1 SQL (all windows at or before t0) on the corrected cohort (N = 2,461)
# Two changes only: (1) cohort keys; (2) the old-replay gate is disabled (the historical generator could not be reproduced). SQL and computations unchanged.
# Patient-level intermediate files stay in private_data; only the aggregate table is copied to tables/.
suppressPackageStartupMessages(library(dplyr))
ROOT <- Sys.getenv("CRRT_REVISION_DIR")
PRIV <- file.path(Sys.getenv("CRRT_REVISION_DIR"), "private_data"); OUTP <- file.path(PRIV, "S14_rich"); dir.create(OUTP, showWarnings = FALSE)
OUTD <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "tables")
coh <- readRDS(file.path(PRIV, "cohort_round1.rds"))
need <- c("stay_id","subject_id","hadm_id","icu_intime","icu_outtime","age","gender","icu_los_hours","kdigo3_time","t0")
miss <- setdiff(need, names(coh))
if (length(miss)) {
  bl <- readRDS(file.path(Sys.getenv("CRRT_INPUT_DIR"), "ch24_q4_corrected_inputs/baseline.rds"))
  coh <- left_join(coh, bl[, c("stay_id", intersect(miss, names(bl)))], by = "stay_id")
}
stopifnot(all(need %in% names(coh)), nrow(coh) == 2461, !anyNA(coh$t0))
keys <- coh[, c(need, intersect("anchor_year_group", names(coh)))]
saveRDS(keys, file.path(OUTP, "keys_2461.rds"))
src <- readLines(file.path(Sys.getenv("CRRT_REPO_DIR"), "revision_2026_11/figures_tables/dependencies/TABLE1_READONLY_ORIGINAL_ADAPTER.R"), warn = FALSE)
g <- grep("^if\\(mode=='new'&&!any", src); stopifnot(length(g) == 1)
src[g] <- sub("^if\\(mode=='new'", "if(FALSE&&mode=='new'", src[g])   # old-replay gate disabled: the historical generator could not be reproduced; the corrected cohort is computed directly with the original SQL
ad <- file.path(PRIV, "S14_adapter_nogate.R"); writeLines(src, ad)
st <- system2("Rscript", c(shQuote(ad), "--root", shQuote(ROOT), "--new", shQuote(file.path(OUTP, "keys_2461.rds")),
                           "--out", shQuote(OUTP), "--run-rich"))
stopifnot(st == 0)
file.copy(file.path(OUTP, "DISPLAYED_MAIN_FIELD_REPLAY.csv"), file.path(OUTD, "S14_rich_new2461.csv"), overwrite = TRUE)
r <- read.csv(file.path(OUTD, "S14_rich_new2461.csv"))
cat("rows:", nrow(r), "| computed:", sum(!is.na(r$calculated) & r$calculated != ""), "\n")
print(r[, c("label", "expected", "calculated")], right = FALSE)
