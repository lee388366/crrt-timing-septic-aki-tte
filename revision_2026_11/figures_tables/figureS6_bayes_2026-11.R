# Supplemental Figure S6 (Bayesian reanalysis, three primary priors): runs dependencies/S12_bayesian_three_priors.R unchanged
# except: (1) input directory -> tables/figS6_bayes_inputs (same formulas as bayes_evalue_qba_2026-11.py); (2) output directory; (3) panel title S12 -> S6
PKG  <- Sys.getenv("CRRT_OUTPUT_DIR")
SD   <- file.path(Sys.getenv("CRRT_REPO_DIR"), "revision_2026_11/figures_tables/dependencies")
SCRIPT_DIR <- SD
source(file.path(SD, "00_style_S6.R"))
NEWIN <- file.path(PKG, "tables/figS6_bayes_inputs"); inp <- function(f) file.path(NEWIN, f)
OUT_DIR <- file.path(PKG, "figures"); dir.create(OUT_DIR, showWarnings = FALSE)
code <- readLines(file.path(SD, "S12_bayesian_three_priors.R"), warn = FALSE)
code <- code[!grepl("^source\\(|^SCRIPT_DIR <- local|^  if \\(length\\(f\\)\\)", code)]
code <- gsub('title = "Figure S12"', 'title = "Figure S6"', code, fixed = TRUE)
code <- gsub('"FigS12_bayesian_three_priors', '"FigS6_bayesian_three_priors', code, fixed = TRUE)
eval(parse(text = code))
cat("Figure S6 saved to", OUT_DIR, "\n")
