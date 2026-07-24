# =============================================================================
# Strict clean-t0 v2 Table 1, aligned to the prior MMR manuscript specification.
# Population: primary SOFA 8-14 equipoise cohort.
# Columns: Characteristic | Overall | Early CRRT within 24 h | No CRRT within 24 h | SMD | p-value
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(flextable)
  library(officer)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
codex_output_dir <- "~/Documents/Codex/2026-07-17/referenced-chatgpt-conversation-this-is-untrusted/outputs"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(codex_output_dir, recursive = TRUE, showWarnings = FALSE)

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "mimiciv31",
  host = "localhost",
  port = 5432,
  user = "postgres",
  password = Sys.getenv("PGPASSWORD")
)
on.exit(dbDisconnect(con), add = TRUE)

baseline <- read_csv(file.path(data_dir, "baseline.csv"), show_col_types = FALSE)

augment <- dbGetQuery(con, "
WITH base AS (
  SELECT stay_id, hadm_id, t0
  FROM data_extract_crrt_strictk3_v2.baseline
  WHERE sofa_t0 BETWEEN 8 AND 14
),
nee AS (
  SELECT
    b.stay_id,
    avg(coalesce(v.norepinephrine, 0)) AS ne_avg,
    avg(coalesce(v.epinephrine, 0)) AS epi_avg,
    avg(coalesce(v.dopamine, 0)) AS dopamine_avg,
    avg(coalesce(v.phenylephrine, 0)) AS phenylephrine_avg,
    avg(coalesce(v.vasopressin, 0)) AS vasopressin_u_per_hour_avg,
    avg(
      coalesce(v.norepinephrine, 0)
      + coalesce(v.epinephrine, 0)
      + coalesce(v.dopamine, 0) / 100.0
      + coalesce(v.phenylephrine, 0) * 0.1
      + coalesce(v.vasopressin, 0) * 2.5
    ) AS nee_old_unscaled_vasopressin,
    avg(
      coalesce(v.norepinephrine, 0)
      + coalesce(v.epinephrine, 0)
      + coalesce(v.dopamine, 0) / 100.0
      + coalesce(v.phenylephrine, 0) * 0.1
      + (coalesce(v.vasopressin, 0) / 60.0) * 2.5
    ) AS nee_at_t0_corrected
  FROM base b
  JOIN mimiciv_derived.vasoactive_agent v
    ON v.stay_id = b.stay_id
   AND v.starttime <= b.t0
   AND v.endtime >= b.t0 - interval '6 hours'
  GROUP BY b.stay_id
)
SELECT
  b.stay_id,
  ad.race,
  fdw.weight AS weight_kg,
  n.ne_avg,
  n.epi_avg,
  n.dopamine_avg,
  n.phenylephrine_avg,
  n.vasopressin_u_per_hour_avg,
  n.nee_old_unscaled_vasopressin,
  n.nee_at_t0_corrected
FROM base b
LEFT JOIN mimiciv_hosp.admissions ad USING (hadm_id)
LEFT JOIN mimiciv_derived.first_day_weight fdw USING (stay_id)
LEFT JOIN nee n USING (stay_id)
")

race_group <- function(x) {
  x <- toupper(ifelse(is.na(x), "", x))
  dplyr::case_when(
    grepl("^WHITE", x) ~ "White",
    grepl("^BLACK", x) ~ "Black/African American",
    grepl("^ASIAN", x) ~ "Asian",
    grepl("HISPANIC|LATINO", x) ~ "Hispanic/Latino",
    TRUE ~ "Other/Unknown"
  )
}

df <- baseline %>%
  filter(sofa_t0 >= 8, sofa_t0 <= 14) %>%
  left_join(augment, by = "stay_id") %>%
  mutate(
    early_crrt = early_crrt %in% TRUE,
    female = gender == "F",
    race_grp = factor(
      race_group(race),
      levels = c("White", "Black/African American", "Asian", "Hispanic/Latino", "Other/Unknown")
    ),
    nee_at_t0_mmr = nee_at_t0_corrected
  )

stopifnot(nrow(df) == 1502)
stopifnot(sum(df$early_crrt) == 146)
stopifnot(sum(!df$early_crrt) == 1356)

fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

fmt_p <- function(p, suffix) {
  if (is.na(p)) return("")
  label <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  paste0(label, suffix)
}

fmt_cont <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    q[2], q[1], q[3]
  )
}

fmt_cat <- function(x) {
  denom <- sum(!is.na(x))
  if (denom == 0) return("")
  n <- sum(x %in% TRUE, na.rm = TRUE)
  sprintf("%s (%.1f)", fmt_n(n), 100 * n / denom)
}

smd_cont <- function(x, g) {
  a <- x[g %in% TRUE]
  b <- x[g %in% FALSE]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  pooled <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  if (!is.finite(pooled) || pooled == 0) return(0)
  abs(mean(a) - mean(b)) / pooled
}

smd_bin <- function(x, g) {
  x <- x %in% TRUE
  a <- x[g %in% TRUE]
  b <- x[g %in% FALSE]
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  p1 <- mean(a, na.rm = TRUE)
  p0 <- mean(b, na.rm = TRUE)
  p <- (p1 + p0) / 2
  denom <- sqrt(p * (1 - p))
  if (!is.finite(denom) || denom == 0) return(0)
  abs(p1 - p0) / denom
}

p_cont <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (length(unique(g[ok])) < 2) return("")
  suppressWarnings(fmt_p(wilcox.test(x[ok] ~ g[ok])$p.value, "ᵃ"))
}

p_cat <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  tab <- table(x[ok], g[ok])
  if (nrow(tab) < 2 || ncol(tab) < 2) return("")
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) {
    p <- if (all(dim(tab) == c(2, 2))) {
      fisher.test(tab)$p.value
    } else {
      set.seed(20260717)
      fisher.test(tab, simulate.p.value = TRUE, B = 20000)$p.value
    }
    fmt_p(p, "ᶜ")
  } else {
    suppressWarnings(fmt_p(chisq.test(tab, correct = FALSE)$p.value, "ᵇ"))
  }
}

overall <- df
early <- df %>% filter(early_crrt)
deferred <- df %>% filter(!early_crrt)

section_row <- function(label) {
  tibble(
    row_type = "section",
    Characteristic = label,
    Overall = "",
    Early = "",
    Deferred = "",
    SMD = "",
    `p-value` = ""
  )
}

continuous_row <- function(label, var, digits) {
  tibble(
    row_type = "data",
    Characteristic = label,
    Overall = fmt_cont(overall[[var]], digits),
    Early = fmt_cont(early[[var]], digits),
    Deferred = fmt_cont(deferred[[var]], digits),
    SMD = sprintf("%.3f", smd_cont(overall[[var]], overall$early_crrt)),
    `p-value` = p_cont(overall[[var]], overall$early_crrt)
  )
}

binary_row <- function(label, expr_all, expr_early, expr_deferred, p_vec) {
  tibble(
    row_type = "data",
    Characteristic = label,
    Overall = fmt_cat(expr_all),
    Early = fmt_cat(expr_early),
    Deferred = fmt_cat(expr_deferred),
    SMD = sprintf("%.3f", smd_bin(expr_all, overall$early_crrt)),
    `p-value` = p_cat(p_vec, overall$early_crrt)
  )
}

race_header <- tibble(
  row_type = "data",
  Characteristic = "Race/ethnicity, n (%)",
  Overall = "",
  Early = "",
  Deferred = "",
  SMD = "",
  `p-value` = p_cat(overall$race_grp, overall$early_crrt)
)

race_rows <- bind_rows(lapply(levels(overall$race_grp), function(level) {
  binary_row(
    paste0("  ", level),
    overall$race_grp == level,
    early$race_grp == level,
    deferred$race_grp == level,
    rep(level, nrow(overall))
  ) %>%
    mutate(`p-value` = "")
}))

table1 <- bind_rows(
  section_row("Demographics"),
  continuous_row("Age, years", "age", 1),
  binary_row("Female sex, n (%)", overall$female, early$female, deferred$female, overall$female),
  race_header,
  race_rows,
  continuous_row("Body weight, kg", "weight_kg", 1),
  section_row("Illness severity"),
  continuous_row("SOFA at t0", "sofa_t0", 0),
  continuous_row("APACHE III (APS III)", "apsiii", 0),
  continuous_row("Charlson comorbidity index", "cci", 0),
  section_row("Shock & metabolic physiology at t0"),
  continuous_row("MAP, mmHg", "map_t0", 1),
  continuous_row("Norepinephrine-equivalent dose, µg/kg/min", "nee_at_t0_mmr", 2),
  continuous_row("Lactate, mmol/L", "lactate_t0", 2),
  section_row("Kidney at t0"),
  continuous_row("Baseline creatinine, mg/dL", "creatinine_baseline", 2),
  continuous_row("Creatinine at t0, mg/dL", "creatinine_t0", 2),
  continuous_row("Urine-output rate, mL/kg/h", "uo_rate_t0", 2),
  binary_row(
    "Stage 3 by creatinine, n (%)",
    overall$stage3_by_creat %in% TRUE,
    early$stage3_by_creat %in% TRUE,
    deferred$stage3_by_creat %in% TRUE,
    overall$stage3_by_creat %in% TRUE
  ),
  binary_row(
    "Stage 3 by urine output, n (%)",
    overall$stage3_by_uo %in% TRUE,
    early$stage3_by_uo %in% TRUE,
    deferred$stage3_by_uo %in% TRUE,
    overall$stage3_by_uo %in% TRUE
  )
)

table1_export <- table1 %>%
  select(-row_type) %>%
  rename(
    `Overall (N = 1,502)` = Overall,
    `Early CRRT within 24 h (n = 146)` = Early,
    `No CRRT within 24 h (n = 1,356)` = Deferred
  )

table1_csv_project <- file.path(results_dir, "strictK3_v2_table1_mmr_spec_sofa8_14.csv")
table1_md_project <- file.path(results_dir, "strictK3_v2_table1_mmr_spec_sofa8_14.md")
table1_docx_project <- file.path(results_dir, "strictK3_v2_table1_mmr_spec_sofa8_14.docx")
table1_csv_canonical <- file.path(results_dir, "strictK3_v2_table1_sofa8_14.csv")
table1_md_canonical <- file.path(results_dir, "strictK3_v2_table1_sofa8_14.md")
table1_docx_canonical <- file.path(results_dir, "strictK3_v2_table1_sofa8_14.docx")
table1_csv_output <- file.path(codex_output_dir, "Table1_clean_t0_v2_MMR_spec.csv")
table1_md_output <- file.path(codex_output_dir, "Table1_clean_t0_v2_MMR_spec.md")
table1_docx_output <- file.path(codex_output_dir, "Table1_clean_t0_v2_MMR_spec.docx")

write_csv(table1_export, table1_csv_project)
write_csv(table1_export, table1_csv_canonical)
write_csv(table1_export, table1_csv_output)

md_lines <- c(
  "# Table 1. Baseline characteristics of the study cohort, by treatment group",
  "",
  "| Characteristic | Overall (N = 1,502) | Early CRRT within 24 h (n = 146) | No CRRT within 24 h (n = 1,356) | SMD | p-value |",
  "|---|---:|---:|---:|---:|---:|",
  apply(table1, 1, function(x) {
    label <- x[["Characteristic"]]
    if (x[["row_type"]] == "section") label <- paste0("**", label, "**")
    sprintf(
      "| %s | %s | %s | %s | %s | %s |",
      label,
      x[["Overall"]],
      x[["Early"]],
      x[["Deferred"]],
      x[["SMD"]],
      x[["p-value"]]
    )
  }),
  "",
  "ᵃ Wilcoxon rank-sum/Mann-Whitney test. ᵇ χ² test. ᶜ Fisher exact test; simulated p-value used for multi-category tables when needed.",
  "SMD, standardized mean difference. APACHE III, Acute Physiology and Chronic Health Evaluation III; CRRT, continuous renal replacement therapy; IQR, interquartile range; MAP, mean arterial pressure; NEE, norepinephrine-equivalent dose; SOFA, Sequential Organ Failure Assessment; t0, index time.",
  "",
  "NEE unit note: vasopressin from `mimiciv_derived.vasoactive_agent` was treated as U/hour and converted to U/min before applying the 2.5 norepinephrine-equivalent multiplier."
)
writeLines(md_lines, table1_md_project)
writeLines(md_lines, table1_md_canonical)
writeLines(md_lines, table1_md_output)

set_flextable_defaults(font.family = "Times New Roman", font.size = 9, padding = 2)
ft <- flextable(table1_export)
ft <- bold(ft, i = which(table1$row_type == "section"), bold = TRUE, part = "body")
ft <- bold(ft, j = 1, bold = TRUE, part = "header")
ft <- align(ft, j = 2:6, align = "center", part = "all")
ft <- align(ft, j = 1, align = "left", part = "all")
ft <- width(ft, j = 1, width = 2.35)
ft <- width(ft, j = 2:4, width = 1.55)
ft <- width(ft, j = 5:6, width = 0.75)
ft <- theme_booktabs(ft)
ft <- autofit(ft)

doc <- read_docx()
doc <- body_add_par(doc, "Table 1. Baseline Characteristics of the Study Cohort, by Treatment Group", style = "heading 1")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(
  doc,
  "ᵃ Wilcoxon rank-sum/Mann-Whitney test. ᵇ χ² test. ᶜ Fisher exact test; simulated p-value used for multi-category tables when needed. SMD, standardized mean difference. APACHE III, Acute Physiology and Chronic Health Evaluation III; CRRT, continuous renal replacement therapy; IQR, interquartile range; MAP, mean arterial pressure; NEE, norepinephrine-equivalent dose; SOFA, Sequential Organ Failure Assessment; t0, index time. NEE unit note: vasopressin from mimiciv_derived.vasoactive_agent was treated as U/hour and converted to U/min before applying the 2.5 norepinephrine-equivalent multiplier.",
  style = "Normal"
)
print(doc, target = table1_docx_project)
print(doc, target = table1_docx_canonical)
print(doc, target = table1_docx_output)

nee_audit <- df %>%
  summarise(
    n = n(),
    early_n = sum(early_crrt),
    old_overall_q1 = quantile(nee_old_unscaled_vasopressin, 0.25, na.rm = TRUE),
    old_overall_median = median(nee_old_unscaled_vasopressin, na.rm = TRUE),
    old_overall_q3 = quantile(nee_old_unscaled_vasopressin, 0.75, na.rm = TRUE),
    corrected_overall_q1 = quantile(nee_at_t0_mmr, 0.25, na.rm = TRUE),
    corrected_overall_median = median(nee_at_t0_mmr, na.rm = TRUE),
    corrected_overall_q3 = quantile(nee_at_t0_mmr, 0.75, na.rm = TRUE),
    corrected_early_q1 = quantile(nee_at_t0_mmr[early_crrt], 0.25, na.rm = TRUE),
    corrected_early_median = median(nee_at_t0_mmr[early_crrt], na.rm = TRUE),
    corrected_early_q3 = quantile(nee_at_t0_mmr[early_crrt], 0.75, na.rm = TRUE),
    corrected_deferred_q1 = quantile(nee_at_t0_mmr[!early_crrt], 0.25, na.rm = TRUE),
    corrected_deferred_median = median(nee_at_t0_mmr[!early_crrt], na.rm = TRUE),
    corrected_deferred_q3 = quantile(nee_at_t0_mmr[!early_crrt], 0.75, na.rm = TRUE),
    vasopressin_q3_u_per_hour = quantile(vasopressin_u_per_hour_avg, 0.75, na.rm = TRUE),
    vasopressin_p90_u_per_hour = quantile(vasopressin_u_per_hour_avg, 0.90, na.rm = TRUE)
  )
write_csv(nee_audit, file.path(results_dir, "strictK3_v2_table1_nee_unit_audit.csv"))
write_csv(nee_audit, file.path(codex_output_dir, "Table1_clean_t0_v2_NEE_unit_audit.csv"))

cat("Wrote Table 1 files:\n")
cat(" - ", table1_csv_project, "\n", sep = "")
cat(" - ", table1_md_project, "\n", sep = "")
cat(" - ", table1_docx_project, "\n", sep = "")
cat("Canonical Table 1 files updated:\n")
cat(" - ", table1_csv_canonical, "\n", sep = "")
cat(" - ", table1_md_canonical, "\n", sep = "")
cat(" - ", table1_docx_canonical, "\n", sep = "")
cat("Copied deliverables to:\n")
cat(" - ", table1_csv_output, "\n", sep = "")
cat(" - ", table1_md_output, "\n", sep = "")
cat(" - ", table1_docx_output, "\n", sep = "")
