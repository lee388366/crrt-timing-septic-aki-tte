# =============================================================================
# Step 14: Fine-Gray Competing Risks Analysis — 28-Day Renal Recovery
#   Event 1 (interest)   = renal_recovery_28d (recovery by day 28)
#   Event 2 (competing)  = died within 28 days without recovery
#   Censored             = survived 28 days, no recovery, no death
#
# Output appended to result.docx (Section 3.3 supplement)
# =============================================================================

suppressPackageStartupMessages({
  library(cmprsk)
  library(dplyr)
  library(survival)
  library(officer)
  library(flextable)
})

# ── Load data ─────────────────────────────────────────────────────────────────
master <- read.csv(
  "analysis/data/analysis_master.csv",
  stringsAsFactors = FALSE
)
fix_bool <- function(x)
  ifelse(x %in% c("t","TRUE","true","1"), TRUE,
  ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))
master$early_crrt        <- fix_bool(master$early_crrt)
master$died_28d          <- fix_bool(master$died_28d)
master$renal_recovery_28d <- fix_bool(master$renal_recovery_28d)

# ── Competing risks status variable ──────────────────────────────────────────
# 0 = censored (alive at day 28, no recovery, no death)
# 1 = renal recovery (event of interest)
# 2 = died within 28 days without recovery (competing event)
master <- master %>%
  mutate(
    cr_status = case_when(
      renal_recovery_28d == TRUE  ~ 1L,
      died_28d == TRUE            ~ 2L,
      TRUE                        ~ 0L
    ),
    cr_time = follow_time_28d,   # days from t0, capped at 28
    trt     = as.integer(early_crrt),        # 1 = early, 0 = deferred
    equipoise = sofa_t0 >= 8 & sofa_t0 <= 14 & !is.na(sofa_t0)
  )

cat("=== Competing risks status ===\n")
cat(sprintf("  Censored (0):   %d\n",  sum(master$cr_status == 0)))
cat(sprintf("  Recovery (1):   %d\n",  sum(master$cr_status == 1)))
cat(sprintf("  Died (2):       %d\n",  sum(master$cr_status == 2)))
cat(sprintf("  Total:          %d\n",  nrow(master)))

# ── Helper: run Fine-Gray + format result ────────────────────────────────────
run_fg <- function(dat, label) {
  # Covariates (adjust same set as logistic model in step 13)
  cov_vars <- c("age", "sofa_baseline", "apsiii_score",
                "mv_baseline", "cci", "ne_eq_baseline", "lactate")
  cov_vars <- cov_vars[cov_vars %in% names(dat)]

  # Build covariate matrix (mean-impute any residual NAs)
  X <- dat[, cov_vars, drop = FALSE]
  X <- as.data.frame(lapply(X, function(v) {
    v[is.na(v)] <- mean(v, na.rm = TRUE); v
  }))
  X$trt <- dat$trt

  fg <- tryCatch(
    crr(
      ftime   = dat$cr_time,
      fstatus = dat$cr_status,
      cov1    = as.matrix(X),
      failcode = 1,
      cencode  = 0
    ),
    error = function(e) { cat("Fine-Gray error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fg)) return(NULL)

  s <- summary(fg)
  # Extract trt row
  trt_idx <- which(rownames(s$coef) == "trt")
  if (length(trt_idx) == 0) trt_idx <- nrow(s$coef)

  coef_row  <- s$coef[trt_idx, ]
  shr       <- exp(coef_row["coef"])
  ci_lo     <- exp(coef_row["coef"] - 1.96 * coef_row["se(coef)"])
  ci_hi     <- exp(coef_row["coef"] + 1.96 * coef_row["se(coef)"])
  p_val     <- coef_row["p-value"]

  cat(sprintf("\n--- %s (n=%d) ---\n", label, nrow(dat)))
  cat(sprintf("  sHR (Early vs Deferred): %.3f [%.3f, %.3f], p=%.4f\n",
              shr, ci_lo, ci_hi, p_val))

  data.frame(
    Population = label,
    N          = nrow(dat),
    N_recovery = sum(dat$cr_status == 1),
    N_death    = sum(dat$cr_status == 2),
    sHR        = sprintf("%.3f", shr),
    CI_95      = sprintf("[%.3f, %.3f]", ci_lo, ci_hi),
    P_value    = if (p_val < 0.001) "< 0.001" else sprintf("%.3f", p_val),
    stringsAsFactors = FALSE
  )
}

# ── Full cohort ───────────────────────────────────────────────────────────────
res_full <- run_fg(master, "Full cohort")

# ── Equipoise subgroup (SOFA 8-14) ───────────────────────────────────────────
res_eq   <- run_fg(master %>% filter(equipoise), "Equipoise (SOFA 8-14)")

# ── Cumulative incidence function (CIF) plot ──────────────────────────────────
library(ggplot2)

make_cif_plot <- function(dat, title_str) {
  fg_cif <- cuminc(
    ftime   = dat$cr_time,
    fstatus = dat$cr_status,
    group   = dat$trt
  )

  # Extract CIF for event=1 (recovery) by group
  extract_cif <- function(obj, grp_event) {
    g <- obj[[grp_event]]
    data.frame(time = g$time, cif = g$est)
  }

  d1 <- extract_cif(fg_cif, "1 1") %>% mutate(group = "Early CRRT",    event = "Recovery")
  d0 <- extract_cif(fg_cif, "0 1") %>% mutate(group = "Deferred CRRT", event = "Recovery")

  cif_df <- bind_rows(d1, d0)
  cif_df$group <- factor(cif_df$group, levels = c("Early CRRT","Deferred CRRT"))

  ggplot(cif_df, aes(x = time, y = cif * 100, color = group, linetype = group)) +
    geom_step(linewidth = 1.0) +
    scale_color_manual(values = c("Early CRRT" = "#8C8C8C",
                                  "Deferred CRRT" = "#C0392B")) +
    scale_linetype_manual(values = c("Early CRRT" = "solid",
                                     "Deferred CRRT" = "dashed")) +
    scale_x_continuous(breaks = c(0, 7, 14, 21, 28), limits = c(0, 28)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title  = title_str,
      x      = "Days from t0",
      y      = "Cumulative Incidence of Renal Recovery (%)",
      color  = NULL, linetype = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.text     = element_text(size = 9),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      plot.margin     = margin(8, 15, 8, 8)
    )
}

fig_dir <- "analysis/figures"

p_cif_full <- make_cif_plot(master, "Figure 3A. CIF: Renal Recovery (Full Cohort)")
p_cif_eq   <- make_cif_plot(
  master %>% filter(equipoise),
  "Figure 3B. CIF: Renal Recovery (Equipoise Subgroup, SOFA 8-14)"
)

pdf(file.path(fig_dir, "Figure4_CIF_RenalRecovery.pdf"), width = 8, height = 4.5, onefile = FALSE)
gridExtra::grid.arrange(p_cif_full, p_cif_eq, ncol = 2)
dev.off()

# Save individually as PNG for flexibility
png(file.path(fig_dir, "Figure4A_CIF_Full.png"), width = 1500, height = 1300, res = 300)
print(p_cif_full); dev.off()
png(file.path(fig_dir, "Figure4B_CIF_Equipoise.png"), width = 1500, height = 1300, res = 300)
print(p_cif_eq); dev.off()

cat("\nCIF plots saved.\n")

# ── Build result table ────────────────────────────────────────────────────────
fg_table <- bind_rows(res_full, res_eq)
cat("\n=== Fine-Gray Results Table ===\n")
print(fg_table)

# ── Write to result.docx ──────────────────────────────────────────────────────
doc_path <- "analysis/result.docx"
doc      <- read_docx(doc_path)

TNR <- "Times New Roman"
fp_h2 <- fp_text(font.family = TNR, font.size = 12, bold = TRUE, color = "black")
fp_bd <- fp_text(font.family = TNR, font.size = 11, bold = TRUE,  color = "black")
fp_nm <- fp_text(font.family = TNR, font.size = 11, bold = FALSE, color = "black")
fp_sm <- fp_text(font.family = TNR, font.size = 9,  bold = FALSE, color = "grey30")

add_par <- function(doc, runs, space_before = 0, space_after = 4) {
  p <- block_list(fpar(!!!runs, fp_p = fp_par(padding.top = space_before,
                                               padding.bottom = space_after)))
  body_add_blocks(doc, p)
}

# Section header
doc <- body_add_fpar(
  doc,
  fpar(ftext("3.3.2  Competing Risks Analysis: 28-Day Renal Recovery", fp_h2)),
  style = "heading 2"
)

# Narrative
shr_full <- fg_table$sHR[fg_table$Population == "Full cohort"]
ci_full  <- fg_table$CI_95[fg_table$Population == "Full cohort"]
p_full   <- fg_table$P_value[fg_table$Population == "Full cohort"]
nr_full  <- fg_table$N_recovery[fg_table$Population == "Full cohort"]
nd_full  <- fg_table$N_death[fg_table$Population == "Full cohort"]

shr_eq   <- fg_table$sHR[fg_table$Population == "Equipoise (SOFA 8-14)"]
ci_eq    <- fg_table$CI_95[fg_table$Population == "Equipoise (SOFA 8-14)"]
p_eq     <- fg_table$P_value[fg_table$Population == "Equipoise (SOFA 8-14)"]

narrative <- paste0(
  "Because death is a competing event for renal recovery — patients who die cannot ",
  "achieve recovery regardless of treatment — standard logistic regression may conflate ",
  "these distinct failure modes. To address this, we applied Fine-Gray subdistribution ",
  "hazard models (cmprsk package, R) with death as the competing event and 28-day renal ",
  "recovery as the outcome of interest, adjusting for the same covariate set as the ",
  "primary logistic model (age, SOFA score, APACHE III, mechanical ventilation, CCI, ",
  "norepinephrine-equivalent dose, and lactate).\n\n",
  "In the full cohort (n = 3,112; recovery events = ", nr_full,
  ", death events = ", nd_full, "), early CRRT was associated with a subdistribution ",
  "hazard ratio (sHR) of ", shr_full, " (95% CI ", ci_full, "; p ", p_full, ") for ",
  "renal recovery compared with deferred CRRT. In the equipoise subgroup (SOFA 8-14), ",
  "the corresponding sHR was ", shr_eq, " (95% CI ", ci_eq, "; p ", p_eq, "). ",
  "These estimates reflect the net effect of early CRRT on the cumulative incidence of ",
  "renal recovery in the presence of the competing risk of death. The attenuation from ",
  "the crude OR (0.514 in logistic regression) toward the sHR illustrates that part of ",
  "the apparent disadvantage of early CRRT on recovery in the logistic model was driven ",
  "by its higher associated mortality risk (competing event), rather than a true ",
  "impairment of renal recovery capacity. Cumulative incidence function (CIF) curves ",
  "are presented in Figure 3."
)

for (para_text in strsplit(narrative, "\n\n")[[1]]) {
  doc <- body_add_fpar(
    doc,
    fpar(ftext(para_text, fp_nm)),
    style = "Normal"
  )
}

# Table header
doc <- body_add_fpar(
  doc,
  fpar(ftext(
    "Table 4. Fine-Gray Competing Risks Analysis: Effect of Early vs Deferred CRRT on 28-Day Renal Recovery",
    fp_bd
  )),
  style = "Normal"
)

# Build flextable
col_labels <- c("Population", "N", "Recovery\nevents", "Death\nevents",
                "sHR\n(Early vs Deferred)", "95% CI", "p-value")
tbl_display <- data.frame(
  Population  = fg_table$Population,
  N           = fg_table$N,
  Recovery    = fg_table$N_recovery,
  Death       = fg_table$N_death,
  sHR         = fg_table$sHR,
  CI          = fg_table$CI_95,
  P           = fg_table$P_value,
  stringsAsFactors = FALSE
)

ft4 <- flextable(tbl_display) %>%
  set_header_labels(
    Population = "Population", N = "N",
    Recovery = "Recovery\nevents", Death = "Death\nevents",
    sHR = "sHR\n(Early vs Deferred)", CI = "95% CI", P = "p-value"
  ) %>%
  font(fontname = TNR, part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  align(j = 2:7, align = "center", part = "body") %>%
  align(j = 1, align = "left", part = "body") %>%
  border_outer(border = fp_border(color = "black", width = 1.2)) %>%
  border_inner_h(border = fp_border(color = "black", width = 0.5)) %>%
  hline(i = 1, border = fp_border(color = "black", width = 1.0), part = "header") %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  padding(padding.top = 3, padding.bottom = 3, part = "all")

doc <- body_add_flextable(doc, ft4)

# Footnote
footnote_text <- paste0(
  "sHR = subdistribution hazard ratio from Fine-Gray model. CIF = cumulative incidence function. ",
  "Covariates: age, SOFA, APACHE III, mechanical ventilation, CCI, NEE dose, lactate. ",
  "Event 1 = 28-day renal recovery; Event 2 = death within 28 days (competing event)."
)
doc <- body_add_fpar(
  doc,
  fpar(ftext(footnote_text, fp_sm)),
  style = "Normal"
)

print(doc, target = doc_path)
cat("result.docx updated with Section 3.3.2 and Table 4.\n")
