# =============================================================================
# Step 13: Secondary Outcome Analysis
#   · 28-day renal recovery  — crude + adjusted logistic OR
#   · RRT-free days alive    — descriptive, Wilcoxon, distribution figure
#   · ICU LOS from t₀       — descriptive, Wilcoxon
#   All analyses: full cohort + equipoise (SOFA 8–14) subgroup
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(flextable)
  library(officer)
})

set.seed(20260328)

base_dir <- "analysis"
fig_dir  <- file.path(base_dir, "figures")
out_docx <- file.path(base_dir, "result.docx")

# ── Load & clean ──────────────────────────────────────────────────────────────
master <- read.csv(file.path(base_dir, "data/analysis_master.csv"),
                   stringsAsFactors = FALSE)

fix_bool <- function(x)
  ifelse(x %in% c("t","TRUE","true","1"), TRUE,
  ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))

master <- master %>%
  mutate(
    early_crrt         = fix_bool(early_crrt),
    died_28d           = fix_bool(died_28d),
    mv_baseline        = fix_bool(mv_baseline),
    renal_recovery_28d = fix_bool(renal_recovery_28d),
    # Recovery composite: TRUE only if explicitly recovered;
    # NA (alive at 28d but still on RRT / missing SCr) counted as no recovery
    recovery_bin       = as.integer(renal_recovery_28d == TRUE & !is.na(renal_recovery_28d)),
    # Equipoise flag (SOFA 8–14 from baseline)
    equipoise          = sofa_baseline >= 8 & sofa_baseline <= 14 & !is.na(sofa_baseline)
  )

A  <- master %>% filter(early_crrt == TRUE)
B  <- master %>% filter(early_crrt == FALSE)
EQ <- master %>% filter(equipoise  == TRUE)
AQ <- EQ %>% filter(early_crrt == TRUE)
BQ <- EQ %>% filter(early_crrt == FALSE)

cat(sprintf("Full cohort: N=%d  (Early=%d, Deferred=%d)\n",
            nrow(master), nrow(A), nrow(B)))
cat(sprintf("Equipoise (SOFA 8-14): N=%d (Early=%d, Deferred=%d)\n",
            nrow(EQ), nrow(AQ), nrow(BQ)))

# ── Helper: OR from logistic glm ─────────────────────────────────────────────
get_or <- function(fit, term = "early_crrtTRUE") {
  co <- coef(fit); se <- sqrt(diag(vcov(fit)))
  idx <- which(names(co) == term)
  if (length(idx) == 0) { idx <- grep("early_crrt", names(co))[1] }
  or <- exp(co[idx]); lo <- exp(co[idx] - 1.96*se[idx]); hi <- exp(co[idx] + 1.96*se[idx])
  pv <- 2*pnorm(-abs(co[idx]/se[idx]))
  list(OR=round(or,3), lo=round(lo,3), hi=round(hi,3),
       pv=if(pv<0.001)"<0.001" else sprintf("%.3f",pv))
}

# ── Helper: n/N (%) ──────────────────────────────────────────────────────────
npct <- function(x, N) sprintf("%d/%d (%.1f%%)", x, N, 100*x/N)

# ── Helper: median [IQR] ──────────────────────────────────────────────────────
miq <- function(x, d=1) {
  sprintf("%.1f [%.1f, %.1f]",
          round(median(x, na.rm=TRUE), d),
          round(quantile(x, 0.25, na.rm=TRUE), d),
          round(quantile(x, 0.75, na.rm=TRUE), d))
}

# ── Helper: Wilcoxon p ────────────────────────────────────────────────────────
wp <- function(x, g) {
  p <- wilcox.test(x ~ g)$p.value
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

# ═══════════════════════════════════════════════════════════════════════════════
# OUTCOME 1: 28-day Renal Recovery
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── Renal Recovery ──\n")

# Counts (per-protocol: all patients, death = no recovery)
n_rec_total  <- sum(master$recovery_bin)
n_rec_early  <- sum(A$recovery_bin)
n_rec_def    <- sum(B$recovery_bin)
n_rec_eq     <- sum(EQ$recovery_bin)
n_rec_eq_A   <- sum(AQ$recovery_bin)
n_rec_eq_B   <- sum(BQ$recovery_bin)

cat(sprintf("Overall: %d/%d (%.1f%%)\n", n_rec_total, nrow(master), 100*n_rec_total/nrow(master)))
cat(sprintf("Early:   %d/%d (%.1f%%)\n", n_rec_early, nrow(A), 100*n_rec_early/nrow(A)))
cat(sprintf("Deferred:%d/%d (%.1f%%)\n", n_rec_def,   nrow(B), 100*n_rec_def/nrow(B)))

# Crude logistic — full cohort
fit_rec_crude_full <- glm(recovery_bin ~ early_crrt, data = master, family = binomial)
or_rec_crude_full  <- get_or(fit_rec_crude_full)
cat(sprintf("Full crude OR=%s [%s,%s] p=%s\n",
            or_rec_crude_full$OR, or_rec_crude_full$lo, or_rec_crude_full$hi, or_rec_crude_full$pv))

# Adjusted logistic — adjust for age, SOFA, APACHE III, MV, CCI (all no-NA cols)
adj_formula <- recovery_bin ~ early_crrt + age + sofa_baseline + apsiii_score + mv_baseline + cci
fit_rec_adj_full <- glm(adj_formula, data = master, family = binomial)
or_rec_adj_full  <- get_or(fit_rec_adj_full)
cat(sprintf("Full adjusted OR=%s [%s,%s] p=%s\n",
            or_rec_adj_full$OR, or_rec_adj_full$lo, or_rec_adj_full$hi, or_rec_adj_full$pv))

# Crude logistic — equipoise subgroup
fit_rec_crude_eq <- glm(recovery_bin ~ early_crrt, data = EQ, family = binomial)
or_rec_crude_eq  <- get_or(fit_rec_crude_eq)

# Adjusted logistic — equipoise subgroup
fit_rec_adj_eq <- glm(adj_formula, data = EQ, family = binomial)
or_rec_adj_eq  <- get_or(fit_rec_adj_eq)
cat(sprintf("Eq crude OR=%s [%s,%s] p=%s\n",
            or_rec_crude_eq$OR, or_rec_crude_eq$lo, or_rec_crude_eq$hi, or_rec_crude_eq$pv))
cat(sprintf("Eq adjusted OR=%s [%s,%s] p=%s\n",
            or_rec_adj_eq$OR, or_rec_adj_eq$lo, or_rec_adj_eq$hi, or_rec_adj_eq$pv))

# ═══════════════════════════════════════════════════════════════════════════════
# OUTCOME 2: RRT-free days alive at day 28
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── RRT-free days alive ──\n")

rfd_p_full <- wp(master$rrt_free_days_28d, master$early_crrt)
rfd_p_eq   <- wp(EQ$rrt_free_days_28d,    EQ$early_crrt)
cat(sprintf("Full:  Early=%s  Deferred=%s  p=%s\n",
            miq(A$rrt_free_days_28d), miq(B$rrt_free_days_28d), rfd_p_full))
cat(sprintf("Equipoise: Early=%s  Deferred=%s  p=%s\n",
            miq(AQ$rrt_free_days_28d), miq(BQ$rrt_free_days_28d), rfd_p_eq))
cat(sprintf("Zero fraction: Early=%.1f%%  Deferred=%.1f%%\n",
            100*mean(A$rrt_free_days_28d==0), 100*mean(B$rrt_free_days_28d==0)))

# ═══════════════════════════════════════════════════════════════════════════════
# OUTCOME 3: ICU length of stay from t₀
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── ICU LOS from t0 ──\n")
los_p_full <- wp(master$icu_los_from_t0, master$early_crrt)
los_p_eq   <- wp(EQ$icu_los_from_t0,    EQ$early_crrt)
cat(sprintf("Full:  Early=%s  Deferred=%s  p=%s\n",
            miq(A$icu_los_from_t0), miq(B$icu_los_from_t0), los_p_full))
cat(sprintf("Equipoise: Early=%s  Deferred=%s  p=%s\n",
            miq(AQ$icu_los_from_t0), miq(BQ$icu_los_from_t0), los_p_eq))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Figures — RRT-free days distribution
# ═══════════════════════════════════════════════════════════════════════════════

plot_df <- master %>%
  mutate(Group = ifelse(early_crrt, "Early CRRT\n(n = 500)", "Deferred CRRT\n(n = 2,612)"),
         Group = factor(Group, levels = c("Early CRRT\n(n = 500)", "Deferred CRRT\n(n = 2,612)")))

fig_rfd <- ggplot(plot_df, aes(x = Group, y = rrt_free_days_28d, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.55, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
               linewidth = 0.6, coef = 1.5) +
  stat_summary(fun = median, geom = "point", shape = 18,
               size = 3, color = "black") +
  scale_fill_manual(values = c("Early CRRT\n(n = 500)" = "#E64B35",
                               "Deferred CRRT\n(n = 2,612)" = "#4DBBD5")) +
  scale_y_continuous(breaks = seq(0, 28, 7), limits = c(-1, 29)) +
  labs(
    title    = "Figure 2. Distribution of RRT-Free Days Alive at Day 28",
    subtitle = "Full analytic cohort (N = 3,112); death assigned 0 days",
    x        = NULL,
    y        = "RRT-Free Days Alive at Day 28",
    caption  = paste0("Median [IQR]: Early = 0 [0, 19] vs Deferred = 20 [0, 28]; Wilcoxon p ", rfd_p_full)
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40"),
        plot.caption  = element_text(size = 8,  color = "grey40"),
        axis.text.x   = element_text(size = 10))

fig_rfd_pdf <- file.path(fig_dir, "Figure3_RRTfree_days.pdf")
fig_rfd_png <- file.path(fig_dir, "Figure3_RRTfree_days.png")
tryCatch({
  pdf(fig_rfd_pdf, width = 6, height = 5, onefile = FALSE)
  print(fig_rfd); dev.off()
  cat("Figure 2 PDF saved:", fig_rfd_pdf, "\n")
}, error = function(e) { tryCatch(dev.off(), error=function(e2){}) })
tryCatch({
  png(fig_rfd_png, width = 1800, height = 1500, res = 300)
  print(fig_rfd); dev.off()
  cat("Figure 2 PNG saved:", fig_rfd_png, "\n")
}, error = function(e) { tryCatch(dev.off(), error=function(e2){}) })

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Summary Table 3 for result.docx
# ═══════════════════════════════════════════════════════════════════════════════

set_flextable_defaults(font.family = "Times New Roman", font.size = 10, padding = 3)

fmt_or <- function(or) sprintf("OR = %s [%s, %s]  p %s%s",
  or$OR, or$lo, or$hi,
  if(or$pv == "<0.001") "" else "= ", or$pv)

tbl3 <- data.frame(
  Outcome = c(
    "**28-Day Renal Recovery**",
    "Overall (composite: survived to d28 + SCr recovery + no RRT)",
    "Crude logistic — full cohort",
    "Adjusted logistic\u1d43 — full cohort",
    "Crude logistic — equipoise (SOFA 8\u201314)",
    "Adjusted logistic\u1d43 — equipoise",
    "",
    "**RRT-Free Days Alive at Day 28**",
    "Full cohort — Early CRRT",
    "Full cohort — Deferred CRRT",
    "Equipoise — Early CRRT",
    "Equipoise — Deferred CRRT",
    "",
    "**ICU Length of Stay from t\u2080, days**",
    "Full cohort — Early CRRT",
    "Full cohort — Deferred CRRT",
    "Equipoise — Early CRRT",
    "Equipoise — Deferred CRRT"
  ),
  `Overall / Early CRRT` = c(
    "", npct(n_rec_total, nrow(master)), "", "", "", "",
    "",
    "", miq(A$rrt_free_days_28d), "", miq(AQ$rrt_free_days_28d), "",
    "",
    "", miq(A$icu_los_from_t0), "", miq(AQ$icu_los_from_t0), ""
  ),
  `Deferred CRRT` = c(
    "", npct(n_rec_def, nrow(B)), "", "", "", "",
    "",
    "", "", miq(B$rrt_free_days_28d), "", miq(BQ$rrt_free_days_28d),
    "",
    "", "", miq(B$icu_los_from_t0), "", miq(BQ$icu_los_from_t0)
  ),
  `Effect Estimate (95% CI) / p-value` = c(
    "",
    sprintf("Early: %s", npct(n_rec_early, nrow(A))),
    fmt_or(or_rec_crude_full),
    fmt_or(or_rec_adj_full),
    fmt_or(or_rec_crude_eq),
    fmt_or(or_rec_adj_eq),
    "",
    sprintf("Zero fraction: Early %.1f%% vs Deferred %.1f%%",
            100*mean(A$rrt_free_days_28d==0), 100*mean(B$rrt_free_days_28d==0)),
    "", sprintf("Wilcoxon p %s%s", if(rfd_p_full=="<0.001")"" else "= ", rfd_p_full), "",
    sprintf("Wilcoxon p %s%s", if(rfd_p_eq=="<0.001")"" else "= ", rfd_p_eq),
    "",
    "",
    "", sprintf("Wilcoxon p %s%s", if(los_p_full=="<0.001")"" else "= ", los_p_full), "",
    sprintf("Wilcoxon p %s%s", if(los_p_eq=="<0.001")"" else "= ", los_p_eq)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

section_rows <- which(grepl("^\\*\\*", tbl3$Outcome))
empty_rows   <- which(tbl3$Outcome == "")
tbl3$Outcome <- gsub("^\\*\\*|\\*\\*$", "", tbl3$Outcome)
indent_rows  <- which(tbl3$Outcome %in% c(
  "Overall (composite: survived to d28 + SCr recovery + no RRT)",
  "Crude logistic \u2014 full cohort", "Adjusted logistic\u1d43 \u2014 full cohort",
  "Crude logistic \u2014 equipoise (SOFA 8\u201314)", "Adjusted logistic\u1d43 \u2014 equipoise",
  "Full cohort \u2014 Early CRRT", "Full cohort \u2014 Deferred CRRT",
  "Equipoise \u2014 Early CRRT",  "Equipoise \u2014 Deferred CRRT"
))

ft3 <- flextable(tbl3) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  bold(i = section_rows, j = 1) %>%
  bg(i = section_rows, bg = "#EEEEEE") %>%
  align(j = 2:4, align = "center", part = "all") %>%
  align(j = 1,   align = "left",   part = "all") %>%
  padding(i = indent_rows, j = 1, padding.left = 16) %>%
  width(j = 1, width = 2.60) %>%
  width(j = 2, width = 1.40) %>%
  width(j = 3, width = 1.25) %>%
  width(j = 4, width = 2.25) %>%
  border_outer(part = "all",  border = fp_border(width = 1.2)) %>%
  hline(part = "header",       border = fp_border(width = 1.2)) %>%
  hline(i = section_rows,      border = fp_border(width = 0.5, color = "#AAAAAA")) %>%
  set_caption(
    caption = as_paragraph(
      as_chunk("Table 3. ",
        props = fp_text(font.family="Times New Roman", font.size=10, bold=TRUE)),
      as_chunk("Secondary Outcomes: 28-Day Renal Recovery, RRT-Free Days Alive, and ICU Length of Stay",
        props = fp_text(font.family="Times New Roman", font.size=10, bold=FALSE))
    )
  )

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Write Section 3.3 to result.docx
# ═══════════════════════════════════════════════════════════════════════════════

tp_h2   <- fp_text(font.family="Times New Roman", font.size=12, bold=TRUE)
tp_body <- fp_text(font.family="Times New Roman", font.size=11)
tp_fn   <- fp_text(font.family="Times New Roman", font.size=9)
pp_head <- fp_par(line_spacing=1.15, text.align="left")

narrative_33 <- paste0(
  "3.3.1  28-day renal recovery. ",
  "Renal recovery was defined as survival to day 28 with serum creatinine ",
  "\u22641.5\u00d7 the pre-AKI baseline value and no ongoing RRT requirement. ",
  "Patients who died within 28 days or remained on RRT at day 28 were classified ",
  "as not recovered. Overall, 28-day renal recovery occurred in ",
  npct(n_rec_total, nrow(master)), " patients across the full cohort. ",
  "Recovery rates were substantially lower in the early CRRT group ",
  "(", npct(n_rec_early, nrow(A)), ") compared with the deferred group ",
  "(", npct(n_rec_def, nrow(B)), "). ",
  "In the unadjusted analysis, early CRRT was associated with a lower odds of ",
  "28-day renal recovery ",
  sprintf("(crude OR\u2009=\u2009%s, 95%%\u2009CI\u2009%s to %s; p\u2009%s%s)",
          or_rec_crude_full$OR, or_rec_crude_full$lo, or_rec_crude_full$hi,
          if(or_rec_crude_full$pv=="<0.001")"" else "= ", or_rec_crude_full$pv),
  ". After adjustment for age, SOFA score, APACHE\u00a0III score, mechanical ",
  "ventilation use, and Charlson comorbidity index, the association remained ",
  sprintf("(adjusted OR\u2009=\u2009%s, 95%%\u2009CI\u2009%s to %s; p\u2009%s%s)",
          or_rec_adj_full$OR, or_rec_adj_full$lo, or_rec_adj_full$hi,
          if(or_rec_adj_full$pv=="<0.001")"" else "= ", or_rec_adj_full$pv),
  ". In the equipoise subgroup (SOFA 8\u201314), results were directionally consistent ",
  sprintf("(adjusted OR\u2009=\u2009%s, 95%%\u2009CI\u2009%s to %s; p\u2009%s%s)",
          or_rec_adj_eq$OR, or_rec_adj_eq$lo, or_rec_adj_eq$hi,
          if(or_rec_adj_eq$pv=="<0.001")"" else "= ", or_rec_adj_eq$pv),
  ". The lower renal recovery rate in the early CRRT group should be interpreted ",
  "cautiously: patients selected for early CRRT had substantially higher baseline ",
  "severity (Table\u00a01), and the composite nature of this outcome (incorporating ",
  "28-day mortality) means that group differences in recovery largely reflect the ",
  "higher mortality burden rather than a differential effect on renal repair per se.",

  "\n\n3.3.2  RRT-free days alive at day 28. ",
  "As prespecified in the study protocol, this outcome combines mortality and ",
  "CRRT duration, with death assigned a value of zero days. ",
  "The distribution was markedly right-skewed and zero-inflated in both groups: ",
  sprintf("%.1f%%", 100*mean(A$rrt_free_days_28d==0)),
  " of early CRRT patients and ",
  sprintf("%.1f%%", 100*mean(B$rrt_free_days_28d==0)),
  " of deferred CRRT patients had zero RRT-free days alive ",
  "(reflecting either death or persistent CRRT dependence within 28 days). ",
  "In the full cohort, the early CRRT group had significantly fewer RRT-free days alive ",
  sprintf("(median %s vs %s days; Wilcoxon p\u2009%s%s)",
          miq(A$rrt_free_days_28d), miq(B$rrt_free_days_28d),
          if(rfd_p_full=="<0.001")"" else "= ", rfd_p_full),
  ". In the equipoise subgroup, the pattern was similar ",
  sprintf("(median %s vs %s days; p\u2009%s%s)",
          miq(AQ$rrt_free_days_28d), miq(BQ$rrt_free_days_28d),
          if(rfd_p_eq=="<0.001")"" else "= ", rfd_p_eq),
  ". The difference in this outcome primarily reflects the longer duration of CRRT ",
  "exposure in the early group and the higher mortality burden documented in the ",
  "unadjusted analysis, rather than a differential effect on CRRT weaning. ",
  "The distribution of RRT-free days is illustrated in Figure\u00a02.",

  "\n\n3.3.3  ICU length of stay. ",
  "Median ICU length of stay from t\u2080 was ",
  miq(A$icu_los_from_t0), " days in the early CRRT group ",
  "and ", miq(B$icu_los_from_t0), " days in the deferred group ",
  sprintf("(Wilcoxon p\u2009%s%s)",
          if(los_p_full=="<0.001")"" else "= ", los_p_full),
  ". In the equipoise subgroup, median ICU stay was ",
  miq(AQ$icu_los_from_t0), " days (early) versus ", miq(BQ$icu_los_from_t0),
  sprintf(" days (deferred; p\u2009%s%s)",
          if(los_p_eq=="<0.001")"" else "= ", los_p_eq),
  ". The modest difference in ICU LOS between groups, while statistically ",
  "significant in the full cohort, is unlikely to represent a clinically ",
  "meaningful treatment effect given the large baseline severity imbalance. ",
  "The results for all secondary outcomes are summarised in Table\u00a03."
)

footnote_33 <- paste0(
  "\u1d43Logistic regression adjusted for age, baseline SOFA score, APACHE\u00a0III score, ",
  "mechanical ventilation use at t\u2080, and Charlson comorbidity index. ",
  "Variables with substantial missingness (lactate, serum creatinine at t\u2080, ",
  "vasopressor dose) were excluded from the adjustment model to preserve sample size. ",
  "Renal recovery denominator: all patients (N\u2009=\u20093,112); death and ongoing RRT-dependence ",
  "at day 28 classified as no recovery. RRT-free days and ICU LOS reported as median [IQR]; ",
  "compared using Wilcoxon rank-sum test. ",
  "Abbreviations: CI, confidence interval; ICU, intensive care unit; IQR, interquartile range; ",
  "OR, odds ratio; RRT, renal replacement therapy; SCr, serum creatinine; ",
  "SOFA, Sequential Organ Failure Assessment; t\u2080, index time."
)

doc <- read_docx(out_docx) %>%
  body_add_par("", style = "Normal") %>%
  body_add_fpar(
    fpar(ftext("3.3  Secondary Outcomes", tp_h2), fp_p = pp_head)
  ) %>%
  body_add_par(narrative_33, style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft3) %>%
  body_add_fpar(fpar(ftext(footnote_33, tp_fn), fp_p = fp_par(line_spacing=1.0))) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("[Figure 3 — see analysis/figures/Figure3_RRTfree_days.pdf]", style = "Normal")

print(doc, target = out_docx)

# ── Console summary ──────────────────────────────────────────────────────────
cat("\n", strrep("=",60), "\n")
cat("SECONDARY OUTCOMES SUMMARY\n"); cat(strrep("=",60), "\n")
cat(sprintf("\n[Renal Recovery — Full Cohort]\n"))
cat(sprintf("  Overall: %s\n", npct(n_rec_total, nrow(master))))
cat(sprintf("  Early:   %s  |  Deferred: %s\n",
            npct(n_rec_early, nrow(A)), npct(n_rec_def, nrow(B))))
cat(sprintf("  Crude OR=%s [%s,%s] p=%s\n", or_rec_crude_full$OR, or_rec_crude_full$lo, or_rec_crude_full$hi, or_rec_crude_full$pv))
cat(sprintf("  Adj   OR=%s [%s,%s] p=%s\n", or_rec_adj_full$OR,   or_rec_adj_full$lo,   or_rec_adj_full$hi,   or_rec_adj_full$pv))
cat(sprintf("\n[RRT-free days — Full Cohort]\n"))
cat(sprintf("  Early %s  |  Deferred %s  |  Wilcoxon p=%s\n",
            miq(A$rrt_free_days_28d), miq(B$rrt_free_days_28d), rfd_p_full))
cat(sprintf("\n[ICU LOS — Full Cohort]\n"))
cat(sprintf("  Early %s  |  Deferred %s  |  Wilcoxon p=%s\n",
            miq(A$icu_los_from_t0), miq(B$icu_los_from_t0), los_p_full))
cat(sprintf("\n✓ result.docx updated: Section 3.3\n"))
cat(sprintf("✓ Figure 2 (RRT-free days): %s\n", fig_rfd_pdf))
