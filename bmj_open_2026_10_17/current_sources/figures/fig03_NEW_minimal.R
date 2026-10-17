# =============================================================================
# strictK3 v2 ICM figure regeneration
#
# Purpose:
#   Rebuild manuscript-facing figures from the clean physiologic-KDIGO3 v2
#   analysis, using the MMR/NPG visual specification. Outputs are written to a
#   figures_icm_v2 directory; same-name files are overwritten on rerun.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(survival)
  library(survminer)
  library(forestploter)
  library(grid)
  library(patchwork)
  library(ragg)
})

set.seed(20260716)


# ── CH-01 改动 1/6：路径（原 data/results 只读；输入/输出目录）────────────────────────
base_dir    <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
ch01_dir <- Sys.getenv("CRRT_FIG3_CH01_DIR")
work_dir    <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out")
forest_csv  <- Sys.getenv("FIG3_FOREST_CSV", file.path(work_dir, "run16", "strictK3_v2_forest_plot_source_ch01.csv"))
qaph_csv    <- Sys.getenv("FIG3_QAPH_CSV",   file.path(work_dir, "run14", "NEW", "strictK3_v2_QA_ph_treatment_results.csv"))
out_dir <- Sys.getenv("CRRT_FIG3_OUT")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
COL_EARLY <- "#4DBBD5"
COL_DEFERRED <- "#E64B35"
COL_INK <- "#2C3E50"
COL_HEADER <- "#EAF0F6"
COL_AFTER <- "#27AE60"
COL_THRESHOLD <- "#C0392B"
COL_GREY <- "grey55"
# Arial is preferred by the style guide. This local R build cannot use Arial
# reliably in vector devices, so Helvetica is used as the pre-specified fallback.
BASE_FAMILY <- "Helvetica"
VECTOR_FAMILY <- "Helvetica"

theme_icm <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      text = element_text(color = COL_INK),
      axis.text = element_text(color = COL_INK),
      axis.title = element_text(color = COL_INK),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, lineheight = 1.05),
      legend.title = element_blank(),
      legend.text = element_text(color = COL_INK),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

fmt_n <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE)

# P values follow the manuscript's AMA convention: capital P, no leading zero, 2 decimals.
# Three decimals are kept where 2 would give .00 or .05, and for the three values the
# manuscript itself quotes to 3 decimals (P = .046, .364, .604), so that the figure and the
# text read identically digit for digit.
LEGACY_P3 <- c(0.046, 0.364, 0.604)
fmt_p1 <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("< .001")
  for (v in LEGACY_P3) if (round(p, 3) == v) return(sprintf(".%03d", round(p * 1000)))
  r <- round(p, 2)
  if (abs(r - 0.05) < 1e-9 || r == 0 || r >= 0.995) return(sprintf(".%03d", round(p * 1000)))
  sprintf(".%02d", round(p * 100))
}
# CH26 mechanical format: retain original digits, add omitted leading zero only.
fmt_p <- function(p) sub(".", "0.", vapply(p, fmt_p1, character(1)), fixed = TRUE)

fmt_hr <- function(hr, lo, hi) {
  ifelse(is.na(hr), "", sprintf("%.3f [%.3f-%.3f]", hr, lo, hi))
}

save_gg_pair <- function(plot, stem, width, height, dpi) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  invisible(dev.off())
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = VECTOR_FAMILY,
                 bg = "white")
  print(plot)
  invisible(dev.off())
  c(png = png_path, pdf = pdf_path)
}
read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(x %in% c(TRUE, "TRUE", "true", "t", "1", 1), TRUE,
         ifelse(x %in% c(FALSE, "FALSE", "false", "f", "0", 0), FALSE, NA))
}

save_grid_pair <- function(draw_fun, stem, width, height, dpi) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  draw_fun()
  invisible(dev.off())
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = VECTOR_FAMILY,
                 bg = "white")
  draw_fun()
  invisible(dev.off())
  c(png = png_path, pdf = pdf_path)
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) { warning("missing input: ", path); return(tibble::tibble()) }
  read_csv_quiet(path)
}

# ── CH-01 改动 2/6：数据源（主分析值来自 CH-01 产物；无源则留空，不编值）──────────────
primary_results <- read_csv_quiet(file.path(ch01_dir, "ch01_primary_dod.csv")) %>%
  transmute(
    analysis    = if_else(population == "Equipoise SOFA 8-14",
                          "Primary SOFA 8-14 equipoise", "Full strictK3_v2 primary cohort"),
    n_patients  = n_patients,
    early_n     = NA_real_,
    weighted_hr = hr, ci_lower = ci_lower, ci_upper = ci_upper,
    p_value     = p_value
  )
forest_source <- read_csv_if_exists(forest_csv)
# 无源分析（lead 2026-10-16 裁决③）：path-adjusted 与 5 条 rescue 行已按 CH-05/CH-06 删除，作者未裁决
empty_path <- tibble::tibble(analysis = character(), n_patients = numeric(), early_n = numeric(), HR = numeric(), CI_lo = numeric(), CI_hi = numeric(), p = numeric())
empty_rescue <- tibble::tibble(analysis = character(), n_patients = numeric(), early_n = numeric(), weighted_hr = numeric(), ci_lower = numeric(), ci_upper = numeric(), p_value = numeric())
# 无源行于 RUN-17 恢复（作者裁决后 data-runner 交付，全为 CH-01 新口径）；文件缺失则维持"无源"行为
path_rescue_csv <- Sys.getenv("FIG3_PATH_RESCUE_CSV", file.path(work_dir, "run17", "行值_forest_接口_15行.csv"))
pr_all <- read_csv_if_exists(path_rescue_csv)
if (nrow(pr_all) == 0) {
  path_adjusted <- empty_path
  rescue_results <- empty_rescue
} else {
  path_adjusted <- pr_all %>% filter(analysis == "SOFA 8-14: IPCW includes stage3_path") %>%
    transmute(analysis, n_patients, early_n, HR = weighted_hr, CI_lo = ci_lower, CI_hi = ci_upper, p = p_value)
  rescue_keys <- c("no_rescue_exclusion: SOFA 8-14", "sensitivity_R1_confirm_R2_R2b_R3: SOFA 8-14",
                   "sensitivity_R1_R2_R2b_R3_R4: SOFA 8-14",
                   "sensitivity_liberal_thresholds_R1_R2_R2b_R3: SOFA 8-14",
                   "sensitivity_conservative_thresholds_R1_R2_R2b_R3: SOFA 8-14")
  rescue_results <- pr_all %>% filter(analysis %in% rescue_keys) %>%
    transmute(analysis, n_patients, early_n, weighted_hr, ci_lower, ci_upper, p_value)
  cat("[FIG-03] run17 恢复行: path", nrow(path_adjusted), "| rescue", nrow(rescue_results), "\n")
}
qa_ph <- read_csv_if_exists(qaph_csv)
l0_results <- forest_source %>% filter(grepl("^L0-only denominator:", analysis)) %>%
  transmute(analysis, n_patients, early_n, HR = weighted_hr, CI_lo = ci_lower, CI_hi = ci_upper, p = p_value)
forest_source <- forest_source %>% filter(!grepl("^L0-only denominator:", analysis))

# ── CH-01 改动 3/6：无源行占位（不编值、不用旧值）────────────────────────────────────
unsourced_row <- function(group, display, key) {
  tibble(group = group, label = paste0(display, "  [NO SOURCE]"), n = NA_real_, early = NA_real_,
         hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, unsourced = TRUE)
}
# -----------------------------------------------------------------------------
# Figure 3: primary and robustness forest (CH-01 = 28-day all-cause mortality)
# -----------------------------------------------------------------------------

make_row <- function(group, label, n, early, hr, lo, hi, p) {
  tibble(group = group, label = label, n = n, early = early, hr = hr, lo = lo, hi = hi, p = p, unsourced = FALSE)
}
header_row <- function(label) {
  tibble(group = label, label = label, n = NA_real_, early = NA_real_, hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, unsourced = FALSE)
}

row_from_forest <- function(analysis_key, display, group = NA_character_) {
  x <- forest_source %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) return(unsourced_row(group, display, analysis_key))
  make_row(
    group = ifelse(is.na(group), x$group, group),
    label = display,
    n = x$n_patients,
    early = x$early_n,
    hr = x$weighted_hr,
    lo = x$ci_lower,
    hi = x$ci_upper,
    p = x$p_value
  )
}

row_from_rescue <- function(analysis_key, display, group = "Urgent-indication criterion") {
  x <- rescue_results %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) return(unsourced_row(group, display, analysis_key))
  make_row(group, display, x$n_patients, x$early_n, x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
}

row_from_ph <- function(model_key, display, group = "Model specification") {
  x <- qa_ph %>% filter(model == model_key) %>% slice(1)
  if (nrow(x) != 1) return(unsourced_row(group, display, model_key))
  # CH-10 改动（作者 2026-10-16 裁决）：原为硬编码 NA_real_ ⇒ 图上该格只印 N、缺 early 数
  #   现要求源件必须带 early_n；缺列即停，不再静默印出空白格（"让错写不出来"）
  if (!"early_n" %in% names(x))
    stop("[FIG-03] ", model_key, " 的源件缺 early_n 列 ⇒ 拒绝出图（不编值、不留空）")
  make_row(group, display, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p_value)
}

row_from_path <- function(analysis_key, display, group = "Model specification") {
  x <- path_adjusted %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) return(unsourced_row(group, display, analysis_key))
  make_row(group, display, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p)
}

row_from_l0 <- function(analysis_key, display, group = "Model specification") {
  x <- l0_results %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) return(unsourced_row(group, display, analysis_key))
  make_row(group, display, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p)
}

forest_rows <- bind_rows(
  header_row("Primary and cohort analyses"),
  row_from_forest("Primary SOFA 8-14 equipoise", "Primary SOFA 8-14 (IPCW-weighted)"),
  row_from_forest("Full strictK3_v2 primary cohort", "Full cohort"),
  header_row("SOFA score strata"),
  row_from_forest("SOFA <=7", "SOFA ≤7"),
  row_from_forest("SOFA 8-14", "SOFA 8-14"),
  row_from_forest("SOFA 9-13", "SOFA 9-13"),
  row_from_forest("SOFA >=15", "SOFA ≥15"),
  # CH-08 改动：L0-only 归入（原被硬编码塞进 rescue 组 ⇒ 分组错误）
  # CH-10 改动（作者 2026-10-16 裁决）：原「Model and t₀ sensitivity」一组拆为两组
  #   依据：正文 L58 自身把「calendar-era models / alternative eligibility anchors /
  #         IPCW-denominator respecification」分三类讲；补材亦单有一节
  #         「Sensitivity analysis — Eligibility anchor」⇒ 一拆二与稿件口径一致
  #   副效果：锚点两行不再各重复一次 "anchor"（组头已说）
  #   风险：锚点两行上移到「Eligibility anchor」前，L0-only 随之移到组一末尾（行序变化，数值不动）
  header_row("Model specification"),
  row_from_path("SOFA 8-14: IPCW includes stage3_path", "KDIGO 3 pathway adjusted"),
  row_from_ph("SOFA 8-14 era-stratified", "Calendar-era stratified"),
  row_from_l0("L0-only denominator: SOFA 8-14", "L0-only IPCW denominator"),
  header_row("Eligibility anchor"),
  row_from_forest("KDIGO3_minus6h: SOFA 8-14", "6 h before KDIGO stage 3",
                  group = "Eligibility anchor"),
  row_from_forest("KDIGO2: SOFA 8-14", "At KDIGO stage 2",
                  group = "Eligibility anchor"),
  # CH-09 改动（作者 2026-10-16 裁决）：组名与 5 行标签一并收紧
  #   ① 组名 41 字符 → 26 字符，且改用正文 L135 自己的名词「urgent-indication criterion」
  #   ② 行内不再重复组名已有的词：urgent-indication 由 4 次降为 1 次（只在组头）
  #   ③ 5 行统一为「名词 + 过去分词」，修掉原第 2 行的祈使句（Hyperkalemia must be reconfirmed）
  #      与第 3 行含糊的 also（Urine-output criterion also excluded）
  header_row("Urgent-indication criterion"),
  row_from_rescue("no_rescue_exclusion: SOFA 8-14", "No exclusion applied"),
  row_from_rescue("sensitivity_R1_confirm_R2_R2b_R3: SOFA 8-14", "Hyperkalemia reconfirmed"),
  row_from_rescue("sensitivity_R1_R2_R2b_R3_R4: SOFA 8-14", "Urine-output criterion added"),
  row_from_rescue("sensitivity_liberal_thresholds_R1_R2_R2b_R3: SOFA 8-14", "Cut-offs lowered"),
  row_from_rescue("sensitivity_conservative_thresholds_R1_R2_R2b_R3: SOFA 8-14", "Cut-offs raised")
) %>%
  mutate(
    is_header = is.na(hr),
    Analysis = if_else(is_header, label, paste0("    ", label)),
    N = if_else(is_header, "", paste0(fmt_n(n), if_else(is.na(early), "", paste0(" / ", fmt_n(early))))),
    `HR [95% CI]` = fmt_hr(hr, lo, hi),
    `P value` = fmt_p(p),
    ` ` = paste(rep(" ", 28), collapse = "")
  )

forest_display <- forest_rows %>%
  transmute(Analysis, `N / early` = N, ` `, `HR [95% CI]`, `P value`)

# Locked display table for Fig3 ALIGN / manuscript (3-dp fmt_hr)
write_csv(
  forest_rows %>%
    filter(!is_header) %>%
    transmute(
      label, n, early,
      weighted_hr = hr, ci_lower = lo, ci_upper = hi, p_value = p,
      hr_display = `HR [95% CI]`, p_display = `P value`, n_display = N
    ),
  file.path(out_dir, "strictK3_v2_fig3_forest_display_raw.csv")
)

fp_theme <- forest_theme(
  base_size = 8.3,
  base_family = BASE_FAMILY,
  ci_pch = 18,
  ci_col = COL_INK,
  ci_lwd = 1.08,
  ci_Theight = 0.16,
  refline_gp = gpar(lwd = 0.9, lty = "dashed", col = COL_GREY),
  summary_col = COL_INK,
  summary_fill = COL_HEADER,
  arrow_type = "closed",
  arrow_label_just = "end",
  arrow_gp = gpar(lwd = 0.7, col = "grey35", fontsize = 7),
  xlab_adjust = "refline",
  xlab_gp = gpar(fontsize = 8, fontface = "plain"),
  xaxis_gp = gpar(fontsize = 7.8)
)

fig3 <- forest(
  data = forest_display,
  est = forest_rows$hr,
  lower = forest_rows$lo,
  upper = forest_rows$hi,
  sizes = ifelse(forest_rows$is_header, NA, 0.42),
  ref_line = 1,
  ci_column = 3,
  xlim = c(0.6, 1.6),
  ticks_at = c(0.75, 1.00, 1.25, 1.50),
  xlab = "Hazard Ratio (Early CRRT vs No CRRT within 24 h)",
  arrow_lab = c("Favours Early CRRT", "Favours No CRRT within 24 h"),
  is_summary = forest_rows$is_header,
  theme = fp_theme
)

# CH-11 改动（作者 2026-10-17 裁决）：
#   ① 底部说明灰字**三版全去**（v2/v3/v4 不再印）
#      原 v2 印 "Note: SOFA 8-14 is repeated in the strata panel …"；v3/v4 印 note_A
#      （其中 v4 那句含内部工单代号 "CH-01 source"，不宜入图）。
#   ② 画布收边 8.3×6.55 → **6.2×5.4 in**（作者：周围留白太多）
#      实测依据（.scratch/e16/ws 十档试跑）：forestploter 表格为**绝对尺寸**，
#      内容恒为 3413×2971 px（5.69×4.95 in），画布只决定外圈白边；
#      6.2×5.4 ⇒ 内容占比 91.7%×91.7%，四边余量 0.23–0.26 in，且不发生裁切
#      （6.2×4.9 已会把上下切掉 ⇒ 高度下限 5.0 in）。
#   最小改动：**只改画布尺寸**，全部数值/列宽/字号/行序均不动。
FIG3_CANVAS_W <- 6.2
FIG3_CANVAS_H <- 5.4
fig3_paths <- save_grid_pair(
  function() {
    plot(fig3)
  },
  "Figure_3_Primary_Robustness_Forest_v2",
  FIG3_CANVAS_W, FIG3_CANVAS_H, 600
)

# ── CH-01 改动 4/6：两种标注变体（lead 裁决③：留白形式不由我自定，两版并列供选）──────
build_rows <- function(drop_unsourced) {
  rows <- forest_rows
  if (drop_unsourced) {
    rows <- rows %>% filter(!unsourced)
    # 组内的表头行按位置保留；某组数据行全空则连表头一起去掉
    hdr <- which(rows$is_header)
    keep <- rep(TRUE, nrow(rows))
    for (i in seq_along(hdr)) {
      s <- hdr[i]
      e <- if (i < length(hdr)) hdr[i + 1] - 1 else nrow(rows)
      if (e <= s) keep[s:e] <- FALSE
    }
    rows <- rows[keep, , drop = FALSE]
  }
  rows %>% mutate(
    is_header = is.na(hr) & !unsourced,
    Analysis  = if_else(is_header, label, paste0("    ", label)),
    N         = if_else(is.na(n), "", paste0(fmt_n(n), if_else(is.na(early), "", paste0(" / ", fmt_n(early))))),
    HR_CI     = if_else(unsourced, "\u2014  [no source]", fmt_hr(hr, lo, hi)),
    p_disp    = if_else(unsourced, "\u2014", fmt_p(p)),
    blank     = paste(rep(" ", 28), collapse = "")
  )
}
make_fig3 <- function(rows) {
  d <- rows %>% transmute(Analysis = Analysis, N_early = N, blank = blank, HR_CI = HR_CI, p_disp = p_disp)
  names(d) <- c("Analysis", "N / early", " ", "HR [95% CI]", "P value")
  forest(
    data = d, est = rows$hr, lower = rows$lo, upper = rows$hi,
    sizes = ifelse(rows$is_header | is.na(rows$hr), NA, 0.42),
    ref_line = 1, ci_column = 3, xlim = c(0.6, 1.6), ticks_at = c(0.75, 1.00, 1.25, 1.50),
    xlab = "Hazard Ratio (Early CRRT vs No CRRT within 24 h)",
    arrow_lab = c("Favours Early CRRT", "Favours No CRRT within 24 h"),
    is_summary = rows$is_header, theme = fp_theme
  )
}

# ── CH-01 改动 5/6：输出（无源行存在时装 v3 仅源行版；恢复后自动出 v4 全行版）──────────
unsourced_labels <- forest_rows %>% filter(unsourced) %>% pull(label)
tag <- if (length(unsourced_labels) == 0) "v4_allRows" else "v3_sourceRowsOnly"
note_A <- if (length(unsourced_labels) == 0) {
  "Note: all 16 analyses with a CH-01 source are shown (15 from the CH-01 forest source plus the calendar-era stratified model). Outcome: 28-day all-cause mortality."
} else {
  "Note: this version shows the 10 analyses with a CH-01 source; 6 further rows (stage-3 path-adjusted IPCW and 5 rescue/eligibility variants) are pending author adjudication and are shown only in the accompanying pending-rows list. Outcome: 28-day all-cause mortality."
}

# CH-11 改动（作者 2026-10-17 裁决）：底部灰字一并去掉（见 v2 段落同款说明）。
#   note_A 的定义保留（其文本仍写进 manifest / 台账作留痕），只是**不再印到图上**。
rows_A <- build_rows(TRUE)
fp <- make_fig3(rows_A)
save_grid_pair(function() {
  plot(fp)
}, paste0("Figure_3_Primary_Robustness_Forest_", tag), FIG3_CANVAS_W, FIG3_CANVAS_H, 600)
readr::write_csv(rows_A, file.path(out_dir, paste0("strictK3_v2_fig3_forest_display_", tag, ".csv")))
cat("\n[FIG-03]", tag, "| 图内行:", sum(!is.na(rows_A$hr)), "| 无源未画:", length(unsourced_labels), "\n")
if (length(unsourced_labels)) cat("[FIG-03] 未画的行:", paste(unsourced_labels, collapse = " | "), "\n")

# ── CH-01 改动 6/6：图外「待补行」清单（仅当仍有未恢复行时才复制）─────────────────────
pending_src <- file.path(work_dir, "run15", "无源行清单.csv")
if (length(unsourced_labels) > 0 && file.exists(pending_src)) {
  readr::write_csv(read_csv_quiet(pending_src), file.path(out_dir, "Figure3_待补行清单_6行.csv"))
  cat("[FIG-03] 待补行清单已复制\n")
}
