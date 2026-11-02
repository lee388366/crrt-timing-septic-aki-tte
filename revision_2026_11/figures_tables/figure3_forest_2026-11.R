# Figure 3: plotting code verbatim from fig03_NEW_minimal.R used for the submitted figure (forestploter theme, colours, fonts, number formats,
# forest() arguments, arrow labels, x-axis range unchanged); only the row data change (tables/fig3_forest_rows.csv from fig3_rows_2026-11.R);
# canvas as in the submitted PDF (576 x 482 pt = 8.0 x 6.69 in).
suppressPackageStartupMessages({ library(dplyr); library(readr); library(forestploter); library(grid); library(ragg) })
ROOT <- Sys.getenv("CRRT_OUTPUT_DIR")
out_dir <- file.path(ROOT, "figures")

# ---- Constants and functions below are verbatim from fig03_NEW_minimal.R ----
COL_INK <- "#2C3E50"
COL_HEADER <- "#EAF0F6"
COL_GREY <- "grey55"
BASE_FAMILY <- "Helvetica"
VECTOR_FAMILY <- "Helvetica"
fmt_n <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE)
LEGACY_P3 <- c(0.046, 0.364, 0.604)
fmt_p1 <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("< .001")
  for (v in LEGACY_P3) if (round(p, 3) == v) return(sprintf(".%03d", round(p * 1000)))
  r <- round(p, 2)
  if (abs(r - 0.05) < 1e-9 || r == 0 || r >= 0.995) return(sprintf(".%03d", round(p * 1000)))
  sprintf(".%02d", round(p * 100))
}
fmt_p <- function(p) sub(".", "0.", vapply(p, fmt_p1, character(1)), fixed = TRUE)
fmt_hr <- function(hr, lo, hi) {
  ifelse(is.na(hr), "", sprintf("%.3f [%.3f-%.3f]", hr, lo, hi))
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
make_row <- function(group, label, n, early, hr, lo, hi, p) {
  tibble(group = group, label = label, n = n, early = early, hr = hr, lo = lo, hi = hi, p = p, unsourced = FALSE)
}
header_row <- function(label) {
  tibble(group = label, label = label, n = NA_real_, early = NA_real_, hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, unsourced = FALSE)
}

# ---- Row data (the only change) ----
src <- read_csv(file.path(ROOT, "tables/fig3_forest_rows.csv"), show_col_types = FALSE)
forest_rows <- bind_rows(lapply(unique(src$group), function(g) {
  x <- src %>% filter(group == g)
  bind_rows(header_row(g), make_row(x$group, x$label, x$n, x$early, x$hr, x$lo, x$hi, x$p))
}))

# ---- Verbatim from fig03_NEW_minimal.R (build_rows / fp_theme / make_fig3) ----
build_rows <- function(drop_unsourced) {
  rows <- forest_rows
  rows %>% mutate(
    is_header = is.na(hr) & !unsourced,
    Analysis  = if_else(is_header, label, paste0("    ", label)),
    N         = if_else(is.na(n), "", paste0(fmt_n(n), if_else(is.na(early), "", paste0(" / ", fmt_n(early))))),
    HR_CI     = if_else(unsourced, "—  [no source]", fmt_hr(hr, lo, hi)),
    p_disp    = if_else(unsourced, "—", fmt_p(p)),
    blank     = paste(rep(" ", 28), collapse = "")
  )
}
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
rows_A <- build_rows(TRUE)
fp <- make_fig3(rows_A)
save_grid_pair(function() {
  plot(fp)
}, "04_Figure_3_Sensitivity_Forest_2026-11", 8.0, 6.69, 600)
cat("Figure 3 rows:", sum(!is.na(rows_A$hr)), "| saved to", out_dir, "\n")
