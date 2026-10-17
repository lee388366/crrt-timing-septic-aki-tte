# =============================================================================
# Figure 1 (v2 clean-t₀) — adapted DIRECTLY from the mother-template source
# MMR/scripts/mmr_generate_figure1_flow_largefont.R
#
# Style is preserved verbatim (Arial, black & white, plain boxes/arrows,
# NO in-figure title — the title lives only in the docx caption).
# Only the funnel content/counts are updated to the v2 clean physiologic
# KDIGO stage 3 cohort, and a dedicated rescue-exclusion side box is added.
#
# v2 funnel (data/strictK3_primary_v2_funnel.csv):
#   first physiologic KDIGO3 stays        14,625
#   septic shock (Sepsis-3) + KDIGO3       3,728
#   pre-eligibility exclusions (n=889)  -> 2,839
#   rescue/urgent indication at t₀ (376)-> 2,463  (final full cohort)
#   primary SOFA-at-t₀ 8-14                1,502  (early 146 / deferred 1,356)
#   observed CRRT groups (full): early 195 (7.9%) / no CRRT 2,268 (92.1%)
# =============================================================================
suppressPackageStartupMessages({ library(grid); library(ragg) })

out_dir <- file.path(Sys.getenv("CRRT_MANUSCRIPT_ROOT"), "11_母板重建_数据回补/图件源件_figure_source")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FONT_SANS <- "Arial"
CANVAS_WIDTH_IN <- 7.48
CANVAS_HEIGHT_IN <- 8.90
CANVAS_DPI <- 1000

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# --- primitives copied verbatim from the mother-template script ---
draw_box <- function(x, y, w, h, lines, bold_mask = NULL,
                     fontsize = 10.5, lineheight = 1.08,
                     align = "centre", pad = 0.018) {
  if (is.null(bold_mask)) bold_mask <- rep(FALSE, length(lines))
  grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "black", lwd = 0.85))
  n <- length(lines)
  lh <- convertHeight(unit(fontsize * lineheight, "points"), "npc", valueOnly = TRUE)
  top <- y + (n - 1) / 2 * lh
  hjust <- if (align == "left") 0 else 0.5
  xpos <- if (align == "left") x - w / 2 + pad else x
  for (i in seq_len(n)) {
    grid.text(lines[i], xpos, top - (i - 1) * lh, hjust = hjust,
              gp = gpar(fontfamily = FONT_SANS, fontsize = fontsize,
                        fontface = if (isTRUE(bold_mask[i])) "bold" else "plain",
                        col = "black"))
  }
}
draw_varrow <- function(x, y1, y2) {
  grid.lines(x = c(x, x), y = c(y1, y2),
             arrow = arrow(length = unit(0.065, "inches"), type = "closed"),
             gp = gpar(fill = "black", col = "black", lwd = 0.80))
}
draw_harrow <- function(x1, y, x2) {
  grid.lines(x = c(x1, x2), y = c(y, y),
             arrow = arrow(length = unit(0.065, "inches"), type = "closed"),
             gp = gpar(fill = "black", col = "black", lwd = 0.80))
}

draw_flow <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(width = 0.985, height = 0.985))

  mx <- 0.28; mw <- 0.42
  ex <- 0.75; ew <- 0.46
  g <- 0.006

  y1 <- 0.930; y2 <- 0.755; y3 <- 0.545; y4 <- 0.335; y5 <- 0.110
  h1 <- 0.070; h2 <- 0.086; h3 <- 0.076; h4 <- 0.066; h5 <- 0.090

  ey1 <- (y1 - h1/2 + y2 + h2/2) / 2
  ey2 <- (y2 - h2/2 + y3 + h3/2) / 2
  ey3 <- (y3 - h3/2 + y4 + h4/2) / 2
  ey4 <- (y4 - h4/2 + y5 + h5/2) / 2

  # B1
  draw_box(mx, y1, mw, h1,
    c("First physiologic KDIGO stage 3 AKI ICU stays",
      "MIMIC-IV v3.1 (2008-2022)", "n = 14,625"),
    bold_mask = c(FALSE, FALSE, TRUE), fontsize = 10.2)
  draw_varrow(mx, y1 - h1/2 - g, y2 + h2/2 + g)

  # side 1: inclusion criteria (t₀ definition)
  draw_harrow(mx, ey1, ex - ew/2 - 0.01)
  draw_box(ex, ey1, ew, 0.150,
    c("Included in target trial emulation (clean t₀):",
      "Adult (>= 18 years)",
      "Septic shock (Sepsis-3) at or before t₀",
      "Index time t₀ = first physiologic KDIGO stage 3",
      "(creatinine and/or urine-output criteria),",
      "bounded by ICU admission; 24-h grace period"),
    bold_mask = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    align = "left", fontsize = 9.0, lineheight = 1.05)

  # B2
  draw_box(mx, y2, mw, h2,
    c("Met eligibility criteria", "septic shock + first physiologic",
      "KDIGO stage 3 AKI", "n = 3,728"),
    bold_mask = c(FALSE, FALSE, FALSE, TRUE), fontsize = 10.2)
  draw_varrow(mx, y2 - h2/2 - g, y3 + h3/2 + g)

  # side 2: pre-eligibility exclusions (n = 889)
  draw_harrow(mx, ey2, ex - ew/2 - 0.01)
  draw_box(ex, ey2, ew, 0.190,
    c("Excluded (n = 889)",
      "(1) RRT prior to t₀ (n = 419)",
      "(2) ESRD or chronic dialysis (n = 350)",
      "(3) ICU stay < 24 h (n = 14)",
      "(4) Missing key baseline covariates (n = 0)",
      "(5) Zero 28-day follow-up (n = 8)",
      "(6) Non-first ICU admission per patient (n = 98)"),
    bold_mask = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    align = "left", fontsize = 9.0, lineheight = 1.05)

  # B3
  draw_box(mx, y3, mw, h3,
    c("Eligible first-ICU cohort", "(before urgent-indication exclusion)", "n = 2,839"),
    bold_mask = c(FALSE, FALSE, TRUE), fontsize = 10.2)
  draw_varrow(mx, y3 - h3/2 - g, y4 + h4/2 + g)

  # side 3: rescue / urgent indication exclusion (n = 376)  <-- NEW
  draw_harrow(mx, ey3, ex - ew/2 - 0.01)
  draw_box(ex, ey3, ew, 0.150,
    c("Excluded - predefined biochemical proxies",
      "for urgent dialysis indication at t₀ (n = 376)",
      "R1 hyperkalemia (K+ >= 6.5)",
      "R2 acidosis (pH <= 7.15 or HCO3- <= 12)",
      "R3 azotemia (BUN >= 100)"),
    bold_mask = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    align = "left", fontsize = 9.0, lineheight = 1.05)

  # B4
  draw_box(mx, y4, mw, h4,
    c("Final full analytic cohort", "n = 2,463"),
    bold_mask = c(FALSE, TRUE), fontsize = 10.2)
  draw_varrow(mx, y4 - h4/2 - g, y5 + h5/2 + g)

  # side 4: observed CRRT groups (descriptive)
  draw_harrow(mx, ey4, ex - ew/2 - 0.01)
  draw_box(ex, ey4, ew, 0.130,
    c("Observed CRRT initiation groups", "(descriptive only)",
      "Early CRRT (<= 24 h of t₀): n = 195 (7.9%)",
      "Deferred (no CRRT within 24 h): n = 2,268 (92.1%)"),
    bold_mask = c(TRUE, FALSE, FALSE, FALSE),
    align = "left", fontsize = 9.0, lineheight = 1.06)

  # B5
  draw_box(mx, y5, mw, h5,
    c("Primary SOFA-at-t₀ subgroup", "SOFA at t₀ 8-14", "n = 1,502",
      "Early CRRT n = 146; deferred n = 1,356"),
    bold_mask = c(FALSE, FALSE, TRUE, FALSE), fontsize = 10.2)

  popViewport()
}

png_path <- file.path(out_dir, "Figure_1_Cohort_Flow_v2_MMRstyle.png")
agg_png(png_path, width = CANVAS_WIDTH_IN, height = CANVAS_HEIGHT_IN,
        units = "in", res = CANVAS_DPI)
draw_flow(); dev.off()

pdf_path <- file.path(out_dir, "Figure_1_Cohort_Flow_v2_MMRstyle.pdf")
FONT_SANS <- "Helvetica"
pdf(pdf_path, width = CANVAS_WIDTH_IN, height = CANVAS_HEIGHT_IN,
    family = "Helvetica", useDingbats = FALSE)
draw_flow(); dev.off()

message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
