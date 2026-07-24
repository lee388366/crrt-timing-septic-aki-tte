# =============================================================================
# Figure 1A: STROBE Cohort Flow Diagram
# Visual style: solid white background, sans-serif, bold key counts, thin arrows
# =============================================================================

suppressPackageStartupMessages({
  library(grid)
})

fig_dir <- "analysis/figures"

# Sans-serif for PDF/PNG (matches reference flowchart typography)
FONT_SANS <- "sans"

# lines: character vector; bold_mask: logical same length — bold where TRUE
draw_box <- function(x, y, w, h, lines, bold_mask = NULL,
                     cex = 0.62, align = "centre") {
  if (is.null(bold_mask)) {
    bold_mask <- rep(FALSE, length(lines))
  }
  grid.rect(x, y, w, h,
            gp = gpar(fill = "white", col = "black", lwd = 0.9))
  n <- length(lines)
  lh <- convertHeight(unit(1, "lines"), "npc", valueOnly = TRUE) * cex * 1.08
  top <- y + (n - 1) / 2 * lh
  hjust <- if (align == "left") 0 else 0.5
  xpos <- if (align == "left") x - w / 2 + 0.018 else x
  for (i in seq_len(n)) {
    ff <- if (isTRUE(bold_mask[i])) "bold" else "plain"
    grid.text(
      lines[i], xpos, top - (i - 1) * lh, hjust = hjust,
      gp = gpar(
        fontfamily = FONT_SANS, fontsize = 9.2, cex = cex,
        fontface = ff, col = "black"
      )
    )
  }
}

varrow <- function(x, y1, y2) {
  grid.lines(
    x = c(x, x), y = c(y1, y2),
    arrow = arrow(length = unit(0.055, "inches"), type = "closed"),
    gp = gpar(fill = "black", col = "black", lwd = 0.75)
  )
}

harrow <- function(x1, y, x2) {
  grid.lines(
    x = c(x1, x2), y = c(y, y),
    arrow = arrow(length = unit(0.055, "inches"), type = "closed"),
    gp = gpar(fill = "black", col = "black", lwd = 0.75)
  )
}

render_flow <- function() {
  grid.newpage()
  pushViewport(viewport(width = 1, height = 1))
  grid.rect(gp = gpar(fill = "white", col = NA))
  popViewport()

  pushViewport(viewport(width = 0.94, height = 0.94))

  grid.text(
    "Figure 1A. STROBE cohort selection flow",
    x = 0.5, y = 0.99, default.units = "npc", just = "center",
    gp = gpar(fontfamily = FONT_SANS, fontsize = 10.5, fontface = "bold", col = "black")
  )

  # Centered main column + right-side exclusion boxes
  mx <- 0.38
  mw <- 0.30
  ex <- 0.76
  ew <- 0.42
  g <- 0.006

  y1 <- 0.855
  y2 <- 0.72
  y3 <- 0.50
  y5 <- 0.23

  h1 <- 0.055
  h2 <- 0.05
  h3 <- 0.058

  ey1 <- (y1 - h1 / 2 + y2 + h2 / 2) / 2
  ey2 <- (y2 - h2 / 2 + y3 + h3 / 2) / 2

  # Box 1
  draw_box(
    mx, y1, mw, h1,
    c(
      "ICU admission records in MIMIC-IV v3.1 (2008-2022)"
    ),
    bold_mask = c(FALSE)
  )

  varrow(mx, y1 - h1 / 2 - g, y2 + h2 / 2 + g)

  eh1 <- 0.095
  harrow(mx, ey1, ex - ew / 2 - 0.01)
  draw_box(
    ex, ey1, ew, eh1,
    c(
      "Included in target trial emulation (t0 definition):",
      "Adult (>= 18 years)",
      "Development of KDIGO Stage 3 AKI",
      "AKI onset during continuous norepinephrine infusion",
      "AKI onset in [-24 h, +48 h] of Sepsis-3 diagnosis"
    ),
    bold_mask = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    align = "left", cex = 0.56
  )

  draw_box(
    mx, y2, mw, h2,
    c(
      "Met eligibility criteria",
      "(septic shock + KDIGO Stage 3 AKI)",
      "n = 3,731"
    ),
    bold_mask = c(FALSE, FALSE, TRUE)
  )

  varrow(mx, y2 - h2 / 2 - g, y3 + h3 / 2 + g)

  eh2 <- 0.14
  harrow(mx, ey2, ex - ew / 2 - 0.01)
  draw_box(
    ex, ey2, ew, eh2,
    c(
      "Excluded (n = 619)",
      "(1) RRT prior to t0 (n = 21)",
      "(2) ESRD or chronic dialysis (n = 451)",
      "(3) ICU stay < 24 h (n = 21)",
      "(4) Missing key baseline covariates (n = 0)",
      "(5) Zero follow-up time (n = 5)",
      "(6) Non-first ICU admission per patient (n = 121)"
    ),
    bold_mask = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    align = "left", cex = 0.54
  )

  draw_box(
    mx, y3, mw, h3,
    c(
      "Final analytic cohort",
      "(first ICU stay per patient)",
      "n = 3,112"
    ),
    bold_mask = c(FALSE, FALSE, TRUE)
  )

  tw <- 0.225
  th <- 0.05
  tlx <- mx - 0.125
  trx <- mx + 0.125

  draw_box(
    tlx, y5, tw, th,
    c("Early CRRT (<= 24 h of t0)", "n = 500 (16.1%)"),
    bold_mask = c(FALSE, TRUE)
  )
  draw_box(
    trx, y5, tw, th,
    c("Deferred CRRT (no CRRT within 24 h)", "n = 2,612 (83.9%)"),
    bold_mask = c(FALSE, TRUE)
  )

  sy <- y3 - h3 / 2 - g
  fy <- y5 + th / 2 + g
  my <- (sy + fy) / 2
  grid.lines(x = c(mx, mx), y = c(sy, my),
             gp = gpar(col = "black", lwd = 0.75))
  grid.lines(x = c(tlx, trx), y = c(my, my),
             gp = gpar(col = "black", lwd = 0.75))
  varrow(tlx, my, fy)
  varrow(trx, my, fy)

  popViewport()
}

# ── Export ────────────────────────────────────────────────────────────────────
pdf_path <- file.path(fig_dir, "Figure1A_STROBE_Flow.pdf")
pdf(file = pdf_path, width = 8.5, height = 10)
render_flow()
dev.off()
cat("PDF saved:", pdf_path, "\n")

png_path <- file.path(fig_dir, "Figure1A_STROBE_Flow.png")
png(filename = png_path, width = 8.5, height = 10, units = "in", res = 300,
    family = FONT_SANS)
render_flow()
dev.off()
cat("PNG saved:", png_path, "\n")
