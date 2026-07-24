# =============================================================================
# Supplementary Figure S7: CCW/IPCW methodology schematic
#
# Conceptual figure only. No patient-level data are read.
# Output:
#   figures_icm_v2/Supplement_Figure_S7_CCW_IPCW_v2.png
#   figures_icm_v2/Supplement_Figure_S7_CCW_IPCW_v2.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(ragg)
})

base_dir <- "analysis/strictK3_primary_v2"
out_dir <- file.path(base_dir, "figures_icm_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL_EARLY <- "#4DBBD5"
COL_DEFERRED <- "#E64B35"
COL_INK <- "#2C3E50"
COL_HEADER <- "#EAF0F6"
COL_CENSOR <- "#C0392B"
COL_EVENT <- "#8B1A1A"
COL_EVENT_DARK <- "#3A3A3A"
COL_GREY <- "grey65"
COL_LIGHT <- "#F7F9FB"
BASE_FAMILY <- "Helvetica"
VECTOR_FAMILY <- "Helvetica"

stem <- "Supplement_Figure_S7_CCW_IPCW_v2"

gp_text <- function(size = 8.5, face = "plain", col = COL_INK) {
  gpar(fontfamily = BASE_FAMILY, fontsize = size, fontface = face, col = col)
}

u <- function(x) unit(x, "npc")

panel <- function(x, y, w, h, code, title, subtitle = NULL) {
  grid.roundrect(
    x = u(x), y = u(y), width = u(w), height = u(h),
    r = unit(0.025, "npc"),
    gp = gpar(fill = "white", col = COL_INK, lwd = 1.1)
  )
  tag_w <- 0.052
  tag_h <- 0.052
  grid.roundrect(
    x = u(x - w / 2 + tag_w / 2 + 0.012),
    y = u(y + h / 2 - tag_h / 2 - 0.014),
    width = u(tag_w), height = u(tag_h), r = unit(0.010, "npc"),
    gp = gpar(fill = COL_INK, col = COL_INK, lwd = 0)
  )
  grid.text(code,
            x = u(x - w / 2 + tag_w / 2 + 0.012),
            y = u(y + h / 2 - tag_h / 2 - 0.014),
            gp = gp_text(8.5, "bold", "white"))
  grid.text(title,
            x = u(x - w / 2 + tag_w + 0.026),
            y = u(y + h / 2 - 0.027),
            just = c("left", "center"),
            gp = gp_text(10.2, "bold"))
  if (!is.null(subtitle)) {
    grid.text(subtitle,
              x = u(x - w / 2 + 0.020),
              y = u(y + h / 2 - 0.082),
              just = c("left", "center"),
              gp = gpar(fontfamily = BASE_FAMILY, fontsize = 6.5, col = COL_INK, lineheight = 0.93))
  }
}

draw_arrow <- function(x0, y0, x1, y1, col = COL_INK, lwd = 1.0, lty = 1) {
  grid.lines(
    x = u(c(x0, x1)), y = u(c(y0, y1)),
    arrow = arrow(length = unit(0.07, "inches"), type = "closed"),
    gp = gpar(col = col, lwd = lwd, lty = lty)
  )
}

draw_line <- function(x0, y0, x1, y1, col = COL_INK, lwd = 2.0, lty = 1) {
  grid.lines(x = u(c(x0, x1)), y = u(c(y0, y1)),
             gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round"))
}

draw_tick <- function(x, y, h = 0.016, col = COL_INK) {
  grid.lines(x = u(c(x, x)), y = u(c(y - h / 2, y + h / 2)),
             gp = gpar(col = col, lwd = 0.7))
}

draw_square <- function(x, y, size = 0.011) {
  grid.rect(x = u(x), y = u(y), width = u(size), height = u(size),
            gp = gpar(fill = COL_EVENT_DARK, col = COL_EVENT_DARK, lwd = 0.7))
}

draw_death <- function(x, y, r = 0.0062) {
  grid.circle(x = u(x), y = u(y), r = u(r),
              gp = gpar(fill = COL_EVENT, col = COL_EVENT, lwd = 0.7))
}

draw_censor_x <- function(x, y, size = 0.014, lwd = 1.55) {
  grid.lines(x = u(c(x - size / 2, x + size / 2)), y = u(c(y - size / 2, y + size / 2)),
             gp = gpar(col = COL_CENSOR, lwd = lwd, lineend = "round"))
  grid.lines(x = u(c(x - size / 2, x + size / 2)), y = u(c(y + size / 2, y - size / 2)),
             gp = gpar(col = COL_CENSOR, lwd = lwd, lineend = "round"))
}

draw_person <- function(x, y, col = COL_INK) {
  grid.circle(x = u(x), y = u(y + 0.018), r = u(0.010),
              gp = gpar(fill = col, col = col))
  grid.roundrect(x = u(x), y = u(y - 0.002), width = u(0.032), height = u(0.030),
                 r = unit(0.006, "npc"), gp = gpar(fill = col, col = col, lwd = 0))
}

time_x <- function(x0, x1, t_hours) {
  # Nonlinear visual map: t0=0, 24h=0.35 of axis, day28=1.0.
  # Within the 24-h grace period, preserve hour-level alignment.
  ifelse(t_hours <= 24,
         x0 + (x1 - x0) * 0.35 * (t_hours / 24),
         x0 + (x1 - x0) * (0.35 + 0.65 * ((t_hours - 24) / (28 * 24 - 24))))
}

draw_time_axis <- function(x0, x1, y, label_grace = TRUE) {
  draw_arrow(x0, y, x1, y, col = COL_INK, lwd = 0.8)
  x24 <- time_x(x0, x1, 24)
  draw_tick(x0, y)
  draw_tick(x24, y)
  draw_tick(x1, y)
  grid.text("t0", x = u(x0), y = u(y - 0.030), gp = gp_text(6.9))
  grid.text("24 h", x = u(x24), y = u(y - 0.030), gp = gp_text(6.9))
  grid.text("Day 28", x = u(x1), y = u(y - 0.030), gp = gp_text(6.9))
  if (label_grace) {
    yb <- y + 0.035
    grid.lines(x = u(c(x0, x24)), y = u(c(yb, yb)), gp = gpar(col = COL_INK, lwd = 0.75))
    grid.lines(x = u(c(x0, x0)), y = u(c(yb - 0.007, yb + 0.007)), gp = gpar(col = COL_INK, lwd = 0.75))
    grid.lines(x = u(c(x24, x24)), y = u(c(yb - 0.007, yb + 0.007)), gp = gpar(col = COL_INK, lwd = 0.75))
    grid.text("24-h grace period", x = u((x0 + x24) / 2), y = u(yb + 0.018), gp = gp_text(6.6))
    grid.text("follow-up", x = u((x24 + x1) / 2), y = u(yb + 0.018), gp = gp_text(6.6))
  }
}

draw_panel_a <- function() {
  panel(0.18, 0.545, 0.29, 0.70, "A", "CLONE",
        "Each eligible patient is copied into both strategy arms.")
  x_left <- 0.070
  x_t0 <- 0.132
  x_mid <- 0.185
  x_right <- 0.302
  y_patient <- 0.705
  y_early <- 0.620
  y_defer <- 0.525

  draw_person(x_left, y_patient, COL_INK)
  grid.text("Eligible patient\nat clean t0",
            x = u(x_left), y = u(y_patient - 0.055),
            gp = gp_text(6.8), just = "center")
  draw_arrow(x_left + 0.025, y_patient, x_t0 - 0.015, y_patient, lwd = 0.8)
  grid.circle(x = u(x_t0), y = u(y_patient), r = u(0.014),
              gp = gpar(fill = COL_HEADER, col = COL_INK, lwd = 0.9))
  grid.text("t0", x = u(x_t0), y = u(y_patient), gp = gp_text(6.8, "bold"))
  draw_arrow(x_t0 + 0.015, y_patient + 0.002, x_mid - 0.018, y_early, col = COL_EARLY, lwd = 1.1)
  draw_arrow(x_t0 + 0.015, y_patient - 0.002, x_mid - 0.018, y_defer, col = COL_DEFERRED, lwd = 1.1)

  grid.roundrect(x = u(x_mid + 0.055), y = u(y_early), width = u(0.135), height = u(0.055),
                 r = unit(0.012, "npc"), gp = gpar(fill = "#E9F8FC", col = COL_EARLY, lwd = 1.0))
  grid.text("Early-strategy clone\nCRRT within 24 h",
            x = u(x_mid + 0.055), y = u(y_early), gp = gp_text(6.8, "bold", COL_INK))
  grid.roundrect(x = u(x_mid + 0.055), y = u(y_defer), width = u(0.135), height = u(0.055),
                 r = unit(0.012, "npc"), gp = gpar(fill = "#FDEDEA", col = COL_DEFERRED, lwd = 1.0))
  grid.text("Deferred-strategy clone\nNo CRRT within 24 h",
            x = u(x_mid + 0.055), y = u(y_defer), gp = gp_text(6.8, "bold", COL_INK))

  x0 <- 0.065
  x1 <- 0.307
  y_axis <- 0.342
  draw_time_axis(x0, x1, y_axis, TRUE)
  grid.text("Design premise: t0 = first physiologic KDIGO3\n(creatinine and/or urine-output criteria only)",
            x = u(0.18), y = u(0.255), gp = gp_text(6.9), just = "center")
}

draw_protocol_example <- function(y_base, row_label, scenario = c("within", "after", "none"),
                                  show_death = FALSE) {
  scenario <- match.arg(scenario)
  x0 <- 0.410
  x1 <- 0.625
  x24 <- time_x(x0, x1, 24)
  x_within <- (x0 + x24) / 2
  x_after <- x24 + (x1 - x24) * 0.34
  y_blue <- y_base + 0.023
  y_orange <- y_base - 0.023

  grid.text(row_label, x = u(0.369), y = u(y_base + 0.058),
            just = c("left", "center"), gp = gp_text(6.45, "bold"))
  grid.text("Early", x = u(0.366), y = u(y_blue), just = c("left", "center"),
            gp = gp_text(6.15, "bold", COL_EARLY))
  grid.text("Deferred", x = u(0.366), y = u(y_orange), just = c("left", "center"),
            gp = gp_text(6.15, "bold", COL_DEFERRED))

  # 24 h vertical guide shared by the row.
  grid.lines(x = u(c(x24, x24)), y = u(c(y_base - 0.052, y_base + 0.050)),
             gp = gpar(col = COL_GREY, lwd = 0.7, lty = "dashed"))

  if (scenario == "within") {
    draw_line(x0, y_blue, x1, y_blue, COL_EARLY, 2.0, 1)
    draw_square(x_within, y_blue)
    if (show_death) {
      draw_death(x1, y_blue, 0.006)
    }
    draw_line(x0, y_orange, x_within, y_orange, COL_DEFERRED, 2.0, 1)
    draw_censor_x(x_within, y_orange)
    draw_line(x_within, y_orange, x1, y_orange, COL_DEFERRED, 1.25, "dashed")
  } else if (scenario == "after") {
    draw_line(x0, y_blue, x24, y_blue, COL_EARLY, 2.0, 1)
    draw_censor_x(x24, y_blue)
    draw_line(x24, y_blue, x1, y_blue, COL_EARLY, 1.25, "dashed")
    draw_line(x0, y_orange, x1, y_orange, COL_DEFERRED, 2.0, 1)
    draw_square(x_after, y_orange)
  } else if (scenario == "none") {
    draw_line(x0, y_blue, x24, y_blue, COL_EARLY, 2.0, 1)
    draw_censor_x(x24, y_blue)
    draw_line(x24, y_blue, x1, y_blue, COL_EARLY, 1.25, "dashed")
    draw_line(x0, y_orange, x1, y_orange, COL_DEFERRED, 2.0, 1)
  }
}

draw_panel_b <- function() {
  panel(0.50, 0.545, 0.29, 0.70, "B", "PROTOCOL CENSOR",
        "Artificial censoring occurs when observed care\nviolates the assigned strategy.")

  x0 <- 0.410
  x1 <- 0.625
  x24 <- time_x(x0, x1, 24)
  # Shared visual anchor: the 24-h grace window.
  grid.rect(
    x = u((x0 + x24) / 2), y = u(0.522),
    width = u(x24 - x0), height = u(0.392),
    gp = gpar(fill = "#E9F8FC", col = NA)
  )
  grid.text("24-h grace window",
            x = u((x0 + x24) / 2), y = u(0.735),
            gp = gp_text(6.35, "bold", COL_EARLY))
  grid.text("t0", x = u(x0), y = u(0.781), gp = gp_text(6.5))
  grid.text("24 h", x = u(x24), y = u(0.781), gp = gp_text(6.5))
  grid.text("Day 28", x = u(x1), y = u(0.781), gp = gp_text(6.5))
  grid.lines(x = u(c(x0, x1)), y = u(c(0.762, 0.762)), gp = gpar(col = COL_INK, lwd = 0.65))
  draw_tick(x0, 0.762, 0.013)
  draw_tick(x24, 0.762, 0.013)
  draw_tick(x1, 0.762, 0.013)

  draw_protocol_example(
    y_base = 0.650,
    row_label = "1. CRRT initiated within the 24-h grace window",
    scenario = "within",
    show_death = TRUE
  )
  draw_protocol_example(
    y_base = 0.515,
    row_label = "2. CRRT initiated after the 24-h grace window",
    scenario = "after"
  )
  draw_protocol_example(
    y_base = 0.380,
    row_label = "3. No CRRT initiation",
    scenario = "none"
  )

  grid.text(
    "Rules: Deferred clone censored at the first CRRT within the grace window.\nEarly clone censored at 24 h if CRRT has not been initiated.\nCRRT initiated after 24 h is compatible with the deferred strategy.",
    x = u(0.500), y = u(0.295), just = "center",
    gp = gpar(fontfamily = BASE_FAMILY, fontsize = 5.45, col = COL_INK, lineheight = 0.92)
  )
  grid.text("Death is an outcome event, not protocol censoring.",
            x = u(0.505), y = u(0.238), gp = gp_text(6.55, "italic"), just = "center")
  draw_death(0.392, 0.238, 0.006)
}

draw_weight_row <- function(y, label, color, weight_label, width_scale, note) {
  x0 <- 0.755
  x1 <- 0.925
  x_end <- x0 + (x1 - x0) * width_scale
  draw_line(x0, y, x_end, y, color, 4.0, 1)
  grid.circle(x = u(x0), y = u(y), r = u(0.006), gp = gpar(fill = color, col = color))
  grid.text(label, x = u(0.735), y = u(y), just = c("right", "center"), gp = gp_text(6.7, "bold", color))
  grid.text(weight_label, x = u(0.950), y = u(y), just = c("right", "center"), gp = gp_text(8.6, "bold", COL_INK))
  grid.text(note, x = u(0.755), y = u(y - 0.030), just = c("left", "center"), gp = gp_text(6.0))
}

draw_panel_c <- function() {
  panel(0.82, 0.545, 0.29, 0.70, "C", "WEIGHT",
        "Remaining uncensored clone-time receives\nstabilized cumulative IPCW.")

  grid.text("Fewer comparable clones remain -> larger weight",
            x = u(0.82), y = u(0.730), gp = gp_text(7.3, "bold"), just = "center")
  draw_weight_row(0.650, "Clone-time 1", COL_EARLY, "W1", 0.43, "highest censoring risk")
  draw_weight_row(0.555, "Clone-time 2", COL_DEFERRED, "W2", 0.62, "intermediate censoring risk")
  draw_weight_row(0.460, "Clone-time 3", COL_INK, "W3", 0.84, "lowest censoring risk")
  grid.text("W1 > W2 > W3",
            x = u(0.82), y = u(0.365), gp = gp_text(12, "bold"), just = "center")
  grid.roundrect(x = u(0.82), y = u(0.275), width = u(0.250), height = u(0.088),
                 r = unit(0.012, "npc"),
                 gp = gpar(fill = COL_LIGHT, col = COL_GREY, lwd = 0.8))
  grid.text("Stabilized cumulative IPCW up-weights\nremaining uncensored clone-time to correct\nfor informative artificial censoring.",
            x = u(0.82), y = u(0.275),
            gp = gpar(fontfamily = BASE_FAMILY, fontsize = 6.7, col = COL_INK, lineheight = 0.95),
            just = "center")
}

draw_legend <- function() {
  y <- 0.103
  grid.roundrect(x = u(0.50), y = u(y), width = u(0.92), height = u(0.116),
                 r = unit(0.014, "npc"), gp = gpar(fill = "white", col = COL_GREY, lwd = 0.7))
  grid.text("Legend", x = u(0.060), y = u(y + 0.040), just = c("left", "center"),
            gp = gp_text(7.4, "bold"))

  xs <- c(0.135, 0.310, 0.565)
  draw_square(xs[1], y + 0.040, 0.010)
  grid.text("CRRT initiation", x = u(xs[1] + 0.020), y = u(y + 0.040),
            just = c("left", "center"), gp = gp_text(6.7))
  draw_censor_x(xs[2], y + 0.040, 0.013, 1.4)
  grid.text("protocol (artificial) censoring", x = u(xs[2] + 0.020), y = u(y + 0.040),
            just = c("left", "center"), gp = gp_text(6.7))
  draw_death(xs[3], y + 0.040, 0.006)
  grid.text("death = outcome event, not censoring", x = u(xs[3] + 0.020), y = u(y + 0.040),
            just = c("left", "center"), gp = gp_text(6.7))

  xs2 <- c(0.135, 0.360, 0.585, 0.770)
  draw_line(xs2[1] - 0.040, y - 0.010, xs2[1] - 0.010, y - 0.010, COL_EARLY, 2.2, 1)
  grid.text("Early-strategy clone", x = u(xs2[1]), y = u(y - 0.010), just = c("left", "center"),
            gp = gp_text(6.8, "bold", COL_EARLY))
  draw_line(xs2[2] - 0.040, y - 0.010, xs2[2] - 0.010, y - 0.010, COL_DEFERRED, 2.2, 1)
  grid.text("Deferred-strategy clone", x = u(xs2[2]), y = u(y - 0.010), just = c("left", "center"),
            gp = gp_text(6.8, "bold", COL_DEFERRED))
  draw_line(xs2[3] - 0.040, y - 0.010, xs2[3] - 0.010, y - 0.010, COL_INK, 1.8, 1)
  grid.text("solid line = on-strategy follow-up", x = u(xs2[3]), y = u(y - 0.010),
            just = c("left", "center"), gp = gp_text(6.7))
  draw_line(xs2[4] - 0.040, y - 0.010, xs2[4] - 0.010, y - 0.010, COL_INK, 1.4, "dashed")
  grid.text("dashed line = after artificial censoring", x = u(xs2[4]), y = u(y - 0.010),
            just = c("left", "center"), gp = gp_text(6.7))
}

draw_s7 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("Supplementary Figure S7. Clone-censor-weighting and IPCW schematic",
            x = u(0.035), y = u(0.960), just = c("left", "center"),
            gp = gp_text(12, "bold"))
  grid.text("Clean t0 = first physiologic KDIGO stage 3 (creatinine and/or urine-output criteria only); early CRRT = initiation within 0-24 h; deferred = no CRRT within 24 h.",
            x = u(0.035), y = u(0.925), just = c("left", "center"),
            gp = gp_text(7.5))
  draw_panel_a()
  draw_panel_b()
  draw_panel_c()
  draw_legend()
}

save_grid_pair <- function(draw_fun, stem, width = 11.0, height = 6.6, dpi = 600) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ragg::agg_png(png_path, width = width, height = height, units = "in",
                res = dpi, background = "white")
  draw_fun()
  invisible(dev.off())
  grDevices::pdf(pdf_path, width = width, height = height, family = VECTOR_FAMILY,
                 bg = "white", useDingbats = FALSE)
  draw_fun()
  invisible(dev.off())
  c(png = png_path, pdf = pdf_path)
}

paths <- save_grid_pair(draw_s7, stem)
message("Wrote: ", paths["png"])
message("Wrote: ", paths["pdf"])
