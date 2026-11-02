# Shared settings for Supplemental Figures S2-S4: palette, theme, export helpers
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(survival)
  library(patchwork); library(scales); library(ragg)
})
ROOT <- "/workspace/tte_bmj"
OUT  <- file.path(ROOT, "out"); dir.create(OUT, showWarnings = FALSE)
# Okabe-Ito colour-blind-safe palette (restrained): Early = vermillion, Deferred = blue
COL_EARLY    <- "#D55E00"
COL_DEFERRED <- "#0072B2"
COL_INK      <- "#222222"
COL_GRID     <- "grey90"
COL_GRACE    <- "grey94"
ARM_COLS <- c("Deferred CRRT" = COL_DEFERRED, "Early CRRT" = COL_EARLY)
BASE_FAMILY <- "sans"

theme_bmj <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(text = element_text(colour = COL_INK),
          axis.text = element_text(colour = COL_INK),
          axis.line = element_line(colour = COL_INK, linewidth = 0.4),
          axis.ticks = element_line(colour = COL_INK, linewidth = 0.4),
          plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
          plot.subtitle = element_text(size = base_size - 0.5, hjust = 0, lineheight = 1.1),
          plot.caption = element_text(size = base_size - 2, hjust = 0, colour = "grey30"),
          legend.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

save_fig <- function(plot, stem, width, height, dpi = 300) {
  png <- file.path(OUT, paste0(stem, ".png")); pdf <- file.path(OUT, paste0(stem, ".pdf"))
  ggsave(png, plot, device = ragg::agg_png, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  ggsave(pdf, plot, device = grDevices::cairo_pdf, width = width, height = height, units = "in", bg = "white")
  message("saved: ", png, " | ", pdf); invisible(c(png, pdf))
}
fmt3 <- function(x) formatC(x, format = "f", digits = 3)
