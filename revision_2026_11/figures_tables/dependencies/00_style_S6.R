# CH55 shared style for Figures S10-S14 (layout only; no analysis)
# Liberation Sans is metric-compatible with Arial.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(scales); library(patchwork); library(ragg)
})

FONT       <- "Liberation Sans"
BASE_SIZE  <- 8
COL_VERM   <- "#D55E00"   # HR > 1 side / harm
COL_BLUE   <- "#0072B2"   # HR < 1 side / benefit; after-IPCW
COL_INK    <- "grey15"
COL_MID    <- "grey45"
COL_LIGHT  <- "grey85"

# Input/output directories (override with env vars)
IN_DIR  <- Sys.getenv("CH55_IN_DIR",  normalizePath(file.path(SCRIPT_DIR, "..", "inputs"), mustWork = FALSE))
OUT_DIR <- Sys.getenv("CH55_OUT_DIR", normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
inp <- function(f) {
  p <- file.path(IN_DIR, f)
  if (!file.exists(p)) stop("Missing designated input: ", p, " (do not substitute)")
  p
}

theme_ch55 <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      text              = element_text(colour = COL_INK),
      axis.text         = element_text(colour = COL_INK, size = base_size - 0.5),
      axis.title        = element_text(colour = COL_INK, size = base_size + 0.5),
      axis.line         = element_line(colour = COL_INK, linewidth = 0.35),
      axis.ticks        = element_line(colour = COL_INK, linewidth = 0.35),
      legend.title      = element_text(size = base_size),
      legend.text       = element_text(size = base_size - 0.5),
      legend.key.size   = unit(9, "pt"),
      legend.background = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base_size),
      plot.title        = element_text(face = "bold", size = base_size + 2, hjust = 0,
                                       margin = margin(0, 0, 4, 0)),
      plot.title.position = "plot",
      plot.tag          = element_text(face = "bold", size = base_size + 2),
      plot.margin       = margin(4, 6, 4, 4)
    )
}

# Figure number only (no titles/subtitles/captions inside the image)
fig_title <- function(n) labs(title = paste0("Figure ", n))

save_fig <- function(plot, stem, width, height) {
  png <- file.path(OUT_DIR, paste0(stem, ".png"))
  pdf <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  ragg::agg_png(png, width = width, height = height, units = "in", res = 600, background = "white")
  print(plot); invisible(dev.off())
  grDevices::cairo_pdf(pdf, width = width, height = height, family = FONT, bg = "white")
  print(plot); invisible(dev.off())
  message("Wrote ", png, " and ", pdf)
  invisible(c(png = png, pdf = pdf))
}
