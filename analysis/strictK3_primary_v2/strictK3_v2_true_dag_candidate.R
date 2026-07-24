#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- "analysis/strictK3_primary_v2/figures_true_dag_candidate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nodes <- data.frame(
  id = c("L0", "Lk", "Ak", "Ck", "Lkp1", "Y", "O", "U"),
  label = c(
    "Baseline covariates\nL0",
    "Clinical status\nLk",
    "CRRT initiation\nAk",
    "Protocol censoring\nCk",
    "Subsequent clinical\nstatus L(k+1)",
    "Hospital-recorded\n28-day death Y",
    "Observed analysis\ncontribution",
    "Unmeasured clinical\njudgement U"
  ),
  x = c(0.5, 3.2, 3.2, 5.7, 6.8, 9.4, 9.4, 3.2),
  y = c(3.2, 3.2, 2.0, 1.0, 3.2, 3.2, 1.6, 0.0),
  type = c("measured", "measured", "treatment", "censoring",
           "measured", "outcome", "observed", "unmeasured"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c(
    "L0", "L0", "L0", "L0",
    "Lk", "Lk", "Lk", "Lk",
    "Ak", "Ak", "Ak",
    "Ck", "Y",
    "U", "U", "U"
  ),
  to = c(
    "Lk", "Ak", "Ck", "Y",
    "Ak", "Ck", "Lkp1", "Y",
    "Ck", "Lkp1", "Y",
    "O", "O",
    "Ak", "Ck", "Y"
  ),
  kind = c(
    rep("measured", 13),
    rep("unmeasured", 3)
  ),
  curve = c(
    0, -0.15, -0.25, 0.28,
    0, -0.12, 0, 0.18,
    0, -0.08, 0.15,
    0, 0,
    0, 0.12, 0.32
  ),
  stringsAsFactors = FALSE
)

assert_acyclic <- function(nodes, edges) {
  children <- split(edges$to, edges$from)
  state <- setNames(rep(0L, nrow(nodes)), nodes$id)

  visit <- function(v, stack = character()) {
    if (state[[v]] == 1L) {
      stop("Cycle detected: ", paste(c(stack, v), collapse = " -> "), call. = FALSE)
    }
    if (state[[v]] == 2L) return(invisible(NULL))
    state[[v]] <<- 1L
    for (w in children[[v]]) {
      if (!is.null(w) && !is.na(w)) visit(w, c(stack, v))
    }
    state[[v]] <<- 2L
    invisible(NULL)
  }

  for (v in nodes$id) visit(v)
  invisible(TRUE)
}

assert_acyclic(nodes, edges)

edge_plot <- edges
edge_plot$x <- nodes$x[match(edge_plot$from, nodes$id)]
edge_plot$y <- nodes$y[match(edge_plot$from, nodes$id)]
edge_plot$xend <- nodes$x[match(edge_plot$to, nodes$id)]
edge_plot$yend <- nodes$y[match(edge_plot$to, nodes$id)]

node_fill <- c(
  measured = "#F8FAFC",
  treatment = "#E8F4FA",
  censoring = "#FFF7E6",
  outcome = "#F7F2FF",
  observed = "#F4F7F0",
  unmeasured = "#FFFFFF"
)

node_border <- c(
  measured = "#34495E",
  treatment = "#2C7FB8",
  censoring = "#B7791F",
  outcome = "#6B46C1",
  observed = "#4A5568",
  unmeasured = "#4A5568"
)

p <- ggplot() +
  annotate("text", x = 0.5, y = 4.05, label = "Baseline", fontface = "bold",
           size = 5.2, color = "#0B2545") +
  annotate("text", x = 4.1, y = 4.05, label = "During follow-up, hour k", fontface = "bold",
           size = 5.2, color = "#0B2545") +
  annotate("text", x = 6.8, y = 4.05, label = "Subsequent hour, k + 1", fontface = "bold",
           size = 5.2, color = "#0B2545") +
  annotate("text", x = 9.4, y = 4.05, label = "Outcome / observation", fontface = "bold",
           size = 5.2, color = "#0B2545") +
  geom_vline(xintercept = c(1.75, 4.7, 8.05), linetype = "dashed",
             color = "#A0AEC0", linewidth = 0.6) +
  geom_segment(
    data = edge_plot,
    aes(x = x, y = y, xend = xend, yend = yend,
        linetype = kind, linewidth = kind),
    color = "#102A43",
    arrow = arrow(length = unit(0.20, "inches"), type = "closed"),
    lineend = "round"
  ) +
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = label, fill = type, color = type,
        linetype = ifelse(type == "unmeasured", "unmeasured", "measured")),
    linewidth = 0.75,
    label.r = unit(0.18, "lines"),
    size = 4.4,
    fontface = "plain",
    lineheight = 0.95,
    label.padding = unit(0.35, "lines"),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = node_fill, guide = "none") +
  scale_color_manual(values = node_border, guide = "none") +
  scale_linetype_manual(
    values = c(measured = "solid", unmeasured = "dashed"),
    breaks = c("measured", "unmeasured"),
    labels = c("Measured/design pathway", "Residual unmeasured pathway"),
    name = NULL
  ) +
  scale_linewidth_manual(values = c(measured = 0.95, unmeasured = 0.9), guide = "none") +
  coord_cartesian(xlim = c(-0.5, 10.3), ylim = c(-0.55, 4.35), clip = "off") +
  labs(
    title = "Directed acyclic graph for CRRT timing target-trial emulation",
    subtitle = "Treatment, time-varying clinical status, protocol censoring, and observed outcome contribution",
    caption = "Ck affects observation/analysis contribution, not the underlying mortality event Y. Dashed arrows denote residual unmeasured clinical judgement."
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 22, color = "#0B2545", hjust = 0.5),
    plot.subtitle = element_text(size = 12.5, color = "#334E68", hjust = 0.5,
                                 margin = margin(b = 18)),
    plot.caption = element_text(size = 10.5, color = "#52606D", hjust = 0.5,
                                margin = margin(t = 14)),
    legend.position = "bottom",
    legend.text = element_text(size = 11, color = "#243B53"),
    plot.margin = margin(18, 24, 18, 24)
  )

png_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_candidate.png")
pdf_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_candidate.pdf")
csv_nodes <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_nodes.csv")
csv_edges <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_edges.csv")
dagitty_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_candidate.dagitty")

ggsave(png_file, p, width = 12.5, height = 6.6, dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 12.5, height = 6.6, device = "pdf", bg = "white")
write.csv(nodes, csv_nodes, row.names = FALSE)
write.csv(edges[, c("from", "to", "kind")], csv_edges, row.names = FALSE)

dagitty_lines <- c(
  "dag {",
  "  # Nodes",
  paste0("  ", nodes$id),
  "  # Edges",
  paste0("  ", edges$from, " -> ", edges$to),
  "}"
)
writeLines(dagitty_lines, dagitty_file)

message("Created true DAG candidate files:")
message("  ", png_file)
message("  ", pdf_file)
message("  ", csv_nodes)
message("  ", csv_edges)
message("  ", dagitty_file)
