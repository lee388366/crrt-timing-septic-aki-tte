#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rsvg)
})

out_dir <- "analysis/strictK3_primary_v2/figures_true_dag_candidate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nodes <- data.frame(
  id = c("T0", "L0", "Lk", "Ak", "Lkp1", "Y", "Ck", "O", "U"),
  label = c(
    "Time zero / eligibility\nphysiologic KDIGO stage 3\nseptic shock; no prior RRT",
    "Baseline covariates\nL0",
    "Time-varying clinical\nstatus Lk",
    "CRRT initiation\nAk",
    "Subsequent clinical\nstatus L(k+1)",
    "Hospital-recorded\n28-day death Y",
    "Protocol censoring\nCk",
    "Observed analysis\ncontribution",
    "Unmeasured clinical\njudgement U"
  ),
  x = c(230, 540, 920, 920, 1340, 1780, 1260, 1780, 920),
  y = c(255, 420, 330, 610, 330, 420, 760, 760, 960),
  w = c(320, 280, 300, 280, 330, 320, 300, 320, 340),
  h = c(130, 110, 110, 110, 110, 110, 110, 110, 110),
  type = c("anchor", "measured", "measured", "treatment",
           "measured", "outcome", "censoring", "observed", "unmeasured"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c(
    "L0", "L0", "L0", "L0",
    "Lk", "Lk", "Lk", "Lk",
    "Ak", "Ak", "Ak",
    "Lkp1", "Ck", "Y",
    "U", "U", "U"
  ),
  to = c(
    "Lk", "Ak", "Ck", "Y",
    "Ak", "Ck", "Lkp1", "Y",
    "Lkp1", "Y", "Ck",
    "Y", "O", "O",
    "Ak", "Ck", "Y"
  ),
  kind = c(rep("measured", 14), rep("unmeasured", 3)),
  bend = c(
    0, 0, 0.10, -0.24,
    0, 0, 0, -0.14,
    0, -0.10, 0,
    0, 0, 0,
    0, 0.10, -0.20
  ),
  stringsAsFactors = FALSE
)

assert_acyclic <- function(nodes, edges) {
  causal_nodes <- nodes$id[nodes$id != "T0"]
  causal_edges <- edges[edges$from %in% causal_nodes & edges$to %in% causal_nodes, ]
  children <- split(causal_edges$to, causal_edges$from)
  state <- setNames(rep(0L, length(causal_nodes)), causal_nodes)

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

  for (v in causal_nodes) visit(v)
  invisible(TRUE)
}

assert_acyclic(nodes, edges)

node <- function(id) nodes[match(id, nodes$id), ]

edge_point <- function(from_id, to_id, side = c("start", "end")) {
  side <- match.arg(side)
  a <- node(from_id)
  b <- node(to_id)
  dx <- b$x - a$x
  dy <- b$y - a$y
  if (dx == 0 && dy == 0) stop("Self edge not supported", call. = FALSE)
  scale_x <- if (dx == 0) Inf else (a$w / 2) / abs(dx)
  scale_y <- if (dy == 0) Inf else (a$h / 2) / abs(dy)
  s <- min(scale_x, scale_y) * 1.04
  if (side == "start") {
    return(c(x = a$x + dx * s, y = a$y + dy * s))
  }
  scale_x <- if (dx == 0) Inf else (b$w / 2) / abs(dx)
  scale_y <- if (dy == 0) Inf else (b$h / 2) / abs(dy)
  s <- min(scale_x, scale_y) * 1.08
  c(x = b$x - dx * s, y = b$y - dy * s)
}

escape_xml <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

text_svg <- function(x, y, label, size = 34, weight = 400, fill = "#102A43", anchor = "middle") {
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  line_h <- size * 1.18
  start_y <- y - (length(lines) - 1) * line_h / 2
  spans <- vapply(seq_along(lines), function(i) {
    sprintf("<tspan x=\"%.1f\" y=\"%.1f\">%s</tspan>",
            x, start_y + (i - 1) * line_h, escape_xml(lines[i]))
  }, character(1))
  sprintf(
    "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"%s\" font-size=\"%s\" font-weight=\"%s\" fill=\"%s\">%s</text>",
    x, y, anchor, size, weight, fill, paste(spans, collapse = "")
  )
}

node_style <- function(type) {
  switch(
    type,
    anchor = c(fill = "#F7FAFC", stroke = "#718096", dash = "8 8", text = "#334E68"),
    measured = c(fill = "#F8FAFC", stroke = "#34495E", dash = "", text = "#102A43"),
    treatment = c(fill = "#E8F4FA", stroke = "#2C7FB8", dash = "", text = "#1B75BB"),
    censoring = c(fill = "#FFF7E6", stroke = "#B7791F", dash = "", text = "#A16207"),
    outcome = c(fill = "#F7F2FF", stroke = "#6B46C1", dash = "", text = "#6B46C1"),
    observed = c(fill = "#F8FAF4", stroke = "#4A5568", dash = "", text = "#3F4B5B"),
    unmeasured = c(fill = "#FFFFFF", stroke = "#4A5568", dash = "10 8", text = "#3F4B5B")
  )
}

node_svg <- function(row) {
  st <- node_style(row$type)
  dash <- if (nzchar(st[["dash"]])) sprintf(" stroke-dasharray=\"%s\"", st[["dash"]]) else ""
  paste0(
    sprintf(
      "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" rx=\"18\" fill=\"%s\" stroke=\"%s\" stroke-width=\"4\"%s/>",
      row$x - row$w / 2, row$y - row$h / 2, row$w, row$h, st[["fill"]], st[["stroke"]], dash
    ),
    text_svg(row$x, row$y + 6, row$label, size = ifelse(row$id == "T0", 24, 30),
             weight = 400, fill = st[["text"]])
  )
}

edge_svg <- function(row) {
  start <- edge_point(row$from, row$to, "start")
  end <- edge_point(row$from, row$to, "end")
  stroke <- if (row$kind == "unmeasured") "#243B53" else "#102A43"
  dash <- if (row$kind == "unmeasured") " stroke-dasharray=\"14 12\"" else ""
  width <- if (row$kind == "unmeasured") 4.2 else 4.0
  if (abs(row$bend) < 1e-9) {
    return(sprintf(
      "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"%.1f\" fill=\"none\" marker-end=\"url(#arrow)\"%s/>",
      start[["x"]], start[["y"]], end[["x"]], end[["y"]], stroke, width, dash
    ))
  }
  mx <- (start[["x"]] + end[["x"]]) / 2
  my <- (start[["y"]] + end[["y"]]) / 2
  dx <- end[["x"]] - start[["x"]]
  dy <- end[["y"]] - start[["y"]]
  norm <- sqrt(dx^2 + dy^2)
  nx <- -dy / norm
  ny <- dx / norm
  cx <- mx + nx * row$bend * norm
  cy <- my + ny * row$bend * norm
  sprintf(
    "<path d=\"M %.1f %.1f Q %.1f %.1f %.1f %.1f\" stroke=\"%s\" stroke-width=\"%.1f\" fill=\"none\" marker-end=\"url(#arrow)\"%s/>",
    start[["x"]], start[["y"]], cx, cy, end[["x"]], end[["y"]], stroke, width, dash
  )
}

svg_width <- 2200
svg_height <- 1250

svg <- paste0(
  "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", svg_width, "\" height=\"", svg_height,
  "\" viewBox=\"0 0 ", svg_width, " ", svg_height, "\">",
  "<defs>",
  "<style><![CDATA[",
  "text{font-family:Arial, Helvetica, sans-serif;dominant-baseline:middle}",
  ".small{font-size:24px;fill:#52606D}",
  ".section{font-size:36px;font-weight:700;fill:#0B2545}",
  "]]></style>",
  "<marker id=\"arrow\" markerWidth=\"18\" markerHeight=\"18\" refX=\"17\" refY=\"9\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">",
  "<path d=\"M 0 0 L 18 9 L 0 18 z\" fill=\"#102A43\"/>",
  "</marker>",
  "</defs>",
  "<rect x=\"0\" y=\"0\" width=\"2200\" height=\"1250\" fill=\"#FFFFFF\"/>",
  text_svg(1100, 78, "Supplementary Figure S1. Directed acyclic graph", size = 56, weight = 700, fill = "#0B2545"),
  text_svg(1100, 138, "Causal structure for CRRT timing target-trial emulation with clone-censor-weighting and IPCW", size = 28, fill = "#334E68"),
  "<line x1=\"390\" y1=\"210\" x2=\"390\" y2=\"1060\" stroke=\"#A0AEC0\" stroke-width=\"3\" stroke-dasharray=\"9 12\"/>",
  "<line x1=\"1110\" y1=\"210\" x2=\"1110\" y2=\"1060\" stroke=\"#CBD5E0\" stroke-width=\"3\" stroke-dasharray=\"9 12\"/>",
  "<line x1=\"1585\" y1=\"210\" x2=\"1585\" y2=\"1060\" stroke=\"#CBD5E0\" stroke-width=\"3\" stroke-dasharray=\"9 12\"/>",
  "<text x=\"220\" y=\"200\" class=\"section\" text-anchor=\"middle\">Design anchor</text>",
  "<text x=\"535\" y=\"200\" class=\"section\" text-anchor=\"middle\">Baseline at t0</text>",
  "<text x=\"930\" y=\"200\" class=\"section\" text-anchor=\"middle\">Hour k</text>",
  "<text x=\"1350\" y=\"200\" class=\"section\" text-anchor=\"middle\">Hour k + 1</text>",
  "<text x=\"1780\" y=\"200\" class=\"section\" text-anchor=\"middle\">Outcome / observation</text>",
  paste(vapply(seq_len(nrow(edges)), function(i) edge_svg(edges[i, ]), character(1)), collapse = "\n"),
  paste(vapply(seq_len(nrow(nodes)), function(i) node_svg(nodes[i, ]), character(1)), collapse = "\n"),
  "<line x1=\"585\" y1=\"1115\" x2=\"660\" y2=\"1115\" stroke=\"#102A43\" stroke-width=\"4\" marker-end=\"url(#arrow)\"/>",
  "<text x=\"690\" y=\"1115\" class=\"small\" text-anchor=\"start\">Measured or design pathway</text>",
  "<line x1=\"1045\" y1=\"1115\" x2=\"1120\" y2=\"1115\" stroke=\"#243B53\" stroke-width=\"4\" stroke-dasharray=\"14 12\" marker-end=\"url(#arrow)\"/>",
  "<text x=\"1150\" y=\"1115\" class=\"small\" text-anchor=\"start\">Residual unmeasured pathway</text>",
  text_svg(1100, 1180, "T0 is a design/selection anchor, not a treatment variable. Ck affects observation and analysis contribution; it does not directly cause the underlying death event Y.", size = 24, fill = "#52606D"),
  "</svg>"
)

svg_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_redrawn_v2.svg")
png_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_redrawn_v2.png")
pdf_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_redrawn_v2.pdf")
csv_nodes <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_redrawn_v2_nodes.csv")
csv_edges <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_redrawn_v2_edges.csv")
dagitty_file <- file.path(out_dir, "Supplement_Figure_S1_true_DAG_redrawn_v2.dagitty")

writeLines(svg, svg_file)
rsvg_png(svg_file, png_file, width = svg_width * 2)
rsvg_pdf(svg_file, pdf_file)
write.csv(nodes, csv_nodes, row.names = FALSE)
write.csv(edges[, c("from", "to", "kind")], csv_edges, row.names = FALSE)

dagitty_lines <- c(
  "dag {",
  "  # T0 is a design/selection anchor and is intentionally excluded from the causal edge list.",
  paste0("  ", nodes$id[nodes$id != "T0"]),
  "  # Causal edges",
  paste0("  ", edges$from, " -> ", edges$to),
  "}"
)
writeLines(dagitty_lines, dagitty_file)

message("Created redrawn true DAG files:")
message("  ", svg_file)
message("  ", png_file)
message("  ", pdf_file)
message("  ", csv_nodes)
message("  ", csv_edges)
message("  ", dagitty_file)
