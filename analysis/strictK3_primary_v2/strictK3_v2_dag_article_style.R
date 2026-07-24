#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rsvg)
})

out_dir <- "analysis/strictK3_primary_v2/figures_true_dag_candidate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nodes <- data.frame(
  id = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "U1", "U2", "X", "Y"),
  label = c(
    "Age, sex, chronic\ncomorbidity",
    "Sepsis and shock\nseverity",
    "AKI severity and\nKDIGO3 pathway",
    "Organ support and\nhemodynamic burden",
    "Biochemical derangement\nand rescue risk",
    "Care context, calendar\nperiod, ICU practice",
    "Clinician judgement,\ngoals of care",
    "Unrecorded volume status,\ntrajectory, reversibility",
    "CRRT timing strategy\nwithin 24 h",
    "Hospital-recorded\n28-day death"
  ),
  x = c(260, 520, 840, 1160, 1160, 840, 520, 840, 1460, 1780),
  y = c(520, 220, 180, 220, 520, 840, 840, 1030, 520, 520),
  type = c(rep("Z", 6), rep("U", 2), "X", "Y"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c(
    "Z1", "Z1", "Z1", "Z1", "Z1", "Z1", "Z1",
    "Z2", "Z2", "Z2", "Z2", "Z2",
    "Z3", "Z3", "Z3", "Z3",
    "Z4", "Z4", "Z4",
    "Z5", "Z5",
    "Z6", "Z6",
    "U1", "U1",
    "U2", "U2", "U2",
    "X"
  ),
  to = c(
    "Z2", "Z3", "Z4", "Z5", "Z6", "X", "Y",
    "Z3", "Z4", "Z5", "X", "Y",
    "Z4", "Z5", "X", "Y",
    "Z5", "X", "Y",
    "X", "Y",
    "X", "Y",
    "X", "Y",
    "Z5", "X", "Y",
    "Y"
  ),
  kind = c(rep("measured", 23), rep("unmeasured", 5), "effect"),
  curve = c(
    0.00, -0.15, -0.24, 0.00, 0.10, 0.05, -0.18,
    0.00, -0.10, 0.00, 0.02, -0.10,
    0.00, 0.00, 0.00, -0.06,
    0.00, 0.00, -0.04,
    0.00, 0.02,
    0.00, 0.10,
    0.00, 0.18,
    0.00, 0.00, 0.14,
    0.00
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

node <- function(id) nodes[match(id, nodes$id), ]

escape_xml <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

text_svg <- function(x, y, label, size = 24, weight = 400, fill = "#111827", anchor = "middle") {
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  line_h <- size * 1.15
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

node_color <- function(type) {
  switch(
    type,
    Z = c(fill = "#00A77E", stroke = "#008A68", label = "#111827", opacity = "1"),
    U = c(fill = "#5DD8C8", stroke = "#34BBA9", label = "#111827", opacity = "0.58"),
    X = c(fill = "#E6A000", stroke = "#C78300", label = "#111827", opacity = "1"),
    Y = c(fill = "#5BB9E8", stroke = "#3B9DCE", label = "#111827", opacity = "1")
  )
}

edge_point <- function(from_id, to_id, side = c("start", "end"), radius = 54) {
  side <- match.arg(side)
  a <- node(from_id)
  b <- node(to_id)
  dx <- b$x - a$x
  dy <- b$y - a$y
  d <- sqrt(dx^2 + dy^2)
  if (d == 0) stop("Self edge not supported", call. = FALSE)
  if (side == "start") return(c(x = a$x + dx / d * radius, y = a$y + dy / d * radius))
  c(x = b$x - dx / d * radius, y = b$y - dy / d * radius)
}

edge_svg <- function(row) {
  start <- edge_point(row$from, row$to, "start")
  end <- edge_point(row$from, row$to, "end")
  stroke <- if (row$kind == "unmeasured") "#4B5563" else "#111827"
  dash <- if (row$kind == "unmeasured") " stroke-dasharray=\"12 10\"" else ""
  width <- if (row$kind == "effect") 4.4 else 3.2
  if (abs(row$curve) < 1e-9) {
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
  cx <- mx + nx * row$curve * norm
  cy <- my + ny * row$curve * norm
  sprintf(
    "<path d=\"M %.1f %.1f Q %.1f %.1f %.1f %.1f\" stroke=\"%s\" stroke-width=\"%.1f\" fill=\"none\" marker-end=\"url(#arrow)\"%s/>",
    start[["x"]], start[["y"]], cx, cy, end[["x"]], end[["y"]], stroke, width, dash
  )
}

node_svg <- function(row) {
  col <- node_color(row$type)
  radius <- 58
  number_size <- if (row$type %in% c("X", "Y")) 42 else 34
  paste0(
    sprintf(
      "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"%s\" fill=\"%s\" stroke=\"%s\" stroke-width=\"3\" opacity=\"%s\"/>",
      row$x, row$y, radius, col[["fill"]], col[["stroke"]], col[["opacity"]]
    ),
    text_svg(row$x, row$y - 6, row$id, size = number_size, weight = 700, fill = "#FFFFFF"),
    sprintf(
      "<rect x=\"%.1f\" y=\"%.1f\" width=\"260\" height=\"58\" rx=\"7\" fill=\"#FFFFFF\" stroke=\"#6B7280\" stroke-width=\"1.6\"/>",
      row$x - 130, row$y + 42
    ),
    text_svg(row$x, row$y + 72, row$label, size = 18, fill = col[["label"]])
  )
}

svg_width <- 2040
svg_height <- 1220

svg <- paste0(
  "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", svg_width, "\" height=\"", svg_height,
  "\" viewBox=\"0 0 ", svg_width, " ", svg_height, "\">",
  "<defs>",
  "<style><![CDATA[text{font-family:Arial, Helvetica, sans-serif;dominant-baseline:middle}.title{font-size:44px;font-weight:700;fill:#111827}.subtitle{font-size:24px;fill:#4B5563}.note{font-size:22px;fill:#4B5563}]]></style>",
  "<marker id=\"arrow\" markerWidth=\"15\" markerHeight=\"15\" refX=\"14\" refY=\"7.5\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">",
  "<path d=\"M 0 0 L 15 7.5 L 0 15 z\" fill=\"#111827\"/>",
  "</marker>",
  "</defs>",
  "<rect width=\"2040\" height=\"1220\" fill=\"#FFFFFF\"/>",
  "<text x=\"1020\" y=\"70\" text-anchor=\"middle\" class=\"title\">DAG for early CRRT initiation and hospital-recorded 28-day death</text>",
  "<text x=\"1020\" y=\"118\" text-anchor=\"middle\" class=\"subtitle\">Article-style causal diagram: measured covariates Z, treatment strategy X, outcome Y, and unmeasured factors U</text>",
  paste(vapply(seq_len(nrow(edges)), function(i) edge_svg(edges[i, ]), character(1)), collapse = "\n"),
  paste(vapply(seq_len(nrow(nodes)), function(i) node_svg(nodes[i, ]), character(1)), collapse = "\n"),
  "<circle cx=\"560\" cy=\"1130\" r=\"20\" fill=\"#00A77E\"/><text x=\"590\" y=\"1130\" class=\"note\" text-anchor=\"start\">Measured covariates / design variables</text>",
  "<circle cx=\"1035\" cy=\"1130\" r=\"20\" fill=\"#5DD8C8\" opacity=\"0.58\"/><text x=\"1065\" y=\"1130\" class=\"note\" text-anchor=\"start\">Unmeasured factors</text>",
  "<circle cx=\"1335\" cy=\"1130\" r=\"20\" fill=\"#E6A000\"/><text x=\"1365\" y=\"1130\" class=\"note\" text-anchor=\"start\">Treatment strategy X</text>",
  "<circle cx=\"1605\" cy=\"1130\" r=\"20\" fill=\"#5BB9E8\"/><text x=\"1635\" y=\"1130\" class=\"note\" text-anchor=\"start\">Outcome Y</text>",
  "</svg>"
)

svg_file <- file.path(out_dir, "Supplement_Figure_S1_DAG_article_style_CRRT_v1.svg")
png_file <- file.path(out_dir, "Supplement_Figure_S1_DAG_article_style_CRRT_v1.png")
pdf_file <- file.path(out_dir, "Supplement_Figure_S1_DAG_article_style_CRRT_v1.pdf")
csv_nodes <- file.path(out_dir, "Supplement_Figure_S1_DAG_article_style_CRRT_v1_nodes.csv")
csv_edges <- file.path(out_dir, "Supplement_Figure_S1_DAG_article_style_CRRT_v1_edges.csv")
dagitty_file <- file.path(out_dir, "Supplement_Figure_S1_DAG_article_style_CRRT_v1.dagitty")

writeLines(svg, svg_file)
rsvg_png(svg_file, png_file, width = svg_width * 2)
rsvg_pdf(svg_file, pdf_file)
write.csv(nodes, csv_nodes, row.names = FALSE)
write.csv(edges[, c("from", "to", "kind")], csv_edges, row.names = FALSE)

dagitty_lines <- c(
  "dag {",
  paste0("  ", nodes$id),
  "  # Edges",
  paste0("  ", edges$from, " -> ", edges$to),
  "}"
)
writeLines(dagitty_lines, dagitty_file)

message("Created article-style CRRT DAG files:")
message("  ", svg_file)
message("  ", png_file)
message("  ", pdf_file)
message("  ", csv_nodes)
message("  ", csv_edges)
message("  ", dagitty_file)
