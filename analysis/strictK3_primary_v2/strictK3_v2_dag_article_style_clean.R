#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rsvg)
})

out_dir <- "analysis/strictK3_primary_v2/figures_true_dag_candidate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

basename <- "Supplement_Figure_S1_DAG_article_style_CRRT_v2_clean"

nodes <- data.frame(
  id = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "U1", "U2", "X", "Y"),
  label = c(
    "Age, sex, chronic\ncomorbidity",
    "Sepsis and shock\nseverity",
    "AKI severity and\nKDIGO stage 3 pathway",
    "Organ support and\nhemodynamic burden",
    "Biochemical derangement\nand rescue risk",
    "Care context, calendar\nperiod, ICU practice",
    "Clinician judgement,\ngoals of care",
    "Unrecorded volume status,\ntrajectory, reversibility",
    "CRRT timing strategy\nwithin 24 h",
    "Hospital-recorded\n28-day death"
  ),
  x = c(230, 560, 900, 560, 900, 230, 560, 900, 1420, 1880),
  y = c(290, 170, 290, 650, 650, 850, 1030, 1030, 520, 520),
  type = c(rep("Z", 6), rep("U", 2), "X", "Y"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c(
    # Key measured causes among baseline and time-varying covariates.
    "Z1", "Z1", "Z1",
    "Z2", "Z2", "Z2", "Z2",
    "Z3", "Z4",
    # Measured common causes of treatment choice and outcome.
    "Z1", "Z2", "Z3", "Z4", "Z5", "Z6",
    "Z1", "Z2", "Z3", "Z4", "Z5", "Z6",
    # Unmeasured common causes.
    "U1", "U1",
    "U2", "U2", "U2",
    # Target causal contrast.
    "X"
  ),
  to = c(
    "Z2", "Z3", "Z6",
    "Z3", "Z4", "Z5", "Z6",
    "Z5", "Z5",
    "X", "X", "X", "X", "X", "X",
    "Y", "Y", "Y", "Y", "Y", "Y",
    "X", "Y",
    "Z5", "X", "Y",
    "Y"
  ),
  kind = c(
    rep("covariate_path", 9),
    rep("measured_common_cause", 12),
    rep("unmeasured_common_cause", 5),
    "target_effect"
  ),
  curve = c(
    0.00, -0.12, 0.05,
    0.00, 0.08, -0.04, 0.20,
    0.00, 0.00,
    0.30, 0.16, 0.06, -0.12, -0.04, -0.30,
    0.34, 0.20, 0.10, -0.20, -0.10, -0.36,
    -0.16, -0.30,
    0.00, -0.10, -0.24,
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

text_svg <- function(x, y, label, size = 24, weight = 400,
                     fill = "#111827", anchor = "middle") {
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  line_h <- size * 1.14
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
    Z = c(fill = "#00A77E", stroke = "#008B6C", label = "#111827", opacity = "1"),
    U = c(fill = "#7FDED2", stroke = "#4FC7BA", label = "#111827", opacity = "0.52"),
    X = c(fill = "#EDA900", stroke = "#C98C00", label = "#111827", opacity = "1"),
    Y = c(fill = "#59BCEB", stroke = "#3197C7", label = "#111827", opacity = "1")
  )
}

edge_style <- function(kind) {
  switch(
    kind,
    covariate_path = c(stroke = "#718096", width = "2.0", dash = "", marker = "arrowMeasured"),
    measured_common_cause = c(stroke = "#111827", width = "2.4", dash = "", marker = "arrowMeasured"),
    unmeasured_common_cause = c(stroke = "#4B5563", width = "2.6", dash = " stroke-dasharray=\"13 10\"", marker = "arrowUnmeasured"),
    target_effect = c(stroke = "#111827", width = "4.8", dash = "", marker = "arrowEffect")
  )
}

edge_point <- function(from_id, to_id, side = c("start", "end"), radius = 60) {
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
  st <- edge_style(row$kind)
  if (abs(row$curve) < 1e-9) {
    return(sprintf(
      "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"%s\" fill=\"none\" marker-end=\"url(#%s)\"%s/>",
      start[["x"]], start[["y"]], end[["x"]], end[["y"]],
      st[["stroke"]], st[["width"]], st[["marker"]], st[["dash"]]
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
    "<path d=\"M %.1f %.1f Q %.1f %.1f %.1f %.1f\" stroke=\"%s\" stroke-width=\"%s\" fill=\"none\" marker-end=\"url(#%s)\"%s/>",
    start[["x"]], start[["y"]], cx, cy, end[["x"]], end[["y"]],
    st[["stroke"]], st[["width"]], st[["marker"]], st[["dash"]]
  )
}

node_svg <- function(row) {
  col <- node_color(row$type)
  radius <- 60
  number_size <- if (row$type %in% c("X", "Y")) 44 else 35
  paste0(
    sprintf(
      "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"%s\" fill=\"%s\" stroke=\"%s\" stroke-width=\"3\" opacity=\"%s\"/>",
      row$x, row$y, radius, col[["fill"]], col[["stroke"]], col[["opacity"]]
    ),
    text_svg(row$x, row$y - 5, row$id, size = number_size, weight = 700, fill = "#FFFFFF"),
    sprintf(
      "<rect x=\"%.1f\" y=\"%.1f\" width=\"285\" height=\"64\" rx=\"7\" fill=\"#FFFFFF\" stroke=\"#6B7280\" stroke-width=\"1.4\"/>",
      row$x - 142.5, row$y + 48
    ),
    text_svg(row$x, row$y + 80, row$label, size = 18, fill = col[["label"]])
  )
}

svg_width <- 2160
svg_height <- 1340

marker <- function(id, fill) {
  paste0(
    "<marker id=\"", id, "\" markerWidth=\"15\" markerHeight=\"15\" refX=\"14\" refY=\"7.5\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">",
    "<path d=\"M 0 0 L 15 7.5 L 0 15 z\" fill=\"", fill, "\"/>",
    "</marker>"
  )
}

legend_item <- function(x, y, fill, label, opacity = "1") {
  paste0(
    sprintf("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"20\" fill=\"%s\" opacity=\"%s\"/>", x, y, fill, opacity),
    text_svg(x + 34, y, label, size = 22, fill = "#4B5563", anchor = "start")
  )
}

svg <- paste0(
  "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", svg_width, "\" height=\"", svg_height,
  "\" viewBox=\"0 0 ", svg_width, " ", svg_height, "\">",
  "<defs>",
  "<style><![CDATA[text{font-family:Arial, Helvetica, sans-serif;dominant-baseline:middle}.note{font-size:22px;fill:#4B5563}]]></style>",
  marker("arrowMeasured", "#111827"),
  marker("arrowUnmeasured", "#4B5563"),
  marker("arrowEffect", "#111827"),
  "</defs>",
  "<rect width=\"2160\" height=\"1340\" fill=\"#FFFFFF\"/>",
  paste(vapply(seq_len(nrow(edges)), function(i) edge_svg(edges[i, ]), character(1)), collapse = "\n"),
  paste(vapply(seq_len(nrow(nodes)), function(i) node_svg(nodes[i, ]), character(1)), collapse = "\n"),
  legend_item(520, 1245, "#00A77E", "Measured covariates Z"),
  legend_item(850, 1245, "#7FDED2", "Unmeasured factors U", opacity = "0.52"),
  legend_item(1160, 1245, "#EDA900", "Treatment strategy X"),
  legend_item(1485, 1245, "#59BCEB", "Outcome Y"),
  "</svg>"
)

svg_file <- file.path(out_dir, paste0(basename, ".svg"))
png_file <- file.path(out_dir, paste0(basename, ".png"))
pdf_file <- file.path(out_dir, paste0(basename, ".pdf"))
csv_nodes <- file.path(out_dir, paste0(basename, "_nodes.csv"))
csv_edges <- file.path(out_dir, paste0(basename, "_edges.csv"))
dagitty_file <- file.path(out_dir, paste0(basename, ".dagitty"))

writeLines(svg, svg_file)
rsvg_png(svg_file, png_file, width = svg_width * 2)
rsvg_pdf(svg_file, pdf_file)
write.csv(nodes, csv_nodes, row.names = FALSE)
write.csv(edges[, c("from", "to", "kind")], csv_edges, row.names = FALSE)

dagitty_lines <- c(
  "dag {",
  "  X [exposure]",
  "  Y [outcome]",
  paste0("  ", setdiff(nodes$id, c("X", "Y"))),
  "  # Edges",
  paste0("  ", edges$from, " -> ", edges$to),
  "}"
)
writeLines(dagitty_lines, dagitty_file)

message("Created clean article-style CRRT DAG files:")
message("  ", svg_file)
message("  ", png_file)
message("  ", pdf_file)
message("  ", csv_nodes)
message("  ", csv_edges)
message("  ", dagitty_file)
