#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rsvg)
})

out_dir <- "analysis/strictK3_primary_v2/figures_true_dag_candidate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

basename <- "Supplement_Figure_S1_DAG_journal_minimal_CRRT_v3"

nodes <- data.frame(
  id = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "U1", "U2", "X", "Y"),
  label = c(
    "Age, sex,\ncomorbidity",
    "Sepsis and\nshock severity",
    "AKI severity and\nKDIGO3 pathway",
    "Organ support and\nhaemodynamic burden",
    "Biochemical derangement\nand rescue risk",
    "Care context and\nICU practice",
    "Clinician judgement,\ngoals of care",
    "Unrecorded volume status,\ntrajectory, reversibility",
    "CRRT initiation\nwithin 24 h",
    "Hospital-recorded\n28-day death"
  ),
  x = c(170, 455, 755, 455, 755, 170, 455, 755, 1125, 1500),
  y = c(315, 160, 300, 535, 575, 750, 880, 880, 445, 445),
  type = c(rep("Z", 6), rep("U", 2), "X", "Y"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c(
    # Structure among measured covariate domains.
    "Z1", "Z1", "Z1",
    "Z2", "Z2", "Z2",
    "Z3", "Z4",
    # Main measured backdoor paths into treatment.
    "Z2", "Z3", "Z4", "Z5", "Z6",
    # Main measured causes of outcome.
    "Z1", "Z2", "Z3", "Z4", "Z5", "Z6",
    # Residual/unmeasured confounding.
    "U1", "U1", "U2", "U2", "U2",
    # Causal contrast of interest.
    "X"
  ),
  to = c(
    "Z2", "Z3", "Z6",
    "Z3", "Z4", "Z5",
    "Z5", "Z5",
    "X", "X", "X", "X", "X",
    "Y", "Y", "Y", "Y", "Y", "Y",
    "X", "Y", "Z5", "X", "Y",
    "Y"
  ),
  kind = c(
    rep("background", 8),
    rep("measured", 11),
    rep("unmeasured", 5),
    "effect"
  ),
  curve = c(
    0.00, -0.05, 0.00,
    0.00, 0.10, -0.04,
    0.00, 0.00,
    0.12, 0.03, -0.03, -0.10, -0.18,
    0.24, 0.16, 0.08, -0.10, -0.05, -0.22,
    -0.10, -0.26, 0.00, -0.08, -0.24,
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

text_svg <- function(x, y, label, size = 19, weight = 400,
                     fill = "#111111", anchor = "middle") {
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  line_h <- size * 1.16
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
    Z = c(fill = "#F5F7F6", stroke = "#2F6F5E", text = "#111111", opacity = "1"),
    U = c(fill = "#FFFFFF", stroke = "#8B8B8B", text = "#555555", opacity = "1"),
    X = c(fill = "#FFF4D6", stroke = "#9A6A00", text = "#111111", opacity = "1"),
    Y = c(fill = "#EAF4FB", stroke = "#2C6F9E", text = "#111111", opacity = "1")
  )
}

edge_style <- function(kind) {
  switch(
    kind,
    background = c(stroke = "#8C96A3", width = "1.6", dash = "", marker = "arrowBg"),
    measured = c(stroke = "#222222", width = "2.0", dash = "", marker = "arrowMain"),
    unmeasured = c(stroke = "#707070", width = "1.8", dash = " stroke-dasharray=\"8 7\"", marker = "arrowU"),
    effect = c(stroke = "#111111", width = "3.0", dash = "", marker = "arrowMain")
  )
}

edge_point <- function(from_id, to_id, side = c("start", "end"), radius = 42) {
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
  st <- node_style(row$type)
  dash <- if (row$type == "U") " stroke-dasharray=\"6 5\"" else ""
  paste0(
    sprintf(
      "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"42\" fill=\"%s\" stroke=\"%s\" stroke-width=\"2.2\" opacity=\"%s\"%s/>",
      row$x, row$y, st[["fill"]], st[["stroke"]], st[["opacity"]], dash
    ),
    text_svg(row$x, row$y, row$id, size = if (row$type %in% c("X", "Y")) 31 else 26,
             weight = 500, fill = st[["text"]]),
    text_svg(row$x, row$y + 72, row$label, size = 18, fill = st[["text"]])
  )
}

marker <- function(id, fill) {
  paste0(
    "<marker id=\"", id, "\" markerWidth=\"12\" markerHeight=\"12\" refX=\"11\" refY=\"6\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">",
    "<path d=\"M 0 0 L 12 6 L 0 12 z\" fill=\"", fill, "\"/>",
    "</marker>"
  )
}

svg_width <- 1650
svg_height <- 1040

svg <- paste0(
  "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", svg_width, "\" height=\"", svg_height,
  "\" viewBox=\"0 0 ", svg_width, " ", svg_height, "\">",
  "<defs>",
  "<style><![CDATA[text{font-family:Arial, Helvetica, sans-serif;dominant-baseline:middle}]]></style>",
  marker("arrowBg", "#8C96A3"),
  marker("arrowMain", "#222222"),
  marker("arrowU", "#707070"),
  "</defs>",
  "<rect width=\"1650\" height=\"1040\" fill=\"#FFFFFF\"/>",
  paste(vapply(seq_len(nrow(edges)), function(i) edge_svg(edges[i, ]), character(1)), collapse = "\n"),
  paste(vapply(seq_len(nrow(nodes)), function(i) node_svg(nodes[i, ]), character(1)), collapse = "\n"),
  "<text x=\"105\" y=\"995\" text-anchor=\"start\" font-size=\"17\" fill=\"#555555\">Solid arrows indicate assumed measured causal paths; dashed arrows indicate unmeasured common causes.</text>",
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

message("Created journal-minimal CRRT DAG files:")
message("  ", svg_file)
message("  ", png_file)
message("  ", pdf_file)
message("  ", csv_nodes)
message("  ", csv_edges)
message("  ", dagitty_file)
