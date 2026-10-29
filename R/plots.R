# ==============================================================================
# plots.R
# Looking at the inputs and results of the scar-pattern indices:
#   - plot_core()            : one core in 3-D, with its outline and plane
#   - export_outline_html()  : every core with an outline, in one HTML page
#   - plot_random_interval() : raw values against the random interval
# ==============================================================================

#' Plot one core in 3-D with its outline
#'
#' Draws the scars of one core as arrows from start to end. If an outline is
#' given, also draws the outline, the fitted plane of intersection and its
#' centre, and colours each scar by its Centripetality Index (CPI): blue for
#' scars running towards the centre (180 degrees), grey along the edge (90),
#' red away from the centre (0). Use it to check that scars and outline share
#' a coordinate frame (the starts should lie on the outline) before reading
#' CPI and PRA.
#'
#' @param data A scar table. If it holds several cores, `id` selects one.
#' @param id The ID of the core to draw. May be omitted when `data` holds a
#'   single core.
#' @param outline Optional outline of the core (see
#'   [fit_intersection_plane()]) or a plane fitted with
#'   [fit_intersection_plane()].
#' @inheritParams spi_ei_batch
#'
#' @return A `plotly` figure.
#'
#' @examples
#' t <- seq(0, 2 * pi, length.out = 13)[-13]
#' start <- cbind(30 * cos(t), 30 * sin(t), 0)
#' end <- start * 0.4
#' end[, 3] <- -5
#' one <- data.frame(ID = "disc", Start_X = start[, 1], Start_Y = start[, 2],
#'                   Start_Z = start[, 3], End_X = end[, 1], End_Y = end[, 2],
#'                   End_Z = end[, 3])
#' u <- seq(0, 2 * pi, length.out = 200)
#' plot_core(one, outline = cbind(30 * cos(u), 30 * sin(u), 0))
#'
#' @seealso [export_outline_html()], [cpi_pra_batch()],
#'   [export_alignment_html_svd()]
#' @export
plot_core <- function(data, id = NULL, outline = NULL, id_col = "ID",
                      start_cols = c("Start_X", "Start_Y", "Start_Z"),
                      end_cols = c("End_X", "End_Y", "End_Z"), min_length = 0) {
  .check_scar_table(data, id_col, start_cols, end_cols, need_id = FALSE,
                    min_length = min_length)
  if (id_col %in% names(data)) {
    ids <- unique(data[[id_col]])
    if (is.null(id)) {
      if (length(ids) > 1L) {
        stop("`data` holds ", length(ids), " cores; choose one with `id`.", call. = FALSE)
      }
      id <- ids[1]
    }
    sel <- as.character(data[[id_col]]) == as.character(id)
    if (!any(sel)) stop("No core with ID ", id, " in `data`.", call. = FALSE)
    data <- data[sel, , drop = FALSE]
  }
  .core_figure(data, id, outline, start_cols, end_cols, min_length)$fig
}

#' Export every core and its outline to an interactive HTML page
#'
#' Writes a self-contained HTML page with one 3-D view per core that has an
#' outline (see [plot_core()]) and a drop-down menu to switch between them.
#' Each view states the core's frame check and flatness (see
#' [cpi_pra_batch()]), so that cores whose scars and outline do not line up
#' can be found quickly.
#'
#' @param data A scar table with all cores.
#' @param outlines A named list of outlines, as for [cpi_pra_batch()].
#' @param out_path Path of the HTML file to write.
#' @inheritParams spi_ei_batch
#'
#' @return Invisibly, `out_path`.
#'
#' @examples
#' \dontrun{
#' outlines <- read_outline_iges(list.files("outlines", full.names = TRUE))
#' export_outline_html(scars, outlines, "outline_check.html")
#' }
#'
#' @seealso [plot_core()], [cpi_pra_batch()]
#' @export
export_outline_html <- function(data, outlines, out_path, id_col = "ID",
                                start_cols = c("Start_X", "Start_Y", "Start_Z"),
                                end_cols = c("End_X", "End_Y", "End_Z"),
                                min_length = 0) {
  .check_scar_table(data, id_col, start_cols, end_cols, min_length = min_length)
  outlines <- .check_outlines(outlines)
  key <- as.character(data[[id_col]])
  ids <- intersect(unique(key), names(outlines))
  if (!length(ids)) stop("No outline matches a core ID.", call. = FALSE)
  panels <- vapply(ids, function(k) {
    f <- .core_figure(data[key == k, , drop = FALSE], k, outlines[[k]],
                      start_cols, end_cols, min_length)
    get_panel_json(f$fig)
  }, "")
  js <- paste0('allPanels[', vapply(ids, function(k) jsonlite::toJSON(k, auto_unbox = TRUE), ""),
               '] = ', panels, ';', collapse = "\n")
  ids_json <- jsonlite::toJSON(ids, auto_unbox = FALSE)

  page <- htmltools::browsable(htmltools::tagList(
    htmltools::tags$script(src = "https://cdn.plot.ly/plotly-2.27.0.min.js"),
    htmltools::tags$h3(style = "font-family:sans-serif; text-align:center; margin:16px 0 4px;",
                       "Cores, outlines and planes of intersection"),
    htmltools::tags$div(
      style = "text-align:center; margin-bottom:10px;",
      htmltools::tags$label("Select core: ", style = "font-family:sans-serif; font-size:13px;"),
      htmltools::tags$select(id = "coreSelect", style = "font-size:13px; padding:3px 8px;",
                             lapply(ids, function(k) htmltools::tags$option(value = k, k)))),
    htmltools::tags$p(
      style = "font-family:sans-serif; text-align:center; color:#666; margin:0 0 12px; font-size:13px;",
      htmltools::HTML(paste(
        "Scars coloured by CPI: <b style='color:#2166ac'>blue</b> towards the centre,",
        "<b style='color:#999999'>grey</b> along the edge, <b style='color:#b2182b'>red</b>",
        "away from it &nbsp;|&nbsp; black line = outline &nbsp;|&nbsp; gray plane = plane of",
        "intersection &nbsp;|&nbsp; black point = centre. Scar starts should lie on the outline."))),
    htmltools::tags$div(id = "plot0", style = "height:640px; padding:0 12px 12px;"),
    htmltools::tags$script(htmltools::HTML(paste0(
      "var allPanels = {};\n", "var allIds = ", ids_json, ";\n", js, "\n",
      "function renderCore(id) {\n",
      "  var p = allPanels[id];\n",
      "  Plotly.react('plot0', p.data, p.layout);\n",
      "}\n",
      "renderCore(allIds[0]);\n",
      "document.getElementById('coreSelect').addEventListener('change',\n",
      "  function() { renderCore(this.value); });\n")))))
  htmltools::save_html(page, out_path)
  invisible(out_path)
}

#' Plot raw index values against their random interval
#'
#' Plots the raw SPI, ER or IR of each core against its scar count, over the
#' random interval (the central 95% of values from cores with random scar
#' directions; see [random_interval()]). Cores outside the band depart from
#' what random directions produce at their scar count.
#'
#' @param x The output of [spi_ei_batch()], with the columns `n_scar` and
#'   `<parameter>_raw`.
#' @param parameter `"SPI"`, `"ER"` or `"IR"`.
#' @param group Optional name of a column of `x` to colour the cores by.
#'
#' @return A `ggplot` object (requires the ggplot2 package).
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   scars <- do.call(rbind, lapply(1:6, function(i) {
#'     n <- 4 + 3 * i
#'     s <- matrix(rnorm(3 * n, sd = 10), n, 3)
#'     d <- matrix(rnorm(3 * n), n, 3)
#'     d[, 1] <- d[, 1] + (i %% 2)
#'     data.frame(ID = paste0("core", i), Start_X = s[, 1], Start_Y = s[, 2],
#'                Start_Z = s[, 3], End_X = s[, 1] + d[, 1],
#'                End_Y = s[, 2] + d[, 2], End_Z = s[, 3] + d[, 3])
#'   }))
#'   fab <- spi_ei_batch(scars, rarefied = FALSE)
#'   plot_random_interval(fab, "SPI")
#' }
#'
#' @seealso [random_interval()], [spi_ei_batch()]
#' @export
plot_random_interval <- function(x, parameter = c("SPI", "ER", "IR"), group = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_random_interval() needs the ggplot2 package: install.packages(\"ggplot2\").",
         call. = FALSE)
  }
  parameter <- .match_ri_param(match.arg(toupper(parameter), .RI_PARAMS))
  col <- paste0(parameter, "_raw")
  if (!is.data.frame(x) || !all(c("n_scar", col) %in% names(x))) {
    stop("`x` must be the output of spi_ei_batch() with columns n_scar and ", col, ".",
         call. = FALSE)
  }
  if (!is.null(group) && !group %in% names(x)) {
    stop("`group` names no column of `x`: ", group, ".", call. = FALSE)
  }
  pts <- data.frame(n = x$n_scar, value = x[[col]])
  pts$group <- if (is.null(group)) "core" else as.character(x[[group]])
  pts <- pts[!is.na(pts$value) & pts$n >= 3, , drop = FALSE]
  if (!nrow(pts)) stop("No core has a raw ", parameter, " value.", call. = FALSE)
  band <- suppressMessages(random_interval(parameter, seq(3, max(pts$n))))
  band <- band[!is.na(band$lower), , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band,
                         ggplot2::aes(x = .data$n, ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::labs(x = "Number of scars", y = paste("Raw", parameter),
                  caption = "Grey band: central 95% of values from random scar directions")
  if (is.null(group)) {
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$n, y = .data$value))
  } else {
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$n, y = .data$value,
                                                          colour = .data$group)) +
      ggplot2::labs(colour = group)
  }
  p + ggplot2::theme_bw()
}

# ---- building the figure of one core -------------------------------------------------

# Scars coloured from red (CPI 0) through grey (90) to blue (180).
.cpi_colour <- function(cpi) {
  ramp <- grDevices::colorRamp(c("#b2182b", "#999999", "#2166ac"))
  out <- rep("steelblue", length(cpi))
  ok <- !is.na(cpi)
  if (any(ok)) {
    rgb <- ramp(pmin(pmax(cpi[ok], 0), 180) / 180)
    out[ok] <- grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
  }
  out
}

# The figure of one core (the rows of `data` are that core's scars).
.core_figure <- function(data, id, outline, start_cols, end_cols, min_length) {
  x <- .scar_vectors(as.matrix(data[, start_cols, drop = FALSE]),
                     as.matrix(data[, end_cols, drop = FALSE]),
                     min_length = min_length)
  if (!x$n) stop("Core ", id, " has no usable scars.", call. = FALSE)
  S <- x$S
  E <- S + x$v_hat * x$L
  fig <- plotly::plot_ly()
  title <- paste0("<b>", if (is.null(id)) "Core" else id, "</b>")
  if (is.null(outline)) {
    fig <- add_scars_3d(fig, S[, 1], S[, 2], S[, 3], E[, 1], E[, 2], E[, 3])
  } else {
    plane <- .plane_for(outline, if (is.null(id)) "" else id)
    ang <- .ring_angles(S, x$v_hat, plane$centre, plane$normal)
    fig <- .add_coloured_scars_3d(fig, S, E, .cpi_colour(ang$cpi),
                                  sprintf("scar %d<br>CPI %.1f<br>PRA %.1f",
                                          seq_len(x$n), ang$cpi, ang$pra))
    fig <- .add_outline_3d(fig, plane$loop, plane$centre)
    fig <- add_tilted_plane_3d(fig, plane$centre, plane$normal, plane$diagonal * 0.4)
    sd <- .start_dist(S, plane)
    title <- paste0(title, sprintf(
      " &nbsp; start_dist %.3f%s &nbsp; flatness %.3f &nbsp; mean CPI %.1f, mean PRA %.1f",
      sd, if (sd <= .FRAME_OK_FRAC) "" else " (frame check failed)",
      plane$flatness, mean(ang$cpi, na.rm = TRUE), mean(ang$pra)))
  }
  list(fig = apply_panel_layout(fig, title))
}

# Scars as lines with cones at their ends, each in its own colour.
.add_coloured_scars_3d <- function(fig, S, E, colours, text) {
  for (i in seq_len(nrow(S))) {
    d <- E[i, ] - S[i, ]
    fig <- fig %>% plotly::add_trace(
      type = "scatter3d", mode = "lines+markers",
      x = c(S[i, 1], E[i, 1]), y = c(S[i, 2], E[i, 2]), z = c(S[i, 3], E[i, 3]),
      line = list(color = colours[i], width = 4),
      marker = list(color = colours[i], size = c(3, 0.1)),
      text = text[i], hoverinfo = "text", showlegend = FALSE)
    fig <- fig %>% plotly::add_trace(
      type = "cone", x = list(E[i, 1]), y = list(E[i, 2]), z = list(E[i, 3]),
      u = list(d[1]), v = list(d[2]), w = list(d[3]),
      sizemode = "scaled", sizeref = 0.08,
      colorscale = list(list(0, colours[i]), list(1, colours[i])), cmin = 0, cmax = 1,
      showscale = FALSE, showlegend = FALSE, anchor = "tip", hoverinfo = "skip",
      lighting = list(ambient = 0.9, diffuse = 0.5))
  }
  fig
}

# The outline of one core and its centre after each step of an alignment:
# `steps` are functions of an n x 3 matrix applied in turn. Returns a list of
# list(loop, centre), raw first; NULL without an outline.
.outline_steps <- function(outline, steps) {
  if (is.null(outline)) return(NULL)
  plane <- fit_intersection_plane(outline)
  cur <- list(loop = plane$loop, centre = matrix(plane$centre, 1L))
  out <- list(cur)
  for (f in steps) {
    cur <- list(loop = f(cur$loop), centre = f(cur$centre))
    out[[length(out) + 1L]] <- cur
  }
  lapply(out, function(o) list(loop = o$loop, centre = as.vector(o$centre)))
}

# Draw step `k` of .outline_steps() (nothing without an outline).
.add_outline_step <- function(fig, steps, k) {
  if (is.null(steps)) return(fig)
  .add_outline_3d(fig, steps[[k]]$loop, steps[[k]]$centre)
}

# The outline as a black line and its centre as a black point.
.add_outline_3d <- function(fig, loop, centre = NULL) {
  fig <- fig %>% plotly::add_trace(
    type = "scatter3d", mode = "lines",
    x = loop[, 1], y = loop[, 2], z = loop[, 3],
    line = list(color = "black", width = 3), hoverinfo = "skip", showlegend = FALSE)
  if (!is.null(centre)) {
    fig <- fig %>% plotly::add_trace(
      type = "scatter3d", mode = "markers",
      x = centre[1], y = centre[2], z = centre[3],
      marker = list(color = "black", size = 5), text = "centre", hoverinfo = "text",
      showlegend = FALSE)
  }
  fig
}
