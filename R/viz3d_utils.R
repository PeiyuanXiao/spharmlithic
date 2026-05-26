# ==============================================================================
# viz3d_utils.R
# Plotly 3-D visualisation helpers
# ==============================================================================

#' Add scar line segments and arrow cones to a Plotly 3-D figure
#'
#' Draws each scar as a line segment with a directional cone at its tip.
#' An optional index allows one scar to be rendered in a highlight style
#' (thicker, pink) — used in the Lin 2024 pipeline to emphasise the longest
#' scar.
#'
#' @param fig A `plotly` figure object.
#' @param sx,sy,sz Numeric vectors of scar start-point coordinates.
#' @param ex,ey,ez Numeric vectors of scar end-point coordinates.
#' @param highlight_idx Integer or `NULL`. Index of the scar to highlight.
#'   When `NULL` (default) all scars are rendered uniformly in steel-blue.
#'
#' @return The modified `plotly` figure object.
#'
#' @details
#' Each scar is rendered as two Plotly traces: a `scatter3d` line for the
#' shaft and a `cone` trace for the arrowhead at the end-point.  The
#' highlighted scar (if any) uses a wider line (`width = 6` vs `2`) and a
#' larger cone (`sizeref = 0.08` vs `0.05`).
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig <- plot_ly()
#' fig <- add_scars_3d(fig, sx = c(0, 1), sy = c(0, 0), sz = c(0, 0),
#'                          ex = c(1, 2), ey = c(0, 1), ez = c(0, 0))
#' fig
#'
#' # Highlight the second scar (e.g. the longest)
#' fig2 <- plot_ly()
#' fig2 <- add_scars_3d(fig2, sx = c(0, 1), sy = c(0, 0), sz = c(0, 0),
#'                            ex = c(1, 2), ey = c(0, 1), ez = c(0, 0),
#'                            highlight_idx = 2)
#' fig2
#' }
#'
#' @importFrom plotly add_trace
#' @noRd
add_scars_3d <- function(fig, sx, sy, sz, ex, ey, ez,
                         highlight_idx = NULL) {
  for (i in seq_along(sx)) {
    is_hl <- !is.null(highlight_idx) && i == highlight_idx
    clr   <- if (is_hl) "pink"  else "steelblue"
    lwd   <- if (is_hl) 6      else 2
    csz   <- if (is_hl) 0.08   else 0.05
    
    dx <- ex[i] - sx[i]
    dy <- ey[i] - sy[i]
    dz <- ez[i] - sz[i]
    
    fig <- fig %>% plotly::add_trace(
      type = "scatter3d", mode = "lines",
      x = c(sx[i], ex[i]), y = c(sy[i], ey[i]), z = c(sz[i], ez[i]),
      line = list(color = clr, width = lwd),
      showlegend = FALSE
    )
    fig <- fig %>% plotly::add_trace(
      type       = "cone",
      x = list(ex[i]), y = list(ey[i]), z = list(ez[i]),
      u = list(dx),    v = list(dy),    w = list(dz),
      sizemode   = "scaled", sizeref = csz,
      colorscale = list(list(0, clr), list(1, clr)),
      cmin = 0, cmax = 1,
      showscale = FALSE, showlegend = FALSE,
      anchor    = "tip",
      lighting  = list(ambient = 0.9, diffuse = 0.5)
    )
  }
  fig
}


#' Add a directional arrow to a Plotly 3-D figure
#'
#' Draws a line from `origin` in `direction` (scaled to `scale` length) with
#' a cone arrowhead at the tip.
#'
#' @param fig A `plotly` figure object.
#' @param origin Numeric vector of length 3. Arrow base coordinates.
#' @param direction Numeric vector of length 3. Arrow direction (need not be
#'   normalised).
#' @param scale Numeric scalar. Arrow length.
#' @param color Character. Arrow colour (default `"red"`).
#'
#' @return The modified `plotly` figure object.
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig <- plot_ly()
#' fig <- add_arrow_3d(fig, origin = c(0, 0, 0), direction = c(0, 0, 1),
#'                     scale = 1, color = "red")
#' fig
#' }
#'
#' @importFrom plotly add_trace
#' @noRd
add_arrow_3d <- function(fig, origin, direction, scale, color = "red") {
  d   <- direction / sqrt(sum(direction^2)) * scale
  tip <- origin + d
  
  fig <- fig %>% plotly::add_trace(
    type = "scatter3d", mode = "lines",
    x = c(origin[1], tip[1]),
    y = c(origin[2], tip[2]),
    z = c(origin[3], tip[3]),
    line = list(color = color, width = 6),
    showlegend = FALSE
  )
  fig %>% plotly::add_trace(
    type       = "cone",
    x = list(tip[1]), y = list(tip[2]), z = list(tip[3]),
    u = list(d[1]),   v = list(d[2]),   w = list(d[3]),
    sizemode   = "scaled", sizeref = 0.2,
    colorscale = list(list(0, color), list(1, color)),
    cmin = 0, cmax = 1,
    showscale = FALSE, showlegend = FALSE,
    anchor    = "tail",
    lighting  = list(ambient = 0.9, diffuse = 0.5)
  )
}


#' Add a tilted plane mesh to a Plotly 3-D figure
#'
#' Draws a semi-transparent quadrilateral mesh centred at `center` with the
#' given `normal` direction. Used to visualise the best-fit plane of scar
#' orientation data before the rotation step.
#'
#' @param fig A `plotly` figure object.
#' @param center Numeric vector of length 3. Plane centre coordinates.
#' @param normal Numeric vector of length 3. Plane normal (need not be
#'   normalised).
#' @param half_size Numeric scalar. Half-side-length of the displayed square.
#'
#' @return The modified `plotly` figure object.
#'
#' @details
#' Two tangent vectors \eqn{\mathbf{u}} and \eqn{\mathbf{v}} spanning the
#' plane are constructed via Gram-Schmidt orthogonalisation against the
#' normalised `normal`.  The reference vector used for \eqn{\mathbf{u}} is
#' `c(1, 0, 0)` unless `normal` is nearly parallel to the X-axis, in which
#' case `c(0, 1, 0)` is used instead.  The four corners of the square are
#' then `center ± half_size * u ± half_size * v`.
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig <- plot_ly()
#' fig <- add_tilted_plane_3d(fig, center = c(0, 0, 0),
#'                            normal = c(0.5, 0.5, 1), half_size = 2)
#' fig
#' }
#'
#' @importFrom plotly add_trace
#' @noRd
add_tilted_plane_3d <- function(fig, center, normal, half_size) {
  n   <- normal / sqrt(sum(normal^2))
  ref <- if (abs(n[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
  u   <- ref - sum(ref * n) * n
  u   <- u / sqrt(sum(u^2))
  v   <- c(
    n[2] * u[3] - n[3] * u[2],
    n[3] * u[1] - n[1] * u[3],
    n[1] * u[2] - n[2] * u[1]
  )
  corners <- rbind(
    center + half_size * u + half_size * v,
    center - half_size * u + half_size * v,
    center - half_size * u - half_size * v,
    center + half_size * u - half_size * v
  )
  fig %>% plotly::add_trace(
    type = "mesh3d",
    x = corners[, 1], y = corners[, 2], z = corners[, 3],
    i = c(0, 0), j = c(1, 2), k = c(2, 3),
    intensity  = c(0, 0, 0, 0),
    colorscale = list(list(0, "lightgray"), list(1, "lightgray")),
    cmin = 0, cmax = 1,
    opacity = 0.35, flatshading = TRUE,
    showscale = FALSE, showlegend = FALSE
  )
}


#' Add a horizontal (XY-parallel) plane mesh to a Plotly 3-D figure
#'
#' Draws a semi-transparent flat plane at height `z0`, centred at (`cx`, `cy`).
#' Used to visualise the XY reference plane after the rotation step has aligned
#' the scar pattern normal to the Z-axis.
#'
#' @param fig A `plotly` figure object.
#' @param cx,cy Numeric scalars. Centre X and Y coordinates.
#' @param z0 Numeric scalar. Plane height (Z coordinate).
#' @param half_size Numeric scalar. Half-side-length of the square.
#'
#' @return The modified `plotly` figure object.
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig <- plot_ly()
#' fig <- add_plane_3d(fig, cx = 0, cy = 0, z0 = 0, half_size = 2)
#' fig
#' }
#'
#' @importFrom plotly add_trace
#' @noRd
add_plane_3d <- function(fig, cx, cy, z0, half_size) {
  h  <- half_size
  xs <- c(cx - h, cx + h, cx + h, cx - h)
  ys <- c(cy - h, cy - h, cy + h, cy + h)
  zs <- rep(z0, 4)
  fig %>% plotly::add_trace(
    type = "mesh3d",
    x = xs, y = ys, z = zs,
    i = c(0, 0), j = c(1, 2), k = c(2, 3),
    intensity  = c(0, 0, 0, 0),
    colorscale = list(list(0, "lightgray"), list(1, "lightgray")),
    cmin = 0, cmax = 1,
    opacity = 0.35, flatshading = TRUE,
    showscale = FALSE, showlegend = FALSE
  )
}


#' Default Plotly 3-D scene layout parameters
#'
#' Returns a list suitable for the `scene` element of `plotly::layout()`,
#' with sensible camera angle, axis titles, and equal aspect ratio.
#'
#' @return A named list for use in `plotly::layout(scene = make_scene())`.
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig <- plot_ly() %>%
#'   plotly::layout(scene = make_scene())
#' fig
#' }
#'
#' @seealso `panel_layout()`
#' @noRd
make_scene <- function() {
  list(
    camera     = list(eye = list(x = 1.6, y = 1.6, z = 1.1)),
    xaxis      = list(title = "X", showgrid = TRUE, zeroline = TRUE),
    yaxis      = list(title = "Y", showgrid = TRUE, zeroline = TRUE),
    zaxis      = list(title = "Z", showgrid = TRUE, zeroline = TRUE),
    aspectmode = "data"
  )
}


#' Plotly panel layout parameters
#'
#' Returns a `layout` list for a single Plotly panel, including a centred
#' title, tight margins, and a light background colour.
#'
#' @param title_text Character string. Panel title (HTML tags are supported,
#'   e.g. `"<b>Step 0</b>: Raw data"`).
#'
#' @return A named list for use in `plotly::layout()`.
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig <- plot_ly() %>%
#'   plotly::layout(panel_layout("<b>Step 0</b>: Raw data"))
#' fig
#' }
#'
#' @seealso `make_scene()`, `build_panel_scar()`, `build_panel_morph()`
#' @noRd
panel_layout <- function(title_text) {
  list(
    title         = list(text = title_text, font = list(size = 13),
                         x = 0.5, xanchor = "center"),
    margin        = list(t = 50, b = 5, l = 5, r = 5),
    paper_bgcolor = "#f5f7fa",
    scene         = make_scene()
  )
}


#' Serialise a Plotly figure to a JSON string for HTML embedding
#'
#' Builds the Plotly figure object and returns a JSON string containing
#' both the `data` and `layout` fields, ready to be embedded inside a
#' `<script>` block in a standalone HTML page.
#'
#' @param p A `plotly` figure object.
#'
#' @return A character string of JSON.
#'
#' @details
#' Standalone HTML export (see [export_alignment_html_svd()] and
#' [export_alignment_html_lin2024()]) serialises all specimen panels to
#' JSON at build time so that the resulting file requires no R session or
#' server to render.  This function handles the per-panel serialisation,
#' producing a compact `{"data": ..., "layout": ...}` string that is
#' assigned to a JavaScript variable and later consumed by `Plotly.react()`.
#'
#' @examples
#' \dontrun{
#' library(plotly)
#' fig  <- plot_ly() %>% add_scars_3d(0, 0, 0, 1, 0, 0)
#' json <- get_panel_json(fig)
#' cat(substr(json, 1, 120))
#' }
#'
#' @seealso [export_alignment_html_svd()], [export_alignment_html_lin2024()]
#' @importFrom plotly plotly_build
#' @importFrom jsonlite toJSON
#' @noRd
get_panel_json <- function(p) {
  built       <- plotly::plotly_build(p)
  data_json   <- jsonlite::toJSON(built$x$data,   auto_unbox = TRUE,
                                  null = "null", force = TRUE)
  layout_json <- jsonlite::toJSON(built$x$layout, auto_unbox = TRUE,
                                  null = "null", force = TRUE)
  paste0('{"data":', data_json, ',"layout":', layout_json, '}')
}