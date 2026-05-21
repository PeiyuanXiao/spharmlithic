# ==============================================================================
# align_lin2024.R
# Scar orientation alignment — Lin 2024 pipeline
# ==============================================================================

#' Get scar lengths from a data frame
#'
#' Returns the length of each scar.  If a `Length` column already exists it
#' is used directly; otherwise lengths are computed from the start- and
#' end-point coordinates.
#'
#' @param df A data frame containing either a `Length` column or columns
#'   `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`.
#'
#' @return A numeric vector of scar lengths, one per row.
#'
#' @details
#' The `Length` column (if present) takes precedence over coordinate-derived
#' lengths.  This allows pre-computed lengths from external software to be
#' used directly without recalculation.
#'
#' @examples
#' \dontrun{
#' lens <- get_scar_length(df_one)
#' longest <- df_one[which.max(lens), ]
#' }
#'
#' @seealso [align_morph()]
#' @export
get_scar_length <- function(df) {
  if ("Length" %in% names(df)) return(df$Length)
  dx <- df$End_X - df$Start_X
  dy <- df$End_Y - df$Start_Y
  dz <- df$End_Z - df$Start_Z
  sqrt(dx^2 + dy^2 + dz^2)
}


#' Align scar orientation data using the Lin 2024 pipeline
#'
#' Implements the two-step alignment procedure described in Lin et al. (2024):
#' \enumerate{
#'   \item **Rotate** — uses the pre-measured morphological plane normal
#'     (`Norm_X/Y/Z`) to rotate the scar data so that the normal points along
#'     the Z-axis.
#'   \item **Translate** — shifts all coordinates so that the start-point of
#'     the longest scar lies at the XY origin.
#' }
#'
#' @param df_group A data frame for a **single** specimen.  Must contain
#'   columns `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`,
#'   `Norm_X`, `Norm_Y`, `Norm_Z`, and optionally `Length`.
#'
#' @return The input data frame with nine additional columns:
#' \describe{
#'   \item{s_x, s_y, s_z}{Aligned start-point coordinates.}
#'   \item{e_x, e_y, e_z}{Aligned end-point coordinates.}
#'   \item{d_x, d_y, d_z}{Aligned unit direction vectors.}
#' }
#'
#' @details
#' This is the **morphology-driven** alignment method: the plane normal is
#' supplied externally (e.g. measured in Geomagic) rather than estimated from
#' the scar data.  This makes the rotation step independent of scar
#' configuration, which can be advantageous for cores with few or poorly
#' distributed scars.
#'
#' For the data-driven alternative that estimates the plane normal from the
#' scars themselves, see [align_scar()].
#'
#' @references
#' Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., & Sutikna,
#' T. (2024). A new method for quantifying flake scar organisation on cores
#' using orientation statistics. \emph{Journal of Archaeological Science},
#' \strong{167}, 105998.
#'
#' @examples
#' \dontrun{
#' # Single specimen
#' df_one <- subset(raw_data, ID == "S001")
#' aligned_one <- align_morph(df_one)
#'
#' # All specimens via the batch wrapper
#' aligned_all <- align_morph_batch(raw_data)
#' }
#'
#' @seealso [align_scar()], [get_scar_length()], [align_morph_batch()]
#'
#' @importFrom dplyr mutate
#' @export
align_morph <- function(df_group) {
  
  # --- Compute unit direction vectors ---
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10
  
  df_group$Direct_X <- ifelse(valid, dx / len, 0)
  df_group$Direct_Y <- ifelse(valid, dy / len, 0)
  df_group$Direct_Z <- ifelse(valid, dz / len, 0)
  
  # --- Step 1: Rotate — align morphological plane normal to the Z-axis ---
  normal <- as.numeric(df_group[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal <- normal / sqrt(sum(normal^2))
  
  R1 <- get_rot_matrix(normal, c(0, 0, 1))
  
  S <- as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)
  
  df_group$s_x <- S[, 1]; df_group$s_y <- S[, 2]; df_group$s_z <- S[, 3]
  df_group$e_x <- E[, 1]; df_group$e_y <- E[, 2]; df_group$e_z <- E[, 3]
  df_group$d_x <- D[, 1]; df_group$d_y <- D[, 2]; df_group$d_z <- D[, 3]
  
  # --- Step 2: Translate — shift the start-point of the longest scar to (0, 0) ---
  longest_idx <- which.max(get_scar_length(df_group))
  shift_x     <- df_group$s_x[longest_idx]
  shift_y     <- df_group$s_y[longest_idx]
  
  df_group <- df_group %>%
    dplyr::mutate(
      s_x = s_x - shift_x,
      s_y = s_y - shift_y,
      e_x = e_x - shift_x,
      e_y = e_y - shift_y
    )
  
  return(df_group)
}


#' Build the three-panel Lin 2024 alignment visualisation for one specimen
#'
#' Reconstructs the two Lin 2024 alignment steps for specimen `demo_id` and
#' returns three Plotly figures (Step 0 – Step 2) for HTML embedding.  The
#' longest scar is highlighted in pink in all panels.
#'
#' @param demo_id A character or numeric specimen identifier present in
#'   `raw_data$ID`.
#' @param raw_data A data frame containing all specimens.  Must contain
#'   columns `ID`, `Start_X/Y/Z`, `End_X/Y/Z`, `Norm_X/Y/Z`, `Pos_X/Y/Z`,
#'   and optionally `Length`.
#'
#' @return A named list with elements `p0`, `p1`, `p2` (Plotly figures).
#'
#' @details
#' This function is primarily called by [export_alignment_html_lin2024()] to
#' assemble the full multi-specimen HTML report.  It can also be called
#' directly to inspect a single specimen's alignment interactively.
#'
#' @examples
#' \dontrun{
#' panels <- build_panel_morph("S001", raw_data)
#' panels$p0  # raw data
#' panels$p2  # fully aligned
#' }
#'
#' @seealso [align_morph()], [export_alignment_html_lin2024()]
#'
#' @importFrom plotly plot_ly
#' @importFrom dplyr filter
#' @export
build_panel_morph <- function(demo_id, raw_data) {
  df <- dplyr::filter(raw_data, .data$ID == demo_id)
  
  s0         <- as.matrix(df[, c("Start_X", "Start_Y", "Start_Z")])
  e0         <- as.matrix(df[, c("End_X",   "End_Y",   "End_Z"  )])
  normal_raw <- as.numeric(df[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  normal_raw <- normal_raw / sqrt(sum(normal_raw^2))
  center_raw <- as.numeric(df[1, c("Pos_X", "Pos_Y", "Pos_Z")])
  
  longest_idx <- which.max(get_scar_length(df))
  arr_scale   <- max(dist(s0)) * 0.25
  half_sz     <- max(dist(s0)) * 0.55
  
  R1        <- get_rot_matrix(normal_raw, c(0, 0, 1))
  s1        <- s0 %*% t(R1)
  e1        <- e0 %*% t(R1)
  center_r1 <- as.numeric(R1 %*% center_raw)
  
  shift_x   <- s1[longest_idx, 1]
  shift_y   <- s1[longest_idx, 2]
  s2        <- s1; e2 <- e1
  s2[, 1]   <- s1[, 1] - shift_x; s2[, 2] <- s1[, 2] - shift_y
  e2[, 1]   <- e1[, 1] - shift_x; e2[, 2] <- e1[, 2] - shift_y
  z_longest <- s2[longest_idx, 3]
  
  p0 <- plotly::plot_ly() %>%
    add_scars_3d(s0[,1], s0[,2], s0[,3],
                 e0[,1], e0[,2], e0[,3],
                 highlight_idx = longest_idx) %>%
    add_arrow_3d(center_raw, normal_raw, arr_scale) %>%
    add_tilted_plane_3d(center_raw, normal_raw, half_sz) %>%
    plotly::layout(panel_layout("<b>Step 0</b>: Raw data \u2014 arbitrary orientation"))
  
  p1 <- plotly::plot_ly() %>%
    add_scars_3d(s1[,1], s1[,2], s1[,3],
                 e1[,1], e1[,2], e1[,3],
                 highlight_idx = longest_idx) %>%
    add_arrow_3d(center_r1, c(0, 0, 1), arr_scale) %>%
    add_plane_3d(center_r1[1], center_r1[2], center_r1[3], half_sz) %>%
    plotly::layout(panel_layout("<b>Step 1</b>: Rotate \u2014 normal aligned to Z-axis"))
  
  p2 <- plotly::plot_ly() %>%
    add_scars_3d(s2[,1], s2[,2], s2[,3],
                 e2[,1], e2[,2], e2[,3],
                 highlight_idx = longest_idx) %>%
    add_arrow_3d(c(0, 0, z_longest), c(0, 0, 1), arr_scale) %>%
    add_plane_3d(0, 0, z_longest, half_sz) %>%
    plotly::layout(panel_layout(
      "<b>Step 2 (Lin 2024)</b>: Translate \u2014 longest scar start to (0,0)"))
  
  list(p0 = p0, p1 = p1, p2 = p2)
}


#' Export an interactive Lin 2024 alignment HTML page
#'
#' Builds three-panel Lin 2024 alignment visualisations for every specimen in
#' `raw_data` and writes a standalone interactive HTML file to `out_path`.
#'
#' @param raw_data A data frame with the required columns (see
#'   [build_panel_morph()]).
#' @param out_path Character. File path for the output HTML.
#'
#' @return Invisibly, the `out_path` string.
#'
#' @details
#' The exported HTML is fully self-contained (Plotly loaded from CDN) and
#' requires no R session to view.  All specimen panels are serialised to
#' JSON at export time; switching specimens in the browser is instant.
#'
#' For the SVD-based equivalent, see [export_alignment_html_svd()].
#'
#' @examples
#' \dontrun{
#' export_alignment_html_lin2024(raw_data, "output/alignment_lin2024.html")
#' }
#'
#' @seealso [build_panel_morph()], [export_alignment_html_svd()]
#'
#' @importFrom htmltools tagList tags browsable save_html HTML
#' @importFrom jsonlite toJSON
#' @export
export_alignment_html_lin2024 <- function(raw_data, out_path) {
  all_ids     <- unique(raw_data$ID)
  panels_list <- lapply(as.character(all_ids), build_panel_morph,
                        raw_data = raw_data)
  names(panels_list) <- as.character(all_ids)
  
  js_data_lines <- sapply(as.character(all_ids), function(id) {
    ps <- panels_list[[id]]
    paste0(
      'allPanels["', id, '"] = {',
      '"p0":', get_panel_json(ps$p0), ',',
      '"p1":', get_panel_json(ps$p1), ',',
      '"p2":', get_panel_json(ps$p2),
      '};'
    )
  })
  js_data_block <- paste(js_data_lines, collapse = "\n")
  ids_json      <- jsonlite::toJSON(as.character(all_ids), auto_unbox = FALSE)
  
  grid <- htmltools::browsable(
    htmltools::tagList(
      htmltools::tags$script(
        src = "https://cdn.plot.ly/plotly-2.27.0.min.js"),
      
      htmltools::tags$h3(
        style = "font-family:sans-serif; text-align:center; margin:16px 0 4px;",
        "Core Alignment Pipeline (Lin 2024)"
      ),
      
      htmltools::tags$div(
        style = "text-align:center; margin-bottom:10px;",
        htmltools::tags$label("Select specimen: ",
                              style = "font-family:sans-serif; font-size:13px;"),
        htmltools::tags$select(
          id    = "specimenSelect",
          style = "font-size:13px; padding:3px 8px;",
          lapply(as.character(all_ids), function(id)
            htmltools::tags$option(value = id, id))
        )
      ),
      
      htmltools::tags$p(
        style = "font-family:sans-serif; text-align:center; color:#666;
                 margin:0 0 12px; font-size:13px;",
        htmltools::HTML(
          "&#9642; <b style='color:steelblue'>Blue</b> = Flaking scars &nbsp;|&nbsp;
           <b style='color:pink'>Pink</b> = Longest scar &nbsp;|&nbsp;
           <b style='color:red'>Red arrow</b> = Plane normal &nbsp;|&nbsp;
           <b style='color:lightgray'>Gray plane</b> = Best fitting plane")
      ),
      
      htmltools::tags$div(
        style = "display:grid; grid-template-columns:1fr 1fr 1fr;
                 gap:8px; padding:0 12px 12px;",
        
        htmltools::tags$div(
          htmltools::tags$div(id = "plot0", style = "height:500px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "Compute unit direction vectors from scar endpoints.<br>
               The morphological plane normal (Norm_X/Y/Z) is measured<br>
               externally via Geomagic as the <b>reference orientation</b>.")
          )
        ),
        
        htmltools::tags$div(
          htmltools::tags$div(id = "plot1", style = "height:500px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "Construct rotation matrix R&#8321; via Rodrigues' formula.<br>
               Rotate all coordinates so the morphological normal aligns with <b>Z-axis</b>.<br>
               The best-fit plane now lies on the XY plane.")
          )
        ),
        
        htmltools::tags$div(
          htmltools::tags$div(id = "plot2", style = "height:500px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "Identify the longest scar and extract its rotated start-point.<br>
               Subtract its XY coordinates from all points to <b>anchor the origin</b>.<br>
               Preserves absolute spatial relationships between scars.")
          )
        )
      ),
      
      htmltools::tags$script(htmltools::HTML(paste0(
        "var allPanels = {};\n",
        "var allIds = ", ids_json, ";\n",
        js_data_block, "\n",
        "
        function renderSpecimen(id) {
          var ps = allPanels[id];
          ['p0','p1','p2'].forEach(function(k, i) {
            Plotly.react('plot' + i, ps[k].data, ps[k].layout);
          });
        }
        renderSpecimen(allIds[0]);
        document.getElementById('specimenSelect').addEventListener('change',
          function() { renderSpecimen(this.value); });
        "
      )))
    )
  )
  
  htmltools::save_html(grid, out_path)
  invisible(out_path)
}
