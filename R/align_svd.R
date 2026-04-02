# ==============================================================================
# align_svd.R
# 片疤方向对齐 — SVD 法
# ==============================================================================

#' Align scar orientation data using SVD (three-step pipeline)
#'
#' Performs the three-step SVD alignment pipeline on a single-specimen data
#' frame:
#' \enumerate{
#'   \item **Rotate** — SVD of the unit-direction matrix extracts the
#'     best-fit plane normal; rotate it onto the Z-axis via
#'     [get_rot_matrix()].
#'   \item **Translate** — shift the point-cloud centroid to the origin.
#'   \item **In-plane rotate** — SVD on the XY projections finds the main
#'     scar direction and rotates it onto the X-axis.
#' }
#'
#' @param df_group A data frame for a **single** specimen.  Must contain
#'   columns `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`.
#'
#' @return The input data frame with nine additional columns:
#' \describe{
#'   \item{s_x, s_y, s_z}{Aligned start-point coordinates.}
#'   \item{e_x, e_y, e_z}{Aligned end-point coordinates.}
#'   \item{d_x, d_y, d_z}{Aligned unit direction vectors.}
#' }
#'
#' @seealso [get_rot_matrix()], [align_morph()]
#'
#' @importFrom dplyr mutate
#' @export
align_scar <- function(df_group) {

  # --- Step 1: 旋转 — 将 SVD 法线对齐到 Z 轴 ---
  dx  <- df_group$End_X - df_group$Start_X
  dy  <- df_group$End_Y - df_group$Start_Y
  dz  <- df_group$End_Z - df_group$Start_Z
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10

  df_group$Direct_X <- ifelse(valid, dx / len, 0)
  df_group$Direct_Y <- ifelse(valid, dy / len, 0)
  df_group$Direct_Z <- ifelse(valid, dz / len, 0)

  if (sum(valid) >= 3) {
    U      <- cbind(dx[valid] / len[valid],
                    dy[valid] / len[valid],
                    dz[valid] / len[valid])
    normal <- svd(U)$v[, 3]
  }

  normal <- normal / sqrt(sum(normal^2))
  if (normal[3] < 0) normal <- -normal

  R1 <- get_rot_matrix(normal, c(0, 0, 1))

  S <- as.matrix(df_group[, c("Start_X", "Start_Y", "Start_Z")]) %*% t(R1)
  E <- as.matrix(df_group[, c("End_X",   "End_Y",   "End_Z"  )]) %*% t(R1)
  D <- as.matrix(df_group[, c("Direct_X","Direct_Y","Direct_Z")]) %*% t(R1)

  df_group$s_x <- S[, 1]; df_group$s_y <- S[, 2]; df_group$s_z <- S[, 3]
  df_group$e_x <- E[, 1]; df_group$e_y <- E[, 2]; df_group$e_z <- E[, 3]
  df_group$d_x <- D[, 1]; df_group$d_y <- D[, 2]; df_group$d_z <- D[, 3]

  # --- Step 2: 平移 — 将点云质心移至原点 ---
  global_center_x <- mean(c(df_group$s_x, df_group$e_x))
  global_center_y <- mean(c(df_group$s_y, df_group$e_y))
  global_center_z <- mean(c(df_group$s_z, df_group$e_z))

  df_group <- df_group %>%
    dplyr::mutate(
      s_x = s_x - global_center_x,
      s_y = s_y - global_center_y,
      s_z = s_z - global_center_z,
      e_x = e_x - global_center_x,
      e_y = e_y - global_center_y,
      e_z = e_z - global_center_z
    )

  # --- Step 3: XY 内旋转 — 将刮痕主方向对齐到 X 轴 ---
  xy_dirs  <- cbind(df_group$d_x, df_group$d_y)
  main_dir <- svd(xy_dirs)$v[, 1]
  mean_dir <- colMeans(xy_dirs)
  if (sum(main_dir * mean_dir) < 0) main_dir <- -main_dir

  theta <- atan2(main_dir[2], main_dir[1])
  R2    <- matrix(c( cos(-theta), -sin(-theta), 0,
                     sin(-theta),  cos(-theta), 0,
                     0,            0,           1), 3, 3, byrow = TRUE)

  S2 <- as.matrix(df_group[, c("s_x", "s_y", "s_z")]) %*% t(R2)
  E2 <- as.matrix(df_group[, c("e_x", "e_y", "e_z")]) %*% t(R2)
  D2 <- as.matrix(df_group[, c("d_x", "d_y", "d_z")]) %*% t(R2)

  df_group$s_x <- S2[, 1]; df_group$s_y <- S2[, 2]; df_group$s_z <- S2[, 3]
  df_group$e_x <- E2[, 1]; df_group$e_y <- E2[, 2]; df_group$e_z <- E2[, 3]
  df_group$d_x <- D2[, 1]; df_group$d_y <- D2[, 2]; df_group$d_z <- D2[, 3]

  return(df_group)
}


#' Build the four-panel SVD alignment visualisation for one specimen
#'
#' Reconstructs the three alignment steps for specimen `demo_id` and returns
#' four Plotly figures (Step 0 – Step 3) suitable for embedding in an
#' interactive HTML page.
#'
#' @param demo_id A character or numeric specimen identifier present in
#'   `raw_data$ID`.
#' @param raw_data A data frame containing all specimens.  Must contain
#'   columns `ID`, `Start_X/Y/Z`, `End_X/Y/Z`, and optionally `Norm_X/Y/Z`.
#'
#' @return A named list with elements `p0`, `p1`, `p2`, `p3` (Plotly figures).
#'
#' @seealso [export_alignment_html_svd()]
#'
#' @importFrom plotly plot_ly
#' @importFrom dplyr filter
#' @export
build_panel_scar <- function(demo_id, raw_data) {
  df <- dplyr::filter(raw_data, .data$ID == demo_id)

  s0 <- as.matrix(df[, c("Start_X", "Start_Y", "Start_Z")])
  e0 <- as.matrix(df[, c("End_X",   "End_Y",   "End_Z"  )])

  center_raw <- c(mean(c(s0[, 1], e0[, 1])),
                  mean(c(s0[, 2], e0[, 2])),
                  mean(c(s0[, 3], e0[, 3])))
  arr_scale  <- max(dist(s0)) * 0.25
  half_sz    <- max(dist(s0)) * 0.55

  dx  <- e0[, 1] - s0[, 1]
  dy  <- e0[, 2] - s0[, 2]
  dz  <- e0[, 3] - s0[, 3]
  len <- sqrt(dx^2 + dy^2 + dz^2)
  valid <- len > 1e-10

  if (sum(valid) >= 3) {
    U          <- cbind(dx[valid] / len[valid],
                        dy[valid] / len[valid],
                        dz[valid] / len[valid])
    normal_svd <- svd(U)$v[, 3]
  } else {
    normal_svd <- as.numeric(df[1, c("Norm_X", "Norm_Y", "Norm_Z")])
  }
  normal_svd <- normal_svd / sqrt(sum(normal_svd^2))
  if (normal_svd[3] < 0) normal_svd <- -normal_svd

  R1        <- get_rot_matrix(normal_svd, c(0, 0, 1))
  s1        <- s0 %*% t(R1)
  e1        <- e0 %*% t(R1)
  center_r1 <- c(mean(c(s1[, 1], e1[, 1])),
                 mean(c(s1[, 2], e1[, 2])),
                 mean(c(s1[, 3], e1[, 3])))

  s2 <- sweep(s1, 2, center_r1, "-")
  e2 <- sweep(e1, 2, center_r1, "-")

  d2       <- e2 - s2
  len2     <- sqrt(rowSums(d2^2))
  valid2   <- len2 > 1e-10
  xy_dirs  <- cbind(d2[valid2, 1], d2[valid2, 2])
  main_dir <- svd(xy_dirs)$v[, 1]
  mean_dir <- colMeans(xy_dirs)
  if (sum(main_dir * mean_dir) < 0) main_dir <- -main_dir
  theta    <- atan2(main_dir[2], main_dir[1])
  R2       <- matrix(c( cos(-theta), -sin(-theta), 0,
                        sin(-theta),  cos(-theta), 0,
                        0,            0,           1), 3, 3, byrow = TRUE)
  s3 <- s2 %*% t(R2)
  e3 <- e2 %*% t(R2)

  p0 <- plotly::plot_ly() %>%
    add_scars_3d(s0[,1], s0[,2], s0[,3], e0[,1], e0[,2], e0[,3]) %>%
    add_arrow_3d(center_raw, normal_svd, arr_scale) %>%
    add_tilted_plane_3d(center_raw, normal_svd, half_sz) %>%
    plotly::layout(panel_layout("<b>Step 0</b>: Raw data — SVD normal shown"))

  p1 <- plotly::plot_ly() %>%
    add_scars_3d(s1[,1], s1[,2], s1[,3], e1[,1], e1[,2], e1[,3]) %>%
    add_arrow_3d(center_r1, c(0, 0, 1), arr_scale) %>%
    add_plane_3d(center_r1[1], center_r1[2], center_r1[3], half_sz) %>%
    plotly::layout(panel_layout("<b>Step 1</b>: Rotate — SVD normal aligned to Z-axis"))

  p2 <- plotly::plot_ly() %>%
    add_scars_3d(s2[,1], s2[,2], s2[,3], e2[,1], e2[,2], e2[,3]) %>%
    add_arrow_3d(c(0, 0, 0), c(0, 0, 1), arr_scale) %>%
    add_plane_3d(0, 0, 0, half_sz) %>%
    plotly::layout(panel_layout("<b>Step 2</b>: Translate — center moved to origin"))

  p3 <- plotly::plot_ly() %>%
    add_scars_3d(s3[,1], s3[,2], s3[,3], e3[,1], e3[,2], e3[,3]) %>%
    add_arrow_3d(c(0, 0, 0), c(0, 0, 1), arr_scale) %>%
    add_arrow_3d(c(0, 0, 0), c(1, 0, 0), arr_scale, color = "orange") %>%
    add_plane_3d(0, 0, 0, half_sz) %>%
    plotly::layout(panel_layout("<b>Step 3</b>: Rotate XY — PCA main axis aligned to X-axis"))

  list(p0 = p0, p1 = p1, p2 = p2, p3 = p3)
}


#' Export an interactive SVD alignment HTML page
#'
#' Builds four-panel SVD alignment visualisations for every specimen in
#' `raw_data` and writes a standalone interactive HTML file to `out_path`.
#' The page contains a drop-down menu to switch between specimens.
#'
#' @param raw_data A data frame with columns `ID`, `Start_X/Y/Z`, `End_X/Y/Z`.
#' @param out_path Character. File path for the output HTML (e.g.
#'   `"output/scar_alignment_svd.html"`).
#'
#' @return Invisibly, the `out_path` string.
#'
#' @seealso [build_panel_scar()]
#'
#' @importFrom htmltools tagList tags browsable save_html HTML
#' @importFrom jsonlite toJSON
#' @export
export_alignment_html_svd <- function(raw_data, out_path) {
  all_ids     <- unique(raw_data$ID)
  panels_list <- lapply(as.character(all_ids), build_panel_scar,
                        raw_data = raw_data)
  names(panels_list) <- as.character(all_ids)

  js_data_lines <- sapply(as.character(all_ids), function(id) {
    ps <- panels_list[[id]]
    paste0(
      'allPanels["', id, '"] = {',
      '"p0":', get_panel_json(ps$p0), ',',
      '"p1":', get_panel_json(ps$p1), ',',
      '"p2":', get_panel_json(ps$p2), ',',
      '"p3":', get_panel_json(ps$p3),
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
        "Core Alignment Pipeline (SVD normal)"
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
           <b style='color:red'>Red arrow</b> = SVD normal &nbsp;|&nbsp;
           <b style='color:orange'>Orange arrow</b> = PCA main axis (X) &nbsp;|&nbsp;
           <b style='color:lightgray'>Gray plane</b> = SVD best-fit plane")
      ),

      htmltools::tags$div(
        style = "display:grid; grid-template-columns:1fr 1fr 1fr 1fr;
                 gap:8px; padding:0 12px 12px;",
        htmltools::tags$div(
          htmltools::tags$div(id = "plot0", style = "height:480px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "Compute unit direction vectors from scar endpoints.<br>
               SVD decomposition extracts the 3rd right singular vector<br>
               as the <b>best-fit plane normal</b> of the scar pattern.")
          )
        ),
        htmltools::tags$div(
          htmltools::tags$div(id = "plot1", style = "height:480px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "Construct rotation matrix R&#8321; via Rodrigues' formula.<br>
               Rotate all coordinates so the SVD normal aligns with <b>Z-axis</b>.<br>
               The best-fit plane now lies on the XY plane.")
          )
        ),
        htmltools::tags$div(
          htmltools::tags$div(id = "plot2", style = "height:480px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "Compute the centroid of all scar endpoints in rotated space.<br>
               Subtract centroid from all coordinates to <b>center the cloud at origin</b>.<br>
               Removes positional bias between specimens.")
          )
        ),
        htmltools::tags$div(
          htmltools::tags$div(id = "plot3", style = "height:480px;"),
          htmltools::tags$p(
            style = "font-family:sans-serif; font-size:12px; color:#444;
                     text-align:center; margin:4px 8px 12px;",
            htmltools::HTML(
              "SVD on XY projections of direction vectors finds the <b>main axis</b>.<br>
               Rotation matrix R&#8322; aligns this axis to the <b>X-axis</b>.<br>
               All specimens now share a standard orientation for comparison.")
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
          ['p0','p1','p2','p3'].forEach(function(k, i) {
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
